#!/usr/bin/env python3
"""Prepare manuscript source-data CSVs and an auditable workbook manifest."""

from __future__ import annotations

import hashlib
import json
import shutil
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--project-root", required=True, type=Path, help="Root containing analysis results and input tables")
parser.add_argument("--output-dir", required=True, type=Path, help="Directory for machine-readable Source Data CSV files")
args = parser.parse_args()

ROOT = args.project_root.resolve()
OUT = args.output_dir.resolve()
FINAL_OUTPUT_DIR = OUT.parent
RESULTS = ROOT / "ros_nonadditive_analysis" / "results"
FIG1 = ROOT / "project_archive_work" / "source_inputs" / "fig1b"
FIG7_GO = ROOT / "project_archive_work" / "source_inputs" / "fig7_go"
FIG7_CC = ROOT / "project_archive_work" / "source_inputs" / "fig7_plantcellchat"
KO_ROOT = ROOT / "sctenifold_downstream_analysis"

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def write_csv(data: pd.DataFrame, filename: str) -> Path:
    path = OUT / filename
    data.to_csv(path, index=False, encoding="utf-8-sig")
    return path

def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".txt":
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)

OUT.mkdir(parents=True, exist_ok=True)
catalog: list[dict[str, object]] = []

def register(
    workbook: str,
    sheet: str,
    source: Path,
    output_name: str,
    description: str,
    derivation: str = "Copied without row filtering; UTF-8 CSV export only.",
    data: pd.DataFrame | None = None,
) -> None:
    table = read_table(source) if data is None else data
    output = write_csv(table, output_name)
    catalog.append(
        {
            "workbook": workbook,
            "sheet": sheet,
            "csv_file": output.name,
            "description": description,
            "derivation": derivation,
            "rows": len(table),
            "columns": len(table.columns),
            "source_path": str(source.relative_to(ROOT)),
            "source_sha256": sha256(source),
            "output_sha256": sha256(output),
        }
    )

fig1_book = "Source Data Fig. 1.xlsx"
register(fig1_book, "Fig1B sample annotation", FIG1 / "Fig1B_sample_position_annotation.csv",
         "Fig1B_sample_position_annotation.csv",
         "Sample-to-position, replicate and longitudinal-zone annotation for the 32 LCM RNA-seq libraries.")
register(fig1_book, "Fig1B VST position means", FIG1 / "Fig1B_position_mean_VST_matrix.csv",
         "Fig1B_position_mean_VST_matrix.csv",
         "DESeq2 variance-stabilized expression averaged across two replicates at each of 16 root positions.")
register(fig1_book, "Fig1B heatmap Z matrix", FIG1 / "Fig1B_heatmap_gene_Zscore_matrix.csv",
         "Fig1B_heatmap_gene_Zscore_matrix.csv",
         "Per-gene Z-scored 16-position matrix used for the Fig. 1B heatmap.")
register(fig1_book, "Fig1B gene clusters", FIG1 / "Fig1B_gene_cluster_assignments.csv",
         "Fig1B_gene_cluster_assignments.csv",
         "K-means cluster assignment for each gene displayed in Fig. 1B.")
register(fig1_book, "Fig1B M EZ gene sets", FIG1 / "Fig1B_M_EZ_module_gene_sets.csv",
         "Fig1B_M_EZ_module_gene_sets.csv",
         "Gene lists assigned to meristematic-zone and elongation-zone developmental programs.")
register(fig1_book, "Fig1B rebuild QC", FIG1 / "Fig1B_rebuild_QC.csv",
         "Fig1B_rebuild_QC.csv",
         "Reconstruction dimensions and non-finite-value quality checks.")

fig7_book = "Source Data Fig. 7.xlsx"
for direction in ("up", "down"):
    for tissue, label in (("cortex", "Cortex"), ("epi", "Epidermis")):
        path = FIG7_GO / f"GO_res_{direction}_{tissue}.txt"
        register(fig7_book, f"GO {direction} {label}", path,
                 f"Fig7_GO_{direction}_{tissue}.csv",
                 f"Complete GO over-representation results for {direction}-regulated non-additive genes in {label}.")
for hybrid in ("B73xA554", "B73xOh40B"):
    for subset, label in (("all_CE", "all cortex-epidermis"), ("spatial_fish", "spatial/FISH-supported")):
        path = FIG7_CC / f"{hybrid}_{subset}_communications.csv"
        register(fig7_book, f"{hybrid[:8]} {subset[:10]}", path,
                 f"Fig7_{hybrid}_{subset}_communications.csv",
                 f"PlantCellChat communication table for {hybrid}, {label} interactions.")

fig8_book = "Source Data Fig. 8.xlsx"
for hybrid in ("B73xA554", "B73xOh40B"):
    base = KO_ROOT / f"{hybrid}_sctenifold_KO_analysis" / "combined_comparison"
    for filename, short, description in (
        ("all_KO_sig_genes_long.csv", "sig genes", "Long-format significant downstream perturbation genes for all KO runs."),
        ("KO_sig_gene_counts.csv", "sig counts", "Number of significant downstream genes per KO."),
        ("pairwise_sig_gene_overlap_heatmap.csv", "overlap", "Pairwise counts of shared significant downstream genes used for Fig. 8A."),
        ("sig_gene_presence_matrix.csv", "presence", "Binary KO-by-gene significance presence matrix underlying the overlap counts."),
    ):
        register(fig8_book, f"8A {hybrid[:8]} {short}", base / filename,
                 f"Fig8A_{hybrid}_{filename}", description)

primary_path = RESULTS / "04_ros_primary_tissue_enrichment.csv"
primary = pd.read_csv(primary_path)
register(fig8_book, "8C ROS primary stats", primary_path,
         "Fig8C_ROS_primary_tissue_enrichment_all_sets.csv",
         "Fisher, Wilcoxon/rank-biserial and GSEA statistics for ROS KOs in the predefined primary tissues.", data=primary)
gsea_plot = primary.loc[primary["set_type"].eq("all")].copy()
register(fig8_book, "8C GSEA plotted data", primary_path,
         "Fig8C_GSEA_all_NAG_plotted_data.csv",
         "All-NAG rows supplying NES, adjusted GSEA P value and leading-edge genes for Fig. 8C.",
         "Filtered to set_type == 'all'; no additional statistical recalculation.", gsea_plot)

tiers_path = RESULTS / "05_ros_candidate_evidence_tiers.csv"
tiers = pd.read_csv(tiers_path)
register(fig8_book, "8E tier evidence", tiers_path,
         "Fig8E_ROS_candidate_evidence_tiers.csv",
         "Per-gene cross-hybrid ROS perturbation, non-additivity and leading-edge evidence with Tier assignment.", data=tiers)
tier_lists = tiers.loc[tiers["Tier"].isin(["Tier 1", "Tier 2", "Tier 3"]), ["gene", "Tier"]].copy()
register(fig8_book, "8E Tier gene lists", tiers_path,
         "Fig8E_Tier1_Tier2_Tier3_gene_lists.csv",
         "Explicit gene membership for Tier 1, Tier 2 and Tier 3.",
         "Filtered to Tier 1/2/3 and retained gene and Tier columns.", tier_lists)

def split_leading_edge(value: object) -> list[str]:
    if pd.isna(value) or not str(value).strip():
        return []
    return [gene.strip() for gene in str(value).split(";") if gene.strip()]

le_rows: list[tuple[str, str]] = []
epi_primary = primary.loc[(primary["Hybrid"] == "B73xA554") & (primary["set_type"] == "all")]
for ko, frame in epi_primary.groupby("KO", sort=True):
    for gene in set(sum((split_leading_edge(v) for v in frame["leading_edge"]), [])):
        le_rows.append((str(ko), gene))
le_table = pd.DataFrame(le_rows, columns=["KO", "gene"]).drop_duplicates()
recurrent = le_table.groupby("gene")["KO"].nunique()
recurrent_genes = set(recurrent.loc[recurrent >= 3].index)

module_rules = {
    "ROS_Tier1_shared_core": set(tiers.loc[tiers["Tier"].eq("Tier 1"), "gene"]),
    "ROS_Tier1_Tier2_consensus": set(tiers.loc[tiers["Tier"].isin(["Tier 1", "Tier 2"]), "gene"]),
    "ROS_Epidermis_recurrent_LE": recurrent_genes,
    "ROS_sensitive_positive_MPH": set(tiers.loc[tiers["Tier"].isin(["Tier 1", "Tier 2"]) & (tiers["n_positive_nonadd_B73xA554"] > 0), "gene"]),
    "ROS_sensitive_negative_MPH": set(tiers.loc[tiers["Tier"].isin(["Tier 1", "Tier 2"]) & (tiers["n_negative_nonadd_B73xA554"] > 0), "gene"]),
}
membership = pd.DataFrame(
    [(module, gene) for module, genes in module_rules.items() for gene in sorted(genes)],
    columns=["module", "gene"],
)
register(fig8_book, "8F module gene lists", tiers_path,
         "Fig8F_spatial_NAG_module_gene_lists.csv",
         "Flat membership table for the five ROS–non-additive gene modules scored in spatial transcriptomics.",
         "Reconstructed exactly from the published tier and recurrent leading-edge rules in 03_spatial_ros_nonadditive.R.", membership)
for filename, sheet, description in (
    ("16_spatial_module_gene_manifest.csv", "8F module coverage", "Requested and detected genes for each UCell module in the spatial object."),
    ("19_spatial_layer_delta_by_slice.csv", "8F slice layer delta", "Per-slice layer score summaries and delta from the within-slice median."),
    ("20_spatial_layer_consistency_descriptive.csv", "8F heatmap data", "Across-slice median delta, direction counts and sign-test statistics used for the Fig. 8F heatmap."),
):
    register(fig8_book, sheet, RESULTS / filename, f"Fig8F_{filename}", description)

go_path = RESULTS / "23_candidate_GO_enrichment.csv"
go = pd.read_csv(go_path)
register(fig8_book, "8G GO all modules", go_path,
         "Fig8G_candidate_GO_enrichment_all_modules.csv",
         "Complete GO enrichment results for all five ROS–non-additive candidate modules.", data=go)
go_fig = go.loc[go["module"].eq("ROS_Epidermis_recurrent_LE")].copy()
register(fig8_book, "8G GO plotted module", go_path,
         "Fig8G_ROS_Epidermis_recurrent_LE_GO.csv",
         "GO enrichment results for the recurrent epidermal leading-edge module displayed in Fig. 8G.",
         "Filtered to module == 'ROS_Epidermis_recurrent_LE'.", go_fig)
register(fig8_book, "8G annotation coverage", RESULTS / "25_functional_enrichment_coverage.csv",
         "Fig8G_functional_enrichment_annotation_coverage.csv",
         "Gene-set sizes, annotation coverage and enrichment-background definition.")

explore_book = "Source Data Exploratory KO module comparison.xlsx"
register(explore_book, "Per KO plotted values", RESULTS / "08_module_comparison_per_ko.csv",
         "Exploratory_KO_module_per_KO_plotted_values.csv",
         "Exact points plotted in 04_exploratory_module_NES_comparison.png; y is median_rank_biserial, not NES.")
register(explore_book, "Module summaries", RESULTS / "09_module_comparison_summary_exploratory.csv",
         "Exploratory_KO_module_summary_crossbars.csv",
         "Module-level descriptive medians and ranges; black crossbars use median_rank_biserial.")
all_enrichment_path = RESULTS / "02_ko_nonadditive_enrichment_all.csv"
all_enrichment = pd.read_csv(all_enrichment_path)
primary_all = all_enrichment.loc[(all_enrichment["analysis_family"] == "primary") & (all_enrichment["set_type"] == "all")].copy()
register(explore_book, "Primary stratum inputs", all_enrichment_path,
         "Exploratory_KO_module_primary_stratum_inputs.csv",
         "Primary-tissue all-NAG rows aggregated to one plotted value per KO.",
         "Filtered to analysis_family == 'primary' and set_type == 'all'.", primary_all)
register(explore_book, "KO input manifest", RESULTS / "01_input_manifest.csv",
         "Exploratory_KO_module_input_manifest.csv",
         "KO module, runner and tested/significant-gene audit table.")

catalog_df = pd.DataFrame(catalog)
catalog_path = write_csv(catalog_df, "Source_Data_catalog.csv")

checks: list[dict[str, object]] = []
def check(name: str, passed: bool, observed: object, expected: object) -> None:
    checks.append({"check": name, "passed": bool(passed), "observed": observed, "expected": expected})

zcols = [f"S{i}" for i in range(1, 17)]
z = pd.read_csv(FIG1 / "Fig1B_heatmap_gene_Zscore_matrix.csv")
check("Fig1B heatmap dimensions", z.shape == (10139, 18), str(z.shape), "(10139, 18 incl gene+cluster)")
check("Fig1B finite Z values", np.isfinite(z[zcols].to_numpy()).all(), int((~np.isfinite(z[zcols].to_numpy())).sum()), 0)
check("Fig8E unique genes", tiers["gene"].is_unique, int(tiers["gene"].duplicated().sum()), 0)
check("Fig8E explicit Tier list count", len(tier_lists) == int(tiers["Tier"].isin(["Tier 1", "Tier 2", "Tier 3"]).sum()), len(tier_lists), "matches evidence table")
manifest = pd.read_csv(RESULTS / "16_spatial_module_gene_manifest.csv")
membership_counts = membership.groupby("module")["gene"].nunique()
expected_counts = manifest.set_index("module")["n_requested"]
check("Fig8F reconstructed module sizes", membership_counts.reindex(expected_counts.index).equals(expected_counts.astype("int64")), membership_counts.to_dict(), expected_counts.to_dict())
heat = pd.read_csv(RESULTS / "20_spatial_layer_consistency_descriptive.csv")
check("Fig8F heatmap complete module-layer grid", len(heat) == 55, len(heat), 55)
check("Exploratory plot one row per Hybrid-Module-KO", not primary_all.empty and pd.read_csv(RESULTS / "08_module_comparison_per_ko.csv").duplicated(["Hybrid", "Module", "KO"]).sum() == 0, "checked", "no duplicate keys")

checks_df = pd.DataFrame(checks)
checks_path = write_csv(checks_df, "Source_Data_validation.csv")
if not checks_df["passed"].all():
    raise RuntimeError("Source-data validation failed:\n" + checks_df.loc[~checks_df["passed"]].to_string(index=False))

config = {
    "output_dir": str(FINAL_OUTPUT_DIR),
    "csv_dir": str(OUT.resolve()),
    "workbooks": [
        {
            "filename": workbook,
            "sheets": [
                {"name": row.sheet, "csv": row.csv_file, "description": row.description,
                 "derivation": row.derivation, "source_path": row.source_path}
                for row in catalog_df.loc[catalog_df["workbook"].eq(workbook)].itertuples()
            ],
        }
        for workbook in (fig1_book, fig7_book, fig8_book, explore_book)
    ],
}
(OUT / "workbook_config.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
print(catalog_df.groupby("workbook").agg(sheets=("sheet", "count"), rows=("rows", "sum")))
print(checks_df.to_string(index=False))
