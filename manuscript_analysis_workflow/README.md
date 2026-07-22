# Manuscript analysis workflow

This top-level folder contains the code used for the maize-root heterosis manuscript, arranged in analysis order. It is independent of the repository's existing `scripts/` directory. No expression matrices, Seurat objects, virtual-knockout results, figures or other data are included.

## Scope

The workflow links longitudinal root programs, cortex and epidermis single-cell states, spatial transcriptomics, non-additive expression, cell-cell communication and scTenifoldKnk virtual perturbation. The final stages focus on the relationship among ROS regulation, non-additive genes and heterosis.

## Run order

1. `workflow/01_longitudinal_reference`: reconstruct meristematic and elongation/maturation programs and the Fig. 1B matrix.
2. `workflow/02_single_cell_cortex_epidermis`: re-cluster cortex and score cortex/epidermis developmental and non-additive programs.
3. `workflow/03_pseudotime_scmmd`: infer Monocle3 trajectories and compare cortex distributions with scMMD.
4. `workflow/04_spatial_validation`: score non-additive programs in cortex spatial layers.
5. `workflow/05_GO_and_cell_communication`: run Fig. 7 GO and PlantCellChat analyses.
6. `workflow/06_virtual_knockout`: prepare matrices, run standard or safe single KOs, run the double KO and summarize scTenifoldKnk output.
7. `workflow/07_ROS_nonadditive_integration`: perform Fisher, Wilcoxon/rank-biserial and GSEA analyses; define ROS-sensitive tiers; score spatial modules; and run candidate GO/KEGG enrichment.
8. `workflow/08_source_data`: assemble figure-linked machine-readable Source Data tables.

Set `MAIZE_PROJECT_ROOT` to the directory containing the `data/`, `results/`, `sctenifold_results/` and `ros_nonadditive_analysis/` inputs before running portable R and shell entry points. Historical scripts were minimally sanitized to remove machine-specific paths; their exact pre-sanitization provenance and checksums are recorded in `SOURCE_PROVENANCE.tsv`.

## Core packages

- R: Seurat, SeuratWrappers, monocle3, DESeq2, apeglm, data.table, fgsea, ComplexHeatmap, UCell, clusterProfiler, AnnotationHub, PlantCellChat, scTenifoldKnk, scTenifoldNet and zellkonverter.
- Python: scanpy, anndata, numpy, pandas and scmmd.

These are established packages for single-cell, spatial, differential-expression, gene-set and virtual-knockout analysis. Third-party package source trees are deliberately not vendored in this folder.

## Figure-to-code map

| Figure | Main source-data or analysis code |
| --- | --- |
| Fig. 1B | `01_longitudinal_reference/01_call_DEGs_and_cluster_programs.R`; `02_rebuild_Fig1B_source_matrix.R` |
| Figs. 1-6 | `02_single_cell_cortex_epidermis`, `03_pseudotime_scmmd`, `04_spatial_validation` |
| Fig. 7 | `05_GO_and_cell_communication` |
| Fig. 8A | `06_virtual_knockout/05_analyze_KO_downstream.R` |
| Fig. 8B-D | `07_ROS_nonadditive_integration/01_core_enrichment_and_tiers.R` |
| Fig. 8E | tier output from `01_core_enrichment_and_tiers.R` |
| Fig. 8F and 8H | `03_spatial_ROS_modules.R` |
| Fig. 8G | `04_candidate_functional_enrichment.R` |
| Source Data | `08_source_data/01_prepare_source_data.py` |

## Important interpretation note

`04_exploratory_module_NES_comparison.png` is a legacy filename. The plotted y value is `median_rank_biserial`, not NES. Each point is one KO summarized across its prespecified primary tissue strata, and the black crossbar is the within-hybrid/module median. Because KO modules were produced in different batches and with different runners, this panel is descriptive and must not be used as a formal module-ranking test.

## Reproducibility

ROS integration was recorded under R 4.5.1 with data.table 1.17.8, ggplot2 4.0.0, fgsea 1.34.2, ComplexHeatmap 2.26.1, Seurat 5.3.0 and UCell 2.12.0. Random seeds and non-default parameters are stated in the scripts. Source Data should be used to audit exact plotted values; large primary data and processed objects should be deposited separately in an appropriate data repository.
