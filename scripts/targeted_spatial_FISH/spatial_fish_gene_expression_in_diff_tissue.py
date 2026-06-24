import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
import math
import os
import gc
from scipy import ndimage as ndi

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--gene_file", type=str)
args = parser.parse_args()

if args.gene_file:
    with open(args.gene_file) as f:
        genes = [x.strip() for x in f if x.strip()]



# ============================================================
# SETTINGS
# ============================================================

gene_to_cell_1 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/B73/05_expression/gene_to_cell.csv"
mask_1 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/B73/04_cell_segment/mask_expand_B73.tif"
gene_to_cell_2 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/B73_Oh40B/05_expression/gene_to_cell.csv"
mask_2 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/B73_Oh40B/04_cell_segment/mask_expand.tif"
gene_to_cell_3 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/Oh40B/05_expression/gene_to_cell.csv"
mask_3 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/Oh40B/04_cell_segment/mask_expand_Oh40B.tif"
gene_to_cell_4 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/B73_A554/05_expression/gene_to_cell.csv"
mask_4 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/B73_A554/04_cell_segment/mask_expand.tif"
gene_to_cell_5 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/A554/05_expression/gene_to_cell.csv"
mask_5 = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/A554/04_cell_segment/mask_expand.tif"


datasets = [
    {"name": "B73", "csv": gene_to_cell_1, "mask": mask_1},
    {"name": "B73_Oh40B", "csv": gene_to_cell_2, "mask": mask_2},
    {"name": "Oh40B", "csv": gene_to_cell_3, "mask": mask_3},
    {"name": "B73_A554", "csv": gene_to_cell_4, "mask": mask_4},
    {"name": "A554", "csv": gene_to_cell_5, "mask": mask_5},
]

output_dir = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/cropped_gene_expression"
os.makedirs(output_dir, exist_ok=True)

# mask meaning:
# 229 = grey cell interiors
# 0   = black boundaries
# 255 = outside tissue
CELL_VAL = 229
LINE_VAL = 0
BG_VAL = 255

# choose aggregation per segmented cell: "median", "mean", or "max"
aggregation = "mean"

# shared scale mode
# recommended:
global_vmin = 0.0
global_quantile_for_vmax = 0.99

# appearance
default_unmapped_cell_color = [0.85, 0.85, 0.85]
boundary_color = [0.0, 0.0, 0.0]
background_color = [1.0, 1.0, 1.0]

# crop tissue automatically
auto_crop = False

# colormap
cmap = LinearSegmentedColormap.from_list(
    "exprmap",
    ["#ffffff", "#ffe5e5", "#fcae91", "#fb6a4a", "#cb181d"]
)



# ============================================================
# HELPERS
# ============================================================
def aggregate_by_label(df: pd.DataFrame, method: str) -> pd.Series:
    grouped = df.groupby("label")["expr"]
    if method == "median":
        return grouped.median()
    if method == "mean":
        return grouped.mean()
    if method == "max":
        return grouped.max()
    raise ValueError(f"Unsupported aggregation method: {method}")


def process_dataset(csv_file: str, mask_file: str, gene: str) -> dict:
    df = pd.read_csv(csv_file)
    mask = tifffile.imread(mask_file)

    # Parse coordinates from "(y, x)"
    coords = df["cell"].str.extract(r"\(([^,]+),([^)]+)\)").astype(float)
    df["y"] = coords[0]
    df["x"] = coords[1]

    cell_mask = mask == CELL_VAL
    line_mask = mask == LINE_VAL
    bg_mask = mask == BG_VAL

    # label segmented cells
    cell_labels, n_cells = ndi.label(cell_mask)

    # map coordinates to label ids
    rows = np.clip(np.round(df["y"]).astype(int), 0, mask.shape[0] - 1)
    cols = np.clip(np.round(df["x"]).astype(int), 0, mask.shape[1] - 1)

    labels = cell_labels[rows, cols]

    # snap points on boundaries/background to nearest cell
    nearest_idx = ndi.distance_transform_edt(
        cell_labels == 0,
        return_distances=False,
        return_indices=True
    )
    nearest_label_map = cell_labels[tuple(nearest_idx)]

    zero_mask = labels == 0
    labels[zero_mask] = nearest_label_map[rows[zero_mask], cols[zero_mask]]

    df["label"] = labels
    df = df[df["label"] > 0].copy()

    # gene expression
    df["expr"] = np.log1p(pd.to_numeric(df[gene], errors="coerce"))
    df = df.dropna(subset=["expr"])

    # aggregate per segmented cell
    label_expr = aggregate_by_label(df, aggregation)

    # paint image at pixel level
    paint = np.full(mask.shape, np.nan, dtype=float)
    for lab, val in label_expr.items():
        paint[cell_labels == lab] = val

    return {
        "df": df,
        "mask": mask,
        "cell_mask": cell_mask,
        "line_mask": line_mask,
        "bg_mask": bg_mask,
        "cell_labels": cell_labels,
        "paint": paint,
        "n_cells": n_cells,
    }


def crop_to_tissue(rgb: np.ndarray, cell_mask: np.ndarray, line_mask: np.ndarray):
    tissue = cell_mask | line_mask
    coords = np.where(tissue)
    if len(coords[0]) == 0:
        return rgb
    r0, r1 = coords[0].min(), coords[0].max()
    c0, c1 = coords[1].min(), coords[1].max()
    return rgb[r0:r1 + 1, c0:c1 + 1]


def build_rgb(result: dict, norm: Normalize, cmap) -> np.ndarray:
    mask = result["mask"]
    cell_mask = result["cell_mask"]
    line_mask = result["line_mask"]
    bg_mask = result["bg_mask"]
    paint = result["paint"]

    rgb = np.ones((*mask.shape, 3), dtype=float)
    rgb[bg_mask] = background_color
    rgb[line_mask] = boundary_color
    rgb[cell_mask] = default_unmapped_cell_color

    valid = cell_mask & ~np.isnan(paint)
    clipped = np.clip(paint, norm.vmin, norm.vmax)
    rgb[valid] = cmap(norm(clipped[valid]))[:, :3]

    if auto_crop:
        rgb = crop_to_tissue(rgb, cell_mask, line_mask)

    return rgb


def save_single_plot(rgb: np.ndarray, name: str, gene: str, norm: Normalize, cmap):
    fig, ax = plt.subplots(figsize=(5, 10))
    ax.imshow(rgb, origin="upper")
    ax.set_title(f"{name} | {gene}")
    ax.axis("off")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("log1p(expression)")

    plt.tight_layout()
    fig.savefig(f"{name}_{gene}_shared_scale.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{name}_{gene}_shared_scale.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_multi_panel(results, rgbs, gene, norm, cmap, output_dir):
    n = len(rgbs)
    ncols = min(5, n)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4.2 * ncols + 1.2, 8 * nrows),
        constrained_layout=False
    )

    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = [axes]

    for ax, res, rgb in zip(axes, results, rgbs):
        ax.imshow(rgb, origin="upper")
        ax.set_title(res["name"], fontsize=14)
        ax.axis("off")

    for ax in axes[len(rgbs):]:
        ax.axis("off")

    # Manually reserve space on the right for the colorbar
    fig.subplots_adjust(
        left=0.03,
        right=0.88,
        top=0.93,
        bottom=0.03,
        wspace=0.18,
        hspace=0.10
    )

    # Dedicated colorbar axis: [left, bottom, width, height]
    cax = fig.add_axes([0.90, 0.22, 0.018, 0.56])

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("log1p(expression)", fontsize=12)

    fig.suptitle(f"{gene} (shared color scale)", y=0.97, fontsize=16)

    pdf_path = os.path.join(output_dir, f"{gene}_multi_panel_shared_scale.pdf")
    print("Saving PDF to:", pdf_path)
    fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print("Saved exists:", os.path.exists(pdf_path))
    plt.close(fig)

    

def process_gene_from_preprocessed(pre, gene, aggregation="mean"):
    df = pre["df"].copy()

    if gene not in df.columns:
        raise KeyError(f"{gene} not found in dataset {pre['name']}")

    df["expr"] = np.log1p(pd.to_numeric(df[gene], errors="coerce"))
    df = df.dropna(subset=["expr"])

    label_expr = aggregate_by_label(df, aggregation)

    paint = np.full(pre["mask"].shape, np.nan, dtype=float)
    for lab, val in label_expr.items():
        paint[pre["cell_labels"] == lab] = val

    return {
        "name": pre["name"],
        "mask": pre["mask"],
        "cell_mask": pre["cell_mask"],
        "line_mask": pre["line_mask"],
        "bg_mask": pre["bg_mask"],
        "cell_labels": pre["cell_labels"],
        "paint": paint,
    }


# ============================================================
# 1) FAST PRECOMPUTATION PER SAMPLE
# ============================================================
def precompute_geometry(pre):
    """
    Adds to pre:
      - geom: per-cell geometry table
      - outer_labels: set of outer-boundary cell labels
      - external_bg: external background mask
      - marker_cache: empty dict for marker-expression caching
    """
    cell_labels = pre["cell_labels"]
    cell_mask = pre["cell_mask"]
    line_mask = pre["line_mask"]
    bg_mask = pre["bg_mask"]

    labels = np.unique(cell_labels)
    labels = labels[labels > 0]

    # Distance to outside
    dist_out = ndi.distance_transform_edt(cell_mask)

    # Tissue centerline: midpoint of tissue width for each row
    tissue = cell_mask | line_mask
    centers = np.full(tissue.shape[0], np.nan, dtype=np.float32)
    for r in range(tissue.shape[0]):
        cols = np.where(tissue[r])[0]
        if len(cols) > 0:
            centers[r] = (cols.min() + cols.max()) / 2.0

    # Fast centroids via center_of_mass
    centroids = ndi.center_of_mass(
        np.ones_like(cell_labels, dtype=np.uint8),
        labels=cell_labels,
        index=labels
    )

    row_centroid = np.array([c[0] for c in centroids], dtype=np.float32)
    col_centroid = np.array([c[1] for c in centroids], dtype=np.float32)

    # Fast area via ndi.sum
    area = ndi.sum(
        np.ones_like(cell_labels, dtype=np.uint8),
        labels=cell_labels,
        index=labels
    ).astype(np.int32)

    # Mean distance to outside via label-wise sum / area
    dist_sum = ndi.sum(dist_out, labels=cell_labels, index=labels)
    mean_dist_out = (dist_sum / area).astype(np.float32)

    # Distance to centerline approximated from centroid position
    center_at_centroid = centers[
        np.clip(np.round(row_centroid).astype(int), 0, len(centers) - 1)
    ]
    mean_dist_center = np.abs(col_centroid - center_at_centroid).astype(np.float32)

    geom = pd.DataFrame({
        "label": labels.astype(np.int32),
        "row_centroid": row_centroid,
        "col_centroid": col_centroid,
        "area": area,
        "mean_dist_out": mean_dist_out,
        "mean_dist_center": mean_dist_center
    })

    # External background only
    bg_labels, _ = ndi.label(bg_mask)
    border_labels = np.unique(np.concatenate([
        bg_labels[0, :],
        bg_labels[-1, :],
        bg_labels[:, 0],
        bg_labels[:, -1]
    ]))
    border_labels = border_labels[border_labels > 0]
    external_bg = np.isin(bg_labels, border_labels)

    # Cells touching the external background
    outer_edge = cell_mask & ndi.binary_dilation(external_bg)
    outer_labels = np.unique(cell_labels[outer_edge])
    outer_labels = outer_labels[outer_labels > 0]

    pre["geom"] = geom
    pre["outer_labels"] = set(outer_labels.tolist())
    pre["external_bg"] = external_bg
    pre["marker_cache"] = {}

    return pre


# ============================================================
# 2) CACHED MARKER EXPRESSION PER SAMPLE
# ============================================================
def get_marker_expression_cached(pre, gene, aggregation="mean"):
    key = (gene, aggregation)
    if key in pre["marker_cache"]:
        return pre["marker_cache"][key]

    if gene not in pre["df"].columns:
        raise KeyError(f"{gene} not found in {pre['name']}")

    df = pre["df"][["label", gene]].copy()
    df["expr"] = np.log1p(pd.to_numeric(df[gene], errors="coerce"))
    df = df.dropna(subset=["expr"])

    if aggregation == "mean":
        out = df.groupby("label", as_index=False)["expr"].mean()
    elif aggregation == "median":
        out = df.groupby("label", as_index=False)["expr"].median()
    else:
        raise ValueError("aggregation must be mean or median")

    out = out.rename(columns={"expr": gene})
    pre["marker_cache"][key] = out
    return out


# ============================================================
# 3) SHARED VALID LONGITUDINAL ZONE
# ============================================================
def get_geom_with_valid_zone(pre, down_px=800, up_px=4000):
    geom = pre["geom"].copy()
    tissue = pre["cell_mask"] | pre["line_mask"]
    tip_row = np.where(tissue)[0].max()

    lower = max(tip_row - up_px, 0)
    upper = max(tip_row - down_px,0)

    geom["in_valid_zone"] = (
        (geom["row_centroid"] >= lower) &
        (geom["row_centroid"] <= upper)
    )
    geom_valid = geom[geom["in_valid_zone"]].copy()
    return geom_valid


# ============================================================
# 4) FAST TISSUE DEFINITIONS
# ============================================================
def define_epidermis_fast(pre, epi_marker, marker_quantile=0.75,
                          aggregation="mean"):
    geom = get_geom_with_valid_zone(pre)
    expr = get_marker_expression_cached(pre, epi_marker, aggregation=aggregation)

    df = geom.merge(expr, on="label", how="left")
    df["is_outer"] = df["label"].isin(pre["outer_labels"])

    outer_only = df[df["is_outer"] & df["in_valid_zone"]].copy()
    thr = float(outer_only[epi_marker].quantile(marker_quantile))

    df["is_epidermis"] = (
        df["is_outer"] &
        df["in_valid_zone"] &
        (df[epi_marker] >= thr)
    )

    return df[[
        "label", "row_centroid", "in_valid_zone",
        epi_marker, "is_outer", "is_epidermis"
    ]], thr


def define_cortex_fast(pre, cortex_marker, epi_df, marker_quantile=0.60,
                       aggregation="mean"):
    geom = get_geom_with_valid_zone(pre)
    expr = get_marker_expression_cached(pre, cortex_marker, aggregation=aggregation)

    df = geom.merge(expr, on="label", how="left")
    df = df.merge(epi_df[["label", "is_epidermis"]], on="label", how="left")
    df["is_epidermis"] = df["is_epidermis"].fillna(False)

    candidates = (~df["is_epidermis"]) & (df["in_valid_zone"])

    q = df.loc[candidates, "mean_dist_center"].quantile([0.4, 0.90])
    low_c, high_c = float(q.iloc[0]), float(q.iloc[1])
    marker_thr = float(df.loc[candidates, cortex_marker].quantile(marker_quantile))

    df["is_cortex"] = (
        candidates &
        (df["mean_dist_center"] >= low_c) &
        (df["mean_dist_center"] <= high_c) &
        (df[cortex_marker] >= marker_thr)
    )

    return df[[
        "label", "row_centroid", "in_valid_zone",
        "mean_dist_center", cortex_marker, "is_cortex"
    ]], marker_thr



def annotate_three_tissues_fast(
    pre,
    epi_marker,
    cortex_marker,
    epi_q=0.75,
    cortex_q=0.60,
    #stele_center_q=0.35,
    aggregation="mean"
):
    epi_df, epi_thr = define_epidermis_fast(
        pre,
        epi_marker=epi_marker,
        marker_quantile=epi_q,
        aggregation=aggregation)

    cortex_df, cortex_thr = define_cortex_fast(
        pre,
        cortex_marker=cortex_marker,
        epi_df=epi_df,
        marker_quantile=cortex_q,
        aggregation=aggregation
    )

    geom = get_geom_with_valid_zone(pre)

    annot = geom[["label", "row_centroid", "in_valid_zone"]].copy()
    annot = annot.merge(epi_df[["label", "is_epidermis"]], on="label", how="left")
    annot = annot.merge(cortex_df[["label", "is_cortex"]], on="label", how="left")
    #annot = annot.merge(stele_df[["label", "is_stele"]], on="label", how="left")
    annot = annot.fillna(False)

    annot["cell_type"] = "other"
    #annot.loc[annot["is_stele"], "cell_type"] = "stele"
    annot.loc[annot["is_cortex"], "cell_type"] = "cortex"
    annot.loc[annot["is_epidermis"], "cell_type"] = "epidermis"

    thresholds = {
        "epi_marker_threshold": epi_thr,
        "cortex_marker_threshold": cortex_thr
    }

    return annot, thresholds


# ============================================================
# 5) PLOTTING
# ============================================================
def plot_three_tissues(pre, annot, length_px=4000, save_path=None):
    cell_labels = pre["cell_labels"]
    cell_mask = pre["cell_mask"]
    line_mask = pre["line_mask"]
    bg_mask = pre["bg_mask"]

    tissue =cell_mask | line_mask
    tip_row = np.where(tissue)[0].max()

    row_start = max(tip_row - length_px, 0)
    row_end = tip_row

    lut = {
        "epidermis": [1.0, 0.2, 0.2],   # red
        "cortex":    [0.2, 0.6, 1.0],   # blue
        #"stele":     [0.2, 0.8, 0.2],   # green
        "other":     [0.85, 0.85, 0.85]
    }

    rgb = np.ones((*cell_labels.shape, 3), dtype=np.float32)
    rgb[bg_mask] = [1, 1, 1]
    rgb[line_mask] = [0, 0, 0]
    rgb[cell_mask] = lut["other"]
    
    rgb_crop = rgb[row_start:row_end, :]
    cell_labels_crop = cell_labels[row_start:row_end, :]

    for _, row in annot.iterrows():
        rgb_crop[cell_labels_crop == row["label"]] = lut[row["cell_type"]]

    fig, ax = plt.subplots(figsize=(4, 10))
    ax.imshow(rgb_crop, origin="upper")
    ax.set_title(f"{pre['name']} tissue annotation")
    ax.axis("off")

    if save_path is not None:
        pdf_path = os.path.join(save_path, f"{pre['name']}_tissue_seperation.pdf")
        print("Saving PDF to:", pdf_path)
        fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
        plt.close(fig)
        

# ============================================================
# 6) TARGET GENE QUANTIFICATION
# ============================================================
def quantify_target_by_tissue_fast(pre, annot, target_gene, aggregation="mean"):
    if target_gene not in pre["df"].columns:
        raise KeyError(f"{target_gene} not found in {pre['name']}")

    df = pre["df"][["label", target_gene]].copy()
    df["expr"] = np.log1p(pd.to_numeric(df[target_gene], errors="coerce"))
    df = df.dropna(subset=["expr"])

    if aggregation == "mean":
        expr_tbl = df.groupby("label", as_index=False)["expr"].mean()
    elif aggregation == "median":
        expr_tbl = df.groupby("label", as_index=False)["expr"].median()
    else:
        raise ValueError("aggregation must be mean or median")

    merged = expr_tbl.merge(annot[["label", "cell_type", "in_valid_zone"]], on="label", how="left")
    merged = merged[merged["in_valid_zone"] == True].copy()

    summary = merged.groupby("cell_type")["expr"].agg(
        n_cells="count",
        mean_expr="mean",
        median_expr="median",
        std_expr="std"
    ).reset_index()

    return merged, summary


# ============================================================
# 7) ONE-SAMPLE WRAPPER
# ============================================================
def annotate_and_quantify_one_sample_fast(
    pre,
    epi_marker,
    cortex_marker,
    target_gene,
    epi_q=0.75,
    cortex_q=0.60,
    #stele_center_q=0.35,
    aggregation="mean",
    plot=True,
    save_plot_path=None
):
    annot, thresholds = annotate_three_tissues_fast(
        pre=pre,
        epi_marker=epi_marker,
        cortex_marker=cortex_marker,
        epi_q=epi_q,
        cortex_q=cortex_q,
        #stele_center_q=stele_center_q,
        aggregation=aggregation
    )

    if plot:
        plot_three_tissues(pre, annot, save_path=save_plot_path)

    merged, summary = quantify_target_by_tissue_fast(
        pre=pre,
        annot=annot,
        target_gene=target_gene,
        aggregation=aggregation
    )

    return {
        "annotation_table": annot,
        "thresholds": thresholds,
        "cell_expression_table": merged,
        "summary_table": summary
    }
    
   
# ============================================================
# MAIN
# ============================================================

preprocessed = []

for ds in datasets:
    print(f"Preprocessing {ds['name']} ...")

    df = pd.read_csv(ds["csv"])
    mask = tifffile.imread(ds["mask"])

    coords = df["cell"].str.extract(r"\(([^,]+),([^)]+)\)").astype(float)
    df["y"] = coords[0]
    df["x"] = coords[1]

    cell_mask = mask == CELL_VAL
    line_mask = mask == LINE_VAL
    bg_mask = mask == BG_VAL

    cell_labels, n_cells = ndi.label(cell_mask)

    rows = np.clip(np.round(df["y"]).astype(int), 0, mask.shape[0] - 1)
    cols = np.clip(np.round(df["x"]).astype(int), 0, mask.shape[1] - 1)

    labels = cell_labels[rows, cols]

    nearest_idx = ndi.distance_transform_edt(
        cell_labels == 0,
        return_distances=False,
        return_indices=True
    )
    nearest_label_map = cell_labels[tuple(nearest_idx)]

    zero_mask = labels == 0
    labels[zero_mask] = nearest_label_map[rows[zero_mask], cols[zero_mask]]

    df["label"] = labels
    df = df[df["label"] > 0].copy()

    preprocessed.append({
        "name": ds["name"],
        "df": df,
        "mask": mask,
        "cell_mask": cell_mask,
        "line_mask": line_mask,
        "bg_mask": bg_mask,
        "cell_labels": cell_labels,
        "n_cells": n_cells,
    })

print("Preprocessing done.")

for pre in preprocessed:
    precompute_geometry(pre)
    
num = 0    

for gene in genes:
    num = num + 1
    if num==1:
        save_plot=True
    else:
        save_plot=False
    
    all_cells = []
    all_summary = []

    for pre in preprocessed:
        res = annotate_and_quantify_one_sample_fast(
            pre=pre,
            epi_marker="Zm00001d037547",
            cortex_marker="Zm00001d019565",
            target_gene=gene,
            epi_q=0.25,
            cortex_q=0.35,
            aggregation="mean",
            plot=save_plot,
            save_plot_path=output_dir
        )
    
        tmp_cells = res["cell_expression_table"].copy()
        tmp_cells["genotype"] = pre["name"]
        all_cells.append(tmp_cells)
    
        tmp_summary = res["summary_table"].copy()
        tmp_summary["genotype"] = pre["name"]
        all_summary.append(tmp_summary)
    
    all_cells = pd.concat(all_cells, ignore_index=True)
    all_summary = pd.concat(all_summary, ignore_index=True)
    
    print(all_summary)
    csv_path = os.path.join(output_dir, f"{gene}_expression_in_tissues.csv")
    all_cells.to_csv(csv_path, index=False)

print("Done.")
