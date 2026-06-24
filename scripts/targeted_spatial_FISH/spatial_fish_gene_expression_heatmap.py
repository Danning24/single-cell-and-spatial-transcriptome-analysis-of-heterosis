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

output_dir = "/home/go67saw/single_cell_RNA_seq/data/spatial_fish/multi_panel_pdfs"
os.makedirs(output_dir, exist_ok=True)

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


for gene in genes:
    pdf_path = os.path.join(output_dir, f"{gene}_multi_panel_shared_scale.pdf")
    if os.path.exists(pdf_path):
        print(f"Skipping {gene}, already exists.")
        saved_pdfs.append(pdf_path)
        continue
    print(f"Processing {gene} ...")

    results = []
    all_vals = []

    skip_gene = False

    for pre in preprocessed:
        try:
            res = process_gene_from_preprocessed(pre, gene, aggregation=aggregation)
            results.append(res)

            vals = res["paint"][res["cell_mask"] & ~np.isnan(res["paint"])]
            if len(vals) > 0:
                all_vals.append(vals)

        except KeyError as e:
            print(f"Skipping {gene}: {e}")
            skip_gene = True
            break

    if skip_gene or len(all_vals) == 0:
        continue

    all_vals = np.concatenate(all_vals)

    shared_vmax = np.max(all_vals)
    norm = Normalize(vmin=global_vmin, vmax=shared_vmax)

    print(f"Shared scale for {gene}: vmin={global_vmin:.4f}, vmax={shared_vmax:.4f}")

    rgbs = [build_rgb(res, norm, cmap) for res in results]

    save_multi_panel(results, rgbs, gene, norm, cmap, output_dir)
    
    del results, all_vals, rgbs
    gc.collect()

print("Done.")
