"""Supplementary Fig. S16C: held-out relative MMD in one shared cortex space."""

import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import pdist
from scipy.stats import norm
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.preprocessing import StandardScaler

# Paths and fixed analysis parameters.
project_dir = Path(os.environ["MAIZE_PROJECT_DIR"])
input_dir = project_dir / "4.Cortex_cells/data_output/cortex_h5ad"
author_code_dir = project_dir / "4.Cortex_cells/analysis_outs/MMD_parent_comparison_20261010"
output_dir = Path(os.getenv("MAIZE_RESULTS_DIR", "results")) / "relative_mmd"
output_dir.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(author_code_dir))
from author_relative_mmd import MMD_3_Sample_Test, MMD_unbiased

seed, n_train, n_test = 20261010, 500, 2000
n_hvg, n_pcs, n_draws, draw_size = 2000, 10, 100, 100
rng = np.random.default_rng(seed)

# Assemble seven unique populations from SCT-corrected count exports.
inputs = [
    ("F1_NSS_vs_B73_SS.h5ad", None),
    ("F1_NSS_vs_P1_NSS.h5ad", "P1_NSS"),
    ("F1_SS_vs_P1_SS.h5ad", None),
    ("F1_TL_vs_P1_TL.h5ad", None),
]
parts = []
for filename, population in inputs:
    obj = sc.read_h5ad(input_dir / filename)
    if population is not None:
        obj = obj[obj.obs["celltype"] == population].copy()
    parts.append(ad.AnnData(X=obj.X.copy(), obs=obj.obs[["celltype", "orig.ident"]].copy(), var=obj.var.copy()))
adata = ad.concat(parts, join="inner", merge="same")

# Fit gene selection, scaling and PCA only on balanced training cells.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
labels = adata.obs["celltype"].to_numpy()
training = np.concatenate([
    rng.choice(np.flatnonzero(labels == label), n_train, replace=False)
    for label in sorted(np.unique(labels))
])
train = adata[training].copy()
sc.pp.highly_variable_genes(train, flavor="seurat", n_top_genes=n_hvg)
genes = train.var_names[train.var["highly_variable"]]
scaler = StandardScaler()
train_scaled = np.minimum(scaler.fit_transform(train[:, genes].X.toarray().astype(np.float64)), 10)
pca = PCA(n_components=n_pcs, svd_solver="arpack", random_state=seed)
train_pca = pca.fit_transform(train_scaled)
coordinates = pca.transform(np.minimum(scaler.transform(adata[:, genes].X.toarray().astype(np.float64)), 10))
sigma = float(np.median(pdist(train_pca)))

# Hold out training cells; use one shared B73 test sample for all comparisons.
is_training = np.zeros(adata.n_obs, dtype=bool)
is_training[training] = True
populations = {
    label: coordinates[rng.choice(np.flatnonzero((labels == label) & ~is_training), n_test, replace=False)]
    for label in sorted(np.unique(labels))
}
tests, draws = [], []
for group in ["NSS", "SS", "TL"]:
    x, y, z = populations[f"F1_{group}"], populations["B73_SS"], populations[f"P1_{group}"]
    _, statistic, _, first, second = MMD_3_Sample_Test(x, y, z, sigma=sigma, computeMMDs=True)
    p_value = float(2 * norm.sf(abs(statistic)))
    tests.append({"group": group, "MMD_F1_B73": first, "MMD_F1_P1": second,
                  "delta": first - second, "z_statistic": statistic,
                  "p_two_sided": p_value, "p_Bonferroni": min(3 * p_value, 1)})

    # Plotting draws are not independent replicates and do not determine P values.
    for iteration in range(1, n_draws + 1):
        a, b, c = (population[rng.integers(n_test, size=draw_size)] for population in (x, y, z))
        gamma = 1 / (2 * sigma**2)
        kaa = rbf_kernel(a, a, gamma=gamma)
        first = MMD_unbiased(kaa, rbf_kernel(b, b, gamma=gamma), rbf_kernel(a, b, gamma=gamma))
        second = MMD_unbiased(kaa, rbf_kernel(c, c, gamma=gamma), rbf_kernel(a, c, gamma=gamma))
        draws.append({"group": group, "iteration": iteration, "MMD_F1_B73": first,
                      "MMD_F1_P1": second, "delta": first - second})

pd.DataFrame(tests).to_csv(output_dir / "relative_mmd_tests.csv", index=False)
pd.DataFrame(draws).to_csv(output_dir / "paired_mmd_draws.csv", index=False)
pd.Series(genes, name="gene").to_csv(output_dir / "training_hvg.csv", index=False)
pd.Series({"seed": seed, "sigma": sigma, "training_per_population": n_train,
           "test_per_population": n_test}).to_csv(output_dir / "parameters.csv", header=["value"])
