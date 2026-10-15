# Maize manuscript analysis workflow

Nine analysis scripts for cortex development, non-additive gene (NAG) scores,
spatial comparisons, and virtual knockout analyses in maize roots. The scripts
use existing processed objects and finalized annotations; raw-data preprocessing
and NAG calling are outside their scope. Data and archived analysis results are
not distributed here.

## Layout and stages

All nine scripts are directly inside `workflow/`. Run them from that directory
so that the bundled `inputs/knockout_targets.csv` is found.

| Script | Analysis | Figure |
| --- | --- | --- |
| `01_cortex_clustering.R` | Cortex re-clustering and NSS/SS/TL UMAP density | Fig. 3A-B |
| `02_developmental_scores.R` | Root-slice M/E programmes and cortex module scores | Fig. 3C |
| `03_monocle3_trajectory.R` | Monocle3 graph and interactive root selection | Fig. 3D |
| `04_relative_mmd.py` | Relative MMD in a shared space with held-out cells | Supplementary Fig. S16C |
| `05_cortex_nag_correlations.R` | Pooled, hybrid-specific, and directional NAG scores | Fig. 3E-G; S15d |
| `06_spatial_nag_scores.R` | Spatial NAG maps and paired section-median comparisons | Fig. 3H-I; S17 |
| `07_virtual_knockout.R` | Three historical knockout profiles across two hybrids | 36 targets; 72 knockouts |
| `08_knockout_summary.R` | Top-decile enrichment, tissue preference, and shared effects | Fig. 4a-b; S18a |
| `09_knockout_enrichment.R` | GSEA and GO for eight C-N candidates | Fig. 4c; S18b |

The target CSV contains `gene`, `symbol`, `module`, `batch`, and `batch_order`.
Stable gene IDs identify knockouts; symbols are display labels. The corrected
GOT1/THP3 display label retains gene ID `Zm00001d043382`.

## Environment and required inputs

From the repository root:

```sh
cd manuscript_analysis_workflow
conda env create -f environment.yml
conda activate maize-heterosis-analysis
cd workflow
```

The Conda file targets **Linux** and records the Python versions used for the
relative-MMD analysis plus an R 4.4/Bioconductor 3.20 reconstruction environment.
It supplies the core dependencies for stages 01-02, 04-06, and 08-09. It is not a
complete historical lock file, and a fresh end-to-end run has not been performed.

For stage 07, install the archived, released scTenifold packages into this
environment from R, in this order:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
remotes::install_version("scTenifoldNet", version = "1.3",
                         dependencies = NA, upgrade = "never")
remotes::install_version("scTenifoldKnk", version = "1.0.1",
                         dependencies = NA, upgrade = "never")
```

These [Net](https://cran.r-project.org/src/contrib/Archive/scTenifoldNet/) and
[Knk](https://cran.r-project.org/src/contrib/Archive/scTenifoldKnk/) releases match
the interfaces used by stage 07. An unconstrained install may select a newer API
and should not be substituted silently.

Stage 03 additionally requires an already configured trajectory R environment
with `monocle3`, `SeuratWrappers 0.4.0`, `Seurat`, `SingleCellExperiment`, and
`ggplot2`. These trajectory extras are **not installed by this Conda file**.
[SeuratWrappers](https://github.com/satijalab/seurat-wrappers) has no tagged stable
release to pin and is unavailable from these Conda channels; no development-branch
installation is provided. Use the existing trajectory environment for this stage
and retain its `sessionInfo()` rather than treating this file as its replacement.
Keep that environment separate when installing the core environment above.

Core packages are Seurat for single-cell/spatial objects and module scores,
DESeq2/apeglm for root-slice differential expression, Monocle3/SeuratWrappers for
trajectories, scTenifoldKnk/scTenifoldNet for knockouts, and fgsea/clusterProfiler
for enrichment. The Python MMD stage uses Scanpy, AnnData, NumPy, pandas, SciPy,
and scikit-learn. These established analysis packages retain the original methods.
Use stable releases and preserve `sessionInfo.txt` files produced by R stages.
An installable environment is not evidence that every historical analysis used
those exact versions; the original knockout environments and RNG states were not
fully recorded.

Set `MAIZE_PROJECT_DIR` to a data root containing the following relative paths:

| Stage | Input below `MAIZE_PROJECT_DIR` |
| --- | --- |
| 01 | `0.data/scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS` |
| 02, 03, 05 | `4.Cortex_cells/data_output/Cortex.rds`, unless overridden |
| 04 | `4.Cortex_cells/data_output/cortex_h5ad/` containing `F1_NSS_vs_B73_SS.h5ad`, `F1_NSS_vs_P1_NSS.h5ad`, `F1_SS_vs_P1_SS.h5ad`, and `F1_TL_vs_P1_TL.h5ad` |
| 04 | `4.Cortex_cells/analysis_outs/MMD_parent_comparison_20261010/author_relative_mmd.py` |
| 05 | `0.data/2026TUM_Danning_Wang_To_HK/To HK/sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv` |
| 06 | In the same `To HK` directory: `all_spatial_merge_filtered_integrated_remove_noise_bin50_annotation.RDS` and `spatial_heterotic_genes_results.csv` |
| 07 | `7.sctenifold_KO/sctenifold20260527_ROS_rel/scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS`, unless overridden |
| 08, 09 | `7.sctenifold_KO/analysis_outs/tables/00_KO_input_manifest.csv`, unless overridden; the result CSVs listed in its `file` column are also required |
| 09 | `7.sctenifold_KO/sctenifold20260527_ROS_rel/sctenifold_downstream_analysis/maize_GO_annotation_maizeGAMER/maizeGAMER_GO_gene_map_used.csv`, unless overridden |

Stage 02 additionally needs the GSE172277 root-slice count table through
`MAIZE_ROOT_SLICE_COUNTS`. Preserve its original 32-library column order:
S1a/S1b through S16a/S16b. Stages 08-09 require the finalized single-cell NAG
table through `MAIZE_NAG_TABLE`.

Stage 04 imports the external `author_relative_mmd.py` implementation; installing
the environment alone does not provide this file. Retain the original copy with
`MMD_3_Sample_Test` and `MMD_unbiased` from the
[author implementation](https://github.com/eugenium/MMD). The recorded SHA256 is
`0d6795d62ae3d854e16bbc67a5b8cfe543d95f75035fc167f9033f7ebb363d87`.
The recorded fixed GO mapping MD5 is `5332140b225b782a5232daf458159150`.

## Configuration and execution

Paths are configured through environment variables. Relative output and target
paths are interpreted from the working directory, `workflow/`.

| Variable | Stages | Setting |
| --- | --- | --- |
| `MAIZE_PROJECT_DIR` | All | Required data/project root |
| `MAIZE_RESULTS_DIR` | All | Output root; default `results`, with one thematic subdirectory per stage |
| `MAIZE_ROOT_SLICE_COUNTS` | 02 | Required path to the GSE172277 count table |
| `MAIZE_CORTEX_RDS` | 02, 03, 05 | Optional cortex object; default path is listed above |
| `MAIZE_MONOCLE_ROOT_FILE` | 03 | Optional newline-separated principal-graph root IDs; valid only for the same learned graph |
| `MAIZE_KO_SEURAT` | 07 | Optional source Seurat RDS; default path is listed above |
| `MAIZE_KO_TARGETS` | 07-09 | Default `inputs/knockout_targets.csv` |
| `MAIZE_CORES` | 07 | Falls back to `SLURM_CPUS_PER_TASK`, then `SLURM_NTASKS`, then 8; historical first two batches use package/detected-core network defaults |
| `MAIZE_KO_MANIFEST` | 08-09 | Optional manifest; defaults to the historical manifest listed above |
| `MAIZE_NAG_TABLE` | 08-09 | Required finalized single-cell NAG table path |
| `MAIZE_GO_MAPPING` | 09 | Optional fixed GO mapping; default path is listed above |
| `MAIZE_GSEA_SEED` | 09 | Optional integer seed for a new GSEA run; historical multilevel RNG state is unavailable |

Example configuration in a POSIX shell; replace the placeholders with your paths:

```sh
export MAIZE_PROJECT_DIR="/path/to/maize_project"
export MAIZE_RESULTS_DIR="/path/to/new_results"
export MAIZE_ROOT_SLICE_COUNTS="/path/to/GSE172277_RootSlices_genecounts.txt.gz"
export MAIZE_NAG_TABLE="$MAIZE_PROJECT_DIR/0.data/2026TUM_Danning_Wang_To_HK/To HK/sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv"
export MAIZE_GSEA_SEED=20261015
```

For the historical processed inputs, keep `MAIZE_CORTEX_RDS` unset and run the
required stages independently:

```sh
Rscript 01_cortex_clustering.R
Rscript 02_developmental_scores.R
python 04_relative_mmd.py
Rscript 05_cortex_nag_correlations.R
Rscript 06_spatial_nag_scores.R
Rscript 08_knockout_summary.R
Rscript 09_knockout_enrichment.R
```

These commands do not automatically chain newly computed objects. To start with
stage 01's new clustering instead, use:

```sh
Rscript 01_cortex_clustering.R
export MAIZE_CORTEX_RDS="$MAIZE_RESULTS_DIR/cortex_clustering/Cortex.rds"
Rscript 02_developmental_scores.R
export MAIZE_CORTEX_RDS="$MAIZE_RESULTS_DIR/developmental_scores/Cortex_M_E_scores.rds"
Rscript 05_cortex_nag_correlations.R
```

Stage 05 requires `M_Scores1` and `E_Scores1`; when using new objects, use the
scored object from stage 02. Stage 04 uses existing H5AD exports rather than
exporting stage 01's new cortex object.

Stage 03 requires an interactive R session unless a valid root file is supplied.
From `workflow/`, with the same environment variables available in R:

```r
source("03_monocle3_trajectory.R")
```

Choose the lower-UMAP, high-M region. The script saves the full trajectory object
and selected root IDs. The original interactive selection was not recorded;
root IDs from one graph cannot generally be reused after learning another graph.

Stage 07 is computationally expensive. It is unnecessary if stages 08-09 will use
the archived 72 knockout results. To compute all batch/hybrid combinations:

```sh
for batch in original6 additional22 completion8; do
  for hybrid in B73xA554 B73xOh40B; do
    Rscript 07_virtual_knockout.R "$batch" "$hybrid"
  done
done
```

Stage 07 writes result RDS/CSV files under
`$MAIZE_RESULTS_DIR/virtual_knockout/<batch>/`. It does not generate the manifest
used by stages 08-09. For new or relocated results, provide a CSV manifest with
columns `Hybrid`, `KO`, and `file`, one row per hybrid/target, and set
`MAIZE_KO_MANIFEST` to it. `KO` must match the stable `gene` IDs in the target CSV;
each `file` must resolve to its corresponding differential-regulation CSV
containing `gene`, `Z`, and `p.adj`. Use absolute file paths or paths relative to
`workflow/`. Preserve each knockout's tested-gene background.

## Interpretation and reproducibility

- Cortex labels retain the manuscript mapping: raw 0 = Cortex1, raw 1 = Cortex3,
  raw 2 = Cortex2. Fresh clustering may permute numeric labels or change the
  embedding; archived objects are needed for the exact original figures.
- Stages 05-06 and 08-09 consume finalized NAG classifications. They do not
  reproduce the upstream NAG-calling model. `MPH.pct` direction is determined
  from its sign. Cell-level Pearson correlations describe associations and do
  not establish biological replication.
- Stage 06 scores the jointly stored cortex reference before selecting ten
  hybrid sections (1-6 and 8-11). Radial tests pair section medians, with BH
  adjustment within each module. Bins are not independent biological replicates;
  ten sections do not automatically represent ten independent plants.
- Stage 07 preserves `original6`, `additional22`, and `completion8` feature sets.
  A554 retains the historical safe one-dimensional alignment workaround; Oh40B
  uses two-dimensional alignment. The completion batch retains rare targets and
  reuses one WT network per hybrid. Seeds are recorded in the script, but package,
  BLAS, and parallel differences can still affect reruns. These historical
  profiles should not be merged into a common gene background.
- Stage 04 uses 500 training and 2,000 held-out cells per population, a shared
  fitted space, and Bonferroni adjustment across three contrasts. Its 100 plotting
  draws do not determine significance or constitute biological replicates.
- Stage 09 applies BH separately to each GSEA category's 40 planned tests. GO
  tests pool BP/MF/CC per hybrid/knockout, use `minGSSize = 5`, and retain terms
  with Count >= 5 and FDR < 0.05. Exact historical GSEA Monte Carlo P values
  require archived outputs or RNG state; setting `MAIZE_GSEA_SEED` makes a new
  run reproducible within a compatible environment.

Outputs are written under `MAIZE_RESULTS_DIR`; R stages also record
`sessionInfo.txt` (stage 07 records one file per hybrid/batch). This compact
workflow has not been validated by a fresh end-to-end rerun of all nine stages.
For this publication update, all eight R scripts passed the official R parser
(R 4.5.1), the Python script passed `ast.parse()`, and the target CSV was checked
for its 36 unique gene IDs and 6/22/8 batch sizes. These are syntax and packaging
checks; they do not establish dependency loading or numerical reproducibility.
The pinned core packages have Linux-compatible builds in the named channels;
the full joint dependency solve and environment installation have not been
completed.
