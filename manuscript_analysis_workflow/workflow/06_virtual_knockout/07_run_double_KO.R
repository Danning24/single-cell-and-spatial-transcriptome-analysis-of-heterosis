project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: run_double_ko.R <group>")
}

group <- args[[1]]
target_genes <- c("Zm00001d007234", "Zm00001d043795")
target_key <- paste(target_genes, collapse = "_plus_")

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("SLURM_NTASKS", "1")))
if (is.na(ncores) || ncores < 1) ncores <- 1L

Sys.setenv(
  OMP_NUM_THREADS = ncores,
  OPENBLAS_NUM_THREADS = ncores,
  MKL_NUM_THREADS = ncores,
  VECLIB_MAXIMUM_THREADS = ncores,
  NUMEXPR_NUM_THREADS = ncores
)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(scTenifoldKnk)
})

input_file <- file.path(
  "new_ko_genes_scripts",
  "double_ko",
  "inputs",
  paste0("countMatrix_", group, ".rds")
)

if (!file.exists(input_file)) {
  stop("input matrix not found: ", input_file)
}

expression_matrix <- readRDS(input_file)
expression_matrix <- as(expression_matrix, "dgCMatrix")

missing_genes <- setdiff(target_genes, rownames(expression_matrix))
if (length(missing_genes) > 0) {
  stop("Double KO genes not found in count matrix: ", paste(missing_genes, collapse = ","))
}

result_dir <- file.path(
  "new_ko_genes_scripts",
  "double_ko",
  "results",
  paste(group, target_key, sep = "__")
)
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

out_prefix <- file.path(result_dir, paste(group, target_key, sep = "__"))

cat("START\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("GROUP\t", group, "\n", sep = "")
cat("DOUBLE_KO\t", paste(target_genes, collapse = "+"), "\n", sep = "")
cat("MATRIX_DIM\t", nrow(expression_matrix), "\t", ncol(expression_matrix), "\n", sep = "")
cat("CORES\t", ncores, "\n", sep = "")
cat("PARAMETERS\tseed=2026; qc_mtThreshold=0.1; qc_minLSize=1000; nc_lambda=0; nc_nNet=10; nc_nCells=500; nc_nComp=3; nc_scaleScores=TRUE; nc_symmetric=FALSE; nc_q=0.9; td_K=3; td_maxIter=1000; td_maxError=1e-05; td_nDecimal=3; ma_nDim=2\n", sep = "")

set.seed(2026)
perturbation_results <- scTenifoldKnk(
  countMatrix = expression_matrix,
  gKO = target_genes,

  qc_mtThreshold = 0.1,
  qc_minLSize = 1000,

  nc_lambda = 0,
  nc_nNet = 10,
  nc_nCells = 500,
  nc_nComp = 3,
  nc_scaleScores = TRUE,
  nc_symmetric = FALSE,
  nc_q = 0.9,

  td_K = 3,
  td_maxIter = 1000,
  td_maxError = 1e-05,
  td_nDecimal = 3,

  ma_nDim = 2
)

saveRDS(perturbation_results, paste0(out_prefix, ".sctenifoldKnk.rds"), compress = "xz")

if (!"diffRegulation" %in% names(perturbation_results)) {
  stop("Result object does not contain diffRegulation")
}

diffreg <- as.data.table(perturbation_results$diffRegulation)
fwrite(diffreg, paste0(out_prefix, ".diffRegulation.csv"))

if ("p.adj" %in% colnames(diffreg) && "Z" %in% colnames(diffreg)) {
  sig <- diffreg[p.adj < 0.05][order(p.adj, -abs(Z))]
} else if ("p.value" %in% colnames(diffreg)) {
  sig <- diffreg[p.value < 0.05]
} else {
  sig <- diffreg
}

fwrite(sig, paste0(out_prefix, ".sig_genes.csv"))
fwrite(
  data.table(
    group = group,
    KO = paste(target_genes, collapse = "+"),
    n_sig = nrow(sig),
    runner = "double_ko_standard_ma2"
  ),
  paste0(out_prefix, ".sig_gene_counts.csv")
)

try({
  png(paste0(out_prefix, ".plotKO.png"), width = 1800, height = 1400, res = 180)
  plotKO(perturbation_results, gKO = target_genes)
  dev.off()
}, silent = TRUE)

cat("END\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
