project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_sctenifold_knk.R <group> <gene>")
}

group <- args[[1]]
gene <- args[[2]]

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("SLURM_NTASKS", "1")))
if (is.na(ncores) || ncores < 1) {
  ncores <- 1L
}

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
  "sctenifold_inputs",
  paste0("countMatrix_", group, ".rds")
)

if (!file.exists(input_file)) {
  stop("input matrix not found: ", input_file)
}

count_matrix <- readRDS(input_file)
count_matrix <- as(count_matrix, "dgCMatrix")

if (!gene %in% rownames(count_matrix)) {
  stop("KO gene not found in count matrix rownames: ", gene)
}

dir.create(
  file.path("new_ko_genes_scripts", "sctenifold_results"),
  showWarnings = FALSE,
  recursive = TRUE
)

out_prefix <- file.path(
  "new_ko_genes_scripts",
  "sctenifold_results",
  paste(group, gene, sep = "__")
)

cat("START\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("GROUP\t", group, "\n", sep = "")
cat("GENE\t", gene, "\n", sep = "")
cat("MATRIX_DIM\t", nrow(count_matrix), "\t", ncol(count_matrix), "\n", sep = "")
cat("CORES\t", ncores, "\n", sep = "")

set.seed(20260521)
result <- scTenifoldKnk(
  countMatrix = count_matrix,
  gKO = gene
)

saveRDS(result, paste0(out_prefix, ".sctenifoldKnk.rds"), compress = "xz")

if (!"diffRegulation" %in% names(result)) {
  stop("Result object does not contain diffRegulation")
}

diffreg <- as.data.table(result$diffRegulation)
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
  data.table(group = group, KO = gene, n_sig = nrow(sig)),
  paste0(out_prefix, ".sig_gene_counts.csv")
)

try({
  png(paste0(out_prefix, ".plotKO.png"), width = 1800, height = 1400, res = 180)
  plotKO(result, gKO = gene)
  dev.off()
}, silent = TRUE)

cat("END\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
