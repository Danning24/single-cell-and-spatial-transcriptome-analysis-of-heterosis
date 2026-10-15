# scTenifoldKnk/scTenifoldNet reproduce the original network workflow.
# SeuratObject and Matrix extract RNA counts and SCT variable features.
# Run separately for each batch and hybrid, using its recorded R environment:
# Rscript 07_virtual_knockout.R original6 B73xA554
# Batches: original6, additional22, completion8. Hybrids: B73xA554, B73xOh40B.

ncores <- as.integer(Sys.getenv("MAIZE_CORES", Sys.getenv(
  "SLURM_CPUS_PER_TASK", Sys.getenv("SLURM_NTASKS", "8")
)))
Sys.setenv(OMP_NUM_THREADS = ncores, OPENBLAS_NUM_THREADS = ncores,
           MKL_NUM_THREADS = ncores, VECLIB_MAXIMUM_THREADS = ncores,
           NUMEXPR_NUM_THREADS = ncores)
library(SeuratObject)
library(Matrix)
library(data.table)
library(scTenifoldKnk)
library(scTenifoldNet)

# Inputs and outputs -----------------------------------------------------------
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
ko_root <- file.path(project_dir, "7.sctenifold_KO")
source_rds <- Sys.getenv("MAIZE_KO_SEURAT", file.path(
  ko_root, "sctenifold20260527_ROS_rel",
  "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"
))
targets <- fread(Sys.getenv("MAIZE_KO_TARGETS", "inputs/knockout_targets.csv"))
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", "results"), "virtual_knockout")
args <- commandArgs(trailingOnly = TRUE)
batch_name <- args[[1]]
hybrid <- args[[2]]
batch_targets <- targets[batch == batch_name][order(batch_order), gene]
output_dir <- file.path(output_dir, batch_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Preserve the three historical feature sets; do not merge their target lists.
object <- readRDS(source_rds)
counts <- LayerData(object, assay = "RNA", layer = "counts")
features <- unique(c(VariableFeatures(object[["SCT"]]), batch_targets))
cells <- rownames(object[[]])[object$orig.ident == hybrid]
if (batch_name == "completion8") cells <- intersect(colnames(counts), cells)
counts <- as(counts[features, cells, drop = FALSE], "dgCMatrix")
fwrite(data.table(batch = batch_name, hybrid, gene = rownames(counts)),
       file.path(output_dir, paste0(hybrid, "_input_genes.csv")))

# Historical A554 workaround ---------------------------------------------------
# The installed manifoldAlignment failed when only one positive dimension
# remained. Retain the original safe_ma1 override rather than changing the model.
safe_alignment <- function(X, Y) {
  shared <- intersect(rownames(X), rownames(Y))
  X <- X[shared, shared]
  Y <- Y[shared, shared]
  L <- diag(length(shared))
  wX <- X + 1
  wY <- Y + 1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- -rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  RhpcBLASctl::omp_set_num_threads(ncores)
  RhpcBLASctl::blas_set_num_threads(ncores)
  E <- suppressWarnings(RSpectra::eigs(W, 2, "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- matrix(as.numeric(E$vectors), nrow = nrow(E$vectors))
  index <- order(E$values)
  values <- E$values[index]
  vectors <- E$vectors[, index, drop = FALSE]
  vectors <- vectors[, values > 1e-8, drop = FALSE]
  aligned <- vectors[, 1, drop = FALSE]
  colnames(aligned) <- "NLMA 1"
  rownames(aligned) <- c(paste0("X_", shared), paste0("Y_", shared))
  aligned
}

# Network reconstruction -------------------------------------------------------
build_network <- function(X, retain_targets = FALSE) {
  X <- getFromNamespace("scQC", "scTenifoldKnk")(
    X, mtThreshold = 0.1, minLSize = 1000
  )
  keep <- Matrix::rowMeans(X != 0) >= 0.05
  if (retain_targets) keep <- keep | rownames(X) %in% batch_targets
  X <- X[keep, , drop = FALSE]
  networks <- makeNetworks(
    X, q = 0.9, nNet = 10, nCells = 500, nComp = 3,
    scaleScores = TRUE, symmetric = FALSE,
    nCores = if (batch_name == "completion8") ncores else parallel::detectCores()
  )
  WT <- tensorDecomposition(
    networks, K = 3, maxError = 1e-5, maxIter = 1000, nDecimal = 3
  )$X
  WT <- as.matrix(getFromNamespace("strictDirection", "scTenifoldKnk")(
    WT, lambda = 0
  ))
  diag(WT) <- 0
  t(WT)
}

# The completion batch reused one WT network per hybrid and retained rare targets.
if (batch_name == "completion8") {
  set.seed(20260521)
  shared_WT <- build_network(counts, retain_targets = TRUE)
}

# Virtual knockouts ------------------------------------------------------------
for (gene in batch_targets) {
  if (batch_name != "completion8") set.seed(20260521)
  if (hybrid == "B73xOh40B" && batch_name != "completion8") {
    result <- scTenifoldKnk(
      countMatrix = counts, gKO = gene,
      qc_mtThreshold = 0.1, qc_minLSize = 1000,
      nc_lambda = 0, nc_nNet = 10, nc_nCells = 500, nc_nComp = 3,
      nc_scaleScores = TRUE, nc_symmetric = FALSE, nc_q = 0.9,
      td_K = 3, td_maxIter = 1000, td_maxError = 1e-5,
      td_nDecimal = 3, ma_nDim = 2
    )
  } else {
    WT <- if (batch_name == "completion8") shared_WT else build_network(counts)
    KO <- WT
    KO[gene, ] <- 0
    alignment <- if (hybrid == "B73xA554") safe_alignment(WT, KO) else
      manifoldAlignment(WT, KO, d = 2, nCores = ncores)
    result <- list(
      tensorNetworks = list(WT = Matrix(WT), KO = Matrix(KO)),
      manifoldAlignment = alignment,
      diffRegulation = getFromNamespace("dRegulation", "scTenifoldKnk")(
        alignment, gene
      )
    )
  }
  prefix <- file.path(output_dir, paste(hybrid, gene, sep = "__"))
  suffix <- if (hybrid == "B73xA554") ".safe_ma1" else ""
  saveRDS(result, paste0(prefix, suffix, ".sctenifoldKnk.rds"), compress = "xz")
  fwrite(as.data.table(result$diffRegulation),
         paste0(prefix, suffix, ".diffRegulation.csv"))
}
capture.output(sessionInfo(), file = file.path(output_dir, paste0(hybrid, "_sessionInfo.txt")))
