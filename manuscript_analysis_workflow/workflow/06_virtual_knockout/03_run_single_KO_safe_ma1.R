project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_sctenifold_knk_safe.R <group> <gene>")
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
  library(scTenifoldNet)
})

safe_manifold_alignment <- function(X, Y, d = 1, nCores = parallel::detectCores()) {
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]

  L <- diag(length(sharedGenes))
  wX <- X + 1
  wY <- Y + 1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L

  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)

  RhpcBLASctl::omp_set_num_threads(nCores)
  RhpcBLASctl::blas_set_num_threads(nCores)

  E <- suppressWarnings(RSpectra::eigs(W, max(2, d * 2), "SR"))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors, 2, as.numeric))

  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[, newOrder, drop = FALSE]

  keep <- E$values > 1e-08
  E$vectors <- E$vectors[, keep, drop = FALSE]

  if (ncol(E$vectors) < d) {
    warning("Only ", ncol(E$vectors), " positive alignment dimensions available; using d = ", ncol(E$vectors))
    d <- ncol(E$vectors)
  }
  if (d < 1) {
    stop("No positive alignment dimensions available")
  }

  alignedNet <- E$vectors[, seq_len(d), drop = FALSE]
  colnames(alignedNet) <- paste0("NLMA ", seq_len(d))
  rownames(alignedNet) <- c(paste0("X_", sharedGenes), paste0("Y_", sharedGenes))
  alignedNet
}

run_knk_safe <- function(countMatrix, gKO) {
  countMatrix <- getFromNamespace("scQC", "scTenifoldKnk")(countMatrix)

  if (ncol(countMatrix) > 500) {
    countMatrix <- countMatrix[rowMeans(countMatrix != 0) >= 0.05, ]
  } else {
    countMatrix <- countMatrix[rowSums(countMatrix != 0) >= 25, ]
  }

  WT <- scTenifoldNet::makeNetworks(
    X = countMatrix,
    q = 0.9,
    nNet = 10,
    nCells = 500,
    scaleScores = TRUE,
    symmetric = FALSE,
    nComp = 3
  )

  WT <- scTenifoldNet::tensorDecomposition(
    xList = WT,
    K = 3,
    maxError = 1e-05,
    maxIter = 1000,
    nDecimal = 3
  )

  WT <- WT$X
  WT <- getFromNamespace("strictDirection", "scTenifoldKnk")(WT, lambda = 0)
  WT <- as.matrix(WT)
  diag(WT) <- 0
  WT <- t(WT)

  KO <- WT
  KO[gKO, ] <- 0

  MA <- safe_manifold_alignment(WT, KO, d = 1, nCores = ncores)
  DR <- getFromNamespace("dRegulation", "scTenifoldKnk")(MA, gKO)

  outputList <- list()
  outputList$tensorNetworks$WT <- Matrix(WT)
  outputList$tensorNetworks$KO <- Matrix(KO)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  outputList
}

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
cat("RUNNER\tSAFE_MA_DIM1\n")

set.seed(20260521)
result <- run_knk_safe(count_matrix, gene)

saveRDS(result, paste0(out_prefix, ".safe_ma1.sctenifoldKnk.rds"), compress = "xz")

diffreg <- as.data.table(result$diffRegulation)
fwrite(diffreg, paste0(out_prefix, ".safe_ma1.diffRegulation.csv"))

if ("p.adj" %in% colnames(diffreg) && "Z" %in% colnames(diffreg)) {
  sig <- diffreg[p.adj < 0.05][order(p.adj, -abs(Z))]
} else if ("p.value" %in% colnames(diffreg)) {
  sig <- diffreg[p.value < 0.05]
} else {
  sig <- diffreg
}

fwrite(sig, paste0(out_prefix, ".safe_ma1.sig_genes.csv"))
fwrite(
  data.table(group = group, KO = gene, n_sig = nrow(sig), runner = "safe_ma1"),
  paste0(out_prefix, ".safe_ma1.sig_gene_counts.csv")
)

try({
  png(paste0(out_prefix, ".safe_ma1.plotKO.png"), width = 1800, height = 1400, res = 180)
  plotKO(result, gKO = gene)
  dev.off()
}, silent = TRUE)

cat("END\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
