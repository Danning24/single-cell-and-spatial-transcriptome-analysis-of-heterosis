# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

# ============================================================
# run_sctenifold_knk_safe.R
#
# 作用：
#   对某一个 group 中的某一个 gene 做 scTenifoldKnk 虚拟敲除。
#
# 为什么需要 safe 版本：
#   之前 B73xA554 用原始 scTenifoldKnk 跑时，manifoldAlignment 出现：
#     Error in E$vectors[, seq_len(d)] : incorrect number of dimensions
#
#   原因是 manifold alignment 后只剩 1 个有效对齐维度，
#   原包代码默认按 2 维矩阵取值，R 自动降维后就报错。
#
# 本脚本做的修正：
#   核心流程仍然和 scTenifoldKnk 一致，
#   只是把 manifoldAlignment 改成 safe_manifold_alignment()，
#   在只剩 1 个有效维度时保留 1 维继续计算。
#
# 适用对象：
#   主要用于 B73xA554。
#
# 输入：
#   new_ko_genes_scripts/sctenifold_inputs/countMatrix_<group>.rds
#
# 输出：
#   new_ko_genes_scripts/sctenifold_results/<group>__<gene>.safe_ma1.*
# ============================================================


# ------------------------------------------------------------
# 1. 读取 group 和 gene 参数
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_sctenifold_knk_safe.R <group> <gene>")
}

group <- args[[1]]
gene <- args[[2]]


# ------------------------------------------------------------
# 2. 读取 Slurm CPU 核数并设置线程
# ------------------------------------------------------------
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


# ------------------------------------------------------------
# 3. 加载包
# ------------------------------------------------------------
# scTenifoldKnk：提供 QC、差异调控等内部函数。
# scTenifoldNet：提供网络构建、张量分解和原始 manifoldAlignment 思路。
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(scTenifoldKnk)
  library(scTenifoldNet)
})


# ------------------------------------------------------------
# 4. safe manifold alignment 函数
# ------------------------------------------------------------
# 这个函数基本照 scTenifoldNet::manifoldAlignment() 的逻辑写，
# 但关键位置加了 drop = FALSE，并允许有效维度不足 2 时使用 1 维。
#
# 输入：
#   X：WT 网络矩阵
#   Y：KO 网络矩阵
#   d：希望保留的维度，这里默认 1
#
# 输出：
#   alignedNet：对齐后的 WT/KO manifold 坐标
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

  # 原包这里如果只剩 1 列，R 可能会降维成向量。
  # 加 drop = FALSE 可以强制保留矩阵结构。
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


# ------------------------------------------------------------
# 5. safe 版本的 scTenifoldKnk 主流程
# ------------------------------------------------------------
# 这一段把 scTenifoldKnk() 的主要流程展开写出来：
#   1. scQC 过滤细胞
#   2. 去掉低检测率基因
#   3. makeNetworks 构建多个网络
#   4. tensorDecomposition 合并网络
#   5. 构建 KO 网络
#   6. safe manifold alignment
#   7. dRegulation 计算差异调控
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


# ------------------------------------------------------------
# 6. 读取输入矩阵
# ------------------------------------------------------------
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


# ------------------------------------------------------------
# 7. 设置输出路径并打印日志
# ------------------------------------------------------------
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


# ------------------------------------------------------------
# 8. 执行 safe 虚拟敲除并保存结果
# ------------------------------------------------------------
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
