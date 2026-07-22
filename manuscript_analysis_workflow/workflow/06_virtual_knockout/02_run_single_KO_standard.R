# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

# ============================================================
# run_sctenifold_knk.R
#
# 作用：
#   对某一个 group 中的某一个 gene 做 scTenifoldKnk 虚拟敲除。
#
# 适用对象：
#   主要用于 B73xOh40B。
#   之前 B73xOh40B 用标准 scTenifoldKnk 可以正常完成。
#
# 输入：
#   new_ko_genes_scripts/sctenifold_inputs/countMatrix_<group>.rds
#
# 输出：
#   new_ko_genes_scripts/sctenifold_results/<group>__<gene>.sctenifoldKnk.rds
#   new_ko_genes_scripts/sctenifold_results/<group>__<gene>.diffRegulation.csv
#   new_ko_genes_scripts/sctenifold_results/<group>__<gene>.sig_genes.csv
#   new_ko_genes_scripts/sctenifold_results/<group>__<gene>.sig_gene_counts.csv
#
# 运行方式示例：
#   R --slave -f new_ko_genes_scripts/run_sctenifold_knk.R --args B73xOh40B Zm00001d011969
# ============================================================


# ------------------------------------------------------------
# 1. 读取命令行参数
# ------------------------------------------------------------
# 这个脚本需要两个参数：
#   group：要分析的组，例如 B73xOh40B
#   gene ：要虚拟敲除的基因，例如 Zm00001d011969
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: run_sctenifold_knk.R <group> <gene>")
}

group <- args[[1]]
gene <- args[[2]]


# ------------------------------------------------------------
# 2. 读取 Slurm 分配的 CPU 核数
# ------------------------------------------------------------
# sbatch 脚本里用 -n 96 申请 96 核。
# 这里从 Slurm 环境变量里读取核数，并传给底层 BLAS/OpenMP。
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
# 3. 加载需要的 R 包
# ------------------------------------------------------------
# data.table：快速写出 csv。
# Matrix：处理单细胞稀疏表达矩阵。
# scTenifoldKnk：虚拟敲除核心包。
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(scTenifoldKnk)
})


# ------------------------------------------------------------
# 4. 读取当前 group 的输入表达矩阵
# ------------------------------------------------------------
# 输入矩阵由 prepare_inputs_annotated.R 提前生成。
# 矩阵格式：
#   行 = genes
#   列 = cells
#
# 例如：
#   new_ko_genes_scripts/sctenifold_inputs/countMatrix_B73xOh40B.rds
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


# ------------------------------------------------------------
# 5. 检查 KO gene 是否在矩阵中
# ------------------------------------------------------------
# scTenifoldKnk 要求 gKO 必须存在于 countMatrix 的 rownames 中。
# 如果不在，说明准备输入时没有把这个基因补进去，不能继续跑。
if (!gene %in% rownames(count_matrix)) {
  stop("KO gene not found in count matrix rownames: ", gene)
}


# ------------------------------------------------------------
# 6. 设置输出路径
# ------------------------------------------------------------
# 所有新基因结果都放在 new_ko_genes_scripts/sctenifold_results/。
# 不会覆盖之前 6 个基因的旧结果。
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


# ------------------------------------------------------------
# 7. 打印日志信息
# ------------------------------------------------------------
# 这些内容会写入 Slurm .out 文件，方便以后检查任务是否正常启动。
cat("START\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("GROUP\t", group, "\n", sep = "")
cat("GENE\t", gene, "\n", sep = "")
cat("MATRIX_DIM\t", nrow(count_matrix), "\t", ncol(count_matrix), "\n", sep = "")
cat("CORES\t", ncores, "\n", sep = "")


# ------------------------------------------------------------
# 8. 执行 scTenifoldKnk 虚拟敲除
# ------------------------------------------------------------
# set.seed() 用于提高可重复性。
# scTenifoldKnk 内部会：
#   1. 对细胞做 QC
#   2. 构建 WT 网络
#   3. 构建 KO 网络
#   4. 做 manifold alignment
#   5. 输出差异调控基因
set.seed(20260521)
result <- scTenifoldKnk(
  countMatrix = count_matrix,
  gKO = gene
)


# ------------------------------------------------------------
# 9. 保存完整结果对象
# ------------------------------------------------------------
# RDS 文件最完整，后续如果要重新画图或提取其他信息，可以读这个。
saveRDS(result, paste0(out_prefix, ".sctenifoldKnk.rds"), compress = "xz")

if (!"diffRegulation" %in% names(result)) {
  stop("Result object does not contain diffRegulation")
}


# ------------------------------------------------------------
# 10. 导出差异调控结果
# ------------------------------------------------------------
# diffRegulation 是主要结果表，通常包含：
#   gene, distance, Z, FC, p.value, p.adj
# 具体列名以 scTenifoldKnk 输出为准。
diffreg <- as.data.table(result$diffRegulation)
fwrite(diffreg, paste0(out_prefix, ".diffRegulation.csv"))


# ------------------------------------------------------------
# 11. 提取显著差异调控基因
# ------------------------------------------------------------
# 优先按 p.adj < 0.05 筛选。
# 如果没有 p.adj，则退而使用 p.value < 0.05。
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


# ------------------------------------------------------------
# 12. 尝试画 plotKO 图
# ------------------------------------------------------------
# 画图失败不影响主要结果，所以用 try 包起来。
try({
  png(paste0(out_prefix, ".plotKO.png"), width = 1800, height = 1400, res = 180)
  plotKO(result, gKO = gene)
  dev.off()
}, silent = TRUE)

cat("END\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
