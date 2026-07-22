# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

# ============================================================
# prepare_inputs_annotated.R
#
# 作用：
#   为新的候选基因准备 scTenifoldKnk 输入矩阵。
#
# 本脚本只做“准备输入数据”，不直接运行 scTenifoldKnk。
# 它会：
#   1. 读取 scRNA Seurat RDS。
#   2. 提取 B73xA554 和 B73xOh40B 两个组。
#   3. 用 SCT 网络基因作为背景。
#   4. 从 RNA counts 里补入新的 KO 基因。
#   5. 输出两个 countMatrix 和一个任务表。
#
# 推荐运行目录：
#   ${MAIZE_PROJECT_ROOT}
#
# 运行方式：
#   setwd("${MAIZE_PROJECT_ROOT}")
#   source("new_ko_genes_scripts/prepare_inputs_annotated.R")
# ============================================================


# ------------------------------------------------------------
# 1. 设置基础输入文件
# ------------------------------------------------------------
# 这是最开始的 scRNA Seurat 对象。
# 它里面包含 RNA counts、SCT 和 meta.data。
scRNA_rds <- "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"


# ------------------------------------------------------------
# 2. 设置要分析的两个组
# ------------------------------------------------------------
# 这两个名字必须和 meta.data$orig.ident 中的名字完全一致。
groups <- c("B73xA554", "B73xOh40B")


# ------------------------------------------------------------
# 3. 设置新的 KO 基因列表
# ------------------------------------------------------------
# 跑 scTenifoldKnk 时真正需要的只有 geneID。
# 功能分类、星级和文字说明不写进脚本，避免运行脚本太复杂。
ko_genes <- c(
  "Zm00001d011969",
  "Zm00001d028307",
  "Zm00001d034109",
  "Zm00001d006753",
  "Zm00001d011208",
  "Zm00001d018973",
  "Zm00001d048021",
  "Zm00001d015852",
  "Zm00001d013410",
  "Zm00001d009858",
  "Zm00001d028392",
  "Zm00001d019565",
  "Zm00001d037547",
  "Zm00001d039221",
  "Zm00001d012221",
  "Zm00001d019283",
  "Zm00001d030725",
  "Zm00001d026577",
  "Zm00001d002275",
  "Zm00001d043551",
  "Zm00001d051804",
  "Zm00001d030557"
)

# 去重，避免同一个基因重复提交。
ko_genes <- unique(ko_genes)
cat("N_KO_GENES_TOTAL\t", length(ko_genes), "\n", sep = "")


# ------------------------------------------------------------
# 4. 创建新的输出目录
# ------------------------------------------------------------
# 为了不覆盖之前 6 个基因的输入矩阵，这里使用新目录。
outdir <- file.path("new_ko_genes_scripts", "sctenifold_inputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 5. 加载 Matrix 包
# ------------------------------------------------------------
# 用于处理稀疏矩阵。
suppressPackageStartupMessages({
  library(Matrix)
})


# ------------------------------------------------------------
# 6. 读取 scRNA Seurat RDS
# ------------------------------------------------------------
if (!file.exists(scRNA_rds)) {
  stop("Cannot find scRNA RDS: ", scRNA_rds)
}

obj <- readRDS(scRNA_rds)


# ------------------------------------------------------------
# 7. 直接从 Seurat 对象内部取 meta.data 和 assays
# ------------------------------------------------------------
# 这里不用 Seurat::GetAssayData()，是为了减少对 SeuratObject 版本的依赖。
obj_attrs <- attributes(obj)
md <- obj_attrs[["meta.data"]]
assays <- obj_attrs[["assays"]]

if (is.null(md)) {
  stop("Cannot find meta.data in Seurat object")
}
if (!"orig.ident" %in% colnames(md)) {
  stop("metadata column orig.ident not found")
}
if (!"RNA" %in% names(assays)) {
  stop("RNA assay not found")
}


# ------------------------------------------------------------
# 8. 从 RNA assay 中取原始 counts layer
# ------------------------------------------------------------
# Seurat v5 常见结构：assays$RNA$layers$counts
rna_attrs <- attributes(assays[["RNA"]])
rna_layers <- rna_attrs[["layers"]]

if (is.null(rna_layers) || !"counts" %in% names(rna_layers)) {
  stop("RNA counts layer not found")
}

mat <- rna_layers[["counts"]]


# ------------------------------------------------------------
# 9. 工具函数：安全获取 dimnames
# ------------------------------------------------------------
get_dimnames <- function(x) {
  dn <- attributes(x)[["Dimnames"]]
  if (is.null(dn)) dn <- attributes(x)[["dimnames"]]
  if (is.null(dn)) dn <- dimnames(x)
  dn
}


# ------------------------------------------------------------
# 10. 工具函数：安全获取矩阵维度
# ------------------------------------------------------------
dims_of <- function(x) {
  d <- attributes(x)[["Dim"]]
  if (is.null(d)) d <- attributes(x)[["dim"]]
  if (is.null(d)) d <- dim(x)
  d
}


# ------------------------------------------------------------
# 11. 取 RNA counts 的基因名和细胞名
# ------------------------------------------------------------
# 这个对象中，RNA counts 的行名/列名存放在 RNA assay 的 features/cells 里。
rna_feature_map <- rna_attrs[["features"]]
rna_cell_map <- rna_attrs[["cells"]]

rna_feature_names <- get_dimnames(rna_feature_map)[[1]]
rna_cell_names <- get_dimnames(rna_cell_map)[[1]]

if (is.null(rna_feature_names) || is.null(rna_cell_names)) {
  stop("RNA assay feature/cell names are missing")
}


# ------------------------------------------------------------
# 12. 给 RNA counts 补上行名/列名，并转成 dgCMatrix
# ------------------------------------------------------------
dimnames(mat) <- list(rna_feature_names, rna_cell_names)
mat <- as(mat, "dgCMatrix")
cat("RNA_COUNTS_DIM\t", paste(dims_of(mat), collapse = "x"), "\n", sep = "")


# ------------------------------------------------------------
# 13. 从 SCT 中获取网络背景基因
# ------------------------------------------------------------
# 优先使用 SCT 的 var.features。
# 这些基因作为网络背景，避免直接使用 27780 个基因导致计算量太大。
feature_set <- character(0)

if ("SCT" %in% names(assays)) {
  sct_attrs <- attributes(assays[["SCT"]])

  if ("var.features" %in% names(sct_attrs)) {
    feature_set <- sct_attrs[["var.features"]]
  }

  for (slot_name in intersect(c("data", "counts", "scale.data"), names(sct_attrs))) {
    if (length(feature_set) > 0) break

    sct_dn <- get_dimnames(sct_attrs[[slot_name]])
    if (!is.null(sct_dn) && length(sct_dn[[1]]) > 0) {
      feature_set <- sct_dn[[1]]
      break
    }
  }
}


# ------------------------------------------------------------
# 14. 如果 SCT 中没有取到基因，则选 RNA counts 中检测细胞最多的前 5000 个基因
# ------------------------------------------------------------
# 这是备用方案，正常情况下不会走到这里。
if (length(feature_set) == 0) {
  detected <- Matrix::rowSums(mat > 0)
  feature_set <- rna_feature_names[order(detected, decreasing = TRUE)][seq_len(min(5000, length(detected)))]
}

cat("N_NETWORK_FEATURES_BEFORE_ADD_KO\t", length(feature_set), "\n", sep = "")


# ------------------------------------------------------------
# 15. 检查新的 KO 基因是否存在于 RNA counts
# ------------------------------------------------------------
ko_present_in_rna <- ko_genes[ko_genes %in% rna_feature_names]
ko_missing_in_rna <- setdiff(ko_genes, rna_feature_names)

cat("N_KO_PRESENT_IN_RNA\t", length(ko_present_in_rna), "\n", sep = "")
cat("N_KO_MISSING_IN_RNA\t", length(ko_missing_in_rna), "\n", sep = "")

if (length(ko_missing_in_rna) > 0) {
  cat("KO_MISSING_IN_RNA\t", paste(ko_missing_in_rna, collapse = ","), "\n", sep = "")
}


# ------------------------------------------------------------
# 16. 强制把存在于 RNA counts 的 KO 基因加入网络基因集合
# ------------------------------------------------------------
# scTenifoldKnk 的 gKO 必须在 countMatrix 的行名中。
feature_set <- unique(c(feature_set, ko_present_in_rna))


# ------------------------------------------------------------
# 17. 把 feature_set 转换成 RNA counts 的行索引
# ------------------------------------------------------------
row_idx <- match(feature_set, rna_feature_names)
row_idx <- row_idx[!is.na(row_idx)]
feature_set <- rna_feature_names[row_idx]

cat("N_FEATURES_FINAL\t", length(feature_set), "\n", sep = "")
cat("KO_IN_FINAL_FEATURE_SET\t", paste(ko_genes[ko_genes %in% feature_set], collapse = ","), "\n", sep = "")


# ------------------------------------------------------------
# 18. 写出 KO gene table
# ------------------------------------------------------------
# 只保留运行需要的信息：gene 和是否存在于输入矩阵。
ko_info <- data.frame(
  gene = ko_genes,
  present_in_RNA = ko_genes %in% rna_feature_names,
  present_in_final_matrix = ko_genes %in% feature_set,
  stringsAsFactors = FALSE
)

write.table(
  ko_info,
  file = file.path(outdir, "ko_gene_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. 准备 manifest 和 job table
# ------------------------------------------------------------
manifest <- data.frame(
  group = character(),
  matrix_file = character(),
  n_genes = integer(),
  n_cells = integer(),
  present_ko_genes = character(),
  stringsAsFactors = FALSE
)

job_table <- data.frame(
  group = character(),
  gene = character(),
  present_in_matrix = logical(),
  suggested_runner = character(),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 20. 对两个 group 分别提取 countMatrix
# ------------------------------------------------------------
for (grp in groups) {
  cells <- rownames(md)[as.character(md[["orig.ident"]]) == grp]
  cells <- intersect(cells, rna_cell_names)

  if (length(cells) == 0) {
    stop("no cells found for group: ", grp)
  }

  col_idx <- match(cells, rna_cell_names)
  count_matrix <- mat[row_idx, col_idx, drop = FALSE]

  out_file <- file.path(outdir, paste0("countMatrix_", grp, ".rds"))
  saveRDS(count_matrix, out_file, compress = "xz")

  present <- ko_genes[ko_genes %in% rownames(count_matrix)]

  manifest <- rbind(
    manifest,
    data.frame(
      group = grp,
      matrix_file = out_file,
      n_genes = nrow(count_matrix),
      n_cells = ncol(count_matrix),
      present_ko_genes = paste(present, collapse = ","),
      stringsAsFactors = FALSE
    )
  )

  # B73xA554 之前出现过 manifoldAlignment 维度退化问题，建议用 safe runner。
  # B73xOh40B 之前标准 runner 正常完成。
  suggested_runner <- if (grp == "B73xA554") {
    "run_sctenifold_knk_safe.R"
  } else {
    "run_sctenifold_knk.R"
  }

  for (g in ko_genes) {
    job_table <- rbind(
      job_table,
      data.frame(
        group = grp,
        gene = g,
        present_in_matrix = g %in% rownames(count_matrix),
        suggested_runner = suggested_runner,
        stringsAsFactors = FALSE
      )
    )
  }

  cat("EXTRACTED\t", grp, "\t", nrow(count_matrix), "\t", ncol(count_matrix), "\t", out_file, "\n", sep = "")
}


# ------------------------------------------------------------
# 21. 写出 manifest.tsv
# ------------------------------------------------------------
write.table(
  manifest,
  file = file.path(outdir, "manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. 写出 job table
# ------------------------------------------------------------
# 这个表就是后续提交任务的总表。
# 理论上应有：2 groups x 22 genes = 44 rows。
write.table(
  job_table,
  file = file.path(outdir, "job_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. 写出实际可运行基因列表
# ------------------------------------------------------------
writeLines(
  ko_genes[ko_genes %in% feature_set],
  file.path(outdir, "ko_genes_present.txt")
)


# ------------------------------------------------------------
# 24. 打印总结
# ------------------------------------------------------------
cat("OUTPUT_DIR\t", outdir, "\n", sep = "")
cat("N_GROUPS\t", length(groups), "\n", sep = "")
cat("N_JOBS_TOTAL\t", nrow(job_table), "\n", sep = "")
cat("N_JOBS_PRESENT_IN_MATRIX\t", sum(job_table$present_in_matrix), "\n", sep = "")
cat("DONE\n")
