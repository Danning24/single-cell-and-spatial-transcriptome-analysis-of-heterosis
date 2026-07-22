# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

setwd(project_root)
list.files()
counts <- read.table(
  "GSE172277_RootSlices_genecounts.txt.gz",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# 看一下结构
dim(counts)
head(counts)
# Dependency only; install outside the analysis script: BiocManager::install("apeglm")
suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
  library(pheatmap)
})


sample_names <- colnames(counts)

## 2) 样本信息：16个位置，每个2个重复 -> 3大区
position <- rep(paste0("S", 1:16), each = 2)
replicate <- rep(c("rep1","rep2"), times = 16)

# !!! 这里按你的真实“位置->区”改
pos2zone <- data.frame(
  position = paste0("S", 1:16),
  zone = c(rep("M", 7), rep("TZ", 1), rep("EZ", 8)),
  stringsAsFactors = FALSE
)
zone <- pos2zone$zone[match(position, pos2zone$position)]

coldata <- data.frame(
  row.names = sample_names,
  position = factor(position),
  replicate = factor(replicate),
  zone = factor(zone, levels = c("EZ","TZ","M"))
)
print(table(coldata$zone))  # 应该是 EZ/TZ/M 各自的样本数

## ---------------------------
## 3) Helper functions
## ---------------------------
pick_deg <- function(res_df, padj_cut=0.05, lfc_cut=1) {
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df$gene[res_df$padj < padj_cut & abs(res_df$log2FoldChange) > lfc_cut]
}

save_res <- function(res, name) {
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  write.csv(res, paste0(name, ".csv"), row.names = FALSE)
  res
}

## ---------------------------
## 4) Run DESeq2 with ref = EZ (get TZ_vs_EZ, M_vs_EZ)
## ---------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ zone)

# filter low expression (strongly recommended with small n)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds$zone <- relevel(dds$zone, ref = "EZ")
dds <- DESeq(dds)

# check coefficient names
print(resultsNames(dds))

res_TZ_vs_EZ <- lfcShrink(dds, coef = "zone_TZ_vs_EZ", type = "apeglm")
res_M_vs_EZ  <- lfcShrink(dds, coef = "zone_M_vs_EZ",  type = "apeglm")

res_TZ_vs_EZ_df <- save_res(res_TZ_vs_EZ, "DESeq2_TZ_vs_EZ_apeglm")
res_M_vs_EZ_df  <- save_res(res_M_vs_EZ,  "DESeq2_M_vs_EZ_apeglm")

## ---------------------------
## 5) Re-run with ref = M (get TZ_vs_M)
## ---------------------------
dds2 <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ zone)
dds2 <- dds2[keep, ]  # use same filtering for consistency

dds2$zone <- relevel(dds2$zone, ref = "M")
dds2 <- DESeq(dds2)

print(resultsNames(dds2))

res_TZ_vs_M <- lfcShrink(dds2, coef = "zone_TZ_vs_M", type = "apeglm")
res_TZ_vs_M_df <- save_res(res_TZ_vs_M, "DESeq2_TZ_vs_M_apeglm")

## ---------------------------
## 6) Union(unique) of DEGs from the three pairwise comparisons
##    Adjust thresholds as you like
## ---------------------------
padj_cut <- 0.05
lfc_cut  <- 1

deg_TZ_EZ <- pick_deg(res_TZ_vs_EZ_df, padj_cut, lfc_cut)
deg_M_EZ  <- pick_deg(res_M_vs_EZ_df,  padj_cut, lfc_cut)
deg_TZ_M  <- pick_deg(res_TZ_vs_M_df,  padj_cut, lfc_cut)

deg_union <- unique(c(deg_TZ_EZ, deg_M_EZ, deg_TZ_M))
cat("DEG union size =", length(deg_union), "\n")
writeLines(deg_union, "DEG_union_unique.txt")

## ---------------------------
## 7) Clustering (k=4) on VST (use dds with all samples)
##    VST from ANY of the fitted objects is fine; use dds here.
## ---------------------------
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

# subset to union DEGs (keep intersection to be safe)
deg_union2 <- intersect(deg_union, rownames(mat))
mat_deg <- mat[deg_union2, , drop = FALSE]

# z-score per gene (focus on pattern)
mat_z <- t(scale(t(mat_deg)))
mat_z[is.na(mat_z)] <- 0

set.seed(123)
km <- kmeans(mat_z, centers = 4, nstart = 50)

cluster_df <- data.frame(
  gene = rownames(mat_z),
  cluster = factor(km$cluster, levels = 1:4)
)
write.csv(cluster_df, "DEG_union_kmeans_k4.csv", row.names = FALSE)

## ---------------------------
## 8) Heatmap (optional)
## ---------------------------
ord <- order(cluster_df$cluster)

ann_col <- data.frame(zone = coldata$zone, position = coldata$position)
rownames(ann_col) <- rownames(coldata)

ann_row <- data.frame(cluster = cluster_df$cluster)
rownames(ann_row) <- cluster_df$gene

## Export genes per cluster
for (k in 1:4) {
  writeLines(cluster_df$gene[cluster_df$cluster == k], paste0("cluster_", k, "_genes.txt"))
}



mat <- assay(vsd)
deg_union2 <- intersect(cluster_df$gene, rownames(mat))
mat_deg <- mat[deg_union2, , drop = FALSE]
# 获取 position 因子
pos <- coldata$position

# 按位置求均值
mat_pos_mean <- sapply(levels(pos), function(p) {
  rowMeans(mat_deg[, pos == p, drop = FALSE])
})

# 转成矩阵
mat_pos_mean <- as.matrix(mat_pos_mean)

dim(mat_pos_mean)  # gene × 16
cluster_df2 <- cluster_df[cluster_df$gene %in% rownames(mat_pos_mean), ]

ord <- order(cluster_df2$cluster)
mat_z_pos <- t(scale(t(mat_pos_mean)))
mat_z_pos[is.na(mat_z_pos)] <- 0
ann_row <- data.frame(cluster = cluster_df2$cluster)
rownames(ann_row) <- cluster_df2$gene
# 假设列名是 "S1"..."S16"
cn <- colnames(mat_z_pos)

# 提取数字部分并排序
ord_col <- order(as.integer(sub("^S", "", cn)))

# 重新排列列顺序
mat_z_pos <- mat_z_pos[, ord_col, drop = FALSE]

# （可选）同时把列名规范一下
colnames(mat_z_pos) <- paste0("S", as.integer(sub("^S", "", colnames(mat_z_pos))))

pheatmap(
  mat_z_pos[ord, ],
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,  # 保持pos顺序
  annotation_row = ann_row,
  fontsize_col = 10,
  main = "Union DEGs (mean per position)"
)
my_colors <- colorRampPalette(c("#3B4CC0", "white", "#7F0000"))(100)

# 设置对称breaks（保证0是白色）
bk <- seq(-3, 3, length.out = 101)
pheatmap(
  mat_z_pos[ord, ],
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_row,
  color = my_colors,
  breaks = bk,
  fontsize_col = 10
)
ann_colors <- list(
  cluster = c(
    "1" = "#4C78A8",
    "2" = "#F58518",
    "3" = "#B279A2",
    "4" = "#54A24B"
  )
)

pheatmap(
  mat_z_pos[ord, ],
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  color = my_colors,
  breaks = bk,
  fontsize_col = 10
)


cluster_df
head(cluster_df)
cluster_df$cluster <- as.character(cluster_df$cluster)
EZ_genes <- cluster_df$gene[cluster_df$cluster == "2"]

M_genes  <- cluster_df$gene[cluster_df$cluster %in% c("3","4")]
gene_group_df <- data.frame(
  gene = c(EZ_genes, M_genes),
  group = c(
    rep("EZ_Genes", length(EZ_genes)),
    rep("M_Genes", length(M_genes))
  ),
  stringsAsFactors = FALSE
)
write.csv(gene_group_df,
          "EZ_M_Genes_from_clusters.csv",
          row.names = FALSE)
gene_group_df2 <- data.frame(
  gene = cluster_df$gene,
  cluster = cluster_df$cluster,
  group = ifelse(cluster_df$cluster == "2", "EZ_Genes",
                 ifelse(cluster_df$cluster %in% c("3","4"), "M_Genes", NA))
)

# 去掉 NA（比如 cluster1 不要）
gene_group_df2 <- gene_group_df2[!is.na(gene_group_df2$group), ]

write.csv(gene_group_df2,
          "EZ_M_Genes_from_clusters_with_clusterinfo.csv",
          row.names = FALSE)
table(gene_group_df$group)
