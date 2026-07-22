# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
message("Arguments: ", paste(args, collapse = " | "))
project_root <- if (length(args) >= 1L) args[[1]] else getwd()
input_dir <- file.path(project_root, "01_Cortex数据")
output_dir <- file.path(project_root, "project_archive_work", "source_inputs", "fig1b")
if (length(args) >= 2L) {
  input_dir <- args[[2]]
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(input_dir, "GSE172277_RootSlices_genecounts.txt.gz")
cluster_file <- file.path(input_dir, "DEG_union_kmeans_k4.csv")
module_file <- file.path(input_dir, "EZ_M_Genes_from_clusters_with_clusterinfo.csv")
stopifnot(file.exists(counts_file), file.exists(cluster_file), file.exists(module_file))

counts_table <- fread(counts_file, header = TRUE, data.table = TRUE, check.names = FALSE)
gene_ids <- counts_table[[1L]]
counts_table[[1L]] <- NULL
counts_matrix <- as.matrix(counts_table)
rownames(counts_matrix) <- gene_ids
storage.mode(counts_matrix) <- "integer"
stopifnot(ncol(counts_matrix) == 32L)

sample_annotation <- data.table(
  sample = colnames(counts_matrix),
  position = rep(paste0("S", 1:16), each = 2),
  replicate = rep(c("rep1", "rep2"), times = 16)
)
sample_annotation[, zone := fifelse(
  position %chin% paste0("S", 1:7), "M",
  fifelse(position == "S8", "TZ", "EZ")
)]
sample_annotation[, position_index := as.integer(sub("S", "", position))]
setcolorder(sample_annotation, c("sample", "position", "position_index", "replicate", "zone"))
fwrite(sample_annotation, file.path(output_dir, "Fig1B_sample_position_annotation.csv"))

coldata <- as.data.frame(sample_annotation[, .(position, replicate, zone)])
rownames(coldata) <- sample_annotation$sample
coldata$position <- factor(coldata$position, levels = paste0("S", 1:16))
coldata$replicate <- factor(coldata$replicate)
coldata$zone <- factor(coldata$zone, levels = c("EZ", "TZ", "M"))

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~zone
)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
dds <- DESeq(dds, quiet = TRUE)
vsd <- vst(dds, blind = TRUE)
vst_matrix <- assay(vsd)

clusters <- fread(cluster_file)
stopifnot(all(c("gene", "cluster") %in% names(clusters)), !anyDuplicated(clusters$gene))
genes <- intersect(clusters$gene, rownames(vst_matrix))
stopifnot(length(genes) == nrow(clusters))
vst_matrix <- vst_matrix[genes, , drop = FALSE]

position_means <- vapply(paste0("S", 1:16), function(position_name) {
  rowMeans(vst_matrix[, sample_annotation$position == position_name, drop = FALSE])
}, FUN.VALUE = numeric(nrow(vst_matrix)))
rownames(position_means) <- genes

z_matrix <- t(scale(t(position_means)))
z_matrix[is.na(z_matrix)] <- 0

mean_table <- data.table(gene = rownames(position_means), position_means)
z_table <- data.table(gene = rownames(z_matrix), z_matrix)
setnames(mean_table, old = names(mean_table)[-1], new = paste0("S", 1:16))
setnames(z_table, old = names(z_table)[-1], new = paste0("S", 1:16))
mean_table <- clusters[mean_table, on = "gene"]
z_table <- clusters[z_table, on = "gene"]
setorder(mean_table, cluster, gene)
setorder(z_table, cluster, gene)

fwrite(mean_table, file.path(output_dir, "Fig1B_position_mean_VST_matrix.csv"))
fwrite(z_table, file.path(output_dir, "Fig1B_heatmap_gene_Zscore_matrix.csv"))
fwrite(clusters[order(cluster, gene)], file.path(output_dir, "Fig1B_gene_cluster_assignments.csv"))
fwrite(fread(module_file)[order(group, cluster, gene)],
       file.path(output_dir, "Fig1B_M_EZ_module_gene_sets.csv"))

quality <- data.table(
  metric = c(
    "raw_gene_count", "raw_sample_count", "genes_after_low_expression_filter",
    "heatmap_gene_count", "heatmap_position_count", "nonfinite_z_values"
  ),
  value = c(
    nrow(counts_matrix), ncol(counts_matrix), nrow(dds), nrow(z_matrix), ncol(z_matrix),
    sum(!is.finite(z_matrix))
  )
)
fwrite(quality, file.path(output_dir, "Fig1B_rebuild_QC.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "Fig1B_R_sessionInfo.txt"))
print(quality)
