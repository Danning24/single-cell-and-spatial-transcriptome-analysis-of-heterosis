#!/usr/bin/env Rscript
# Figure 3C: derive root developmental programmes and score cortex cells.

library(DESeq2)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Inputs and outputs
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
root_slice_file <- Sys.getenv("MAIZE_ROOT_SLICE_COUNTS")
cortex_file <- Sys.getenv("MAIZE_CORTEX_RDS",
                          unset = file.path(project_dir, "4.Cortex_cells/data_output/Cortex.rds"))
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", unset = "results"), "developmental_scores")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# GSE172277: 16 positions, two libraries per position.
count_matrix <- as.matrix(read.delim(root_slice_file, row.names = 1, check.names = FALSE))
coldata <- data.frame(row.names = colnames(count_matrix),
                      position = rep(paste0("S", 1:16), each = 2),
                      replicate = rep(c("rep1", "rep2"), 16),
                      zone = factor(rep(c(rep("M", 7), "TZ", rep("EZ", 8)), each = 2),
                                    levels = c("EZ", "TZ", "M")))
dds <- DESeqDataSetFromMatrix(count_matrix, coldata, design = ~ zone)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- DESeq(dds[keep, ])
contrasts <- list(TZ_vs_EZ = lfcShrink(dds, coef = "zone_TZ_vs_EZ", type = "apeglm"),
                  M_vs_EZ = lfcShrink(dds, coef = "zone_M_vs_EZ", type = "apeglm"))
dds_m <- DESeqDataSetFromMatrix(count_matrix, coldata, design = ~ zone)
dds_m <- dds_m[keep, ]
dds_m$zone <- relevel(dds_m$zone, ref = "M")
dds_m <- DESeq(dds_m)
contrasts$TZ_vs_M <- lfcShrink(dds_m, coef = "zone_TZ_vs_M", type = "apeglm")

# Use the union of significant genes, then cluster across all 32 libraries.
deg_union <- unique(unlist(lapply(contrasts, function(result) {
  rownames(result)[!is.na(result$padj) & result$padj < 0.05 & abs(result$log2FoldChange) > 1]
})))
vst_matrix <- assay(vst(dds, blind = TRUE))[deg_union, , drop = FALSE]
gene_z <- t(scale(t(vst_matrix)))
gene_z[is.na(gene_z)] <- 0
set.seed(123)
gene_clusters <- kmeans(gene_z, centers = 4, nstart = 50)
programmes <- data.frame(gene = rownames(gene_z), cluster = gene_clusters$cluster) |>
  mutate(group = case_when(cluster == 2 ~ "EZ_Genes",
                            cluster %in% c(3, 4) ~ "M_Genes",
                            TRUE ~ NA_character_))
write.csv(programmes, file.path(output_dir, "developmental_gene_clusters.csv"), row.names = FALSE)
write.csv(filter(programmes, !is.na(group)),
          file.path(output_dir, "EZ_M_Genes_from_clusters_with_clusterinfo.csv"), row.names = FALSE)

# Score the two programmes in separate calls, as in the original analysis.
cortex <- readRDS(cortex_file)
DefaultAssay(cortex) <- "SCT"
m_genes <- intersect(filter(programmes, group == "M_Genes")$gene, rownames(cortex[["SCT"]]))
e_genes <- intersect(filter(programmes, group == "EZ_Genes")$gene, rownames(cortex[["SCT"]]))
cortex <- AddModuleScore(cortex, features = list(m_genes), assay = "SCT", slot = "data",
                         name = "M_Scores", nbin = 24, ctrl = 100, seed = 1)
cortex <- AddModuleScore(cortex, features = list(e_genes), assay = "SCT", slot = "data",
                         name = "E_Scores", nbin = 24, ctrl = 100, seed = 1)
score_plot <- FeaturePlot(cortex, features = c("M_Scores1", "E_Scores1"), combine = FALSE)
score_plot <- lapply(score_plot, function(plot) {
  plot + scale_colour_gradientn(colours = c("#5E3C99", "white", "#B22222")) + NoAxes()
})
ggsave(file.path(output_dir, "cortex_M_E_scores.pdf"),
       wrap_plots(score_plot, ncol = 2), width = 10, height = 4)
saveRDS(cortex, file.path(output_dir, "Cortex_M_E_scores.rds"))
write.csv(data.frame(cell_id = colnames(cortex), M_Score = cortex$M_Scores1,
                      E_Score = cortex$E_Scores1),
          file.path(output_dir, "cortex_M_E_scores.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
