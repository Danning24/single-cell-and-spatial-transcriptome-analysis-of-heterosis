#!/usr/bin/env Rscript
# Figure 3A-B: cortex re-clustering and group-specific UMAP density.

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(viridis)
library(patchwork)

# Inputs and outputs
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
input_file <- file.path(project_dir, "0.data/scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS")
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", unset = "results"), "cortex_clustering")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Re-cluster the jointly integrated cortex population.
sce <- readRDS(input_file)
cortex <- subset(sce, subset = Annotation == "Cortex")
rm(sce)
cortex <- SCTransform(cortex, assay = "RNA", vst.flavor = "v2", variable.features.n = 3000)
cortex <- RunPCA(cortex, npcs = 30)
cortex <- FindNeighbors(cortex, dims = 1:20, k.param = 20, prune.SNN = 1 / 15)
cortex <- FindClusters(cortex, algorithm = 1, resolution = 0.04)
cortex <- RunUMAP(cortex, dims = 1:15, metric = "cosine", 
                  n.neighbors = 30, min.dist = 0.3, spread = 1)

# Labels follow the current manuscript: raw 1 is the lower branch.
cluster_labels <- c("0" = "Cortex1", "1" = "Cortex3", "2" = "Cortex2")
cortex$Cortex_subcluster <- factor(unname(cluster_labels[as.character(cortex$seurat_clusters)]),
                                    levels = c("Cortex1", "Cortex2", "Cortex3"))
cortex$orig.ident <- sub("_Rep[123]$", "", as.character(cortex$orig.ident))
parent_group <- c(A554 = "NSS", Oh40B = "NSS", Oh43 = "NSS",
                  A680 = "SS", B37 = "SS", NC328 = "SS",
                  NC298 = "TL", NC320 = "TL", NC352 = "TL")
parent <- sub("^B73x", "", cortex$orig.ident)
cortex$group <- unname(parent_group[parent])
cortex$group[cortex$orig.ident == "B73"] <- "SS"
cortex$group_FP <- ifelse(grepl("x", cortex$orig.ident), "F1",
                           ifelse(cortex$orig.ident == "B73", "B73", "P1"))

# Clustering and cell-density panels use the same UMAP embedding.
cluster_plot <- DimPlot(cortex, group.by = "Cortex_subcluster", label = TRUE) +
  scale_colour_manual(values = c(Cortex1 = "#E64B35", Cortex2 = "#00A087", Cortex3 = "#4DBBD5")) +
  NoAxes()
ggsave(file.path(output_dir, "cortex_subclusters.pdf"), cluster_plot, width = 6, height = 5)

umap <- as.data.frame(Embeddings(cortex, "umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")
umap$group <- cortex$group
density_plots <- lapply(c("NSS", "SS", "TL"), function(group_name) {
  ggplot(filter(umap, group == group_name), aes(UMAP_1, UMAP_2)) +
    geom_pointdensity(size = 0.01, adjust = 1) +
    scale_colour_viridis(option = "magma") +
    coord_fixed() +
    labs(title = group_name) +
    theme_void() +
    theme(plot.background = element_rect(fill = "black", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          plot.title = element_text(colour = "white"), legend.position = "none")
})
ggsave(file.path(output_dir, "cortex_group_density.pdf"),
       wrap_plots(density_plots, ncol = 3), width = 12, height = 4)

saveRDS(cortex, file.path(output_dir, "Cortex.rds"))
write.csv(data.frame(cell_id = colnames(cortex), orig.ident = cortex$orig.ident,
                      Cortex_subcluster = cortex$Cortex_subcluster,
                      group = cortex$group, group_FP = cortex$group_FP, umap[, 1:2]),
          file.path(output_dir, "cortex_metadata_umap.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
