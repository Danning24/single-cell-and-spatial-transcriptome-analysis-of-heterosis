project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

setwd(project_root)
library(Seurat)
Spatial <- readRDS(file.path(project_root, "data", "all_spatial_merge_filtered_integrated_remove_noise_bin50_annotation.RDS"))
Spatail <- UpdateSeuratObject(Spatail)
DimPlot(Spatail)
table(Spatail$sample)
Spatail@images

library("dplyr")
library(stringr)

spatial_heterotic_genes_results$cluster<-stringr::str_to_upper(spatial_heterotic_genes_results$cluster)
table(spatial_heterotic_genes_results$cluster)
nonadd_cortex <- spatial_heterotic_genes_results %>%
  filter(pattern == "Non-additive") %>%
  group_by(Trait) %>%
  filter(
    all(str_detect(cluster, "CORTEX"))
  ) %>%
  ungroup()

unique(nonadd_cortex$Trait)

result <- nonadd_cortex %>%
  group_by(Trait) %>%
  filter(all(MPH.pct > 0)) %>%
  ungroup()

upnonadd <- unique(result$Trait)

result <- nonadd_cortex %>%
  group_by(Trait) %>%
  filter(all(MPH.pct < 0)) %>%
  ungroup()

downnonadd <- unique(result$Trait)

Spatial_cortex <- subset(Spatail,idents = c("Cortex","Outer_cortex","Inner_cortex"))
Spatial_cortex <- AddModuleScore(Spatial_cortex,
                                 features = list(upnonadd,downnonadd,c(unique(nonadd_cortex$Trait))),
                                 name = c("upNAG_Scores",
                                          "downNAG_Scores",
                                          "NAG_Scores"))

Spatial_cortex <- subset(Spatial_cortex,subset = sample %in% c("hybrid_slice_1_bin50",
                                                               "hybrid_slice_2_bin50","hybrid_slice_3_bin50","hybrid_slice_4_bin50",
                                                               "hybrid_slice_5_bin50","hybrid_slice_6_bin50","hybrid_slice_7_bin50",
                                                               "hybrid_slice_8_bin50","hybrid_slice_9_bin50","hybrid_slice_10_bin50",
                                                               "hybrid_slice_11_bin50"))

DimPlot(Spatial_cortex)
png("A554xB73_Cortex_SpatialDimplot.png",height = 1700,width = 2500,res = 400)
SpatialDimPlot(Spatial_cortex,pt.size.factor = 10,image.alpha = 0.8,ncol = 5)&NoLegend()
dev.off()

library(ggplot2)
library(ggplot2)
png("A554xB73_Cortex_upNonaddFeatureplot.png",height = 1800,width = 4000,res = 400)
SpatialFeaturePlot(Spatial_cortex,features = "upNAG_Scores1",pt.size.factor = 10,
                   image.alpha = 0.8,ncol = 5)&theme(legend.key.size = unit(0.2,"cm"),
                                                     legend.text = element_text(size = 4),
                                                     legend.title = element_text(size = 6),
                                                     legend.direction = "vertical" )&
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(1.2, "cm"),
      barheight = unit(0.2, "cm")
    )
  )
dev.off()

png("A554xB73_Cortex_downNonaddFeatureplot.png",height = 1800,width = 4000,res = 400)
SpatialFeaturePlot(Spatial_cortex,features = "downNAG_Scores2",pt.size.factor = 10,
                   image.alpha = 0.8,ncol = 5)&theme(legend.key.size = unit(0.2,"cm"),
                                                     legend.text = element_text(size = 4),
                                                     legend.title = element_text(size = 6),
                                                     legend.direction = "vertical" )&
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(1.2, "cm"),
      barheight = unit(0.2, "cm")
    )
  )
dev.off()

png("A554xB73_Cortex_NonaddFeatureplot.png",height = 1800,width = 4000,res = 400)
SpatialFeaturePlot(Spatial_cortex,features = "NAG_Scores3",pt.size.factor = 10,
                   image.alpha = 0.8,ncol = 5)&theme(legend.key.size = unit(0.2,"cm"),
                                                     legend.text = element_text(size = 4),
                                                     legend.title = element_text(size = 6),
                                                     legend.direction = "vertical" )&
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(1.2, "cm"),
      barheight = unit(0.2, "cm")
    )
  )
dev.off()

p1 <- VlnPlot(Spatial_cortex,features = "upNAG_Scores1")
p2 <- VlnPlot(Spatial_cortex,features = "downNAG_Scores2")
p3 <- VlnPlot(Spatial_cortex,features = "NAG_Scores3")

p1+NoLegend()+theme(axis.title.x = element_blank())|p2+NoLegend()+theme(axis.title.x = element_blank())|p3+theme(axis.title.x = element_blank())
