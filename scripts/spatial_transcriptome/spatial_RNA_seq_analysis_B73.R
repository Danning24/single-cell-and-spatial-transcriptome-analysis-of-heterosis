library(Seurat)
library(clustree)
library(dplyr)
library(ggplot2)



################################################################################
#   processing, clustering and finding markers for B73 spatial transcriptome
################################################################################


# List files 
setwd('~/spatial_RNA/intermediate_data/B73/bin50/')
files <- list.files(path="./", pattern="bin50.rds")

# Create a list of count matrices
rds_read <- lapply(paste0("./", files), readRDS)

# Assign names
names(rds_read) <- sub("\\.[^.]*$", "", files)

# merge multiple seurat objects to one object with multiple layers
sp_merge <- merge(rds_read[[1]], y = rds_read[2:length(rds_read)], 
                  add.cell.ids = names(rds_read), project="Heterosis_spatial_RNA_bin50")

head(sp_merge@meta.data)
Assays(sp_merge)
Layers(sp_merge)

sp_merge$slice <- sub('_[^_]*$', '', colnames(sp_merge))

# re-join layers to combine counts matrix
sp_merge[["Spatial"]] <- JoinLayers(sp_merge[["Spatial"]])

# visulization of filtering threshold
metadata <- sp_merge@meta.data

setwd('~/spatial_RNA/results/bin50/B73/new_joined_layers/')

pdf(file = './B73_bin50_nfeature_ncount_plot.pdf', width = 6, height = 6)
ggplot(metadata) +
  geom_point(aes(x = nCount_Spatial, y = nFeature_Spatial), shape=21, alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  geom_vline(xintercept = 150, color="red", linetype="dotted")+
  geom_hline(yintercept = 200, color="red", linetype="dotted")
dev.off()

# use the same parameters across all samples
sp_merge_filtered <- subset(sp_merge, nFeature_Spatial > 200)

# slice G1 and G9 contains lower number of cells, so remove them.
sp_merge_filtered <- subset(sp_merge_filtered, !slice %in% c('A04337G3_G1.bin50', 'A04337G3_G9.bin50'))

# Filter genes that are expressed in at least 3 cells (min.cells) and have at least 10 counts (min.count)
sp_merge_filtered <- subset(sp_merge_filtered, features = rownames(sp_merge_filtered)[
  rowSums(sp_merge_filtered@assays$Spatial$counts > 0) >= 3 & rowSums(sp_merge_filtered@assays$Spatial$counts) >= 10])

# split layers for following integration, this step is very neccessary
sp_merge_filtered[["Spatial"]] <- split(sp_merge_filtered[["Spatial"]], f = sp_merge_filtered$slice)
# normalization step
sp_merge_filtered <- SCTransform(sp_merge_filtered, assay = "Spatial", verbose = FALSE)
# by default, PCA is run only 3000 variable features
sp_merge_filtered <- RunPCA(sp_merge_filtered, verbose = FALSE, assay="SCT")
ElbowPlot(sp_merge_filtered, ndims = 50)

# this is very important for integration!
for (i in 1:10) {
  slot(object = sp_merge_filtered@assays$SCT@SCTModel.list[[i]], name="umi.assay") <- "Spatial"
}

# perform integration across slices
sp_merge_filtered_int <- IntegrateLayers(object = sp_merge_filtered,
                                         method = CCAIntegration,
                                         orig.reduction = "pca",
                                         new.reduction = "cca.integration",
                                         normalization.method = "SCT",
                                         verbose = FALSE)
# clustering spatial bins
sp_merge_filtered_int <- FindNeighbors(sp_merge_filtered_int, dims = 1:30)
sp_merge_filtered_int <- FindClusters(sp_merge_filtered_int, resolution = seq(1.5, 2, 0.1))

Idents(sp_merge_filtered_int) <- sp_merge_filtered_int$SCT_snn_res.1.7
sp_merge_filtered_int <- RunUMAP(sp_merge_filtered_int, dims = 1:30)

# decide the resolution = 1.7
p <- clustree(sp_merge_filtered_int, prefix = "SCT_snn_res.")
ggsave(p, filename = "./QC_clustree2.pdf", height = 10, width = 12)

# UMAP plot
plot1 <- DimPlot(sp_merge_filtered_int, reduction = "umap", label = TRUE) 
ggsave('./B73_spatial_integrated_slices_bin50_UMAP_plot_res_1_7.pdf', 
       plot1, width = 7, height = 6)

tissue.col <- c("#C9A1FA", "#EDC66A", '#B46DA9',"#F9B3AD", '#E87B1E', "#018B38",
                "#D9DEE7", "#8498AB", "blue", "green", "yellow", "royalblue",
                 "red", "#46f0f0", "#bcf60c")
# spatial root section plot
plot2 <- SpatialDimPlot(sp_merge_filtered_int, label = TRUE, label.size = 1, 
                        pt.size.factor = 1, ncol = 3) +
                        scale_fill_manual(values = tissue.col)
# remove cluster 5 since they are mostly outside the edge of image
sp_merge_filtered_int <- subset(sp_merge_filtered_int, seurat_clusters != 5)

ggsave('./B73_spatial_integrated_spatial_plot_res_1_7.pdf', 
       plot2, width = 12, height = 15)

saveRDS(sp_merge_filtered_int, '~/spatial_RNA/intermediate_data/sp_B73_integrated_bin50_new.RDS')

#####   annoate spatial clusters with LCM markers  #####

# read differentially expressed genes from LCM transcriptome
all.DE.genes <- readxl::read_xlsx('~/LCM_RNA/results/DEGenes/glmQLfit/LCM_glmQLfit_all_DEGs.xlsx')
plot.genes <- all.DE.genes[all.DE.genes$geneID %in% rownames(sp_merge_filtered_int), ]

# choose top 50 LCM differential genes for each tissue
top_50_LCM_genes <- plot.genes %>%
  group_by(Tissue) %>%
  slice_max(order_by = logCPM, n = 50, with_ties = FALSE) %>%
  ungroup()

# average expression across clusters
avg_exp <- AverageExpression(sp_merge_filtered_int, 
                             features = top_50_LCM_genes$geneID, 
                             return.seurat = FALSE)
# Get expression matrix
avg_matrix <- avg_exp$SCT 
rownames(avg_matrix) <- top_50_LCM_genes$Tissue

# visulize gene expression in each cluster
p3 = pheatmap(avg_matrix, 
              scale = "row",              # scale gene-wise
              cluster_rows = F, cluster_cols = F,
              color = colorRampPalette(c("blue", "white", "red"))(100))

ggsave(paste0('./sp_B73_slice_integrated_anno_by_LCM_markers_res_1_7.pdf'), 
       p3, width = 10, height = 15)

# find all markers of spatial transcriptome
sp_merge_filtered_int <- PrepSCTFindMarkers(sp_merge_filtered_int)
sp.all.markers <- c()

for (i in 0:13) {
  print(i)
  all.markers <- FindMarkers(sp_merge_filtered_int, 
                             ident.1 = i, 
                             assay = 'SCT',
                             only.pos = TRUE,
                             min.pct = 0.15, 
                             logfc.threshold = 0.25,
                             min.diff.pct = 0.1,
                             verbose = FALSE)
  sp.all.markers <- rbind.data.frame(sp.all.markers, cbind(all.markers, cluster=i, markerID=rownames(all.markers)))
}
sp.all.markers <- sp.all.markers[sp.all.markers$p_val_adj <0.05, ]

write_xlsx(sp.all.markers, './cluster_annotation/sp_B73_all_markers_res1_7.xlsx')

