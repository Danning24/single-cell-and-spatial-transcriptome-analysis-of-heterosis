library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)



################################################################################
#####      Integration of all samples for spatial transcriptome            #####
################################################################################


setwd('~/spatial_RNA/intermediate_data/')

files_mat <- list.files(path="./B73/bin50", pattern="bin50.rds")
files_pat <- list.files(path="./A554/slices", pattern="bin50.RDS")
files_hyb <- list.files(path="./hybird/slices", pattern="bin50.RDS")

# Create a list of count matrices
rds_read_mat <- lapply(paste0("./B73/bin50/", files_mat), readRDS)
rds_read_pat <- lapply(paste0("./A554/slices/", files_pat), readRDS)
rds_read_hyb <- lapply(paste0("./hybird/slices/", files_hyb), readRDS)

# assign names 
names(rds_read_mat)<- paste0('B73_bin50_', sub("^[^_]*_([^.]*)\\..*$", "\\1", files_mat))
names(rds_read_pat)<- paste0('A554_', sub("^[^_]*_([^.]*)\\..*$", "\\1", files_pat))
names(rds_read_hyb)<- paste0('B73_A554_', sub("^[^_]*_([^.]*)\\..*$", "\\1", files_hyb))

# merge all rds files
rds_read <- c(rds_read_mat, rds_read_pat, rds_read_hyb)

# merge multiple seurat objects to one object with multiple layers
all_sp_merge <- merge(rds_read[[1]], y = rds_read[2:length(rds_read)], 
                      add.cell.ids = names(rds_read), project="Heterosis_spatial_RNA_bin50")

head(all_sp_merge@meta.data)
Assays(all_sp_merge)
Layers(all_sp_merge)

# re-join layers to combine counts matrix
all_sp_merge[["Spatial"]] <- JoinLayers(all_sp_merge[["Spatial"]])

# visulization of fitering threshold

metadata <- all_sp_merge@meta.data

pdf(file = './all_sp_bin50_nfeature_ncount_plot.pdf', width = 6, height = 6)
ggplot(metadata) +
  geom_point(aes(x=nCount_Spatial, y=nFeature_Spatial), shape=21, alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  geom_vline(xintercept = 40,color="red",linetype="dotted")+
  geom_hline(yintercept=30,color="red", linetype="dotted")

dev.off()

# OR use the same parameters across all samples
sp_merge_filtered <- subset(all_sp_merge, nFeature_Spatial > 30 & nCount_Spatial > 40)

# Filter genes that are expressed in at least 3 cells (min.cells) and have at least 10 counts (min.count)
sp_merge_filtered <- subset(sp_merge_filtered, features = rownames(sp_merge_filtered)[
  rowSums(sp_merge_filtered@assays$Spatial$counts > 0) >= 3 & rowSums(sp_merge_filtered@assays$Spatial$counts) >= 10])

saveRDS(sp_merge_filtered, '~/spatial_RNA/intermediate_data/all_samples/sp_all_merge_filtered_bin50.RDS')

# split layers for following integration, this step is very neccessary!!
sp_merge_filtered$sample <- sub('_[^_]*$', '', colnames(sp_merge_filtered))
sp_merge_filtered[["Spatial"]] <- split(sp_merge_filtered[["Spatial"]], f = sp_merge_filtered$sample)

sp_merge_filtered <- SCTransform(sp_merge_filtered, assay = "Spatial", verbose = FALSE)
# by default, PCA is run only 3000 variable features
sp_merge_filtered <- RunPCA(sp_merge_filtered, verbose = FALSE, assay="SCT")
ElbowPlot(sp_merge_filtered, ndims = 50)

### perform integration across samples ###
sp_merge_filtered_int <- IntegrateLayers(object = sp_merge_filtered,
                                         method = CCAIntegration,
                                         orig.reduction = "pca",
                                         new.reduction = "cca.integration",
                                         normalization.method = "SCT",
                                         verbose = FALSE)

# clustering cells
sp_merge_filtered_int <- FindNeighbors(sp_merge_filtered_int, dims = 1:30, reduction = 'cca.integration')
sp_merge_filtered_int <- FindClusters(sp_merge_filtered_int, resolution = seq(1, 1.5, 0.2))

sp_merge_filtered_int <- RunUMAP(sp_merge_filtered_int, dims = 1:30, reduction = 'cca.integration')

plot1 <- DimPlot(sp_merge_filtered_int, reduction = "umap", label = TRUE) 
ggsave('./spatial_plots/spatial_all_integrated_bin50_UMAP_plot_res_1_2.pdf', plot1, width = 7, height = 6)

# find coluster 19 is noise and remove those cells for the following
Idents(sp_merge_filtered_int) <- sp_merge_filtered_int$SCT_snn_res.1.2
sp_merge_filtered_int2 <- subset(sp_merge_filtered_int, SCT_snn_res.1.2 != 19)
sp_merge_filtered_int2 <- subset(sp_merge_filtered_int2, !sample %in% c( "hybrid_slice_7_bin50" , "hybrid_slice_12_bin50"))
# filtering cells again
sp_merge_filtered_int2 <- subset(sp_merge_filtered_int2, nFeature_Spatial >200)

saveRDS(sp_merge_filtered_int2, '~/spatial_RNA/intermediate_data/all_sp_merge_filtered_integrated_remove_noise_bin50.RDS')


tissue.col <- c("#C9A1FA", "#EDC66A", '#B46DA9',"#F9B3AD", '#E87B1E', "#018B38",
                "#D9DEE7", "#8498AB", 
                "blue", "green", "yellow", "royalblue", "red", "#46f0f0", "#bcf60c", "darkred", 
                "gold3", "pink3", "brown", "darkgreen")



plot2 <- SpatialDimPlot(sp_merge_filtered_int2, label = TRUE, label.size = 1, 
                        pt.size.factor = 1, ncol = 6) + 
                  scale_fill_manual(values = tissue.col)
ggsave('./spatial_plots/spatial_all_integrated_spatial_plot_res_1_2.pdf', 
       plot2, width = 35, height = 30)


# annoate spatial clusters with markers from LCM and previous literature

all.DE.genes <- readxl::read_xlsx('~/LCM_RNA/results/DEGenes/glmQLfit/LCM_glmQLfit_all_DEGs.xlsx')
plot.genes <- all.DE.genes[all.DE.genes$geneID %in% rownames(sp_merge_filtered_int2), ]
# choose top 50 LCM differential genes for each tissue
top_50_LCM_genes <- plot.genes %>%
  group_by(Tissue) %>%
  slice_max(order_by = logCPM, n = 50, with_ties = FALSE) %>%
  ungroup()

# average expression across clusters
avg_exp <- AverageExpression(sp_merge_filtered_int2, 
                             features = top_50_LCM_genes$geneID, 
                             return.seurat = FALSE)
# Get expression matrix
avg_matrix <- avg_exp$SCT  

p3 = pheatmap(avg_matrix, 
              scale = "row",              # scale gene-wise
              cluster_rows = F, cluster_cols = F,
              color = colorRampPalette(c("blue", "white", "red"))(100))

ggsave(paste0('./cluster_annotation/sp_all_integrated_anno_by_LCM_markers_res_1_2.pdf'), 
       p3, width = 15, height = 30)

markers_valid <- data.frame(geneID = c('Zm00001d032822',  'Zm00001d046186',
                                       'Zm00001d012081', 'Zm00001d052380', 'Zm00001d022180',
                                       'Zm00001d021706', 'Zm00001d035689', 'Zm00001d049540',
                                       'Zm00001d021192', 'Zm00001d005472', 'Zm00001d050168', 
                                       'Zm00001d040390'),
                            cellType = c('Epidermis-hs',  'Cortex Mature',
                                         'Cortex-hs', 'Endodermis_ZmSCR',  'Pericycle',
                                         'Pericycle', 'Xylem_ZmXCP1', 'Phloem_NAC77',
                                         'Stele-hs', 'Pericycle-hs', 'Endodermis-hs',
                                         'Initials'))

# pericycle <- c('Zm00001d002543,Zm00001d002601,Zm00001d003493,Zm00001d003725,Zm00001d003730,Zm00001d005791,Zm00001d006213,Zm00001d006547,Zm00001d006548,Zm00001d008669,Zm00001d010575,Zm00001d011236,Zm00001d012837,Zm00001d013066,Zm00001d013067,Zm00001d013392,Zm00001d016723,Zm00001d018981,Zm00001d019045,Zm00001d020584,Zm00001d020585,Zm00001d021300,Zm00001d021433,Zm00001d021477,Zm00001d021706,Zm00001d021707,Zm00001d023901,Zm00001d025381,Zm00001d025406,Zm00001d025913,Zm00001d026015,Zm00001d029394,Zm00001d032070,Zm00001d034479,Zm00001d036036,Zm00001d036250,Zm00001d038381,Zm00001d039733,Zm00001d039790,Zm00001d039821,Zm00001d041672,Zm00001d042730,Zm00001d042930,Zm00001d042954,Zm00001d044246,Zm00001d044247,Zm00001d045268,Zm00001d045475,Zm00001d047787,Zm00001d047788,Zm00001d048491,Zm00001d050100,Zm00001d050697,Zm00001d051478,Zm00001d051591')
# pericycle <- strsplit(pericycle, split = ',')

setwd('~/spatial_RNA/results/integrated_3_samples/')
g1 <- DotPlot(object = sp_merge_filtered_int2, 
              features = markers_valid$geneID,
              cols = c("blue", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = markers_valid$cellType)

ggsave('./cluster_annotation/sp_all_annotated_literature_markers.pdf', 
       g1, width = 10, height = 8)








