library(Seurat)
library(scDblFinder)
library(tidyverse)
library(BPCells)
library(devtools)
library(ggplot2)


#devtools::install_github("immunogenomics/presto")
#install.packages('BPCells', repos = c('https://bnprks.r-universe.dev', 'https://cloud.r-project.org'))
#BiocManager::install('glmGamPoi')



################################################################################
#####   Single cell transcriptome data analysis                            #####
################################################################################

setwd('~/single_cell_RNA_seq/')

#List files 
files <- list.files(path="./data/", recursive=T, pattern="filtered_feature_bc_matrix.h5")

# Create a list of count matrices
h5_read <- lapply(paste0("./data/",
                         files), Read10X_h5)
# Assign names manually
names(h5_read)<-c("B73xA554","B73xOh40B","B73xOh43", "B73xNC328", "B73xB37",
                  "B73xA680", "B73xNC298", "B73xNC352", "B73xNC320",
                  "B73_Rep1", "B73_Rep2", "B73_Rep3", "A554", "Oh40B", "Oh43",
                  "NC328", "B37", "A680", "NC298", "NC352", "NC320")
# Create seurat objects
all_scRNA <- mapply(CreateSeuratObject, counts=h5_read,  
              project=names(h5_read),
              MoreArgs = list(min.cells = 3, min.features = 100))
# merge multiple seurat objects to one object with multiple layers
all_scRNA_merge <- merge(all_scRNA[[1]], y = all_scRNA[2:length(all_scRNA)], 
                 add.cell.ids = names(all_scRNA), project="Heterosis")

head(all_scRNA_merge@meta.data)
Assays(all_scRNA_merge)
Layers(all_scRNA_merge)

# Visualize the number of cell counts per sample
all_scRNA_merge@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = ..count..),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")

# Add percent of mitochondrial
all_scRNA_merge[["percent.mt"]] <- PercentageFeatureSet(all_scRNA_merge, pattern = "^MT-")

metadata <- all_scRNA_merge@meta.data

# density of counts
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 500,color="black",linetype="dotted")+
  geom_vline(xintercept = 200000,color="black",linetype="dotted")

# density of unique genes
ggplot(metadata, aes(x=nFeature_RNA, fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 400,color="black",linetype="dotted")+
  geom_vline(xintercept = 12000,color="black",linetype="dotted")

# mitochondrial percent
ggplot(metadata, aes(x=percent.mt, fill=orig.ident)) +
  geom_density(alpha = 0.2) + 
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 5,color="black",linetype="dotted")

pdf(file = './results/pre_filter_plot3.pdf', width = 20, height = 20)
ggplot(metadata) +
  geom_point(aes(x=nCount_RNA, y=nFeature_RNA, fill=percent.mt >5), shape=21, alpha=0.4) + 
  theme_classic() +
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(.~orig.ident, ncol = 5) +
  geom_vline(xintercept = 1000,color="red",linetype="dotted")+
  geom_hline(yintercept=500,color="red", linetype="dotted")
dev.off()

# use the same parameters across all samples
all_scRNA_merge_filtered <- subset(all_scRNA_merge, nFeature_RNA >500 & nCount_RNA >1000 & percent.mt <5 )
saveRDS(all_scRNA_merge_filtered, './intermediate_data/all_scRNA_merge_filtered.RDS')

# set 
options(future.globals.maxSize = 7e+10)

# Normalize and scale RNA assay (independently of SCT)
all_scRNA_merge_norm <- NormalizeData(all_scRNA_merge_filtered)
all_scRNA_merge_norm <- ScaleData(all_scRNA_merge_norm)

scRNA_sce <- as.SingleCellExperiment(all_scRNA_merge_norm)
sce <- scDblFinder(scRNA_sce, samples = "orig.ident") 
table(sce$scDblFinder.class)

# Explore results and add to seurat object
meta_scdblfinder <- sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) 

head(meta_scdblfinder)
rownames(meta_scdblfinder) <- sce@colData@rownames

all_scRNA_merge_filtered <- AddMetaData(object = all_scRNA_merge_filtered, 
                                    metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))
# Remove doublets
scRNA_merge_filtered_clean <- subset(all_scRNA_merge_filtered, scDblFinder.class == 'singlet')

# Normalization
options(future.globals.maxSize = 7e+10)

scRNA_merge_filtered_clean <- SCTransform(scRNA_merge_filtered_clean, vars.to.regress = "percent.mt", verbose = FALSE)
# by default, PCA is run only 3000 variable features
scRNA_merge_filtered_clean <- RunPCA(scRNA_merge_filtered_clean, verbose = FALSE, assay="SCT")
ElbowPlot(scRNA_merge_filtered_clean, ndims = 50)

### perform integration across samples ###
scRNA_merge_filtered_clean <- IntegrateLayers(object = scRNA_merge_filtered_clean, 
                                              method = CCAIntegration, 
                                              orig.reduction = "pca", 
                                              new.reduction = "cca.integration",
                                              normalization.method = "SCT",
                                              verbose = FALSE)
# clustering cells 
scRNA_merge_filtered_clean <- FindNeighbors(scRNA_merge_filtered_clean, reduction = "cca.integration", dims = 1:30)
scRNA_merge_filtered_clean <- FindClusters(scRNA_merge_filtered_clean, resolution = 0.52)
# UMAP plot
scRNA_merge_filtered_clean <- RunUMAP(scRNA_merge_filtered_clean, dims = 1:30, reduction = "cca.integration")

pdf(file = './results/umap/clusters_after_integration.pdf', width = 18, height = 10)
DimPlot(scRNA_merge_filtered_clean, reduction = "umap", label = T, 
        group.by = c('orig.ident', 'seurat_clusters'),
        alpha=0.4, ncol=2)
dev.off()

### find marker genes for each cluster ###

scRNA_merge_filtered_clean <- PrepSCTFindMarkers(scRNA_merge_filtered_clean)
sc.all.markers <- c()

for (i in 0:20) {
  print(i)
  all.markers <- FindMarkers(scRNA_merge_filtered_clean, 
                             ident.1 = i, 
                             assay = 'SCT',
                             only.pos = TRUE,
                             min.pct = 0.15, 
                             logfc.threshold = 0.25,
                             min.diff.pct = 0.1,
                             verbose = FALSE)
  head(all.markers)
  sc.all.markers <- rbind.data.frame(sc.all.markers, cbind(all.markers, cluster=i, markerID=rownames(all.markers)))
}


sc.all.markers <- sc.all.markers[sc.all.markers$p_val_adj <0.05, ]

write_xlsx(sc.all.markers, './results/sc_RNA_all_markers_res_0_52.xlsx')

# 'Single-cell transcriptomes reveal spatiotemporal heat stress response in maize roots' Through in situ hybridization, we found
# Zm00001d032822 is speciﬁc to epidermal cells, Zm00001d012081
# (PLT29)to cortexcells, Zm00001d050168 to the endodermis,
# Zm00001d021192 (umc2686b)to the stele, Zm00001d032672 to the
# xylem, Zm00001d037032 to the phloem, Zm00001d005472 to the
# pericycle, Zm00001d004089 (PRP18) to the columella, and
# Zm00001d040390 to the initials


markers_from_literature = readxl::read_xlsx('./intermediate_data/Maize root markers.xlsx', sheet = 3)
table(markers_from_literature$clusterName)

DotPlot(object = scRNA_merge_filtered_clean, 
        features = markers_from_literature[markers_from_literature$clusterName == 'Root cortex', ]$gene)

# these genes from DE genes of LCM RNA
DotPlot(object = scRNA_merge_filtered_clean, 
        features = c('Zm00001d032822', 'Zm00001d048954', 
                                                        'Zm00001d006571', 'Zm00001d036066', 
                                                        'Zm00001d048759', 'Zm00001d036069'))






