# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

rm(list = ls());gc()

setwd(project_root)
sce <- readRDS(file.path(project_root, "data", "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"))
marker <- read_excel("marker.xlsx")

#BiocManager::install("clustree")
library(Seurat)
library(clustree)
library(readxl)
library(monocle3)
library(patchwork)
library(SeuratWrappers)

table(sce$Annotation)
Idents(sce)<-sce$Annotation
Idents(sce)<-sce$orig.ident
newname<-c("A554","A680","B37","B73","B73","B73","B73xA554","B73xA680","B73xB37","B73xNC298",
  "B73xNC320","B73xNC328","B73xNC352","B73xOh40B","B73xOh43","NC298","NC320","NC328",
  "NC352","Oh40B","Oh43")
names(newname)<-c("A554","A680","B37","B73_Rep1","B73_Rep2","B73_Rep3","B73xA554","B73xA680","B73xB37","B73xNC298",
                  "B73xNC320","B73xNC328","B73xNC352","B73xOh40B","B73xOh43","NC298","NC320","NC328",
                  "NC352","Oh40B","Oh43")
sce<-RenameIdents(sce,newname)

DimPlot(sce,label = T)
head(sce@assays)
table(sce$seurat_clusters)

#分离出相应的细胞
Cortex<-subset(sce,idents="Cortex")
Cortex <- Cortex %>%
  SCTransform() %>%
  RunPCA(npcs = 30)

Cortex <- Cortex %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.04) %>%
  RunUMAP(dims = 1:15)

Idents(Cortex)<-Cortex$SCT_snn_res.0.04
DimPlot(Cortex,label = T)

Markers<-FindAllMarkers(Cortex,logfc.threshold = 0.58,only.pos = T)
FeaturePlot(Cortex,features = c("Zm00001d022421"))
FeaturePlot(Cortex,features = c("Zm00001d001951"))
FeaturePlot(Cortex,features = c("Zm00001d040390"))
VlnPlot(Cortex,features = c("Zm00001d022421"))+NoLegend()
VlnPlot(Cortex,features = c("Zm00001d001951"))+NoLegend()
VlnPlot(Cortex,features = c("Zm00001d040390"))+NoLegend()



Idents(Cortex)<-Cortex$orig.ident
newname<-c("A554","A680","B37","B37","B73xA554","B73xA680","B73xB37","B73xNC298",
           "B73xNC320","B73xNC328","B73xNC352","B73xOh40B","B73xOh43","NC298","NC320","NC328",
           "NC352","Oh40B","Oh43")
names(newname)<-levels(Cortex$orig.ident)
Idents(Cortex)<-factor(Cortex$orig.ident,levels = c("A554","A680","B37","NC298","NC320","NC328",
                                                    "NC352","Oh40B","Oh43","B73xA554","B73xA680","B73xB37","B73xNC298",
                                                    "B73xNC320","B73xNC328","B73xNC352","B73xOh40B","B73xOh43"))

Cortex<-RenameIdents(Cortex,newname)
Cortex$orig.ident<-Idents(Cortex)

save_feature_plot <- function(gene, seurat_obj) {
  p_list <- FeaturePlot(seurat_obj,
                        features = gene,
                        combine = FALSE,
                        split.by = "orig.ident",
                        pt.size =  0.01)
  p_combined <- wrap_plots(p_list, ncol = 6)
  filename <- paste0(".", "/", gene, "_featureplot.png")
  ggsave(filename, plot = p_combined, width = 27, height = 10)
}
save_feature_plot("Zm00001d022421", Cortex)
save_feature_plot("Zm00001d001951", Cortex)
save_feature_plot("Zm00001d040390", Cortex)

marker <-marker[marker$Gene %in% rownames(Cortex),]
marker_list <- split(marker$Gene, marker$Cluster)
Cortex <- AddModuleScore(Cortex, features = marker_list, name = "CortexType")
VlnPlot(Cortex, features = c("CortexType1", "CortexType2", "CortexType3",
      "CortexType4", "CortexType5", "CortexType6","CortexType7"),
        group.by = "SCT_snn_res.0.1",pt.size=0)

marker <- marker %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%
  ungroup()



Epidermis<-subset(sce,idents="Epidermis")
Epidermis <- Epidermis %>%
  SCTransform() %>%
  RunPCA(npcs = 30)

Epidermis <- Epidermis %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0) %>%
  RunUMAP(dims = 1:15)

Idents(Epidermis)<-Epidermis$SCT_snn_res.0.04
DimPlot(Epidermis,label = T)

##Cortex的monocle3分析

cds <- as.cell_data_set(Cortex)
name<-c("Cortex_1","Cortex_2","Cortex_3")
names(name)<-levels(Cortex)
Cortex<-RenameIdents(Cortex,name)
Cortex$Annotation<-Idents(Cortex)
cds@clusters$UMAP$clusters <- Idents(Cortex)
reducedDims(cds)$UMAP <- Cortex@reductions$umap@cell.embeddings

cds <- cluster_cells(cds, resolution = 1e-3)
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster")

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = F,
           label_branch_points = F)


DimPlot(Epidermis,label = T)
cds_Epidermis <- as.cell_data_set(Epidermis)
name<-c("Epidermis")
names(name)<-levels(Epidermis)
Epidermis<-RenameIdents(Epidermis,name)
cds_Epidermis@clusters$UMAP$clusters <- Idents(Epidermis)
reducedDims(cds_Epidermis)$UMAP <- Epidermis@reductions$umap@cell.embeddings

cds_Epidermis <- cluster_cells(cds_Epidermis, resolution = 1e-3)
cds_Epidermis <- learn_graph(cds_Epidermis)
plot_cells(cds_Epidermis, color_cells_by = "cluster")
FeaturePlot(Epidermis,features = c("Zm00001d052952"))
FeaturePlot(Epidermis,features = c("Zm00001d002191"))
FeaturePlot(Epidermis,features = c("Zm00001d041611"))

cds_Epidermis <- order_cells(cds_Epidermis)
plot_cells(cds_Epidermis,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = F,
           label_branch_points = F)
