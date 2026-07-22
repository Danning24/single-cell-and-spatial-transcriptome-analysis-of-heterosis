# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

setwd(project_root)
sce <- readRDS(file.path(project_root, "data", "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"))
library(Seurat)
library(dplyr)
DimPlot(scRNA_all_annotated_filtered_500gene_1000umi_CCA)
Cortex<-subset(scRNA_all_annotated_filtered_500gene_1000umi_CCA,ident = "Cortex")
Cortex <- Cortex %>%
  SCTransform() %>%
  RunPCA(npcs = 30)

Cortex <- Cortex %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.04) %>%
  RunUMAP(dims = 1:15)


Idents(Cortex)<-Cortex$seurat_clusters

library(ggsci)
DimPlot(Cortex) +
  scale_color_npg()+
  NoLegend()+
  NoAxes()

Idents(Cortex)<-Cortex$orig.ident
newname<-c("A554","A680","B37","B73","B73","B73","B73xA554","B73xA680","B73xB37","B73xNC298",
           "B73xNC320","B73xNC328","B73xNC352","B73xOh40B","B73xOh43","NC298","NC320","NC328",
           "NC352","Oh40B","Oh43")
names(newname)<-c("A554","A680","B37","B73_Rep1","B73_Rep2","B73_Rep3","B73xA554","B73xA680","B73xB37","B73xNC298",
                  "B73xNC320","B73xNC328","B73xNC352","B73xOh40B","B73xOh43","NC298","NC320","NC328",
                  "NC352","Oh40B","Oh43")
Cortex<-RenameIdents(Cortex,newname)
Cortex$orig.ident<-Idents(Cortex)

table(Cortex$orig.ident)
newname<-c("NSS","NSS","NSS","SS","SS","SS","TL","TL","TL","SS",
           "NSS","NSS","NSS","SS","SS","SS","TL","TL","TL")
names(newname)<-unique(Cortex$orig.ident)
Cortex <- RenameIdents(Cortex,newname)
Cortex$group<-Idents(Cortex)


DimPlot(Cortex,split.by = "group") +
  scale_color_npg()+
  NoLegend()+
  NoAxes()


library(ggpointdensity)
library(viridis)
library(ggfx)

# 提取umap坐标数据
Cortex$group_FP <- ifelse(grepl("x", Cortex$orig.ident),
                       "F1", ifelse(Cortex$orig.ident == "B73",
                                    "B73",
                                    "P1"))

table(Cortex$group_FP)

dat <- data.frame(Embeddings(subset(Cortex,subset = group == "TL"), reduction = "umap"))
ggplot(data = dat, mapping = aes(x = umap_1, y = umap_2)) +
  geom_pointdensity(size = 0.01) +
  scale_color_viridis(option = "magma")+
  theme_void()+
  theme(
    plot.background = element_rect(fill = "black", color = NA), # 图表背景
    panel.background = element_rect(fill = "black", color = NA), # 绘图区域背景
    legend.position = "none" # 移除图例
  )


plot_umap_density_one <- function(
    obj,
    ident_value,
    label_color = "white",
    label_x = -Inf,
    label_y = Inf
) {
  # 子集：同时限制 group 与 orig.ident
  sub <- subset(obj, subset = orig.ident == ident_value)

  dat <- as.data.frame(Embeddings(sub, reduction = "umap"))
  ggplot(dat, aes(x = umap_1, y = umap_2)) +
    geom_pointdensity(size = 0.01) +
    scale_color_viridis(option = "magma") +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.position  = "none"
    ) +
    # 标名字：默认左上角
    annotate(
      "text",
      x = label_x, y = label_y,
      label = ident_value,
      hjust = -0.05, vjust = 1.1,
      color = label_color, size = 6, fontface = "bold"
    )
}

idents <- unique(Cortex$orig.ident)

for (id in idents) {
  p <- plot_umap_density_one(Cortex, ident_value = id)
  print(p)
}


dat <- data.frame(Embeddings(subset(Cortex,subset = group == "TL" &
                                      group_FP == "F1"),
                             reduction = "umap"))
ggplot(data = dat, mapping = aes(x = umap_1, y = umap_2)) +
  geom_pointdensity(size = 0.01) +
  scale_color_viridis(option = "magma")+
  theme_void()+
  theme(
    plot.background = element_rect(fill = "black", color = NA), # 图表背景
    panel.background = element_rect(fill = "black", color = NA), # 绘图区域背景
    legend.position = "none" # 移除图例
  )+
  # 标名字：默认左上角
  annotate(
    "text",
    x=-Inf,y=Inf,
    label = "F1_TL",
    hjust = -0.05, vjust = 1.1,
    color = "white", size = 6, fontface = "bold"
  )

table(Cortex$group,Cortex$orig.ident)




library(readr)
genesets <- read_csv("EZ_M_Genes_from_clusters_with_clusterinfo.csv")

Cortex <- AddModuleScore(Cortex,
                       features = list(dplyr::filter(genesets,group == "M_Genes")$gene),
                       name = "M_Scores")

Cortex <- AddModuleScore(Cortex,
                       features = list(dplyr::filter(genesets,group == "EZ_Genes")$gene),
                       name = "E_Scores")




FeaturePlot(Cortex, features = "E_Scores1") +
  scale_color_gradientn(
    colors = c("#5E3C99", "white", "#B22222")
  ) +
  labs(title = " ") +
  theme(
    plot.title = element_text(
      hjust = 0.5,      # 居中
      size = 16,
      face = "bold"
    )
  )+
  NoAxes()
saveRDS(Cortex,"Cortex.rds")


plot_umap_Mscores_one <- function(
    obj,
    ident_value,
    label_color = "black",
    label_x = -Inf,
    label_y = Inf
) {
  # 子集：同时限制 group 与 orig.ident
  sub <- subset(obj, subset = orig.ident == ident_value)
  FeaturePlot(sub, features = "M_Scores1") +
    scale_color_gradientn(
      colors = c("#5E3C99", "white", "#B22222")
    ) +
    labs(title = " ") +
    theme(
      plot.title = element_text(
        hjust = 0.5,      # 居中
        size = 16,
        face = "bold"
      )
    )+
    NoAxes()+
    # 标名字：默认左上角
    annotate(
      "text",
      x = label_x, y = label_y,
      label = ident_value,
      hjust = -0.05, vjust = 1.1,
      color = label_color, size = 6, fontface = "bold"
    )


}

idents <- unique(Cortex$orig.ident)

for (id in idents) {
  p <- plot_umap_Mscores_one(Cortex, ident_value = id)
  print(p)
}


FeaturePlot(subset(Cortex,subset = group == "NSS" &
                     group_FP == "F1"), features = "M_Scores1") +
  scale_color_gradientn(
    colors = c("#5E3C99", "white", "#B22222")
  ) +
  labs(title = " ") +
  theme(
    plot.title = element_text(
      hjust = 0.5,      # 居中
      size = 16,
      face = "bold"
    )
  )+
  NoAxes()+
  # 标名字：默认左上角
  annotate(
    "text",
    x=-Inf,y=Inf,
    label = "F1_NSS",
    hjust = -0.05, vjust = 1.1,
    color = "black", size = 6, fontface = "bold"
  )

cor(Cortex$E_Scores1, Cortex$M_Scores1, method = "pearson")
cor.test(Cortex$E_Scores1, Cortex$M_Scores1)
plot(Cortex$E_Scores1, Cortex$M_Scores1,
     xlab = "M_Scores",
     ylab = "E_Scores", pch = 16, col = "black",cex=0.01)
abline(lm(Cortex$E_Scores1 ~ Cortex$M_Scores1), col = "red", lwd = 2)





##检测Cortex的拟时序和非加性基因的关系
setwd(project_root)
library(readr)
nonaddgenes <- read_csv("sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv")#读取单细胞非加基因数据

library(dplyr)
library(stringr)

#所有材料的Cortex特异性基因的叠加
nonadd_cortex <- nonaddgenes %>%
  filter(pattern == "Non-additive") %>%        # 先筛出所有非加性
  group_by(Trait) %>%                           # 按基因分组
  filter(
    all(str_detect(cluster_anno, "Cortex"))    # 该基因的所有 Non-additive 都在 Cortex
  ) %>%
  ungroup()


length(unique(nonadd_cortex$Trait))#查看筛选后的基因数目

Cortex <- AddModuleScore(Cortex,
                         features = list(unique(nonadd_cortex$Trait)),
                         name = "NonAddgene_Scores")

FeaturePlot(Cortex, features = "NonAddgene_Scores1") +
  scale_color_gradientn(
    colors = c("#000000", "#DDDDDD", "#E69F00")
  ) +
  labs(title = " ") +
  theme(
    plot.title = element_text(
      hjust = 0.5,      # 居中
      size = 16,
      face = "bold"
    )
  )+
  NoAxes()

Cortex$NonAddgene_Scores1
Cortex$E_Scores1

cor(Cortex$NonAddgene_Scores1, Cortex$M_Scores1, method = "pearson")
cor.test(Cortex$NonAddgene_Scores1, Cortex$M_Scores1)



library(ggplot2)

plot(Cortex$NonAddgene_Scores1, Cortex$E_Scores1,
     xlab = "NonAddgene_Scores",
     ylab = "E_Scores", pch = 16, col = "black",cex=0.01)
abline(lm(Cortex$E_Scores1 ~ Cortex$NonAddgene_Scores1), col = "red", lwd = 2)


plot(Cortex$NonAddgene_Scores1, Cortex$M_Scores1,
     xlab = "NonAddgene_Scores",
     ylab = "M_Scores", pch = 16, col = "black",cex=0.01)
abline(lm(Cortex$M_Scores1 ~ Cortex$NonAddgene_Scores1), col = "red", lwd = 2)




##是否是因为混合分析掩饰了杂交品种中的强相关性？
##首先用B73xoh4ob来进行初步验证
table(Cortex@meta.data$orig.ident)


nonadd_cortex <- nonaddgenes %>%
  filter(
    pattern == "Non-additive",
    Hybrid == "B73xOh40B"        # 新增条件
  ) %>%
  group_by(Trait) %>%
  filter(
    all(str_detect(cluster_anno, "Cortex"))
  ) %>%
  ungroup()
length(unique(nonadd_cortex$Trait))

subcortex<-subset(Cortex,subset = orig.ident == "B73xOh40B")
subcortex <- AddModuleScore(subcortex,
                         features = list(unique(nonadd_cortex$Trait)),
                         name = "NonAddgene_Scores")
cor.test(subcortex$NonAddgene_Scores1, subcortex$E_Scores1)
cor.test(subcortex$NonAddgene_Scores1, subcortex$M_Scores1)



##函数模块化
run_nonadd_cortex_cor <- function(nonaddgenes, Cortex, hyb,
                                  pattern_value = "Non-additive",
                                  cortex_keyword = "Cortex",
                                  method = "pearson") {
  # 1) 筛选 nonadd_cortex
  nonadd_cortex <- nonaddgenes %>%
    dplyr::filter(
      pattern == pattern_value,
      Hybrid == hyb
    ) %>%
    dplyr::group_by(Trait) %>%
    dplyr::filter(
      all(stringr::str_detect(cluster_anno, cortex_keyword))
    ) %>%
    dplyr::ungroup()

  traits <- unique(nonadd_cortex$Trait)


  # 2) subset cortex
  subcortex <- subset(Cortex, subset = orig.ident == hyb)

  # 3) AddModuleScore
  subcortex <- Seurat::AddModuleScore(
    subcortex,
    features = list(traits),
    name = "NonAddgene_Scores"
  )

  x <- subcortex$NonAddgene_Scores1
  e <- subcortex$E_Scores1
  m <- subcortex$M_Scores1

  # 4) 两次 cor.test（use="complete.obs" 等价于自动忽略 NA 配对）
  ct_e <- stats::cor.test(x, e, method = method, use = "complete.obs")
  ct_m <- stats::cor.test(x, m, method = method, use = "complete.obs")

  # 5) 输出：r 和 p
  out <- data.frame(
    hyb = hyb,
    n_traits = length(traits),
    test = c("NonAdd_vs_E", "NonAdd_vs_M"),
    r = c(unname(ct_e$estimate), unname(ct_m$estimate)),
    p = c(ct_e$p.value, ct_m$p.value)
  )

  return(out)
}

hyb_list <- c(
  "B73xA554",
  "B73xA680",
  "B73xB37",
  "B73xNC298",
  "B73xNC320",
  "B73xNC328",
  "B73xNC352",
  "B73xOh40B",
  "B73xOh43"
)

# 批量运行并合并结果
all_results <- do.call(
  rbind,
  lapply(hyb_list, function(h) {
    run_nonadd_cortex_cor(
      nonaddgenes = nonaddgenes,
      Cortex = Cortex,
      hyb = h
    )
  })
)

all_results

#找出在所有材料中都保守的cortex特异的基因
traits_all_hybrids <- nonadd_cortex %>%
  group_by(Trait) %>%
  summarise(n_hybrid = n_distinct(Hybrid)) %>%
  filter(n_hybrid == length(unique(nonadd_cortex$Hybrid))) %>%
  pull(Trait)


traits_all_hybrids

FeaturePlot(Cortex, features = "Zm00001d004340",pt.size = 0.01) +
  labs(title = " ") +
  theme(
    plot.title = element_text(
      hjust = 0.5,      # 居中
      size = 16,
      face = "bold"
    )
  )+
  NoAxes()


gene <- "Zm00001d004340"
gene_expr <- FetchData(Cortex, vars = c(gene, "E_Scores1"))

# 做相关
cor_result <- cor.test(gene_expr$Zm00001d004340, gene_expr$E_Scores1, method = "spearman")

cor_result

ggplot(gene_expr, aes(x = Zm00001d004340, y = E_Scores1)) +
  geom_jitter(width = 0.05, height = 0, alpha = 0.35, size = 1) +
  theme_classic() +
  labs(
    x = "log1p(Zm00001d004340 expression)",
    y = "E_Scores1",
    title = paste0(
      "Spearman rho = ", round(cor_result$estimate, 3),
      ", p = ", signif(cor_result$p.value, 3)
    )
  )


#B73xOh40B上调和下调的非加基因分析
nonadd_cortex <- nonaddgenes %>%
  dplyr::filter(
    pattern == "Non-additive",
    Hybrid == "B73xOh40B"
  ) %>%
  dplyr::group_by(Trait) %>%
  dplyr::filter(
    all(stringr::str_detect(cluster_anno, "Cortex"))
  ) %>%
  dplyr::ungroup()

traits <- unique(nonadd_cortex$Trait)

nonadd_cortex_up<-filter(nonadd_cortex,MPH.pct>0)$Trait
nonadd_cortex_down<-filter(nonadd_cortex,MPH.pct<0)$Trait
nonadd_cortex_upf<-setdiff(nonadd_cortex_up, nonadd_cortex_down)
nonadd_cortex_downf<-setdiff(nonadd_cortex_down, nonadd_cortex_up)

# 2) subset cortex
subcortex <- subset(Cortex, subset = orig.ident == "B73xOh40B")

# 3) AddModuleScore
subcortex <- Seurat::AddModuleScore(
  subcortex,
  features = list(nonadd_cortex_upf,nonadd_cortex_downf),
  name = c("up_score","down_score")
)

u <- subcortex$up_score1
d <- subcortex$down_score2
e <- subcortex$E_Scores1

# 4) 两次 cor.test（use="complete.obs" 等价于自动忽略 NA 配对）
ct_u <- stats::cor.test(u, e, method = "pearson", use = "complete.obs")
ct_d <- stats::cor.test(d, e, method = "pearson", use = "complete.obs")
plot(u, e,
     xlab = "NonAddgene_UP_Scores",
     ylab = "E_Scores", pch = 16, col = "black",cex=0.1)
abline(lm(e ~ u), col = "red", lwd = 2)

plot(d, e,
     xlab = "NonAddgene_DOWN_Scores",
     ylab = "E_Scores", pch = 16, col = "black",cex=0.1)
abline(lm(e ~ d), col = "red", lwd = 2)

"Zm00001d004340"%in%nonadd_cortex_upf
