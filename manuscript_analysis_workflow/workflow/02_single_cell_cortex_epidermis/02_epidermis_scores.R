# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

##Epidermis后面的非加基因表达是否与根的垂直发育有关

library(readr)
genesets <- read_csv(file.path(project_root, "data", "EZ_M_Genes_from_clusters_with_clusterinfo.csv"))
Epidermis <- AddModuleScore(Epidermis,
                         features = list(dplyr::filter(genesets,group == "M_Genes")$gene),
                         name = "M_Scores")

Epidermis <- AddModuleScore(Epidermis,
                         features = list(dplyr::filter(genesets,group == "EZ_Genes")$gene),
                         name = "E_Scores")




FeaturePlot(Epidermis, features = "M_Scores1") +
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


##检测Epidermis的拟时序和非加性基因的关系
setwd(project_root)
library(readr)
nonaddgenes <- read_csv("sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv")#读取单细胞非加基因数据

library(dplyr)
library(stringr)

#所有材料的Cortex特异性基因的叠加
nonadd_epi <- nonaddgenes %>%
  filter(pattern == "Non-additive") %>%        # 先筛出所有非加性
  group_by(Trait) %>%                           # 按基因分组
  filter(
    all(str_detect(cluster_anno, "Epidermis"))    # 该基因的所有 Non-additive 都在Epidermis
  ) %>%
  ungroup()


nonadd_epi_consistent <- nonadd_epi %>%
  group_by(Trait) %>%
  filter(
    # 该基因所有非加性表达的 mph.pct 都 > 0
    all(MPH.pct > 0) |
      # 或者都 < 0
      all(MPH.pct < 0)
  ) %>%
  ungroup()

# 可选：添加一个辅助列，标记方向类型
nonadd_epi_consistent <- nonadd_epi_consistent %>%
  group_by(Trait) %>%
  mutate(
    direction = case_when(
      all(MPH.pct > 0) ~ "Positive",
      all(MPH.pct < 0) ~ "Negative",
      TRUE ~ "Mixed"  # 这行实际不会出现，因为已被 filter 排除
    )
  ) %>%
  ungroup()

# 查看结果
head(nonadd_epi_consistent)


##上调的基因
up_nonadd_epi<-unique(filter(nonadd_epi_consistent,direction=="Positive")$Trait)
down_nonadd_epi<-unique(filter(nonadd_epi_consistent,direction=="Negative")$Trait)

Epidermis <- AddModuleScore(Epidermis,
                         features = list(up_nonadd_epi,down_nonadd_epi),
                         name = c("NonAddgene_Scores_UP","NonAddgene_Scores_DOWN"))

FeaturePlot(Epidermis, features = "NonAddgene_Scores_DOWN2") +
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


cor(Epidermis$NonAddgene_Scores_UP1, Epidermis$M_Scores1, method = "pearson")
cor.test(Epidermis$NonAddgene_Scores_UP1, Epidermis$E_Scores1)

cor(Epidermis$NonAddgene_Scores_DOWN2, Epidermis$M_Scores1, method = "pearson")
cor.test(Epidermis$NonAddgene_Scores_DOWN2, Epidermis$E_Scores1)



library(ggplot2)

plot(Epidermis$NonAddgene_Scores_UP1, Epidermis$E_Scores1,
     xlab = "NonAddgene_Scores_UP",
     ylab = "E_Scores", pch = 16, col = "black",cex=0.01)
abline(lm(Epidermis$E_Scores1 ~ Epidermis$NonAddgene_Scores_UP1),
       col = "red", lwd = 2)


plot(Epidermis$NonAddgene_Scores_DOWN2, Epidermis$E_Scores1,
     xlab = "NonAddgene_Scores_DOWN",
     ylab = "E_Scores", pch = 16, col = "black",cex=0.01)
abline(lm(Epidermis$NonAddgene_Scores_DOWN2 ~ Epidermis$E_Scores1), col = "red", lwd = 2)
