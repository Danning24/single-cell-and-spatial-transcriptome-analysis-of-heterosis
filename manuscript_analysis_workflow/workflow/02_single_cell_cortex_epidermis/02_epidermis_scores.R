project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

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
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  )+
  NoAxes()
saveRDS(Cortex,"Cortex.rds")

setwd(project_root)
library(readr)
nonaddgenes <- read_csv("sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv")

library(dplyr)
library(stringr)

nonadd_epi <- nonaddgenes %>%
  filter(pattern == "Non-additive") %>%
  group_by(Trait) %>%
  filter(
    all(str_detect(cluster_anno, "Epidermis"))
  ) %>%
  ungroup()

nonadd_epi_consistent <- nonadd_epi %>%
  group_by(Trait) %>%
  filter(
    all(MPH.pct > 0) |
      all(MPH.pct < 0)
  ) %>%
  ungroup()

nonadd_epi_consistent <- nonadd_epi_consistent %>%
  group_by(Trait) %>%
  mutate(
    direction = case_when(
      all(MPH.pct > 0) ~ "Positive",
      all(MPH.pct < 0) ~ "Negative",
      TRUE ~ "Mixed"
    )
  ) %>%
  ungroup()

head(nonadd_epi_consistent)

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
      hjust = 0.5,
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
