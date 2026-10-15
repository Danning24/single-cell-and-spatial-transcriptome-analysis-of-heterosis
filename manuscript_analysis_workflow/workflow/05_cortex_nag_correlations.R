#!/usr/bin/env Rscript
# Figure 3E-G and Supplementary Figure S15d: cortex NAG scores and correlations.

library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# Inputs and outputs
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
cortex_file <- Sys.getenv("MAIZE_CORTEX_RDS",
  unset = file.path(project_dir, "4.Cortex_cells/data_output/Cortex.rds"))
nag_file <- file.path(project_dir,
  "0.data/2026TUM_Danning_Wang_To_HK/To HK/sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv")
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", unset = "results"), "cortex_nag")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cortex <- readRDS(cortex_file)
nags <- read_csv(nag_file, show_col_types = FALSE) |>
  filter(pattern == "Non-additive")

# Pool all hybrids: retain genes whose NAG records occur only in cortex states.
pooled_nags <- nags |>
  group_by(Trait) |>
  filter(all(str_detect(cluster_anno, "Cortex"))) |>
  ungroup()
pooled_genes <- unique(pooled_nags$Trait)
cortex <- AddModuleScore(cortex, features = list(pooled_genes),
  name = "NonAddgene_Scores", assay = "SCT", slot = "data",
  seed = 1, nbin = 24, ctrl = 100)
pooled_scores <- data.frame(cell_id = colnames(cortex), material = cortex$orig.ident,
  NAG_score = cortex$NonAddgene_Scores1,
  E_Score = cortex$E_Scores1, M_Score = cortex$M_Scores1)
write_csv(pooled_scores, file.path(output_dir, "pooled_cortex_scores.csv"))
write_csv(pooled_nags, file.path(output_dir, "pooled_cortex_NAG_records.csv"))

nag_map <- FeaturePlot(cortex, features = "NonAddgene_Scores1") +
  scale_colour_gradientn(colours = c("#000000", "#DDDDDD", "#E69F00")) +
  labs(title = "Non-additive gene score") + NoAxes()
ggsave(file.path(output_dir, "pooled_NAG_UMAP.pdf"), nag_map, width = 6, height = 5)

# Pearson tests describe cell-level associations, not biological replication.
pearson_result <- function(x, y) {
  test <- cor.test(x, y, method = "pearson", alternative = "two.sided")
  data.frame(n_cells = sum(complete.cases(x, y)),
    Pearson_r = unname(test$estimate), P_value = test$p.value)
}
pooled_tests <- bind_rows(lapply(c("E_Score", "M_Score"), function(programme) {
  cbind(data.frame(programme), pearson_result(pooled_scores$NAG_score,
                                             pooled_scores[[programme]]))
}))
write_csv(pooled_tests, file.path(output_dir, "pooled_correlations.csv"))
pooled_long <- pooled_scores |>
  pivot_longer(c(E_Score, M_Score), names_to = "programme", values_to = "score")
pooled_plot <- ggplot(pooled_long, aes(NAG_score, score)) +
  geom_point(size = 0.1, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "red") +
  facet_wrap(~ programme, scales = "free_y") +
  labs(x = "Non-additive gene score", y = "Developmental programme score") +
  theme_classic()
ggsave(file.path(output_dir, "pooled_NAG_programme_correlations.pdf"),
       pooled_plot, width = 9, height = 4)

# Re-select the cortex-exclusive set and score cells within each hybrid.
hybrids <- c("B73xA554", "B73xA680", "B73xB37", "B73xNC298", "B73xNC320",
             "B73xNC328", "B73xNC352", "B73xOh40B", "B73xOh43")
hybrid_tests <- bind_rows(lapply(hybrids, function(hybrid) {
  genes <- nags |>
    filter(Hybrid == hybrid) |>
    group_by(Trait) |>
    filter(all(str_detect(cluster_anno, "Cortex"))) |>
    ungroup() |>
    pull(Trait) |> unique()
  cells <- subset(cortex, subset = orig.ident == hybrid)
  cells <- AddModuleScore(cells, features = list(genes), name = "NonAddgene_Scores",
    assay = "SCT", slot = "data", seed = 1, nbin = 24, ctrl = 100)
  bind_rows(lapply(c("E_Scores1", "M_Scores1"), function(programme) {
    cbind(data.frame(hybrid, n_genes = length(genes), programme),
      pearson_result(cells$NonAddgene_Scores1, cells[[programme, drop = TRUE]]))
  }))
}))
write_csv(hybrid_tests, file.path(output_dir, "hybrid_correlations.csv"))
hybrid_summary <- hybrid_tests |>
  select(hybrid, n_genes, programme, Pearson_r) |>
  pivot_wider(names_from = programme, values_from = Pearson_r)
hybrid_plot <- ggplot(hybrid_summary, aes(E_Scores1, M_Scores1, size = n_genes)) +
  geom_point(shape = 21, fill = "red", colour = "black") +
  geom_text(aes(label = hybrid), size = 3, vjust = -0.8, show.legend = FALSE) +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "Pearson r: NAG score versus E_Score",
       y = "Pearson r: NAG score versus M_Score", size = "NAGs") +
  theme_classic()
ggsave(file.path(output_dir, "hybrid_correlation_summary.pdf"),
       hybrid_plot, width = 7, height = 5)

# Oh40B directional sets: remove genes represented in both MPH directions.
oh40b_nags <- nags |>
  filter(Hybrid == "B73xOh40B") |>
  group_by(Trait) |>
  filter(all(str_detect(cluster_anno, "Cortex"))) |>
  ungroup()
positive <- filter(oh40b_nags, MPH.pct > 0)$Trait
negative <- filter(oh40b_nags, MPH.pct < 0)$Trait
gene_sets <- list(Positive = setdiff(positive, negative), Negative = setdiff(negative, positive))
oh40b <- subset(cortex, subset = orig.ident == "B73xOh40B")
oh40b <- AddModuleScore(oh40b, features = unname(gene_sets),
  name = c("up_score", "down_score"), assay = "SCT", slot = "data",
  seed = 1, nbin = 24, ctrl = 100)
directional_scores <- data.frame(cell_id = colnames(oh40b),
  Positive = oh40b$up_score1, Negative = oh40b$down_score2, E_Score = oh40b$E_Scores1)
directional_tests <- bind_rows(lapply(names(gene_sets), function(direction) {
  cbind(data.frame(direction, n_genes = length(gene_sets[[direction]])),
    pearson_result(directional_scores[[direction]], directional_scores$E_Score))
}))
write_csv(directional_scores, file.path(output_dir, "Oh40B_directional_scores.csv"))
write_csv(directional_tests, file.path(output_dir, "Oh40B_directional_correlations.csv"))
write_csv(bind_rows(lapply(names(gene_sets), function(direction) {
  data.frame(direction, gene = gene_sets[[direction]])
})), file.path(output_dir, "Oh40B_directional_genes.csv"))
directional_long <- directional_scores |>
  pivot_longer(c(Positive, Negative), names_to = "direction", values_to = "NAG_score")
directional_plot <- ggplot(directional_long, aes(NAG_score, E_Score)) +
  geom_point(size = 0.1, alpha = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "red") +
  facet_wrap(~ direction, scales = "free_x") +
  labs(x = "Directional non-additive gene score", y = "E_Score") + theme_classic()
ggsave(file.path(output_dir, "Oh40B_directional_NAG_Escore.pdf"),
       directional_plot, width = 9, height = 4)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
