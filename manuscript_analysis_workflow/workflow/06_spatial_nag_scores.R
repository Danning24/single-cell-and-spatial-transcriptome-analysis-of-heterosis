#!/usr/bin/env Rscript
# Figure 3H-I and Supplementary Figure S17: spatial NAG modules and radial comparisons.

library(Seurat)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Inputs and outputs
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
input_dir <- file.path(project_dir, "0.data/2026TUM_Danning_Wang_To_HK/To HK")
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", unset = "results"), "spatial_nag")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
spatial <- UpdateSeuratObject(readRDS(file.path(input_dir,
  "all_spatial_merge_filtered_integrated_remove_noise_bin50_annotation.RDS")))
nags <- read_csv(file.path(input_dir, "spatial_heterotic_genes_results.csv"),
                 show_col_types = FALSE) |>
  filter(pattern == "Non-additive") |>
  group_by(Trait) |>
  filter(all(str_detect(str_to_upper(cluster), "CORTEX"))) |>
  ungroup()

# Direction must agree across all retained cortex annotations for each gene.
positive <- nags |> group_by(Trait) |> filter(all(MPH.pct > 0)) |> pull(Trait) |> unique()
negative <- nags |> group_by(Trait) |> filter(all(MPH.pct < 0)) |> pull(Trait) |> unique()
gene_sets <- list(Positive = positive, Negative = negative, All = unique(nags$Trait))
write_csv(bind_rows(lapply(names(gene_sets), function(module) {
  data.frame(module, gene = gene_sets[[module]])
})), file.path(output_dir, "spatial_NAG_genes.csv"))

# Preserve the original reference pool: score all cortex bins before selecting hybrids.
layers <- c("Outer_cortex", "Cortex", "Inner_cortex")
cortex <- subset(spatial, idents = layers)
cortex <- AddModuleScore(cortex, features = unname(gene_sets),
  name = c("upNAG_Scores", "downNAG_Scores", "NAG_Scores"),
  assay = "SCT", slot = "data", seed = 1, nbin = 24, ctrl = 100)
features <- c("upNAG_Scores1", "downNAG_Scores2", "NAG_Scores3")
names(features) <- names(gene_sets)
hybrid_samples <- paste0("hybrid_slice_", c(1:6, 8:11), "_bin50")
scores <- cortex[[]] |>
  tibble::rownames_to_column("bin") |>
  filter(sample %in% hybrid_samples) |>
  transmute(bin, sample, layer = cluster_annotation,
    Positive = upNAG_Scores1, Negative = downNAG_Scores2, All = NAG_Scores3) |>
  pivot_longer(c(Positive, Negative, All), names_to = "module", values_to = "score")
write_csv(scores, file.path(output_dir, "spatial_bin_scores.csv"))

# Keep all tissue bins for display; unscored, non-cortex bins are grey.
for (feature in features) {
  values <- setNames(rep(NA_real_, ncol(spatial)), colnames(spatial))
  values[colnames(cortex)] <- cortex[[feature, drop = TRUE]]
  spatial[[feature]] <- values
}
hybrid <- subset(spatial, subset = sample %in% hybrid_samples)
hybrid$cortex_layer <- ifelse(hybrid$cluster_annotation %in% layers,
                             as.character(hybrid$cluster_annotation), "Other")
layer_colours <- c(Outer_cortex = "#F8766D", Cortex = "#4DB64D",
                   Inner_cortex = "#4E83C7", Other = "grey75")
celltype_plot <- SpatialDimPlot(hybrid, group.by = "cluster_annotation",
                                image.alpha = 0, pt.size.factor = 10, ncol = 5)
layer_plot <- SpatialDimPlot(hybrid, group.by = "cortex_layer", cols = layer_colours,
                             image.alpha = 0, pt.size.factor = 10, ncol = 5)
ggsave(file.path(output_dir, "spatial_cell_types.pdf"), celltype_plot, width = 14, height = 6)
ggsave(file.path(output_dir, "spatial_cortex_layers.pdf"), layer_plot, width = 14, height = 6)
for (module in names(features)) {
  map <- SpatialFeaturePlot(hybrid, features = features[[module]],
    image.alpha = 0, pt.size.factor = 10, ncol = 5) &
    scale_fill_gradientn(colours = c("#2166AC", "white", "#B2182B"), na.value = "grey75")
  ggsave(file.path(output_dir, paste0("spatial_", module, "_NAG.pdf")),
         map, width = 14, height = 6)
}

# The paired unit is a section: ten section medians per cortex layer.
medians <- scores |>
  group_by(sample, layer, module) |>
  summarise(n_bins = n(), median_score = median(score), .groups = "drop")
write_csv(medians, file.path(output_dir, "section_layer_medians.csv"))
comparisons <- combn(layers, 2, simplify = FALSE)
tests <- bind_rows(lapply(names(gene_sets), function(module_name) {
  paired <- medians |>
    filter(module == module_name) |>
    select(sample, layer, median_score) |>
    pivot_wider(names_from = layer, values_from = median_score)
  bind_rows(lapply(comparisons, function(pair) {
    test <- wilcox.test(paired[[pair[1]]], paired[[pair[2]]],
      paired = TRUE, alternative = "two.sided", exact = TRUE)
    data.frame(module = module_name, layer_1 = pair[1], layer_2 = pair[2],
      n_sections = nrow(paired), Wilcoxon_V = unname(test$statistic), P_value = test$p.value)
  })) |>
    mutate(P_adj_BH = p.adjust(P_value, method = "BH"))
}))
write_csv(tests, file.path(output_dir, "paired_section_tests.csv"))

# Plot section summaries rather than treating bins as biological replicates.
medians$layer <- factor(medians$layer, levels = layers)
radial_plot <- ggplot(medians, aes(layer, median_score, fill = layer)) +
  geom_violin(trim = FALSE, colour = "black", linewidth = 0.4) +
  geom_boxplot(width = 0.16, fill = "white", colour = "black", outlier.shape = NA) +
  geom_line(aes(group = sample), colour = "grey60", linewidth = 0.25, alpha = 0.5) +
  geom_point(size = 1) +
  scale_fill_manual(values = layer_colours[layers]) +
  scale_x_discrete(labels = c("Outer\ncortex", "Middle\ncortex", "Inner\ncortex")) +
  facet_wrap(~ module, scales = "free_y") +
  labs(x = NULL, y = "Section-median NAG score",
       caption = "Ten paired sections; two-sided Wilcoxon tests; BH adjustment within each module.") +
  theme_classic() + theme(legend.position = "none")
ggsave(file.path(output_dir, "section_median_radial_comparisons.pdf"),
       radial_plot, width = 11, height = 4)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
