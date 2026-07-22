# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
}

set.seed(20260714)
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(project_root, "ros_nonadditive_analysis")
result_dir <- file.path(analysis_dir, "results")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")
spatial_file <- file.path(
  project_root, "2026TUM_Danning_Wang_To_HK", "To HK",
  "all_spatial_merge_filtered_integrated_remove_noise_bin50_annotation.RDS"
)
tier_file <- file.path(result_dir, "05_ros_candidate_evidence_tiers.csv")
primary_file <- file.path(result_dir, "04_ros_primary_tissue_enrichment.csv")
stopifnot(file.exists(spatial_file), file.exists(tier_file), file.exists(primary_file))

split_genes <- function(x) {
  unique(trimws(unlist(strsplit(x[!is.na(x) & nzchar(x)], ";", fixed = TRUE))))
}

tiers <- fread(tier_file)
primary <- fread(primary_file)

recurrent_le <- primary[
  Hybrid == "B73xA554" & set_type == "all",
  .(gene = split_genes(leading_edge)),
  by = KO
][, .(n_ros_kos = uniqueN(KO)), by = gene][n_ros_kos >= 3, gene]

module_genes <- list(
  ROS_Tier1_shared_core = tiers[Tier == "Tier 1", gene],
  ROS_Tier1_Tier2_consensus = tiers[Tier %chin% c("Tier 1", "Tier 2"), gene],
  ROS_Epidermis_recurrent_LE = recurrent_le,
  ROS_sensitive_positive_MPH = tiers[
    Tier %chin% c("Tier 1", "Tier 2") & n_positive_nonadd_B73xA554 > 0, gene
  ],
  ROS_sensitive_negative_MPH = tiers[
    Tier %chin% c("Tier 1", "Tier 2") & n_negative_nonadd_B73xA554 > 0, gene
  ]
)
module_genes <- lapply(module_genes, unique)

message("Reading local spatial object")
obj <- readRDS(spatial_file)
stopifnot(all(c("sample", "Genotype", "cluster_annotation") %in% colnames(obj[[]])))

hybrid_cells <- rownames(obj[[]])[obj$Genotype == "hybrid"]
hyb <- obj
rm(obj)
DefaultAssay(hyb) <- "SCT"

hybrid_samples <- sort(unique(as.character(hyb$sample[hyb$Genotype == "hybrid"])))
expected_samples <- paste0("hybrid_slice_", 1:11, "_bin50")
sample_inventory <- data.table(
  sample = expected_samples,
  present = expected_samples %chin% hybrid_samples,
  n_bins = vapply(expected_samples, function(s) sum(hyb$sample == s & hyb$Genotype == "hybrid"), integer(1))
)
fwrite(sample_inventory, file.path(result_dir, "15_hybrid_spatial_sample_inventory.csv"))

available_features <- rownames(hyb[["SCT"]])
detected_modules <- lapply(module_genes, intersect, y = available_features)
module_manifest <- rbindlist(lapply(names(module_genes), function(nm) {
  requested <- module_genes[[nm]]
  detected <- detected_modules[[nm]]
  data.table(
    module = nm,
    n_requested = length(requested),
    n_detected = length(detected),
    detection_fraction = length(detected) / length(requested),
    requested_genes = paste(requested, collapse = ";"),
    detected_genes = paste(detected, collapse = ";"),
    score_status = ifelse(length(detected) >= 5, "scored", "not_scored_lt5_detected")
  )
}))
fwrite(module_manifest, file.path(result_dir, "16_spatial_module_gene_manifest.csv"))

score_modules <- detected_modules[lengths(detected_modules) >= 5]
stopifnot(length(score_modules) > 0)
message("UCell scoring modules: ", paste(names(score_modules), collapse = ", "))
hyb <- AddModuleScore_UCell(
  hyb,
  features = score_modules,
  assay = "SCT",
  slot = "counts",
  ncores = 1,
  name = "_UCell"
)
ucell_cols <- paste0(names(score_modules), "_UCell")

message("AddModuleScore sensitivity analysis")
ams_cols <- character()
for (nm in names(score_modules)) {
  prefix <- paste0("AMS_", nm, "_")
  hyb <- AddModuleScore(
    hyb,
    features = list(score_modules[[nm]]),
    assay = "SCT",
    slot = "data",
    name = prefix,
    seed = 20260714,
    search = FALSE
  )
  new_col <- paste0(prefix, "1")
  ams_cols[nm] <- new_col
}

md <- as.data.table(hyb[[]], keep.rownames = "bin")
md_hybrid <- md[Genotype == "hybrid"]
score_long <- melt(
  md_hybrid,
  id.vars = c("bin", "sample", "cluster_annotation", "x", "y"),
  measure.vars = ucell_cols,
  variable.name = "score_column",
  value.name = "ucell_score"
)
score_long[, module := sub("_UCell$", "", as.character(score_column))]
score_long[, score_column := NULL]
fwrite(score_long, file.path(result_dir, "17_hybrid_spatial_ucell_scores_long.csv"))

sample_layer <- score_long[, .(
  n_bins = .N,
  mean_ucell = mean(ucell_score),
  median_ucell = median(ucell_score),
  q25_ucell = quantile(ucell_score, 0.25),
  q75_ucell = quantile(ucell_score, 0.75)
), by = .(sample, cluster_annotation, module)]
fwrite(sample_layer, file.path(result_dir, "18_spatial_sample_layer_module_summary.csv"))

sample_overall <- score_long[, .(slice_median = median(ucell_score)), by = .(sample, module)]
layer_delta <- merge(sample_layer, sample_overall, by = c("sample", "module"), all.x = TRUE)
layer_delta[, delta_from_slice_median := median_ucell - slice_median]
fwrite(layer_delta, file.path(result_dir, "19_spatial_layer_delta_by_slice.csv"))

consistency <- layer_delta[, {
  d <- delta_from_slice_median[is.finite(delta_from_slice_median)]
  d_signed <- d[d != 0]
  k <- sum(d_signed > 0)
  n <- length(d)
  n_signed <- length(d_signed)
  p <- if (n_signed >= 5) binom.test(k, n_signed, p = 0.5)$p.value else NA_real_
  list(
    n_slices = n,
    n_nonzero_slices = n_signed,
    median_delta = median(d),
    positive_slices = k,
    negative_slices = sum(d_signed < 0),
    positive_fraction_nonzero = ifelse(n_signed > 0, k / n_signed, NA_real_),
    sign_test_p = p
  )
}, by = .(module, cluster_annotation)]
consistency[, sign_test_fdr := p.adjust(sign_test_p, method = "BH")]
fwrite(consistency, file.path(result_dir, "20_spatial_layer_consistency_descriptive.csv"))

sensitivity <- rbindlist(lapply(names(score_modules), function(nm) {
  u <- md_hybrid[[paste0(nm, "_UCell")]]
  a <- md_hybrid[[ams_cols[[nm]]]]
  data.table(
    module = nm,
    spearman_rho = cor(u, a, method = "spearman", use = "complete.obs"),
    n_bins = sum(complete.cases(u, a))
  )
}))
fwrite(sensitivity, file.path(result_dir, "21_ucell_addmodulescore_sensitivity.csv"))

image_map <- rbindlist(lapply(Images(hyb), function(im) {
  image_cells <- Cells(hyb[[im]])
  samples <- unique(as.character(hyb$sample[match(image_cells, colnames(hyb))]))
  samples <- samples[!is.na(samples)]
  data.table(
    image = im,
    sample = if (length(samples) == 1) samples else paste(samples, collapse = ";"),
    n_image_bins = length(image_cells)
  )
}))
fwrite(image_map, file.path(result_dir, "22_spatial_image_sample_map.csv"))
hybrid_images <- image_map[sample %chin% hybrid_samples, image]

module_order <- names(score_modules)
layer_order <- c("Epidermis", "Outer_cortex", "Cortex", "Inner_cortex",
                 "Endodermis", "Pericycle", "Phloem", "Xylem", "MetaXylem",
                 "Cells_betw_MX", "Pith")

p_heat <- ggplot(
  consistency,
  aes(x = factor(cluster_annotation, levels = layer_order),
      y = factor(module, levels = rev(module_order)), fill = median_delta)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(positive_slices, "/", n_slices)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Median delta",
       title = "ROS–non-additive spatial module localization",
       subtitle = "Text: slices with positive within-slice delta / slices available") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(figure_dir, "06_spatial_layer_module_consistency.png"), p_heat,
       width = 13, height = 5.8, dpi = 300)
ggsave(file.path(figure_dir, "06_spatial_layer_module_consistency.pdf"), p_heat,
       width = 13, height = 5.8)

p_box <- ggplot(
  sample_layer,
  aes(x = factor(cluster_annotation, levels = layer_order), y = median_ucell, color = module)
) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.55, size = 1) +
  facet_wrap(~module, scales = "free_y", ncol = 2) +
  guides(color = "none") +
  labs(x = NULL, y = "Per-slice median UCell score",
       title = "Spatial distribution across annotated root tissues") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(figure_dir, "07_spatial_layer_module_scores.png"), p_box,
       width = 12, height = 10, dpi = 300)
ggsave(file.path(figure_dir, "07_spatial_layer_module_scores.pdf"), p_box,
       width = 12, height = 10)

message("Writing spatial feature maps")
for (nm in names(score_modules)) {
  plot_dt <- score_long[module == nm]
  p_spatial <- ggplot(plot_dt, aes(x = x, y = y, color = ucell_score)) +
    geom_point(size = 0.65, alpha = 0.95) +
    facet_wrap(~sample, ncol = 5, scales = "free") +
    scale_color_gradientn(
      colors = c("#313695", "#4575B4", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
      na.value = "grey85"
    ) +
    scale_y_reverse() +
    labs(title = gsub("_", " ", nm), color = "UCell") +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 15, face = "bold"),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      aspect.ratio = 1
    )
  ggsave(file.path(figure_dir, paste0("08_spatial_map_", nm, ".png")), p_spatial,
         width = 16, height = 8.5, dpi = 220, limitsize = FALSE)
}

p_anno <- SpatialDimPlot(
  hyb,
  group.by = "cluster_annotation",
  images = hybrid_images,
  ncol = 5,
  pt.size.factor = 1.3,
  crop = TRUE,
  combine = TRUE
) + plot_annotation(title = "Hybrid spatial tissue annotations")
ggsave(file.path(figure_dir, "09_hybrid_spatial_tissue_annotations.png"), p_anno,
       width = 16, height = 9, dpi = 220, limitsize = FALSE)

saveRDS(
  md_hybrid[, c("bin", "sample", "cluster_annotation", ucell_cols, unname(ams_cols)), with = FALSE],
  file.path(result_dir, "hybrid_spatial_ROS_nonadditive_scores.rds")
)

capture.output({
  cat("Hybrid samples present:\n"); print(hybrid_samples)
  cat("\nExpected sample inventory:\n"); print(sample_inventory)
  cat("\nModule manifest:\n"); print(module_manifest[, .(module, n_requested, n_detected, score_status)])
  cat("\nUCell vs AddModuleScore sensitivity:\n"); print(sensitivity)
  cat("\nImage map:\n"); print(image_map)
  cat("\nSession info:\n"); print(sessionInfo())
}, file = file.path(log_dir, "03_spatial_ros_nonadditive.txt"))

message("Spatial ROS–non-additive analysis completed")
