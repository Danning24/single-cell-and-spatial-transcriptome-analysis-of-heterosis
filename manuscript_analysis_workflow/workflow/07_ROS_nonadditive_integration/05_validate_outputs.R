project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

suppressPackageStartupMessages(library(data.table))
if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(root, "ros_nonadditive_analysis")
result_dir <- file.path(analysis_dir, "results")

core <- fread(file.path(result_dir, "02_ko_nonadditive_enrichment_all.csv"))
stopifnot(nrow(core) == 2520L)
for (nm in intersect(c("fisher_p", "wilcox_p", "fgsea_p"), names(core))) {
  stopifnot(all(is.na(core[[nm]]) | (core[[nm]] >= 0 & core[[nm]] <= 1)))
}

inventory <- fread(file.path(result_dir, "15_hybrid_spatial_sample_inventory.csv"))
stopifnot(sum(inventory$present) == 10L)
stopifnot(inventory[sample == "hybrid_slice_7_bin50", !present])

manifest <- fread(file.path(result_dir, "16_spatial_module_gene_manifest.csv"))
stopifnot(nrow(manifest) == 5L, all(manifest$score_status == "scored"))

scores <- fread(file.path(result_dir, "17_hybrid_spatial_ucell_scores_long.csv"))
stopifnot(uniqueN(scores$bin) == 11238L, uniqueN(scores$module) == 5L)
stopifnot(nrow(scores) == 11238L * 5L)

sensitivity <- fread(file.path(result_dir, "21_ucell_addmodulescore_sensitivity.csv"))
stopifnot(all(sensitivity$spearman_rho > 0.85))

coverage <- fread(file.path(result_dir, "25_functional_enrichment_coverage.csv"))
stopifnot(coverage[module == "ROS_Epidermis_recurrent_LE", n_background] == 2233L)

required_figures <- c(
  "01_ros_primary_nonadditive_NES.png",
  "02_ros_mph_direction_NES_heatmap.png",
  "03_ros_core_candidate_heatmap.png",
  "05_double_vs_single_rank_comparison.png",
  "06_spatial_layer_module_consistency.png",
  "08_spatial_map_ROS_Tier1_Tier2_consensus.png",
  "10_candidate_GO_enrichment.png"
)
stopifnot(all(file.exists(file.path(analysis_dir, "figures", required_figures))))
stopifnot(file.exists(file.path(
  analysis_dir, "report", "ROS_KO_非加性基因_空间联合分析报告.docx"
)))

writeLines(c(
  "Validation passed",
  paste("Core enrichment rows:", nrow(core)),
  paste("Hybrid spatial bins:", uniqueN(scores$bin)),
  paste("Hybrid slices present:", sum(inventory$present)),
  paste("Spatial modules scored:", uniqueN(scores$module)),
  paste("UCell/AddModuleScore rho range:",
        sprintf("%.3f-%.3f", min(sensitivity$spearman_rho), max(sensitivity$spearman_rho))),
  paste("Required figures present:", length(required_figures))
), file.path(analysis_dir, "logs", "07_validation.txt"))

message("All output validation checks passed")
