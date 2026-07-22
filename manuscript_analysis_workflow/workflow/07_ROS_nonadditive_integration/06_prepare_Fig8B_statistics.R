# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

suppressPackageStartupMessages(library(data.table))

input_file <- "ros_nonadditive_analysis/results/04_ros_primary_tissue_enrichment.csv"
output_dir <- "ros_nonadditive_analysis/story_package/evidence"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file, na.strings = c("", "NA", "NaN"))
dat <- dat[set_type == "all"]
setorder(dat, Hybrid, KO, celltype)

fisher <- dat[, .(
  Hybrid, KO, celltype, overlap, sig_total, set_total, universe_total,
  odds_ratio, ci_low, ci_high, fold_enrichment, fisher_p,
  fisher_fdr_family, fisher_fdr_global,
  significant_FDR05 = fisher_fdr_family < 0.05
)]

wilcoxon <- dat[, .(
  Hybrid, KO, celltype, set_total, universe_total, rank_auc,
  rank_biserial, wilcox_p, wilcox_fdr_family,
  significant_FDR05 = wilcox_fdr_family < 0.05
)]

gsea <- dat[, .(
  Hybrid, KO, celltype, fgsea_size, NES, fgsea_p, fgsea_fdr_family,
  leading_edge, rank_method,
  significant_FDR05 = fgsea_fdr_family < 0.05
)]

fwrite(fisher, file.path(output_dir, "figure_b_fisher_all_NAG.csv"), na = "")
fwrite(wilcoxon, file.path(output_dir, "figure_b_wilcoxon_all_NAG.csv"), na = "")
fwrite(gsea, file.path(output_dir, "figure_b_gsea_all_NAG.csv"), na = "")

summary <- data.table(
  test = c("Fisher exact test", "Wilcoxon rank-sum test", "GSEA"),
  tested_combinations = nrow(dat),
  significant_FDR05 = c(
    sum(dat$fisher_fdr_family < 0.05, na.rm = TRUE),
    sum(dat$wilcox_fdr_family < 0.05, na.rm = TRUE),
    sum(dat$fgsea_fdr_family < 0.05, na.rm = TRUE)
  )
)
fwrite(summary, file.path(output_dir, "figure_b_test_summary.csv"))
