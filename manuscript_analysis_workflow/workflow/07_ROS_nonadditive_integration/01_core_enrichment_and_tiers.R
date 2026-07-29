project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fgsea)
  library(ComplexHeatmap)
  library(circlize)
})

options(stringsAsFactors = FALSE)
set.seed(20260714)
if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
}

project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(project_root, "ros_nonadditive_analysis")
result_dir <- file.path(analysis_dir, "results")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Project root: ", project_root)

heterosis_file <- file.path(
  project_root, "数据汇总", "汇总文件",
  "TbaleS2 sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv"
)
stopifnot(file.exists(heterosis_file))

module_map <- data.table(
  KO = c(
    "Zm00001d007234", "Zm00001d026154", "Zm00001d029170", "Zm00001d042096", "Zm00001d043795", "Zm00001d048020",
    "Zm00001d011969", "Zm00001d028307", "Zm00001d034109", "Zm00001d006753", "Zm00001d011208", "Zm00001d018973", "Zm00001d048021", "Zm00001d015852",
    "Zm00001d013410", "Zm00001d009858", "Zm00001d028392", "Zm00001d019565", "Zm00001d037547",
    "Zm00001d039221", "Zm00001d012221", "Zm00001d019283", "Zm00001d030725",
    "Zm00001d026577", "Zm00001d002275",
    "Zm00001d043551", "Zm00001d051804", "Zm00001d030557"
  ),
  Module = c(
    rep("ROS", 6), rep("Hormone", 8), rep("Cell wall", 5), rep("ER stress", 4),
    rep("Defense", 2), rep("Special", 3)
  )
)

primary_cortex <- c("Cortex early", "Cortex I mature", "Cortex II mature", "Cortex outer")
is_primary <- function(hybrid, celltype) {
  (hybrid == "B73xA554" & celltype == "Epidermis") |
    (hybrid == "B73xOh40B" & celltype %chin% primary_cortex)
}

message("Reading heterosis table")
heterosis <- fread(heterosis_file)
setnames(heterosis, "Trait", "gene")
stopifnot(!anyDuplicated(heterosis[, .(Hybrid, cluster_anno, gene)]))
heterosis[, direction := fifelse(pattern == "Non-additive" & MPH.pct > 0, "positive",
                          fifelse(pattern == "Non-additive" & MPH.pct < 0, "negative", "additive"))]

read_ros_long <- function(hybrid) {
  f <- file.path(
    project_root, "sctenifold_downstream_analysis",
    paste0(hybrid, "_sctenifold_KO_analysis"), "combined_comparison",
    "all_KO_diffRegulation_downstream_no_selfKO_long.csv"
  )
  x <- fread(f)
  x[, `:=`(
    Hybrid = hybrid,
    Module = "ROS",
    Source = "chaosuan_ROS_local_copy",
    Runner = ifelse(hybrid == "B73xA554", "safe_ma1", "standard")
  )]
  x[, gene_label := NULL]
  x
}

ros <- rbindlist(lapply(c("B73xA554", "B73xOh40B"), read_ros_long), use.names = TRUE)

non_ros_dir <- file.path(analysis_dir, "input", "non_ros_wuzh02")
non_ros_files <- list.files(non_ros_dir, pattern = "diffRegulation\\.csv$", full.names = TRUE)
stopifnot(length(non_ros_files) == 44L)

read_non_ros <- function(f) {
  nm <- basename(f)
  m <- regexec("^(B73xA554|B73xOh40B)__([^.]+)(\\.safe_ma1)?\\.diffRegulation\\.csv$", nm)
  z <- regmatches(nm, m)[[1]]
  if (length(z) == 0L) stop("Cannot parse file: ", nm)
  hybrid <- z[2]
  ko <- z[3]
  x <- fread(f)
  x <- x[gene != ko]
  x[, `:=`(
    KO = ko,
    Hybrid = hybrid,
    Source = "wuzh02",
    Runner = ifelse(hybrid == "B73xA554", "safe_ma1", "standard")
  )]
  x[module_map, Module := i.Module, on = "KO"]
  if (anyNA(x$Module)) stop("Missing module mapping for ", ko)
  x
}

non_ros <- rbindlist(lapply(non_ros_files, read_non_ros), use.names = TRUE)
ko_results <- rbindlist(list(ros, non_ros), use.names = TRUE, fill = TRUE)
setcolorder(ko_results, c("Hybrid", "Module", "KO", "gene", "distance", "Z", "FC", "p.value", "p.adj", "Source", "Runner"))

dup_check <- ko_results[, .N, by = .(Hybrid, KO, gene)][N > 1]
if (nrow(dup_check)) stop("Duplicated Hybrid/KO/gene records detected")

manifest <- ko_results[, .(
  n_tested = uniqueN(gene),
  n_sig = sum(p.adj < 0.05, na.rm = TRUE),
  min_padj = min(p.adj, na.rm = TRUE)
), by = .(Hybrid, Module, KO, Source, Runner)]
setorder(manifest, Module, KO, Hybrid)
fwrite(manifest, file.path(result_dir, "01_input_manifest.csv"))

fisher_one <- function(sig, member) {
  a <- sum(sig & member)
  b <- sum(sig & !member)
  c <- sum(!sig & member)
  d <- sum(!sig & !member)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  baseline <- mean(member)
  observed <- if ((a + b) > 0) a / (a + b) else NA_real_
  data.table(
    overlap = a,
    sig_total = a + b,
    set_total = a + c,
    universe_total = a + b + c + d,
    odds_ratio = unname(ft$estimate),
    ci_low = ft$conf.int[1],
    ci_high = ft$conf.int[2],
    fold_enrichment = observed / baseline,
    fisher_p = ft$p.value
  )
}

rank_one <- function(z, member) {
  n1 <- sum(member)
  n0 <- sum(!member)
  if (n1 < 2L || n0 < 2L) {
    return(data.table(rank_auc = NA_real_, rank_biserial = NA_real_, wilcox_p = NA_real_))
  }
  ranks <- rank(z, ties.method = "average")
  u <- sum(ranks[member]) - n1 * (n1 + 1) / 2
  auc <- u / (n1 * n0)
  wt <- suppressWarnings(wilcox.test(z[member], z[!member], alternative = "greater",
                                     exact = FALSE, correct = FALSE))
  data.table(rank_auc = auc, rank_biserial = 2 * auc - 1, wilcox_p = wt$p.value)
}

analyze_one_ko <- function(x, module, run_rank = FALSE, celltypes_override = NULL) {
  hybrid <- x$Hybrid[1]
  ko <- x$KO[1]
  celltypes <- if (is.null(celltypes_override)) {
    sort(unique(heterosis[Hybrid == hybrid, cluster_anno]))
  } else {
    celltypes_override
  }
  ans <- vector("list", length(celltypes))

  for (i in seq_along(celltypes)) {
    ct <- celltypes[i]
    h <- heterosis[Hybrid == hybrid & cluster_anno == ct, .(gene, pattern, direction, MPH.pct, MPH.padj)]
    d <- merge(x[, .(gene, Z, p.adj)], h, by = "gene", all = FALSE)
    d <- d[is.finite(Z) & !is.na(pattern)]
    d[, sig := p.adj < 0.05]

    memberships <- list(
      all = d$pattern == "Non-additive",
      positive = d$direction == "positive",
      negative = d$direction == "negative"
    )
    fr <- rbindlist(lapply(names(memberships), function(nm) {
      cbind(
        data.table(set_type = nm),
        fisher_one(d$sig, memberships[[nm]]),
        rank_one(d$Z, memberships[[nm]])
      )
    }))

    fg <- data.table(set_type = names(memberships), NES = NA_real_, fgsea_p = NA_real_,
                     fgsea_size = NA_integer_, leading_edge = NA_character_, rank_method = NA_character_)
    formal_primary <- module %chin% c("ROS", "ROS double KO") && is_primary(hybrid, ct)
    if (run_rank && formal_primary) {
      ranks <- d$Z
      names(ranks) <- d$gene
      ranks <- sort(ranks, decreasing = TRUE)
      pathways <- lapply(memberships, function(v) d$gene[v])
      pathways <- pathways[lengths(pathways) >= 10L]
      if (length(pathways)) {
        fg0 <- suppressWarnings(fgseaMultilevel(
          pathways = pathways, stats = ranks,
          minSize = 10, maxSize = length(ranks), eps = 1e-10, nproc = 1
        ))
        fg0 <- as.data.table(fg0)
        if (nrow(fg0)) {
          fg0[, leading_edge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
          fg <- fg0[, .(set_type = pathway, NES, fgsea_p = pval, fgsea_size = size, leading_edge)]
          fg[, rank_method := "fgseaMultilevel"]
        }
      }
    }
    ans[[i]] <- merge(fr, fg, by = "set_type", all.x = TRUE)[, `:=`(
      Hybrid = hybrid, Module = module, KO = ko, celltype = ct,
      analysis_family = ifelse(is_primary(hybrid, ct), "primary", "exploratory")
    )]
  }
  rbindlist(ans)
}

message("Running enrichment tests")
split_ko <- split(ko_results, interaction(ko_results$Hybrid, ko_results$KO, drop = TRUE))
enrichment_cache <- file.path(result_dir, "02_ko_nonadditive_enrichment_all.csv")
if (file.exists(enrichment_cache)) {
  cached <- fread(enrichment_cache)
  if (nrow(cached) == 2520L && all(c("rank_auc", "wilcox_p", "fisher_p") %in% names(cached))) {
    message("Using validated enrichment cache")
    enrichment <- cached
  } else {
    stop("Existing enrichment cache is incomplete; remove it before recomputing")
  }
} else {
  enrichment <- rbindlist(lapply(split_ko, function(x) {
    module <- x$Module[1]
    analyze_one_ko(x, module, run_rank = module == "ROS")
  }), fill = TRUE)
}

enrichment[, fisher_fdr_global := p.adjust(fisher_p, method = "BH")]
enrichment[, fisher_fdr_family := p.adjust(fisher_p, method = "BH"),
           by = .(Module, analysis_family, set_type)]
enrichment[, wilcox_fdr_family := p.adjust(wilcox_p, method = "BH"),
           by = .(Module, analysis_family, set_type)]
enrichment[, fgsea_fdr_family := p.adjust(fgsea_p, method = "BH"),
           by = .(Module, analysis_family, set_type)]
setorder(enrichment, Module, Hybrid, KO, celltype, set_type)
fwrite(enrichment, file.path(result_dir, "02_ko_nonadditive_enrichment_all.csv"))
fwrite(enrichment[Module == "ROS"], file.path(result_dir, "03_ros_nonadditive_enrichment.csv"))
fwrite(enrichment[Module == "ROS" & analysis_family == "primary"],
       file.path(result_dir, "04_ros_primary_tissue_enrichment.csv"))

message("Building ROS candidate evidence table")
ros_ranked <- copy(ros)
ros_ranked[, z_percentile := frank(Z, ties.method = "average") / .N, by = .(Hybrid, KO)]
ros_gene <- ros_ranked[, .(
  n_kos_tested = uniqueN(KO),
  n_kos_sig = sum(p.adj < 0.05, na.rm = TRUE),
  n_kos_top5pct = sum(z_percentile >= 0.95, na.rm = TRUE),
  median_z_percentile = median(z_percentile, na.rm = TRUE),
  max_z = max(Z, na.rm = TRUE),
  min_padj = min(p.adj, na.rm = TRUE)
), by = .(Hybrid, gene)]

het_gene <- heterosis[Hybrid %chin% c("B73xA554", "B73xOh40B"), .(
  n_celltypes_tested = .N,
  n_nonadd_celltypes = sum(pattern == "Non-additive"),
  n_positive_nonadd = sum(direction == "positive"),
  n_negative_nonadd = sum(direction == "negative"),
  max_abs_mph = max(abs(MPH.pct), na.rm = TRUE),
  primary_nonadd = any(pattern == "Non-additive" & is_primary(Hybrid, cluster_anno)),
  primary_celltypes = paste(cluster_anno[pattern == "Non-additive" & is_primary(Hybrid, cluster_anno)], collapse = ";")
), by = .(Hybrid, gene)]
ros_gene <- merge(ros_gene, het_gene, by = c("Hybrid", "gene"), all.x = TRUE)

leading <- enrichment[Module == "ROS" & analysis_family == "primary" & set_type == "all" & !is.na(leading_edge),
  .(Hybrid, KO, genes = strsplit(leading_edge, ";", fixed = TRUE))]
if (nrow(leading)) {
  leading <- leading[, .(gene = unlist(genes)), by = .(Hybrid, KO)]
  leading_support <- leading[, .(n_primary_leading_edges = uniqueN(paste(Hybrid, KO))), by = .(Hybrid, gene)]
  ros_gene <- merge(ros_gene, leading_support, by = c("Hybrid", "gene"), all.x = TRUE)
}
ros_gene[is.na(n_primary_leading_edges), n_primary_leading_edges := 0L]

wide <- dcast(
  ros_gene,
  gene ~ Hybrid,
  value.var = c("n_kos_sig", "n_kos_top5pct", "median_z_percentile", "n_nonadd_celltypes",
                "n_positive_nonadd", "n_negative_nonadd", "max_abs_mph", "primary_nonadd",
                "primary_celltypes", "n_primary_leading_edges"),
  fill = NA
)
for (nm in c("primary_nonadd_B73xA554", "primary_nonadd_B73xOh40B")) {
  if (!nm %in% names(wide)) wide[, (nm) := FALSE]
  wide[is.na(get(nm)), (nm) := FALSE]
}
wide[, Tier := fcase(
  n_kos_sig_B73xA554 >= 3 & n_kos_sig_B73xOh40B >= 3 &
    primary_nonadd_B73xA554 & primary_nonadd_B73xOh40B, "Tier 1",
  (n_kos_sig_B73xA554 >= 3 & primary_nonadd_B73xA554) |
    (n_kos_sig_B73xOh40B >= 3 & primary_nonadd_B73xOh40B), "Tier 2",
  (n_kos_sig_B73xA554 >= 1 & primary_nonadd_B73xA554) |
    (n_kos_sig_B73xOh40B >= 1 & primary_nonadd_B73xOh40B), "Tier 3",
  default = "Background"
)]
wide[, total_sig_support := rowSums(.SD, na.rm = TRUE),
     .SDcols = c("n_kos_sig_B73xA554", "n_kos_sig_B73xOh40B")]
wide[, total_leading_support := rowSums(.SD, na.rm = TRUE),
     .SDcols = c("n_primary_leading_edges_B73xA554", "n_primary_leading_edges_B73xOh40B")]
wide[, tier_order := match(Tier, c("Tier 1", "Tier 2", "Tier 3", "Background"))]
setorder(wide, tier_order, -total_sig_support, -total_leading_support)
wide[, tier_order := NULL]
fwrite(wide, file.path(result_dir, "05_ros_candidate_evidence_tiers.csv"))

message("Analyzing double KO")
double_dir <- file.path(analysis_dir, "input", "double_ko_wuzh02")
double_files <- list.files(double_dir, pattern = "diffRegulation\\.csv$", full.names = TRUE)
stopifnot(length(double_files) == 2L)
double_results <- rbindlist(lapply(double_files, function(f) {
  nm <- basename(f)
  hybrid <- sub("__.*$", "", nm)
  x <- fread(f)[!gene %chin% c("Zm00001d007234", "Zm00001d043795")]
  x[, `:=`(
    Hybrid = hybrid,
    KO = "Zm00001d007234+Zm00001d043795",
    Module = "ROS double KO",
    Source = "wuzh02",
    Runner = ifelse(hybrid == "B73xA554", "standard_ma2_mismatch", "standard_ma2")
  )]
  x
}), use.names = TRUE)

double_enrichment <- rbindlist(lapply(split(double_results, double_results$Hybrid), function(x) {
  analyze_one_ko(x, "ROS double KO", run_rank = TRUE)
}))
double_enrichment[, fisher_fdr_global := p.adjust(fisher_p, method = "BH")]
double_enrichment[, fisher_fdr_family := p.adjust(fisher_p, method = "BH"), by = .(analysis_family, set_type)]
double_enrichment[, wilcox_fdr_family := p.adjust(wilcox_p, method = "BH"), by = .(analysis_family, set_type)]
double_enrichment[, fgsea_fdr_family := p.adjust(fgsea_p, method = "BH"), by = .(analysis_family, set_type)]
double_enrichment[, interpretation := ifelse(Hybrid == "B73xOh40B", "comparable_runner", "exploratory_runner_mismatch")]
fwrite(double_enrichment, file.path(result_dir, "06_double_ko_nonadditive_enrichment.csv"))

double_gene <- rbindlist(lapply(c("B73xA554", "B73xOh40B"), function(hybrid) {
  d <- copy(double_results[Hybrid == hybrid])
  s1 <- copy(ros[Hybrid == hybrid & KO == "Zm00001d007234"])
  s2 <- copy(ros[Hybrid == hybrid & KO == "Zm00001d043795"])
  for (x in list(d, s1, s2)) x[, pct := frank(Z, ties.method = "average") / .N]
  d <- d[, .(gene, double_Z = Z, double_padj = p.adj, double_pct = pct)]
  s1 <- s1[, .(gene, ko007234_Z = Z, ko007234_padj = p.adj, ko007234_pct = pct)]
  s2 <- s2[, .(gene, ko043795_Z = Z, ko043795_padj = p.adj, ko043795_pct = pct)]
  z <- Reduce(function(a, b) merge(a, b, by = "gene", all = FALSE), list(d, s1, s2))
  z[, max_single_pct := pmax(ko007234_pct, ko043795_pct)]
  z[, delta_double_vs_max_single := double_pct - max_single_pct]
  z[, class := fcase(
    double_padj < 0.05 & ko007234_padj >= 0.05 & ko043795_padj >= 0.05, "double_specific",
    double_padj < 0.05 & ko007234_padj < 0.05 & ko043795_padj < 0.05, "shared_all_three",
    double_padj < 0.05 & (ko007234_padj < 0.05 | ko043795_padj < 0.05), "double_plus_one_single",
    double_padj >= 0.05 & (ko007234_padj < 0.05 | ko043795_padj < 0.05), "single_only",
    default = "not_significant"
  )]
  z[, `:=`(
    Hybrid = hybrid,
    interpretation = ifelse(hybrid == "B73xOh40B", "comparable_runner", "exploratory_runner_mismatch")
  )]
  z
}))
fwrite(double_gene, file.path(result_dir, "07_double_vs_single_gene_comparison.csv"))

message("Creating module comparison tables")
per_ko_module <- enrichment[analysis_family == "primary" & set_type == "all", .(
  median_log2_or = median(log2(odds_ratio), na.rm = TRUE),
  median_fold_enrichment = median(fold_enrichment, na.rm = TRUE),
  median_rank_auc = median(rank_auc, na.rm = TRUE),
  median_rank_biserial = median(rank_biserial, na.rm = TRUE),
  min_fisher_fdr = min(fisher_fdr_family, na.rm = TRUE),
  min_wilcox_fdr = min(wilcox_fdr_family, na.rm = TRUE),
  min_fgsea_fdr = if (all(is.na(fgsea_fdr_family))) NA_real_ else min(fgsea_fdr_family, na.rm = TRUE),
  n_primary_strata = .N
), by = .(Hybrid, Module, KO)]
module_summary <- per_ko_module[, .(
  n_KO = uniqueN(KO),
  median_log2_or = median(median_log2_or, na.rm = TRUE),
  range_log2_or_low = min(median_log2_or, na.rm = TRUE),
  range_log2_or_high = max(median_log2_or, na.rm = TRUE),
  median_rank_auc = median(median_rank_auc, na.rm = TRUE),
  median_rank_biserial = median(median_rank_biserial, na.rm = TRUE),
  range_rank_biserial_low = min(median_rank_biserial, na.rm = TRUE),
  range_rank_biserial_high = max(median_rank_biserial, na.rm = TRUE)
), by = .(Hybrid, Module)]
fwrite(per_ko_module, file.path(result_dir, "08_module_comparison_per_ko.csv"))
fwrite(module_summary, file.path(result_dir, "09_module_comparison_summary_exploratory.csv"))

message("Drawing figures")
ros_primary <- enrichment[Module == "ROS" & analysis_family == "primary" & set_type == "all"]
ros_primary[, KO_short := sub("Zm00001d", "d", KO)]

p1 <- ggplot(ros_primary, aes(x = KO_short, y = celltype, fill = NES, size = -log10(fgsea_fdr_family))) +
  geom_point(shape = 21, colour = "grey20", na.rm = TRUE) +
  facet_wrap(~Hybrid, scales = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
  scale_size_continuous(range = c(2, 7), limits = c(0, NA)) +
  labs(x = "ROS virtual KO", y = NULL, fill = "NES", size = "-log10 FDR",
       title = "ROS KO sensitivity of non-additive genes in primary tissues") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
ggsave(file.path(figure_dir, "01_ros_primary_nonadditive_NES.png"), p1, width = 10, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "01_ros_primary_nonadditive_NES.pdf"), p1, width = 10, height = 6)

direction_plot <- enrichment[
  Module == "ROS" & analysis_family == "primary" & set_type %chin% c("positive", "negative")
]
direction_plot[, label := paste(sub("Zm00001d", "d", KO), set_type, sep = " | ")]
p2 <- ggplot(direction_plot, aes(x = label, y = celltype, fill = NES)) +
  geom_tile(colour = "white") +
  facet_wrap(~Hybrid, scales = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
  labs(x = NULL, y = NULL, fill = "NES", title = "Positive and negative MPH gene-set sensitivity to ROS KO") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
ggsave(file.path(figure_dir, "02_ros_mph_direction_NES_heatmap.png"), p2, width = 14, height = 6.5, dpi = 300)
ggsave(file.path(figure_dir, "02_ros_mph_direction_NES_heatmap.pdf"), p2, width = 14, height = 6.5)

top_genes <- wide[Tier != "Background"][order(factor(Tier, levels = c("Tier 1", "Tier 2", "Tier 3")),
                                                -total_sig_support, -total_leading_support)][1:min(.N, 50), gene]
if (length(top_genes) >= 2) {
  hm <- ros_ranked[gene %chin% top_genes, .(gene, column = paste(Hybrid, sub("Zm00001d", "d", KO), sep = " | "), z_percentile)]
  mat <- dcast(hm, gene ~ column, value.var = "z_percentile")
  rn <- mat$gene
  mat <- as.matrix(mat[, -"gene"])
  rownames(mat) <- rn
  mat <- mat[match(top_genes, rownames(mat)), , drop = FALSE]
  ht <- Heatmap(
    mat, name = "Z percentile",
    col = colorRamp2(c(0, 0.5, 1), c("#F7FBFF", "#6BAED6", "#08306B")),
    cluster_rows = FALSE, cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 8),
    column_names_rot = 60,
    column_split = sub(" \\|.*$", "", colnames(mat)),
    column_title = "ROS KO network-perturbation percentile"
  )
  png(file.path(figure_dir, "03_ros_core_candidate_heatmap.png"), width = 2600, height = 2600, res = 300)
  draw(ht)
  dev.off()
  pdf(file.path(figure_dir, "03_ros_core_candidate_heatmap.pdf"), width = 10, height = 10)
  draw(ht)
  dev.off()
}

p4 <- ggplot(per_ko_module, aes(x = Module, y = median_rank_biserial, colour = Module)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, colour = "black") +
  facet_wrap(~Hybrid) +
  labs(x = NULL, y = "Per-KO median rank-biserial effect across primary tissue strata",
       title = "Exploratory standardized comparison across KO modules",
       subtitle = "Different KO batches; individual KO values are shown and module ranking is descriptive") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
ggsave(file.path(figure_dir, "04_exploratory_module_NES_comparison.png"), p4, width = 11, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "04_exploratory_module_NES_comparison.pdf"), p4, width = 11, height = 6)

p5 <- ggplot(double_gene, aes(x = max_single_pct, y = double_pct, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.55, size = 1.2) +
  facet_wrap(~Hybrid) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Maximum percentile of either single KO", y = "Double-KO percentile",
       colour = "Significance class", title = "Double KO versus matched single-KO perturbation ranks") +
  theme_bw(base_size = 10)
ggsave(file.path(figure_dir, "05_double_vs_single_rank_comparison.png"), p5, width = 10, height = 5.5, dpi = 300)
ggsave(file.path(figure_dir, "05_double_vs_single_rank_comparison.pdf"), p5, width = 10, height = 5.5)

capture.output(sessionInfo(), file = file.path(log_dir, "01_core_analysis_sessionInfo.txt"))
message("Core analysis completed")
