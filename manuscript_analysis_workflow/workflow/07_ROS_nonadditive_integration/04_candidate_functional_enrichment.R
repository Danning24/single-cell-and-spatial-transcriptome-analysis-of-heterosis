# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(ggplot2)
})

if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
}

project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(project_root, "ros_nonadditive_analysis")
result_dir <- file.path(analysis_dir, "results")
figure_dir <- file.path(analysis_dir, "figures")
annotation_dir <- file.path(
  project_root, "sctenifold_downstream_analysis", "maize_GO_annotation_maizeGAMER"
)

tiers <- fread(file.path(result_dir, "05_ros_candidate_evidence_tiers.csv"))
primary <- fread(file.path(result_dir, "04_ros_primary_tissue_enrichment.csv"))
nonadd <- fread(file.path(
  project_root, "数据汇总", "汇总文件",
  "TbaleS2 sc_all_heterotic_genes_10pct_cells_pattern_cortex.csv"
))
go_map <- fread(file.path(annotation_dir, "maizeGAMER_GO_gene_map_used.csv"))
kegg_map <- fread(file.path(annotation_dir, "maize_KEGG_gene_pathway_map_used.csv"))

split_genes <- function(x) {
  unique(trimws(unlist(strsplit(x[!is.na(x) & nzchar(x)], ";", fixed = TRUE))))
}

recurrent_le <- primary[
  Hybrid == "B73xA554" & set_type == "all",
  .(gene = split_genes(leading_edge)),
  by = KO
][, .(n_ros_kos = uniqueN(KO)), by = gene][n_ros_kos >= 3, gene]

gene_sets <- list(
  ROS_Tier1 = tiers[Tier == "Tier 1", gene],
  ROS_Tier1_Tier2_consensus = tiers[Tier %chin% c("Tier 1", "Tier 2"), gene],
  ROS_sensitive_positive_MPH = tiers[
    Tier %chin% c("Tier 1", "Tier 2") & n_positive_nonadd_B73xA554 > 0, gene
  ],
  ROS_sensitive_negative_MPH = tiers[
    Tier %chin% c("Tier 1", "Tier 2") & n_negative_nonadd_B73xA554 > 0, gene
  ],
  ROS_Epidermis_recurrent_LE = recurrent_le
)

gene_col <- intersect(c("gene", "Gene", "genes", "Gene_ID", "gene_id"), names(nonadd))[1]
if (is.na(gene_col) && "Trait" %in% names(nonadd) &&
    mean(grepl("^Zm00001d", nonadd$Trait)) > 0.9) gene_col <- "Trait"
stopifnot(!is.na(gene_col))
setnames(nonadd, gene_col, "gene")
epidermis_background <- unique(nonadd[
  Hybrid == "B73xA554" & pattern == "Non-additive" & cluster_anno == "Epidermis",
  gene
])

backgrounds <- lapply(names(gene_sets), function(nm) {
  if (nm == "ROS_Epidermis_recurrent_LE") epidermis_background else tiers$gene
})
names(backgrounds) <- names(gene_sets)

go_term2gene <- unique(go_map[, .(term = GOID, gene)])
go_term2name <- unique(go_map[, .(term = GOID, name = Description)])
kegg_term2gene <- unique(kegg_map[, .(term = TERM_ID, gene)])
kegg_term2name <- unique(kegg_map[, .(term = TERM_ID, name = Description)])

run_custom_enrichment <- function(genes, universe, term2gene, term2name) {
  genes <- unique(intersect(genes, universe))
  x <- enricher(
    gene = genes,
    universe = unique(universe),
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500
  )
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(data.table())
  as.data.table(as.data.frame(x))
}

go_results <- rbindlist(lapply(names(gene_sets), function(nm) {
  rbindlist(lapply(c("BP", "CC", "MF"), function(ont) {
    map_ont <- go_map[Ontology == ont]
    ans <- run_custom_enrichment(
      gene_sets[[nm]], backgrounds[[nm]],
      unique(map_ont[, .(term = GOID, gene)]),
      unique(map_ont[, .(term = GOID, name = Description)])
    )
    if (nrow(ans)) ans[, `:=`(module = nm, ontology = ont)]
    ans
  }), fill = TRUE)
}), fill = TRUE)

kegg_results <- rbindlist(lapply(names(gene_sets), function(nm) {
  ans <- run_custom_enrichment(
    gene_sets[[nm]], backgrounds[[nm]], kegg_term2gene, kegg_term2name
  )
  if (nrow(ans)) ans[, module := nm]
  ans
}), fill = TRUE)

if (nrow(go_results)) {
  setcolorder(go_results, c("module", "ontology", setdiff(names(go_results), c("module", "ontology"))))
  setorder(go_results, module, ontology, p.adjust)
}
if (nrow(kegg_results)) {
  setcolorder(kegg_results, c("module", setdiff(names(kegg_results), "module")))
  setorder(kegg_results, module, p.adjust)
}
fwrite(go_results, file.path(result_dir, "23_candidate_GO_enrichment.csv"))
fwrite(kegg_results, file.path(result_dir, "24_candidate_KEGG_enrichment.csv"))

coverage <- rbindlist(lapply(names(gene_sets), function(nm) {
  genes <- unique(gene_sets[[nm]])
  universe <- unique(backgrounds[[nm]])
  data.table(
    module = nm,
    n_genes = length(genes),
    n_background = length(universe),
    n_GO_annotated = sum(genes %chin% go_term2gene$gene),
    n_KEGG_annotated = sum(genes %chin% kegg_term2gene$gene),
    background_definition = ifelse(
      nm == "ROS_Epidermis_recurrent_LE",
      "B73xA554 Epidermis non-additive genes",
      "all genes evaluated in ROS candidate tiering"
    )
  )
}))
fwrite(coverage, file.path(result_dir, "25_functional_enrichment_coverage.csv"))

plot_enrichment <- function(dt, filename, facet_var) {
  sig <- dt[p.adjust < 0.05]
  if (!nrow(sig)) return(invisible(NULL))
  top <- sig[order(p.adjust), head(.SD, 10), by = c("module", facet_var)]
  top[, label := factor(Description, levels = rev(unique(Description)))]
  p <- ggplot(top, aes(x = Count, y = label, color = -log10(p.adjust), size = Count)) +
    geom_point() +
    facet_grid(rows = vars(module), cols = vars(!!rlang::sym(facet_var)), scales = "free_y", space = "free_y") +
    scale_color_viridis_c(option = "C", direction = -1) +
    labs(x = "Candidate genes", y = NULL, color = expression(-log[10](FDR)), size = "Count") +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0), panel.grid.minor = element_blank())
  ggsave(file.path(figure_dir, paste0(filename, ".png")), p, width = 14, height = 12, dpi = 300)
  ggsave(file.path(figure_dir, paste0(filename, ".pdf")), p, width = 14, height = 12)
}

if (nrow(go_results)) plot_enrichment(go_results, "10_candidate_GO_enrichment", "ontology")
if (nrow(kegg_results)) {
  kegg_results[, category := "KEGG"]
  plot_enrichment(kegg_results, "11_candidate_KEGG_enrichment", "category")
}

message("Candidate GO/KEGG enrichment completed")
