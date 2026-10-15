# fgsea implements preranked enrichment; clusterProfiler implements GO ORA.
# Fixed gene sets, the original GO mapping, and separate FDR families are retained.
library(data.table)
library(fgsea)
library(clusterProfiler)
library(ggplot2)

# Inputs and outputs -----------------------------------------------------------
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
analysis_dir <- file.path(project_dir, "7.sctenifold_KO", "analysis_outs")
manifest <- fread(Sys.getenv("MAIZE_KO_MANIFEST", file.path(
  analysis_dir, "tables", "00_KO_input_manifest.csv"
)))
targets <- fread(Sys.getenv("MAIZE_KO_TARGETS", "inputs/knockout_targets.csv"))
cn <- targets[module == "C-N lipid and energy metabolism"]
manifest <- merge(manifest[, .(Hybrid, KO, file)], cn,
                  by.x = "KO", by.y = "gene", sort = FALSE)
nag <- unique(fread(Sys.getenv("MAIZE_NAG_TABLE"), select = c(
  "Hybrid", "Trait", "cluster_anno", "pattern", "MPH.pct"
)))
go_file <- Sys.getenv("MAIZE_GO_MAPPING", file.path(
  project_dir, "7.sctenifold_KO", "sctenifold20260527_ROS_rel",
  "sctenifold_downstream_analysis", "maize_GO_annotation_maizeGAMER",
  "maizeGAMER_GO_gene_map_used.csv"
))
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", "results"), "knockout_enrichment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
contexts <- data.table(
  Hybrid = c("B73xA554", rep("B73xOh40B", 4)),
  state = c("Epidermis", "Cortex outer", "Cortex II mature",
            "Cortex I mature", "Cortex early")
)

# Preranked GSEA: 8 candidates x 5 contexts x 3 gene-set categories ---------------
# The historical Multilevel run did not record an RNG seed. An optional seed
# makes a new run reproducible; use the archived results for exact original P values.
gsea_seed <- Sys.getenv("MAIZE_GSEA_SEED")
if (nzchar(gsea_seed)) set.seed(as.integer(gsea_seed))
gsea <- rbindlist(lapply(cn$gene, function(ko) {
  rbindlist(lapply(seq_len(nrow(contexts)), function(j) {
    hybrid <- contexts$Hybrid[j]
    state <- contexts$state[j]
    dr <- fread(manifest[Hybrid == hybrid & KO == ko, file])
    dr <- dr[gene != ko & is.finite(Z)]
    annotation <- nag[Hybrid == hybrid & cluster_anno == state]
    x <- merge(dr, annotation, by.x = "gene", by.y = "Trait", sort = FALSE)
    ranks <- sort(setNames(x$Z, x$gene), decreasing = TRUE)
    gene_sets <- list(
      all = x[pattern == "Non-additive", gene],
      positive = x[pattern == "Non-additive" & MPH.pct > 0, gene],
      negative = x[pattern == "Non-additive" & MPH.pct < 0, gene]
    )
    gene_sets <- gene_sets[lengths(gene_sets) >= 10]
    enrichment <- as.data.table(fgseaMultilevel(
      pathways = gene_sets, stats = ranks, minSize = 10,
      maxSize = length(ranks), eps = 1e-10, nproc = 1
    ))
    enrichment[, leading_edge := vapply(leadingEdge, paste, collapse = ";",
                                        FUN.VALUE = character(1))]
    out <- merge(data.table(category = c("all", "positive", "negative")),
                 enrichment[, .(category = pathway, NES, P = pval,
                                 gene_set_size = size, leading_edge)], all.x = TRUE)
    out[, `:=`(Hybrid = hybrid, KO = ko, state = state)]
    out
  }))
}))
# Each family contains 40 planned tests; unavailable results remain NA.
gsea[, FDR := p.adjust(P, "BH"), by = category]
fwrite(gsea, file.path(output_dir, "GSEA_120_tests.csv"))

# GO ORA: 8 candidates x 2 hybrids ---------------------------------------------
annotation <- unique(fread(go_file)[Ontology %in% c("BP", "CC", "MF")])
term2gene <- unique(annotation[, .(GOID, gene)])
term2name <- unique(annotation[, .(GOID, Description)])
ontology <- unique(annotation[, .(ID = GOID, Ontology)])
go_inputs <- list()
go <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  ko <- manifest$KO[i]
  dr <- fread(manifest$file[i])
  universe <- union(dr[is.finite(p.adj), unique(gene)], ko)
  query <- union(dr[is.finite(p.adj) & p.adj < 0.05 & gene != ko, unique(gene)], ko)
  go_inputs[[i]] <<- data.table(Hybrid = manifest$Hybrid[i], KO = ko, gene = query)
  # enricher restricts query and universe to genes present in TERM2GENE.
  # Test BP/CC/MF together. There is no 'ont' argument in enricher.
  enrichment <- enricher(
    gene = query, universe = universe, TERM2GENE = term2gene, TERM2NAME = term2name,
    minGSSize = 5, maxGSSize = length(universe), pAdjustMethod = "BH",
    pvalueCutoff = 1, qvalueCutoff = 1
  )
  out <- merge(as.data.table(enrichment@result), ontology, by = "ID", all.x = TRUE)
  out[, `:=`(Hybrid = manifest$Hybrid[i], KO = ko)]
  out
}), fill = TRUE)
setnames(go, "p.adjust", "FDR_KO_all_ontologies")
go[, significant := Count >= 5 & FDR_KO_all_ontologies < 0.05]
fwrite(go, file.path(output_dir, "GO_all_tests.csv"))
fwrite(go[significant == TRUE], file.path(output_dir, "GO_Count5_FDR05.csv"))
fwrite(rbindlist(go_inputs), file.path(output_dir, "GO_query_genes.csv"))
fwrite(data.table(file = go_file, md5 = unname(tools::md5sum(go_file))),
       file.path(output_dir, "GO_annotation_provenance.csv"))

# Fig. 4c and S18b -------------------------------------------------------------
gsea[, context := paste(Hybrid, state, sep = ": ")]
gsea[, context := factor(context, paste(contexts$Hybrid, contexts$state, sep = ": "))]
gsea[, KO := factor(KO, cn$gene)]
p_gsea <- ggplot(gsea[is.finite(NES)], aes(KO, category)) +
  geom_point(aes(fill = NES, size = pmin(-log10(pmax(FDR, 1e-300)), 4)), shape = 21) +
  geom_text(data = gsea[FDR < 0.05], label = "*", size = 3) +
  facet_grid(context ~ .) + scale_x_discrete(labels = setNames(cn$symbol, cn$gene)) +
  scale_fill_gradient2(low = "#3265B4", mid = "white", high = "#ED4646") +
  theme_bw() + labs(x = NULL, y = NULL, size = "-log10 FDR", fill = "NES")
# Preserve the figure's eight interpretable terms; retain GO:0009913 in full results.
shown_terms <- c("GO:0048767", "GO:0010075", "GO:2000026", "GO:0009734",
                 "GO:0048646", "GO:0071555", "GO:0009627", "GO:0005856")
go_plot <- go[Hybrid == "B73xOh40B" & ID %in% shown_terms & Count >= 5]
go_plot[, KO := factor(KO, cn$gene)]
p_go <- ggplot(go_plot, aes(KO, Description)) +
  geom_point(aes(size = Count, color = -log10(FDR_KO_all_ontologies),
                 shape = significant)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  scale_x_discrete(labels = setNames(cn$symbol, cn$gene)) + theme_bw() +
  labs(x = NULL, y = NULL, color = "-log10 FDR") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(output_dir, "S18b_GSEA.pdf"), p_gsea, width = 8, height = 7)
ggsave(file.path(output_dir, "Fig4c_GO.pdf"), p_go, width = 8, height = 4)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
