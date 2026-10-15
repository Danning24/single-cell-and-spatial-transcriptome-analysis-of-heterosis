# data.table handles the original result tables; stats supplies Fisher/Wilcoxon
# tests and BH adjustment; ggplot2 draws the manuscript summaries.
library(data.table)
library(ggplot2)

# Inputs and outputs -----------------------------------------------------------
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
analysis_dir <- file.path(project_dir, "7.sctenifold_KO", "analysis_outs")
manifest <- fread(Sys.getenv("MAIZE_KO_MANIFEST", file.path(
  analysis_dir, "tables", "00_KO_input_manifest.csv"
)))
targets <- fread(Sys.getenv("MAIZE_KO_TARGETS", "inputs/knockout_targets.csv"))
nag <- unique(fread(Sys.getenv("MAIZE_NAG_TABLE"), select = c(
  "Hybrid", "Trait", "cluster_anno", "pattern", "MPH.pct"
)))
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", "results"), "knockout_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
hybrids <- c("B73xA554", "B73xOh40B")
states <- c("Epidermis", "Cortex outer", "Cortex II mature",
            "Cortex I mature", "Cortex early")
modules <- c("C-N lipid and energy metabolism",
             "Stress transport and protein homeostasis",
             "Hormone metabolism and signaling",
             "ROS and redox homeostasis",
             "Cell wall water transport and cytoskeleton")
manifest <- merge(manifest[, .(Hybrid, KO, file)], targets,
                  by.x = "KO", by.y = "gene", sort = FALSE)
results <- lapply(manifest$file, fread)
names(results) <- paste(manifest$Hybrid, manifest$KO, sep = "__")

# Top-decile NAG enrichment: 36 targets x 2 hybrids x 5 cell states --------------
metrics <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  hybrid <- manifest$Hybrid[i]
  ko <- manifest$KO[i]
  dr <- results[[i]][gene != ko & is.finite(Z)]
  rbindlist(lapply(states, function(state) {
    annotation <- nag[Hybrid == hybrid & cluster_anno == state]
    x <- merge(dr, annotation, by.x = "gene", by.y = "Trait")
    x <- x[pattern %in% c("Additive", "Non-additive")]
    top <- frank(x$Z, ties.method = "average") / nrow(x) >= 0.90
    is_nag <- x$pattern == "Non-additive"
    a <- sum(top & is_nag)
    b <- sum(top & !is_nag)
    N_nag <- sum(is_nag)
    N_additive <- sum(!is_nag)
    fisher <- fisher.test(matrix(c(a, b, N_nag - a, N_additive - b), 2),
                          alternative = "greater")
    data.table(
      Hybrid = hybrid, KO = ko, state, module = manifest$module[i],
      N_nag, N_additive, top_nag = a, top_additive = b,
      log2_FE = log2(((a + 0.5) / (N_nag + 1)) /
                      ((b + 0.5) / (N_additive + 1))),
      fisher_P = fisher$p.value
    )
  }))
}))
metrics[, fisher_FDR := p.adjust(fisher_P, "BH"), by = .(Hybrid, state)]
fwrite(metrics, file.path(output_dir, "cell_state_enrichment.csv"))

# Cortex is the median of four state-specific values, not a pooled-cell score.
tissue <- metrics[, .(log2_FE = median(log2_FE)),
                  by = .(Hybrid, KO, module,
                         Tissue = fifelse(state == "Epidermis", "Epidermis", "Cortex"))]
preference <- dcast(tissue, Hybrid + KO + module ~ Tissue, value.var = "log2_FE")
preference[, D := Epidermis - Cortex]
fwrite(tissue, file.path(output_dir, "tissue_enrichment.csv"))
fwrite(preference, file.path(output_dir, "tissue_preference.csv"))

# Match the same knockout between hybrids; adjust the five module comparisons.
paired <- dcast(preference, KO + module ~ Hybrid, value.var = "D")
tests <- paired[, {
  w <- wilcox.test(B73xA554, B73xOh40B, paired = TRUE,
                   alternative = "two.sided", exact = TRUE)
  .(n_pairs = .N, V = unname(w$statistic), P = w$p.value,
    median_A554 = median(B73xA554), median_Oh40B = median(B73xOh40B))
}, by = module]
tests[, FDR := p.adjust(P, "BH")]
fwrite(tests, file.path(output_dir, "paired_module_tests.csv"))

# Shared effects: each knockout keeps its original tested-gene background. ------
cn <- targets[module == "C-N lipid and energy metabolism", gene]
overlap <- rbindlist(lapply(hybrids, function(hybrid) {
  sets <- setNames(lapply(cn, function(ko) {
    results[[paste(hybrid, ko, sep = "__")]][
      gene != ko & is.finite(p.adj) & p.adj < 0.05, unique(gene)
    ]
  }), cn)
  pairs <- CJ(row = seq_along(cn), column = seq_along(cn))
  pairs[, c("n_row", "n_column", "shared", "union", "Jaccard") := {
    a <- sets[[row]]
    b <- sets[[column]]
    shared <- length(intersect(a, b))
    total <- length(union(a, b))
    list(length(a), length(b), shared, total, if (total) shared / total else 0)
  }, by = .(row, column)]
  pairs[, `:=`(Hybrid = hybrid, KO_row = cn[row], KO_column = cn[column])]
  pairs
}))
fwrite(overlap, file.path(output_dir, "within_hybrid_overlap.csv"))

# Figure S18a and Fig. 4a-b -----------------------------------------------------
tissue[, module := factor(module, modules)]
p_tissue <- ggplot(tissue, aes(module, log2_FE, color = module)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) + geom_point(size = 1.3) +
  facet_grid(Tissue ~ Hybrid) + theme_classic() +
  labs(x = NULL, y = "Top-decile NAG enrichment (log2 FE)") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")
preference[, module := factor(module, rev(modules))]
p_preference <- ggplot(preference, aes(D, module, color = Hybrid)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.65), width = 0.55) +
  geom_point(position = position_dodge(0.65), size = 1.3) +
  scale_color_manual(values = c("#377EB8", "#E5B52F")) + theme_classic() +
  labs(x = "Epidermis - Cortex (top-decile NAG log2 FE)", y = NULL)
overlap[, label := ifelse(row < column, as.character(shared),
                         ifelse(row == column, as.character(n_row), sprintf("%.2f", Jaccard)))]
symbol_labels <- setNames(targets$symbol, targets$gene)
p_overlap <- ggplot(overlap, aes(column, -row)) +
  geom_tile(fill = "white", color = "grey65") +
  geom_tile(data = overlap[row == column], fill = "grey85", color = "grey65") +
  geom_tile(data = overlap[row > column], aes(fill = Jaccard), color = "grey65") +
  geom_text(aes(label = label), size = 2.5) + facet_wrap(~Hybrid) +
  scale_fill_gradient(low = "white", high = "#0071A5", limits = c(0, 1)) +
  scale_x_continuous(breaks = seq_along(cn), labels = symbol_labels[cn], position = "top") +
  scale_y_continuous(breaks = -seq_along(cn), labels = symbol_labels[cn]) +
  coord_fixed() + theme_void() + theme(axis.text = element_text(size = 7))
ggsave(file.path(output_dir, "S18a_tissue_enrichment.pdf"), p_tissue, width = 8, height = 6)
ggsave(file.path(output_dir, "Fig4a_tissue_preference.pdf"), p_preference, width = 8, height = 3)
ggsave(file.path(output_dir, "Fig4b_overlap.pdf"), p_overlap, width = 8, height = 4)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
