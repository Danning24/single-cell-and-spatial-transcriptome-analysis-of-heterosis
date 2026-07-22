# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

suppressPackageStartupMessages({
  library(ggplot2)
})

base_dir <- file.path(project_root, "sctenifold_results")
summary_dir <- file.path(base_dir, "_summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

count_files <- list.files(base_dir, pattern = "sig_gene_counts\\.csv$", recursive = TRUE, full.names = TRUE)
count_files <- count_files[grepl("B73x.*__Zm", count_files)]

read_one <- function(f) {
  x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  folder <- basename(dirname(f))
  parts <- strsplit(folder, "__", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    x$group <- parts[1]
    x$KO <- parts[2]
  }
  x$source_file <- f
  x
}

rows <- lapply(count_files, read_one)
all_cols <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(x) {
  missing_cols <- setdiff(all_cols, names(x))
  for (cc in missing_cols) x[[cc]] <- NA
  x[, all_cols, drop = FALSE]
})
count_df <- do.call(rbind, rows)
keep_cols <- intersect(c("group", "KO", "n_sig", "runner", "source_file"), names(count_df))
count_df <- count_df[, keep_cols, drop = FALSE]
count_df$n_sig <- as.integer(count_df$n_sig)
count_df <- count_df[order(count_df$group, count_df$KO), , drop = FALSE]

write.csv(count_df, file.path(summary_dir, "perturbed_gene_counts_by_KO_group.csv"), row.names = FALSE)

# 按两个组的总扰动数量排序，方便横向比较。
ord <- aggregate(n_sig ~ KO, count_df, sum, na.rm = TRUE)
ord <- ord[order(ord$n_sig, decreasing = TRUE), ]
count_df$KO <- factor(count_df$KO, levels = ord$KO)

p <- ggplot(count_df, aes(x = KO, y = n_sig, fill = group)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  geom_text(
    aes(label = n_sig),
    position = position_dodge(width = 0.76),
    hjust = -0.15,
    size = 3
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("B73xA554" = "#4F8FC7", "B73xOh40B" = "#F29B52")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  labs(
    title = "Significant perturbed genes per virtual knockout",
    x = "KO gene",
    y = "Number of significant perturbed genes",
    fill = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    legend.position = "top",
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 30, 10, 10)
  )

ggsave(file.path(summary_dir, "perturbed_gene_counts_by_KO_group.pdf"), p, width = 9, height = 9)
ggsave(file.path(summary_dir, "perturbed_gene_counts_by_KO_group.png"), p, width = 9, height = 9, dpi = 300)

# 再输出一个每组总数的小表，方便直接看总体量。
group_total <- aggregate(n_sig ~ group, count_df, sum, na.rm = TRUE)
write.csv(group_total, file.path(summary_dir, "perturbed_gene_counts_group_total.csv"), row.names = FALSE)

p2 <- ggplot(group_total, aes(x = group, y = n_sig, fill = group)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = n_sig), vjust = -0.35, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("B73xA554" = "#4F8FC7", "B73xOh40B" = "#F29B52")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Total significant perturbed genes across virtual knockouts",
    x = NULL,
    y = "Total count"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(summary_dir, "perturbed_gene_counts_group_total.pdf"), p2, width = 5.5, height = 4.5)
ggsave(file.path(summary_dir, "perturbed_gene_counts_group_total.png"), p2, width = 5.5, height = 4.5, dpi = 300)

cat("N_ROWS\t", nrow(count_df), "\n", sep = "")
cat("OUTPUT\t", summary_dir, "\n", sep = "")
