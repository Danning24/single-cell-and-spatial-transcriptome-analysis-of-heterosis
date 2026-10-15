#!/usr/bin/env Rscript
# Figure 3D: Monocle3 trajectory on the original cortex UMAP.

library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)

# Inputs and outputs
project_dir <- Sys.getenv("MAIZE_PROJECT_DIR")
cortex_file <- Sys.getenv("MAIZE_CORTEX_RDS",
                          unset = file.path(project_dir, "4.Cortex_cells/data_output/Cortex.rds"))
root_file <- Sys.getenv("MAIZE_MONOCLE_ROOT_FILE", unset = "")
output_dir <- file.path(Sys.getenv("MAIZE_RESULTS_DIR", unset = "results"), "monocle3_trajectory")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Preserve the shared Seurat embedding for graph learning.
cortex <- readRDS(cortex_file)
cds <- as.cell_data_set(cortex)
SingleCellExperiment::reducedDims(cds)$UMAP <- Embeddings(cortex, "umap")
cds <- cluster_cells(cds, resolution = 0.001)
cds <- learn_graph(cds)

# The original root was chosen interactively in the lower, high-M_Score region.
# Supply root-node IDs only when reusing the same learned graph.
if (nzchar(root_file)) {
  cds <- order_cells(cds, root_pr_nodes = readLines(root_file))
} else {
  cds <- order_cells(cds)
}
pseudotime_plot <- plot_cells(cds, color_cells_by = "pseudotime",
                              label_groups_by_cluster = FALSE,
                              label_leaves = FALSE, label_branch_points = FALSE)
ggsave(file.path(output_dir, "cortex_pseudotime.pdf"), pseudotime_plot,
       width = 6, height = 5)
writeLines(cds@principal_graph_aux[["UMAP"]]$root_pr_nodes,
           file.path(output_dir, "root_pr_nodes.txt"))
write.csv(data.frame(cell_id = colnames(cds), pseudotime = pseudotime(cds)),
          file.path(output_dir, "cortex_pseudotime.csv"), row.names = FALSE)
saveRDS(cds, file.path(output_dir, "cortex_trajectory.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
