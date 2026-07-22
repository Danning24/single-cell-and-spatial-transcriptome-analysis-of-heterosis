# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

rm(list=ls());gc()
library(Seurat)
library(SeuratDisk)

setwd(project_root)
Cortex <- readRDS(file.path(project_root, "results", "Cortex.rds"))
Cortex$group_FP <- ifelse(grepl("x", Cortex$orig.ident),
                          "F1", ifelse(Cortex$orig.ident == "B73",
                                       "B73",
                                       "P1"))
Cortex$celltype<-paste0(Cortex$group_FP,"_",Cortex$group)
unique(Cortex$group_mmd)
group<-list(c("F1_SS","P1_SS"),c("F1_SS","B73_SS"),
     c("F1_TL","P1_TL"),c("F1_TL","B73_SS"),
     c("F1_NSS","P1_NSS"),c("F1_NSS","B73_SS"))

out_dir <- "./data_output/cortex_h5ad"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


library(SingleCellExperiment)
library(zellkonverter)

for (pair in group) {
  obj_sub <- subset(Cortex, subset = celltype %in% pair)
  obj_sub$group_mmd <- paste0(pair[1], "_vs_", pair[2])

  # Seurat -> SingleCellExperiment
  sce <- as.SingleCellExperiment(obj_sub)

  tag <- gsub("[^A-Za-z0-9_\\-]", "_", paste0(pair[1], "_vs_", pair[2]))
  h5ad_path <- file.path(out_dir, paste0(tag, ".h5ad"))

  zellkonverter::writeH5AD(sce, file = h5ad_path)
  message("Saved: ", h5ad_path)

  rm(obj_sub, sce); gc()
}
