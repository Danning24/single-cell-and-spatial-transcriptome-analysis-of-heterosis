# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
}

project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(project_root, "ros_nonadditive_analysis")
result_dir <- file.path(analysis_dir, "results")
log_dir <- file.path(analysis_dir, "logs")
spatial_file <- file.path(
  project_root, "2026TUM_Danning_Wang_To_HK", "To HK",
  "all_spatial_merge_filtered_integrated_remove_noise_bin50_annotation.RDS"
)
stopifnot(file.exists(spatial_file))

message("Reading spatial Seurat object: ", spatial_file)
obj <- readRDS(spatial_file)
message("Object loaded")

md <- obj[[]]
audit <- data.table(
  item = c("class", "seurat_version", "n_features_default_assay", "n_bins", "default_assay",
           "assays", "reductions", "images", "metadata_columns"),
  value = c(
    paste(class(obj), collapse = ";"),
    as.character(packageVersion("Seurat")),
    nrow(obj),
    ncol(obj),
    DefaultAssay(obj),
    paste(Assays(obj), collapse = ";"),
    paste(Reductions(obj), collapse = ";"),
    paste(Images(obj), collapse = ";"),
    paste(colnames(md), collapse = ";")
  )
)
fwrite(audit, file.path(result_dir, "10_spatial_object_audit.csv"))

sample_field <- intersect(c("sample", "orig.ident", "Sample", "slice", "Slice"), colnames(md))[1]
if (!is.na(sample_field)) {
  sample_counts <- as.data.table(table(sample = md[[sample_field]], useNA = "ifany"))
  setnames(sample_counts, "N", "n_bins")
  sample_counts[, sample_field := sample_field]
  fwrite(sample_counts, file.path(result_dir, "11_spatial_sample_counts.csv"))
}

ident_counts <- as.data.table(table(annotation = as.character(Idents(obj)), useNA = "ifany"))
setnames(ident_counts, "N", "n_bins")
fwrite(ident_counts, file.path(result_dir, "12_spatial_identity_counts.csv"))

categorical <- rbindlist(lapply(colnames(md), function(nm) {
  x <- md[[nm]]
  nu <- uniqueN(x, na.rm = FALSE)
  if (nu > 100L) return(NULL)
  z <- as.data.table(table(value = as.character(x), useNA = "ifany"))
  setnames(z, "N", "n_bins")
  z[, field := nm]
  z
}), fill = TRUE)
fwrite(categorical, file.path(result_dir, "13_spatial_categorical_metadata_summary.csv"))

feature_summary <- rbindlist(lapply(Assays(obj), function(a) {
  data.table(assay = a, n_features = nrow(obj[[a]]), n_bins = ncol(obj[[a]]))
}))
fwrite(feature_summary, file.path(result_dir, "14_spatial_assay_summary.csv"))

capture.output({
  cat("Spatial file:", spatial_file, "\n")
  print(obj)
  cat("\nDefault assay:", DefaultAssay(obj), "\n")
  cat("\nAssays:\n"); print(Assays(obj))
  cat("\nImages:\n"); print(Images(obj))
  cat("\nMetadata columns:\n"); print(colnames(md))
  cat("\nIdentity counts:\n"); print(sort(table(Idents(obj)), decreasing = TRUE))
  cat("\nSample field:", sample_field, "\n")
  if (!is.na(sample_field)) print(sort(table(md[[sample_field]]), decreasing = TRUE))
  cat("\nSession info:\n"); print(sessionInfo())
}, file = file.path(log_dir, "02_spatial_audit.txt"))

message("Spatial audit completed")
