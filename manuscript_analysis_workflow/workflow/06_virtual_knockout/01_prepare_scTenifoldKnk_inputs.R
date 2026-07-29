project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

scRNA_rds <- "scRNA_all_annotated_filtered_500gene_1000umi_CCA.RDS"

groups <- c("B73xA554", "B73xOh40B")

ko_genes <- c(
  "Zm00001d011969",
  "Zm00001d028307",
  "Zm00001d034109",
  "Zm00001d006753",
  "Zm00001d011208",
  "Zm00001d018973",
  "Zm00001d048021",
  "Zm00001d015852",
  "Zm00001d013410",
  "Zm00001d009858",
  "Zm00001d028392",
  "Zm00001d019565",
  "Zm00001d037547",
  "Zm00001d039221",
  "Zm00001d012221",
  "Zm00001d019283",
  "Zm00001d030725",
  "Zm00001d026577",
  "Zm00001d002275",
  "Zm00001d043551",
  "Zm00001d051804",
  "Zm00001d030557"
)

ko_genes <- unique(ko_genes)
cat("N_KO_GENES_TOTAL\t", length(ko_genes), "\n", sep = "")

outdir <- file.path("new_ko_genes_scripts", "sctenifold_inputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(Matrix)
})

if (!file.exists(scRNA_rds)) {
  stop("Cannot find scRNA RDS: ", scRNA_rds)
}

obj <- readRDS(scRNA_rds)

obj_attrs <- attributes(obj)
md <- obj_attrs[["meta.data"]]
assays <- obj_attrs[["assays"]]

if (is.null(md)) {
  stop("Cannot find meta.data in Seurat object")
}
if (!"orig.ident" %in% colnames(md)) {
  stop("metadata column orig.ident not found")
}
if (!"RNA" %in% names(assays)) {
  stop("RNA assay not found")
}

rna_attrs <- attributes(assays[["RNA"]])
rna_layers <- rna_attrs[["layers"]]

if (is.null(rna_layers) || !"counts" %in% names(rna_layers)) {
  stop("RNA counts layer not found")
}

mat <- rna_layers[["counts"]]

get_dimnames <- function(x) {
  dn <- attributes(x)[["Dimnames"]]
  if (is.null(dn)) dn <- attributes(x)[["dimnames"]]
  if (is.null(dn)) dn <- dimnames(x)
  dn
}

dims_of <- function(x) {
  d <- attributes(x)[["Dim"]]
  if (is.null(d)) d <- attributes(x)[["dim"]]
  if (is.null(d)) d <- dim(x)
  d
}

rna_feature_map <- rna_attrs[["features"]]
rna_cell_map <- rna_attrs[["cells"]]

rna_feature_names <- get_dimnames(rna_feature_map)[[1]]
rna_cell_names <- get_dimnames(rna_cell_map)[[1]]

if (is.null(rna_feature_names) || is.null(rna_cell_names)) {
  stop("RNA assay feature/cell names are missing")
}

dimnames(mat) <- list(rna_feature_names, rna_cell_names)
mat <- as(mat, "dgCMatrix")
cat("RNA_COUNTS_DIM\t", paste(dims_of(mat), collapse = "x"), "\n", sep = "")

feature_set <- character(0)

if ("SCT" %in% names(assays)) {
  sct_attrs <- attributes(assays[["SCT"]])

  if ("var.features" %in% names(sct_attrs)) {
    feature_set <- sct_attrs[["var.features"]]
  }

  for (slot_name in intersect(c("data", "counts", "scale.data"), names(sct_attrs))) {
    if (length(feature_set) > 0) break

    sct_dn <- get_dimnames(sct_attrs[[slot_name]])
    if (!is.null(sct_dn) && length(sct_dn[[1]]) > 0) {
      feature_set <- sct_dn[[1]]
      break
    }
  }
}

if (length(feature_set) == 0) {
  detected <- Matrix::rowSums(mat > 0)
  feature_set <- rna_feature_names[order(detected, decreasing = TRUE)][seq_len(min(5000, length(detected)))]
}

cat("N_NETWORK_FEATURES_BEFORE_ADD_KO\t", length(feature_set), "\n", sep = "")

ko_present_in_rna <- ko_genes[ko_genes %in% rna_feature_names]
ko_missing_in_rna <- setdiff(ko_genes, rna_feature_names)

cat("N_KO_PRESENT_IN_RNA\t", length(ko_present_in_rna), "\n", sep = "")
cat("N_KO_MISSING_IN_RNA\t", length(ko_missing_in_rna), "\n", sep = "")

if (length(ko_missing_in_rna) > 0) {
  cat("KO_MISSING_IN_RNA\t", paste(ko_missing_in_rna, collapse = ","), "\n", sep = "")
}

feature_set <- unique(c(feature_set, ko_present_in_rna))

row_idx <- match(feature_set, rna_feature_names)
row_idx <- row_idx[!is.na(row_idx)]
feature_set <- rna_feature_names[row_idx]

cat("N_FEATURES_FINAL\t", length(feature_set), "\n", sep = "")
cat("KO_IN_FINAL_FEATURE_SET\t", paste(ko_genes[ko_genes %in% feature_set], collapse = ","), "\n", sep = "")

ko_info <- data.frame(
  gene = ko_genes,
  present_in_RNA = ko_genes %in% rna_feature_names,
  present_in_final_matrix = ko_genes %in% feature_set,
  stringsAsFactors = FALSE
)

write.table(
  ko_info,
  file = file.path(outdir, "ko_gene_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

manifest <- data.frame(
  group = character(),
  matrix_file = character(),
  n_genes = integer(),
  n_cells = integer(),
  present_ko_genes = character(),
  stringsAsFactors = FALSE
)

job_table <- data.frame(
  group = character(),
  gene = character(),
  present_in_matrix = logical(),
  suggested_runner = character(),
  stringsAsFactors = FALSE
)

for (grp in groups) {
  cells <- rownames(md)[as.character(md[["orig.ident"]]) == grp]
  cells <- intersect(cells, rna_cell_names)

  if (length(cells) == 0) {
    stop("no cells found for group: ", grp)
  }

  col_idx <- match(cells, rna_cell_names)
  count_matrix <- mat[row_idx, col_idx, drop = FALSE]

  out_file <- file.path(outdir, paste0("countMatrix_", grp, ".rds"))
  saveRDS(count_matrix, out_file, compress = "xz")

  present <- ko_genes[ko_genes %in% rownames(count_matrix)]

  manifest <- rbind(
    manifest,
    data.frame(
      group = grp,
      matrix_file = out_file,
      n_genes = nrow(count_matrix),
      n_cells = ncol(count_matrix),
      present_ko_genes = paste(present, collapse = ","),
      stringsAsFactors = FALSE
    )
  )

  suggested_runner <- if (grp == "B73xA554") {
    "run_sctenifold_knk_safe.R"
  } else {
    "run_sctenifold_knk.R"
  }

  for (g in ko_genes) {
    job_table <- rbind(
      job_table,
      data.frame(
        group = grp,
        gene = g,
        present_in_matrix = g %in% rownames(count_matrix),
        suggested_runner = suggested_runner,
        stringsAsFactors = FALSE
      )
    )
  }

  cat("EXTRACTED\t", grp, "\t", nrow(count_matrix), "\t", ncol(count_matrix), "\t", out_file, "\n", sep = "")
}

write.table(
  manifest,
  file = file.path(outdir, "manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  job_table,
  file = file.path(outdir, "job_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  ko_genes[ko_genes %in% feature_set],
  file.path(outdir, "ko_genes_present.txt")
)

cat("OUTPUT_DIR\t", outdir, "\n", sep = "")
cat("N_GROUPS\t", length(groups), "\n", sep = "")
cat("N_JOBS_TOTAL\t", nrow(job_table), "\n", sep = "")
cat("N_JOBS_PRESENT_IN_MATRIX\t", sum(job_table$present_in_matrix), "\n", sep = "")
cat("DONE\n")
