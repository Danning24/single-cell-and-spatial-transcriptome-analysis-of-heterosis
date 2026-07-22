# Portable project root; set MAIZE_PROJECT_ROOT before running.
project_root <- Sys.getenv("MAIZE_PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) stop("Set MAIZE_PROJECT_ROOT to the analysis data root")

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) args[[1]] else project_root
result_dir <- if (length(args) >= 2) args[[2]] else file.path(project_dir, "sctenifold_results")
output_root <- if (length(args) >= 3) args[[3]] else file.path(project_dir, "sctenifold_downstream_analysis")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

has_go_db <- requireNamespace("GO.db", quietly = TRUE)
has_annotation_dbi <- requireNamespace("AnnotationDbi", quietly = TRUE)
has_keggrest <- requireNamespace("KEGGREST", quietly = TRUE)

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

wrap_text <- function(x, width = 55) {
  vapply(x, function(s) paste(strwrap(as.character(s), width = width), collapse = "\n"), character(1))
}

safe_neglog10 <- function(x) {
  x[is.na(x)] <- 1
  x[x <= 0] <- .Machine$double.xmin
  -log10(x)
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

drop_self_ko <- function(df, ko) {
  if (nrow(df) == 0 || !"gene" %in% names(df)) return(df)
  # 鍗曟暡鏃?ko 鏄竴涓熀鍥狅紱鍙屾暡鏃?ko 浼氬啓鎴?
  # Zm00001d007234_plus_Zm00001d043795銆?
  # 瀵岄泦鍒嗘瀽鍙湅琚壈鍔ㄧ殑涓嬫父鍩哄洜锛屾墍浠ヨ繖閲岃鎶婃墍鏈夎鏁插熀鍥犻兘鍘绘帀銆?
  ko_genes <- unlist(strsplit(ko, "_plus_", fixed = TRUE))
  df[!(tolower(df$gene) %in% tolower(ko_genes)), , drop = FALSE]
}

standardize_gene_col <- function(df) {
  if (nrow(df) == 0) return(df)
  if (!"gene" %in% names(df)) {
    names(df)[1] <- "gene"
  }
  df
}

go_class_labels <- c(BP = "biological_process", CC = "cellular_component", MF = "molecular_function")
go_class_colors <- c(BP = "#4F8FC7", CC = "#94B96E", MF = "#F29B52")

maize_gamer_base <- "https://download.maizegdb.org/GeneFunction_and_Expression/Maize-GAMER/B73_RefGen_v4/d.non_red_gaf"
maize_refgen_base <- "https://download.maizegdb.org/Zm-B73-REFERENCE-GRAMENE-4.0"
maize_gff3_file <- "Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz"
maize_go_terms_file <- "Zm-B73-REFERNECE-GRAMENE-4.0-GMs-GOTerms.csv.gz"
maize_full_gene_file <- "Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.full_gene_data.txt.gz"
uniprot_maize_kegg_file <- "uniprot_maize_UP000007305_kegg.tsv"
maize_gamer_files <- c(
  "maize.B73.AGPv4.arabidopsis.gaf.gz",
  "maize.B73.AGPv4.argot2.5.gaf.gz",
  "maize.B73.AGPv4.fanngo.gaf.gz",
  "maize.B73.AGPv4.Interproscan5.gaf.gz",
  "maize.B73.AGPv4.pannzer.gaf.gz",
  "maize.B73.AGPv4.uniprot.plants.gaf.gz"
)

empty_enrichment <- function() {
  data.frame(
    ID = character(),
    Description = character(),
    Ontology = character(),
    Count = integer(),
    InputSize = integer(),
    TermSize = integer(),
    UniverseSize = integer(),
    GeneRatio = numeric(),
    BgRatio = numeric(),
    pvalue = numeric(),
    genes = character(),
    p.adjust = numeric(),
    stringsAsFactors = FALSE
  )
}

download_file_quiet <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(TRUE)
  ok <- tryCatch({
    download.file(url, destfile = dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok && nzchar(Sys.which("curl"))) {
    status <- suppressWarnings(system2("curl", c("-k", "-L", "-o", dest, url), stdout = TRUE, stderr = TRUE))
    ok <- file.exists(dest) && file.info(dest)$size > 0
  }
  if (!ok && file.exists(dest)) unlink(dest)
  ok
}

extract_attr <- function(attr, key) {
  pat <- paste0("(^|;)", key, "=([^;]+)")
  m <- regexec(pat, attr)
  out <- regmatches(attr, m)
  vapply(out, function(x) if (length(x) >= 3) URLdecode(x[[3]]) else NA_character_, character(1))
}

extract_zm_id <- function(x) {
  m <- regexpr("Zm[0-9A-Za-z]+d[0-9]+", x)
  out <- ifelse(m > 0, regmatches(x, m), NA_character_)
  out
}

make_gene_label <- function(genes, gene_annot = NULL) {
  genes <- as.character(genes)
  if (is.null(gene_annot) || nrow(gene_annot) == 0) return(genes)
  idx <- match(genes, gene_annot$gene)
  symbol <- gene_annot$symbol[idx]
  desc <- gene_annot$description[idx]
  label <- genes
  use_symbol <- !is.na(symbol) & nzchar(symbol) & symbol != genes
  label[use_symbol] <- paste0(symbol[use_symbol], " (", genes[use_symbol], ")")
  use_desc <- !use_symbol & !is.na(desc) & nzchar(desc) & desc != genes
  short_desc <- desc
  long_desc <- !is.na(short_desc) & nchar(short_desc) > 36
  short_desc[long_desc] <- paste0(substr(short_desc[long_desc], 1, 33), "...")
  label[use_desc] <- paste0(short_desc[use_desc], " (", genes[use_desc], ")")
  label
}

load_gene_annotation <- function(annotation_dir, universe_genes) {
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  gff_path <- file.path(annotation_dir, maize_gff3_file)
  gff_url <- paste0(maize_refgen_base, "/", maize_gff3_file)
  download_file_quiet(gff_url, gff_path)

  annot <- data.frame(gene = unique(universe_genes), symbol = unique(universe_genes),
                      description = NA_character_, ncbi_geneid = NA_character_, stringsAsFactors = FALSE)
  if (file.exists(gff_path)) {
    con <- gzfile(gff_path, open = "rt")
    on.exit(close(con), add = TRUE)
    gff <- read.delim(con, header = FALSE, comment.char = "#", stringsAsFactors = FALSE, quote = "", fill = TRUE)
    if (nrow(gff) > 0 && ncol(gff) >= 9) {
      genes <- gff[gff[[3]] == "gene", , drop = FALSE]
      attrs <- genes[[9]]
      gff_annot <- data.frame(
        gene = extract_zm_id(attrs),
        symbol = extract_attr(attrs, "Name"),
        description = extract_attr(attrs, "description"),
        dbxref = extract_attr(attrs, "Dbxref"),
        stringsAsFactors = FALSE
      )
      gff_annot$ncbi_geneid <- sub(".*GeneID:([0-9]+).*", "\\1", gff_annot$dbxref)
      gff_annot$ncbi_geneid[!grepl("^[0-9]+$", gff_annot$ncbi_geneid)] <- NA_character_
      gff_annot <- gff_annot[!is.na(gff_annot$gene) & gff_annot$gene %in% universe_genes, c("gene", "symbol", "description", "ncbi_geneid"), drop = FALSE]
      gff_annot <- gff_annot[!duplicated(gff_annot$gene), , drop = FALSE]
      annot <- merge(annot[, "gene", drop = FALSE], gff_annot, by = "gene", all.x = TRUE)
      annot$symbol[is.na(annot$symbol) | !nzchar(annot$symbol)] <- annot$gene[is.na(annot$symbol) | !nzchar(annot$symbol)]
    }
  }
  annot
}

load_uniprot_kegg_crossref <- function(annotation_dir, universe_genes) {
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(annotation_dir, uniprot_maize_kegg_file)
  url <- "https://rest.uniprot.org/uniprotkb/stream"
  if (file.exists(path)) {
    first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
    if (!length(first_line) || !grepl("^Entry\\t", first_line)) {
      unlink(path)
    }
  }
  if (!file.exists(path)) {
    if (nzchar(Sys.which("curl"))) {
      curl_args <- c(
        "-k", "-L", "--retry", "3", "-G", "-o", path, url,
        "--data-urlencode", "compressed=false",
        "--data-urlencode", "format=tsv",
        "--data-urlencode", "query=organism_id:4577",
        "--data-urlencode", "fields=accession,gene_names,xref_kegg,gene_primary,xref_geneid,protein_name"
      )
      suppressWarnings(system2("curl", curl_args, stdout = TRUE, stderr = TRUE))
    } else {
      direct_url <- paste0(
        url, "?compressed=false&format=tsv",
        "&query=%28organism_id%3A4577%29",
        "&fields=accession,gene_names,xref_kegg,gene_primary,xref_geneid,protein_name"
      )
      download_file_quiet(direct_url, path)
    }
  }
  first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  if (!length(first_line) || !grepl("^Entry\\t", first_line)) return(data.frame())
  x <- tryCatch(read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = ""), error = function(e) data.frame())
  if (nrow(x) == 0) return(data.frame())
  kegg_col <- names(x)[grepl("^KEGG$", names(x), ignore.case = TRUE)][1]
  gene_name_col <- names(x)[grepl("^Gene Names$", names(x), ignore.case = TRUE)][1]
  primary_col <- names(x)[grepl("Gene Names \\(primary\\)", names(x), ignore.case = TRUE)][1]
  protein_col <- names(x)[grepl("Protein names", names(x), ignore.case = TRUE)][1]
  if (is.na(kegg_col) || is.na(gene_name_col)) return(data.frame())
  gene <- extract_zm_id(as.character(x[[gene_name_col]]))
  kegg_raw <- as.character(x[[kegg_col]])
  keep <- !is.na(gene) & gene %in% universe_genes & !is.na(kegg_raw) & nzchar(kegg_raw)
  if (!any(keep)) return(data.frame())
  out <- data.frame(
    gene = gene[keep],
    kegg_gene = kegg_raw[keep],
    symbol = if (!is.na(primary_col)) as.character(x[[primary_col]][keep]) else NA_character_,
    description = if (!is.na(protein_col)) as.character(x[[protein_col]][keep]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out <- tidyr::separate_rows(out, kegg_gene, sep = ";")
  out$kegg_gene <- trimws(out$kegg_gene)
  out <- out[grepl("^zma:", out$kegg_gene), , drop = FALSE]
  out <- out[!duplicated(out[, c("gene", "kegg_gene")]), , drop = FALSE]
  out
}

curl_read_lines <- function(url) {
  if (!nzchar(Sys.which("curl"))) return(character())
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  suppressWarnings(system2("curl", c("-k", "-L", "--retry", "3", "-o", tmp, url), stdout = TRUE, stderr = TRUE))
  if (!file.exists(tmp) || file.info(tmp)$size == 0) return(character())
  readLines(tmp, warn = FALSE)
}

kegg_link_pathway_zma <- function() {
  if (has_keggrest) {
    kg <- tryCatch(KEGGREST::keggLink("pathway", "zma"), error = function(e) character())
    if (length(kg) > 0) return(kg)
  }
  lines <- curl_read_lines("https://rest.kegg.jp/link/pathway/zma")
  if (!length(lines)) return(character())
  x <- read.delim(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  out <- x[[2]]
  names(out) <- x[[1]]
  out
}

kegg_pathway_names_zma <- function() {
  if (has_keggrest) {
    pn <- tryCatch(KEGGREST::keggList("pathway", "zma"), error = function(e) character())
    if (length(pn) > 0) return(pn)
  }
  lines <- curl_read_lines("https://rest.kegg.jp/list/pathway/zma")
  if (!length(lines)) return(character())
  x <- read.delim(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  out <- x[[2]]
  names(out) <- x[[1]]
  out
}

kegg_gene_list_zma <- function() {
  if (has_keggrest) {
    gl <- tryCatch(KEGGREST::keggList("zma"), error = function(e) character())
    if (length(gl) > 0) return(gl)
  }
  lines <- curl_read_lines("https://rest.kegg.jp/list/zma")
  if (!length(lines)) return(character())
  x <- read.delim(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  out <- x[[2]]
  names(out) <- x[[1]]
  out
}

ensure_maize_gaf_files <- function(annotation_dir) {
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(annotation_dir, maize_gamer_files)
  missing <- !file.exists(paths)
  if (any(missing)) {
    for (fn in maize_gamer_files[missing]) {
      url <- paste0(maize_gamer_base, "/", fn)
      dest <- file.path(annotation_dir, fn)
      if (!download_file_quiet(url, dest)) {
        message("Could not download GO annotation: ", fn)
      }
    }
  }
  paths[file.exists(paths)]
}

read_gaf_table <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  x <- tryCatch(
    read.delim(con, header = FALSE, comment.char = "!", stringsAsFactors = FALSE, quote = "", fill = TRUE),
    error = function(e) data.frame()
  )
  if (nrow(x) == 0 || ncol(x) < 5) return(data.frame())
  gene <- x[[2]]
  go_id <- x[[5]]
  qualifier <- if (ncol(x) >= 4) x[[4]] else ""
  data.frame(gene = gene, GOID = go_id, qualifier = qualifier, stringsAsFactors = FALSE)
}

read_maize_go_terms_csv <- function(annotation_dir, universe_genes) {
  path <- file.path(annotation_dir, maize_go_terms_file)
  url <- paste0(maize_refgen_base, "/", maize_go_terms_file)
  if (!download_file_quiet(url, path)) return(data.frame())
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  x <- tryCatch(read.csv(con, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (nrow(x) == 0) return(data.frame())
  gene_col <- names(x)[vapply(x, function(v) any(grepl("Zm[0-9A-Za-z]+d[0-9]+", as.character(v))), logical(1))][1]
  go_col <- names(x)[vapply(x, function(v) any(grepl("^GO:[0-9]+$", as.character(v))), logical(1))][1]
  if (is.na(gene_col) || is.na(go_col)) return(data.frame())
  out <- data.frame(
    gene = extract_zm_id(as.character(x[[gene_col]])),
    GOID = as.character(x[[go_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene) & out$gene %in% universe_genes & grepl("^GO:[0-9]+$", out$GOID), , drop = FALSE]
  unique(out)
}

load_go_annotation <- function(annotation_dir, universe_genes) {
  gaf_paths <- ensure_maize_gaf_files(annotation_dir)
  go_parts <- list()
  if (length(gaf_paths) > 0) {
    gaf_map <- bind_rows(lapply(gaf_paths, read_gaf_table))
    gaf_map <- gaf_map[grepl("^Zm[0-9A-Za-z]+d[0-9]+$", gaf_map$gene) & grepl("^GO:[0-9]+$", gaf_map$GOID), , drop = FALSE]
    gaf_map <- gaf_map[!grepl("NOT", gaf_map$qualifier, fixed = TRUE), c("gene", "GOID"), drop = FALSE]
    go_parts[["gaf"]] <- gaf_map
  }
  csv_map <- read_maize_go_terms_csv(annotation_dir, universe_genes)
  if (nrow(csv_map) > 0) go_parts[["goterms_csv"]] <- csv_map
  if (length(go_parts) == 0) return(data.frame())
  go_map <- bind_rows(go_parts)
  go_map <- unique(go_map)
  go_map <- go_map[go_map$gene %in% universe_genes, , drop = FALSE]
  if (nrow(go_map) == 0) return(data.frame())

  go_ids <- unique(go_map$GOID)
  go_terms <- data.frame(GOID = go_ids, Description = go_ids, Ontology = NA_character_, stringsAsFactors = FALSE)
  if (has_go_db && has_annotation_dbi) {
    term_df <- suppressMessages(AnnotationDbi::select(
      GO.db::GO.db,
      keys = go_ids,
      keytype = "GOID",
      columns = c("GOID", "TERM", "ONTOLOGY")
    ))
    term_df <- term_df[!duplicated(term_df$GOID), , drop = FALSE]
    go_terms <- merge(go_terms[, "GOID", drop = FALSE], term_df, by = "GOID", all.x = TRUE)
    names(go_terms) <- c("GOID", "Description", "Ontology")
    go_terms$Description[is.na(go_terms$Description)] <- go_terms$GOID[is.na(go_terms$Description)]
  }
  go_terms$Ontology[is.na(go_terms$Ontology)] <- "GO"
  go_map <- merge(go_map, go_terms, by = "GOID", all.x = TRUE)
  go_map <- go_map[go_map$Ontology %in% c("BP", "CC", "MF"), , drop = FALSE]
  unique(go_map)
}

run_go_enrichment <- function(sig_genes, universe_genes, go_map, ontology = NULL, min_overlap = 1) {
  sig_genes <- unique(sig_genes[!is.na(sig_genes) & nzchar(sig_genes)])
  sig_genes <- intersect(sig_genes, universe_genes)
  if (!is.null(ontology)) {
    go_map <- go_map[go_map$Ontology == ontology, , drop = FALSE]
  }
  go_map <- go_map[go_map$gene %in% universe_genes, , drop = FALSE]
  if (length(sig_genes) == 0 || nrow(go_map) == 0) return(empty_enrichment())
  M <- length(unique(go_map$gene))
  N <- length(sig_genes)
  pieces <- lapply(split(go_map, go_map$GOID), function(tm) {
    term_genes <- unique(tm$gene)
    hit <- intersect(sig_genes, term_genes)
    x <- length(hit)
    if (x < min_overlap) return(NULL)
    K <- length(term_genes)
    data.frame(
      ID = tm$GOID[[1]],
      Description = tm$Description[[1]],
      Ontology = tm$Ontology[[1]],
      Count = x,
      InputSize = N,
      TermSize = K,
      UniverseSize = M,
      GeneRatio = x / N,
      BgRatio = K / M,
      pvalue = stats::phyper(x - 1, K, M - K, N, lower.tail = FALSE),
      genes = paste(sort(hit), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  if (is.null(out) || nrow(out) == 0) return(empty_enrichment())
  out$p.adjust <- p.adjust(out$pvalue, method = "BH")
  out <- out[order(out$p.adjust, out$pvalue, -out$Count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

run_term_enrichment <- function(sig_genes, universe_genes, term_map, min_overlap = 1) {
  sig_genes <- unique(sig_genes[!is.na(sig_genes) & nzchar(sig_genes)])
  sig_genes <- intersect(sig_genes, universe_genes)
  term_map <- term_map[term_map$gene %in% universe_genes, , drop = FALSE]
  if (length(sig_genes) == 0 || nrow(term_map) == 0) return(empty_enrichment())
  M <- length(unique(term_map$gene))
  N <- length(sig_genes)
  pieces <- lapply(split(term_map, term_map$TERM_ID), function(tm) {
    term_genes <- unique(tm$gene)
    hit <- intersect(sig_genes, term_genes)
    x <- length(hit)
    if (x < min_overlap) return(NULL)
    K <- length(term_genes)
    data.frame(
      ID = tm$TERM_ID[[1]],
      Description = tm$Description[[1]],
      Ontology = if ("Ontology" %in% names(tm)) tm$Ontology[[1]] else "KEGG",
      Count = x,
      InputSize = N,
      TermSize = K,
      UniverseSize = M,
      GeneRatio = x / N,
      BgRatio = K / M,
      pvalue = stats::phyper(x - 1, K, M - K, N, lower.tail = FALSE),
      genes = paste(sort(hit), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  if (is.null(out) || nrow(out) == 0) return(empty_enrichment())
  out$p.adjust <- p.adjust(out$pvalue, method = "BH")
  out <- out[order(out$p.adjust, out$pvalue, -out$Count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

build_kegg_map <- function(gene_annot, universe_genes, annotation_dir = NULL) {
  ga <- data.frame()
  if (nrow(gene_annot) > 0 && "ncbi_geneid" %in% names(gene_annot)) {
    ga <- gene_annot[gene_annot$gene %in% universe_genes & !is.na(gene_annot$ncbi_geneid) & nzchar(gene_annot$ncbi_geneid), , drop = FALSE]
    if (nrow(ga) > 0) {
      ga$kegg_gene <- paste0("zma:", ga$ncbi_geneid)
    }
  }
  if (nrow(ga) == 0 && !is.null(annotation_dir)) {
    kegg_genes <- kegg_gene_list_zma()
    if (length(kegg_genes) > 0) {
      kg_desc <- as.character(kegg_genes)
      kg_zm <- extract_zm_id(kg_desc)
      kg_map <- data.frame(
        gene = kg_zm,
        kegg_gene = names(kegg_genes),
        stringsAsFactors = FALSE
      )
      kg_map <- kg_map[!is.na(kg_map$gene) & kg_map$gene %in% universe_genes, , drop = FALSE]
      if (nrow(kg_map) > 0) {
        ga <- unique(kg_map)
      }
    }
  }
  if (nrow(ga) == 0 && !is.null(annotation_dir)) {
    uni <- load_uniprot_kegg_crossref(annotation_dir, universe_genes)
    if (nrow(uni) > 0) {
      ga <- uni[, c("gene", "kegg_gene"), drop = FALSE]
    }
  }
  if (nrow(ga) == 0) return(data.frame())
  message("KEGG gene cross-reference genes overlapping universe: ", length(unique(ga$gene)))
  kg <- kegg_link_pathway_zma()
  if (length(kg) == 0) return(data.frame())
  kegg_df <- data.frame(
    kegg_gene = sub("^zma:", "zma:", names(kg)),
    pathway_id = sub("^path:", "", as.character(kg)),
    stringsAsFactors = FALSE
  )
  kegg_df <- kegg_df[kegg_df$kegg_gene %in% ga$kegg_gene, , drop = FALSE]
  if (nrow(kegg_df) == 0) return(data.frame())
  path_names <- kegg_pathway_names_zma()
  path_df <- data.frame(
    pathway_id = sub("^path:", "", names(path_names)),
    Description = sub(" - Zea mays \\(maize\\)$", "", as.character(path_names)),
    stringsAsFactors = FALSE
  )
  out <- merge(kegg_df, ga[, c("gene", "kegg_gene"), drop = FALSE], by = "kegg_gene", all.x = TRUE)
  out <- merge(out, path_df, by = "pathway_id", all.x = TRUE)
  out$Description[is.na(out$Description) | !nzchar(out$Description)] <- out$pathway_id[is.na(out$Description) | !nzchar(out$Description)]
  out$TERM_ID <- out$pathway_id
  out$Ontology <- "KEGG"
  unique(out[, c("gene", "TERM_ID", "Description", "Ontology"), drop = FALSE])
}

select_go_terms_for_plot <- function(go_df, top_n_per_ontology = 15) {
  go_df %>%
    filter(Ontology %in% c("BP", "CC", "MF"), Count > 0) %>%
    mutate(
      Ontology = factor(Ontology, levels = c("BP", "CC", "MF")),
      GO_class = go_class_labels[as.character(Ontology)],
      RichFactor = Count / TermSize
    ) %>%
    group_by(Ontology) %>%
    arrange(desc(Count), pvalue, .by_group = TRUE) %>%
    slice_head(n = top_n_per_ontology) %>%
    ungroup()
}

plot_go_bar <- function(go_df, title, out_prefix, top_n_per_ontology = 18) {
  d <- select_go_terms_for_plot(go_df, top_n_per_ontology)
  if (nrow(d) == 0) {
    writeLines(paste("No GO terms for", title), paste0(out_prefix, "_NO_TERMS.txt"))
    return(invisible(NULL))
  }
  d <- d %>% arrange(Ontology, desc(Count), pvalue)
  d$GO_term_label <- d$Description
  if (any(duplicated(d$GO_term_label))) {
    d$GO_term_label <- paste0(d$Description, " (", d$GO_class, ")")
  }
  d$GO_term <- factor(d$GO_term_label, levels = unique(d$GO_term_label))
  p <- ggplot(d, aes(x = GO_term, y = Count, fill = Ontology)) +
    geom_col(width = 0.88) +
    scale_fill_manual(values = go_class_colors, labels = go_class_labels, name = NULL) +
    labs(title = wrap_text(title, 78), x = "GO_term", y = "Num of Genes") +
    theme_grey(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "top",
      legend.text = element_text(face = "italic", size = 11),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, face = "bold", color = "grey35"),
      axis.text.y = element_text(face = "bold"),
      plot.margin = margin(16, 16, 36, 16)
    )
  ggsave(paste0(out_prefix, ".png"), p, width = 15, height = 8, dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 15, height = 8)
}

plot_go_richfactor <- function(go_df, title, out_prefix, top_n = 20) {
  d <- go_df %>%
    filter(Ontology %in% c("BP", "CC", "MF"), Count > 0) %>%
    mutate(RichFactor = Count / TermSize) %>%
    arrange(pvalue, desc(RichFactor), desc(Count)) %>%
    slice_head(n = top_n)
  if (nrow(d) == 0) {
    writeLines(paste("No GO terms for", title), paste0(out_prefix, "_NO_TERMS.txt"))
    return(invisible(NULL))
  }
  d$GO_Term <- factor(wrap_text(d$Description, 48), levels = rev(unique(wrap_text(d$Description, 48))))
  p <- ggplot(d, aes(x = RichFactor, y = GO_Term, size = Count, color = pvalue)) +
    geom_point(alpha = 0.95) +
    scale_color_gradient(low = "red", high = "blue", name = "pvalue") +
    scale_size_continuous(range = c(2.5, 9), name = "Gene_number") +
    labs(title = title, x = "Rich factor", y = "GO_Term") +
    theme_grey(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 22),
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 11),
      plot.margin = margin(16, 28, 16, 24)
    )
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = max(6, 0.36 * nrow(d) + 2), dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = max(6, 0.36 * nrow(d) + 2))
}

plot_kegg_dotplot <- function(kegg_df, title, out_prefix, top_n = 20) {
  if (nrow(kegg_df) == 0) {
    writeLines(paste("No KEGG terms for", title), paste0(out_prefix, "_NO_TERMS.txt"))
    return(invisible(NULL))
  }
  d <- kegg_df %>%
    mutate(RichFactor = Count / TermSize) %>%
    arrange(pvalue, desc(RichFactor), desc(Count)) %>%
    slice_head(n = top_n)
  d$Pathway <- factor(wrap_text(d$Description, 48), levels = rev(unique(wrap_text(d$Description, 48))))
  p <- ggplot(d, aes(x = RichFactor, y = Pathway, size = Count, color = pvalue)) +
    geom_point(alpha = 0.95) +
    scale_color_gradient(low = "red", high = "blue", name = "pvalue") +
    scale_size_continuous(range = c(2.5, 9), name = "Gene_number") +
    labs(title = wrap_text(title, 70), x = "Rich factor", y = "KEGG pathway") +
    theme_grey(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 15),
      axis.text.y = element_text(size = 11),
      plot.margin = margin(16, 28, 16, 24)
    )
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = max(6, 0.36 * nrow(d) + 2), dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = max(6, 0.36 * nrow(d) + 2))
}

plot_top_genes <- function(diff_df, ko, out_prefix, top_n = 25, gene_annot = NULL) {
  if (nrow(diff_df) == 0 || !"Z" %in% names(diff_df)) return(invisible(NULL))
  d <- diff_df %>%
    mutate(absZ = abs(Z)) %>%
    arrange(desc(absZ)) %>%
    slice_head(n = top_n)
  d$gene_label <- make_gene_label(d$gene, gene_annot)
  d$gene_label <- factor(d$gene_label, levels = rev(d$gene_label))
  p <- ggplot(d, aes(x = absZ, y = gene_label, fill = Z)) +
    geom_col(width = 0.75) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "Z") +
    labs(title = paste0(ko, ": top downstream perturbed genes"), x = "|Z|", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold", size = 8))
  ggsave(paste0(out_prefix, ".png"), p, width = 7, height = max(5, 0.23 * nrow(d) + 1.8), dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 7, height = max(5, 0.23 * nrow(d) + 1.8))
}

plot_count_bar <- function(count_df, title, out_prefix) {
  p <- ggplot(count_df, aes(x = reorder(KO, n_sig), y = n_sig, fill = KO)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = n_sig), vjust = -0.3, fontface = "bold", size = 3.8) +
    labs(title = title, x = NULL, y = "Significant perturbed genes") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold", angle = 25, hjust = 1))
  ymax <- max(count_df$n_sig, na.rm = TRUE)
  p <- p + coord_cartesian(ylim = c(0, ymax * 1.2 + 1))
  ggsave(paste0(out_prefix, ".png"), p, width = 7.5, height = 4.8, dpi = 300)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 7.5, height = 4.8)
}

plot_pairwise_overlap <- function(presence_df, title, out_prefix) {
  kos <- setdiff(names(presence_df), "gene")
  mat <- matrix(0, nrow = length(kos), ncol = length(kos), dimnames = list(kos, kos))
  for (i in kos) {
    for (j in kos) {
      mat[i, j] <- sum(presence_df[[i]] & presence_df[[j]])
    }
  }
  write.csv(mat, paste0(out_prefix, ".csv"))
  if (length(kos) < 2) {
    writeLines("Only one KO in this group; pairwise overlap heatmap was skipped.", paste0(out_prefix, "_SKIPPED.txt"))
    return(invisible(NULL))
  }
  png(paste0(out_prefix, ".png"), width = 6.5, height = 5.8, units = "in", res = 300)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
           main = title, color = colorRampPalette(c("white", "#3182bd"))(100))
  dev.off()
  pdf(paste0(out_prefix, ".pdf"), width = 6.5, height = 5.8)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
           main = title, color = colorRampPalette(c("white", "#3182bd"))(100))
  dev.off()
}

plot_z_heatmap <- function(diff_long, title, out_prefix, top_n = 80, gene_annot = NULL) {
  if (nrow(diff_long) == 0 || !"Z" %in% names(diff_long)) return(invisible(NULL))
  top_genes <- diff_long %>%
    group_by(gene) %>%
    summarise(max_absZ = max(abs(Z), na.rm = TRUE), n_KO = n_distinct(KO), .groups = "drop") %>%
    arrange(desc(n_KO), desc(max_absZ)) %>%
    slice_head(n = top_n) %>%
    pull(gene)
  mat_df <- diff_long %>%
    filter(gene %in% top_genes) %>%
    select(gene, KO, Z) %>%
    distinct() %>%
    pivot_wider(names_from = KO, values_from = Z, values_fill = 0)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- make_gene_label(mat_df$gene, gene_annot)
  png(paste0(out_prefix, ".png"), width = 8, height = max(6, 0.15 * nrow(mat) + 2), units = "in", res = 300)
  pheatmap(mat, main = title, color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
           fontsize_row = 6.5, fontsize_col = 9, cluster_cols = FALSE, border_color = NA)
  dev.off()
  pdf(paste0(out_prefix, ".pdf"), width = 8, height = max(6, 0.15 * nrow(mat) + 2))
  pheatmap(mat, main = title, color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
           fontsize_row = 6.5, fontsize_col = 9, cluster_cols = FALSE, border_color = NA)
  dev.off()
}

read_universe <- function(group_name, diff_files) {
  input_path <- file.path(project_dir, "sctenifold_inputs", paste0("countMatrix_", group_name, ".rds"))
  if (file.exists(input_path)) {
    x <- readRDS(input_path)
    rn <- rownames(x)
    if (!is.null(rn) && length(rn) > 0) return(unique(rn))
  }
  genes <- unique(unlist(lapply(diff_files, function(f) {
    d <- standardize_gene_col(read_csv_safe(f))
    if ("gene" %in% names(d)) d$gene else character()
  })))
  genes[!is.na(genes) & nzchar(genes)]
}

parse_result_files <- function(result_dir) {
  # 褰撳墠缁撴灉宸茬粡鏁寸悊鎴愨€滄瘡涓?KO 涓€涓枃浠跺す鈥濓紝鎵€浠ヨ繖閲屽繀椤?recursive = TRUE銆?
  # 鍚屾椂鍏煎锛?
  #   鍗曟暡锛欱73xA554__Zm00001d011969.safe_ma1.diffRegulation.csv
  #   鍙屾暡锛欱73xA554__Zm00001d007234_plus_Zm00001d043795.diffRegulation.csv
  files <- list.files(
    result_dir,
    pattern = "\\.(diffRegulation|sig_genes|sig_gene_counts)\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  base <- basename(files)
  m <- regexec(
    "^(.+)__((?:Zm[0-9A-Za-z]+d[0-9]+)(?:_plus_Zm[0-9A-Za-z]+d[0-9]+)*)(\\.safe_ma1)?\\.(diffRegulation|sig_genes|sig_gene_counts)\\.csv$",
    base,
    perl = TRUE
  )
  parts <- regmatches(base, m)
  rows <- lapply(seq_along(parts), function(i) {
    p <- parts[[i]]
    if (length(p) == 0) return(NULL)
    data.frame(group = p[[2]], KO = p[[3]], safe = p[[4]], type = p[[5]], path = files[[i]], stringsAsFactors = FALSE)
  })
  out <- bind_rows(rows)
  if (nrow(out) == 0) {
    stop("No scTenifoldKnk result CSV files were parsed from: ", result_dir)
  }
  out
}

process_group <- function(group_name, file_index, go_map, kegg_map, gene_annot) {
  message("Processing group: ", group_name)
  group_out <- file.path(output_root, paste0(group_name, "_sctenifold_KO_analysis"))
  per_ko_dir <- file.path(group_out, "per_KO_analysis")
  combined_dir <- file.path(group_out, "combined_comparison")
  dir.create(per_ko_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

  idx <- file_index[file_index$group == group_name, , drop = FALSE]
  kos <- sort(unique(idx$KO))
  universe <- read_universe(group_name, idx$path[idx$type == "diffRegulation"])
  write.csv(data.frame(gene = universe), file.path(group_out, "background_genes_used_for_enrichment.csv"), row.names = FALSE)
  group_go_map <- go_map[go_map$gene %in% universe, , drop = FALSE]
  group_kegg_map <- kegg_map[kegg_map$gene %in% universe, , drop = FALSE]
  group_gene_annot <- gene_annot[gene_annot$gene %in% universe, , drop = FALSE]
  write.csv(group_gene_annot, file.path(group_out, "gene_id_to_name_annotation.csv"), row.names = FALSE)

  all_sig <- list()
  all_diff <- list()
  count_rows <- list()
  for (ko in kos) {
    message("  KO: ", ko)
    sig_path <- idx$path[idx$KO == ko & idx$type == "sig_genes"][1]
    diff_path <- idx$path[idx$KO == ko & idx$type == "diffRegulation"][1]

    # GO/KEGG 缁撴灉璺熼殢姣忎釜 KO 鑷繁鐨?scTenifoldKnk 缁撴灉鏂囦欢澶逛繚瀛橈紝鏂逛究鍚庣画鏌ョ湅銆?
    # 渚嬪锛歴ctenifold_results/B73xA554__Zm00001d011969/GO_KEGG/
    source_result_dir <- dirname(diff_path)
    ko_out <- file.path(source_result_dir, "GO_KEGG")
    dir.create(ko_out, recursive = TRUE, showWarnings = FALSE)
    sig_df <- standardize_gene_col(read_csv_safe(sig_path))
    diff_df <- standardize_gene_col(read_csv_safe(diff_path))
    sig_df <- drop_self_ko(sig_df, ko)
    diff_df <- drop_self_ko(diff_df, ko)
    sig_genes <- if (nrow(sig_df) > 0) unique(sig_df$gene) else character()
    sig_genes <- sig_genes[!is.na(sig_genes) & nzchar(sig_genes)]
    write.csv(data.frame(gene = sig_genes), file.path(ko_out, paste0(ko, "_input_sig_genes_for_enrichment.csv")), row.names = FALSE)
    if (nrow(sig_df) > 0) {
      sig_df$KO <- ko
      all_sig[[ko]] <- sig_df
    }
    if (nrow(diff_df) > 0) {
      diff_df$KO <- ko
      diff_df$gene_label <- make_gene_label(diff_df$gene, group_gene_annot)
      all_diff[[ko]] <- diff_df
      write.csv(diff_df, file.path(ko_out, paste0(ko, "_diffRegulation_downstream_no_selfKO.csv")), row.names = FALSE)
      plot_top_genes(diff_df, ko, file.path(ko_out, paste0(ko, "_top_downstream_perturbed_genes_by_absZ")), gene_annot = group_gene_annot)
    }
    count_rows[[ko]] <- data.frame(KO = ko, n_sig = length(sig_genes), stringsAsFactors = FALSE)

    if (nrow(group_go_map) > 0) {
      go_parts <- lapply(c("BP", "CC", "MF"), function(ont) {
        e <- run_go_enrichment(sig_genes, universe, group_go_map, ontology = ont)
        e$analysis <- paste0("GO_", ont)
        e$KO <- ko
        write.csv(e, file.path(ko_out, paste0(ko, "_GO_", ont, "_enrichment.csv")), row.names = FALSE)
        e
      })
      go_all <- bind_rows(go_parts)
      write.csv(go_all, file.path(ko_out, paste0(ko, "_GO_enrichment_all_ontologies.csv")), row.names = FALSE)
      plot_go_bar(go_all, paste0(ko, " GO enrichment"), file.path(ko_out, paste0(ko, "_GO_barplot")))
      plot_go_richfactor(go_all, "Statistics of GO Enrichment", file.path(ko_out, paste0(ko, "_GO_richfactor_dotplot")))
    } else {
      writeLines("GO annotation was not available or did not overlap with the universe genes.", file.path(ko_out, paste0(ko, "_GO_NOT_RUN.txt")))
    }

    if (nrow(group_kegg_map) > 0) {
      unlink(file.path(ko_out, paste0(ko, "_KEGG_NOT_RUN.txt")))
      kegg <- run_term_enrichment(sig_genes, universe, group_kegg_map)
      kegg$analysis <- rep("KEGG", nrow(kegg))
      kegg$KO <- rep(ko, nrow(kegg))
      write.csv(kegg, file.path(ko_out, paste0(ko, "_KEGG_enrichment.csv")), row.names = FALSE)
      plot_kegg_dotplot(kegg, paste0(ko, " KEGG enrichment"), file.path(ko_out, paste0(ko, "_KEGG_dotplot")))
    } else {
      writeLines("KEGG mapping was not available or did not overlap with the universe genes.", file.path(ko_out, paste0(ko, "_KEGG_NOT_RUN.txt")))
    }
  }

  count_df <- bind_rows(count_rows)
  write.csv(count_df, file.path(combined_dir, "KO_sig_gene_counts.csv"), row.names = FALSE)
  plot_count_bar(count_df, paste0(group_name, ": significant perturbed genes per virtual KO"),
                 file.path(combined_dir, "KO_sig_gene_counts_barplot"))

  sig_long <- bind_rows(all_sig)
  diff_long <- bind_rows(all_diff)
  write.csv(sig_long, file.path(combined_dir, "all_KO_sig_genes_long.csv"), row.names = FALSE)
  write.csv(diff_long, file.path(combined_dir, "all_KO_diffRegulation_downstream_no_selfKO_long.csv"), row.names = FALSE)

  all_genes <- sort(unique(sig_long$gene))
  presence <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  for (ko in kos) {
    presence[[ko]] <- presence$gene %in% sig_long$gene[sig_long$KO == ko]
  }
  presence$n_KO <- rowSums(presence[, kos, drop = FALSE])
  presence <- presence[order(-presence$n_KO, presence$gene), , drop = FALSE]
  write.csv(presence, file.path(combined_dir, "sig_gene_presence_matrix.csv"), row.names = FALSE)
  write.csv(presence[presence$n_KO >= 2, , drop = FALSE], file.path(combined_dir, "shared_sig_genes_in_2plus_KOs.csv"), row.names = FALSE)
  write.csv(presence[presence$n_KO == 1, , drop = FALSE], file.path(combined_dir, "KO_specific_sig_genes.csv"), row.names = FALSE)
  plot_pairwise_overlap(presence[, c("gene", kos), drop = FALSE],
                        paste0(group_name, ": pairwise overlap of significant genes"),
                        file.path(combined_dir, "pairwise_sig_gene_overlap_heatmap"))
  plot_z_heatmap(diff_long, paste0(group_name, ": top perturbed genes Z-score across KOs"),
                 file.path(combined_dir, "top_perturbed_genes_Z_heatmap"), gene_annot = group_gene_annot)

  summary <- data.frame(
    group = group_name,
    n_KO = length(kos),
    background_genes = length(universe),
    GO_annotated_background_genes = length(unique(group_go_map$gene)),
    KEGG_annotated_background_genes = length(unique(group_kegg_map$gene)),
    total_unique_sig_genes = length(all_genes),
    shared_sig_genes_in_2plus_KOs = sum(presence$n_KO >= 2),
    output_dir = group_out,
    stringsAsFactors = FALSE
  )
  write.csv(summary, file.path(group_out, "analysis_summary.csv"), row.names = FALSE)
}

file_index <- parse_result_files(result_dir)
write.csv(file_index, file.path(output_root, "input_result_file_index.csv"), row.names = FALSE)

annotation_dir <- file.path(output_root, "maize_GO_annotation_maizeGAMER")
universe_all <- unique(unlist(lapply(split(file_index, file_index$group), function(idx) read_universe(idx$group[[1]], idx$path[idx$type == "diffRegulation"]))))
gene_annot <- load_gene_annotation(annotation_dir, universe_all)
write.csv(gene_annot, file.path(annotation_dir, "maize_gene_id_to_name_annotation_used.csv"), row.names = FALSE)
go_map <- load_go_annotation(annotation_dir, universe_all)
write.csv(go_map, file.path(annotation_dir, "maizeGAMER_GO_gene_map_used.csv"), row.names = FALSE)
kegg_cache_file <- file.path(annotation_dir, "maize_KEGG_gene_pathway_map_used.csv")
if (file.exists(kegg_cache_file) && file.info(kegg_cache_file)$size > 10) {
  kegg_map <- tryCatch(read.csv(kegg_cache_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) data.frame())
  if (nrow(kegg_map) > 0) {
    message("Using cached KEGG map: ", kegg_cache_file)
  } else {
    kegg_map <- build_kegg_map(gene_annot, universe_all, annotation_dir)
    write.csv(kegg_map, kegg_cache_file, row.names = FALSE)
  }
} else {
  kegg_map <- build_kegg_map(gene_annot, universe_all, annotation_dir)
  write.csv(kegg_map, kegg_cache_file, row.names = FALSE)
}

for (group_name in sort(unique(file_index$group))) {
  process_group(group_name, file_index, go_map, kegg_map, gene_annot)
}

readme <- c(
  "Maize scTenifoldKnk downstream analysis",
  "",
  paste0("Input result directory: ", result_dir),
  paste0("Output root: ", output_root),
  "",
  "Groups are analyzed separately:",
  paste0("- ", sort(unique(file_index$group)), "_sctenifold_KO_analysis", collapse = "\n"),
  "",
  "Per-KO outputs:",
  "- <KO>_diffRegulation_downstream_no_selfKO.csv",
  "- <KO>_top_downstream_perturbed_genes_by_absZ.png/pdf",
  "- <KO>_GO_barplot.png/pdf",
  "- <KO>_GO_richfactor_dotplot.png/pdf",
  "- <KO>_GO_BP/CC/MF_enrichment.csv and <KO>_GO_enrichment_all_ontologies.csv",
  "- <KO>_KEGG_enrichment.csv",
  "- <KO>_KEGG_dotplot.png/pdf",
  "",
  "Combined outputs:",
  "- KO_sig_gene_counts.csv and KO_sig_gene_counts_barplot.png/pdf",
  "- sig_gene_presence_matrix.csv",
  "- shared_sig_genes_in_2plus_KOs.csv",
  "- KO_specific_sig_genes.csv",
  "- pairwise_sig_gene_overlap_heatmap.csv/png/pdf",
  "- top_perturbed_genes_Z_heatmap.png/pdf",
  "",
  "Notes:",
  "- The KO gene itself is removed from downstream affected genes and enrichment inputs; for double KO, both KO genes are removed.",
  "- Gene labels use gene symbol/name from B73 RefGen_v4 GFF3 when available, with Zm00001d ID kept in parentheses.",
  "- GO enrichment uses maize-GAMER B73 RefGen_v4 Zm00001d.2 GAF annotations when downloads are available.",
  "- KEGG enrichment first uses B73 RefGen_v4 GFF3 GeneID cross-references and then falls back to UniProt maize KEGG cross-references for Zm00001d-to-zma mapping."
)
writeLines(readme, file.path(output_root, "README.txt"))

message("All done. Output root: ", output_root)
