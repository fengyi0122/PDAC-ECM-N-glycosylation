script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/")) else file.path(getwd(), "scripts")
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

RAW <- existing_dataset_dir("geo_pdac", "GSE154778")
CANDIDATE_GENES <- read_gene_set("external_six_gene_module")
MARKERS <- list(
  CAF = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRB", "TAGLN", "MYL9"),
  Endothelial = c("PECAM1", "VWF", "KDR", "CDH5", "CLDN5", "RAMP2"),
  Myeloid = c("LYZ", "LST1", "CD68", "C1QA", "C1QB", "CSF1R", "AIF1"),
  T_NK = c("CD3D", "CD3E", "TRAC", "CD8A", "NKG7", "GNLY", "PRF1"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MSLN", "CEACAM6", "MUC1"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN")
)

extract_gene_matrix <- function() {
  wanted <- unique(c(CANDIDATE_GENES, unlist(MARKERS, use.names = FALSE)))
  con <- gzfile(file.path(RAW, "GSE154778_dgeMtx.csv.gz"), open = "rt")
  on.exit(close(con))
  header <- strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]]
  cells <- header[-1]
  rows <- list()
  repeat {
    lines <- readLines(con, n = 1000)
    if (length(lines) == 0) {
      break
    }
    for (line in lines) {
      parts <- strsplit(line, ",", fixed = TRUE)[[1]]
      gene <- parts[[1]]
      if (gene %in% wanted) {
        vals <- suppressWarnings(as.numeric(parts[-1]))
        vals[is.na(vals)] <- 0
        rows[[gene]] <- vals
      }
    }
  }
  mat <- do.call(rbind, rows)
  colnames(mat) <- cells
  rownames(mat) <- names(rows)
  mat <- log1p(mat)
  out_dt <- data.table::as.data.table(mat, keep.rownames = "gene")
  data.table::fwrite(out_dt, file.path(PROCESSED, "gse154778_selected_gene_log1p.tsv"), sep = "\t")
  mat
}

score_signature <- function(mat, genes) {
  found <- intersect(genes, rownames(mat))
  if (length(found) == 0) {
    return(list(score = rep(NA_real_, ncol(mat)), found = character()))
  }
  list(score = colMeans(mat[found, , drop = FALSE], na.rm = TRUE), found = found)
}

mat <- extract_gene_matrix()
cell_scores <- data.table::data.table(cell_id = colnames(mat))
coverage <- list()

candidate <- score_signature(mat, CANDIDATE_GENES)
cell_scores$candidate_module <- candidate$score
coverage[[1]] <- data.frame(signature = "candidate_module", n_genes = length(CANDIDATE_GENES), n_found = length(candidate$found), genes_found = paste(candidate$found, collapse = ","))

for (label in names(MARKERS)) {
  scored <- score_signature(mat, MARKERS[[label]])
  cell_scores[[sprintf("%s_marker_score", label)]] <- scored$score
  coverage[[length(coverage) + 1]] <- data.frame(signature = label, n_genes = length(MARKERS[[label]]), n_found = length(scored$found), genes_found = paste(scored$found, collapse = ","))
}

marker_cols <- names(cell_scores)[endsWith(names(cell_scores), "_marker_score")]
zmarkers <- as.data.frame(lapply(cell_scores[, ..marker_cols], zscore_vector, zero_if_constant = TRUE))
cell_scores$inferred_cell_type <- sub("_marker_score$", "", marker_cols[max.col(as.matrix(zmarkers), ties.method = "first")])
cell_scores$patient <- vapply(strsplit(cell_scores$cell_id, ":", fixed = TRUE), `[`, character(1), 1)
cell_scores$candidate_module_z <- zscore_vector(cell_scores$candidate_module, zero_if_constant = TRUE)

summary <- cell_scores[, .(
  n_cells = .N,
  candidate_module_mean = mean(candidate_module, na.rm = TRUE),
  candidate_module_z_mean = mean(candidate_module_z, na.rm = TRUE),
  candidate_module_detected_fraction = mean(candidate_module > 0, na.rm = TRUE)
), by = inferred_cell_type]
data.table::setorder(summary, -candidate_module_z_mean)

gene_records <- list()
for (gene in intersect(CANDIDATE_GENES, rownames(mat))) {
  tmp <- data.table::data.table(expr = as.numeric(mat[gene, ]), inferred_cell_type = cell_scores$inferred_cell_type)
  gene_summary <- tmp[, .(
    mean_log1p_expr = mean(expr, na.rm = TRUE),
    detected_fraction = mean(expr > 0, na.rm = TRUE)
  ), by = inferred_cell_type]
  gene_summary$gene <- gene
  gene_records[[length(gene_records) + 1]] <- gene_summary
}

data.table::fwrite(cell_scores, file.path(PROCESSED, "gse154778_marker_inferred_cell_scores.tsv"), sep = "\t")
data.table::fwrite(summary, file.path(RESULTS, "gse154778_candidate_module_by_inferred_cell_type.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(gene_records, fill = TRUE), file.path(RESULTS, "gse154778_candidate_gene_by_inferred_cell_type.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(coverage, fill = TRUE), file.path(RESULTS, "gse154778_marker_gene_coverage.tsv"), sep = "\t")
cat("Wrote", file.path(RESULTS, "gse154778_candidate_module_by_inferred_cell_type.tsv"), "\n")
