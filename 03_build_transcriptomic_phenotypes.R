script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/")) else file.path(getwd(), "scripts")
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

RAW <- existing_dataset_dir("cptac", "CPTAC-PDAC")

SIGNATURES <- list(
  cd8_t_score = c("CD8A", "CD8B", "TRAC", "CD3D", "CD3E"),
  cytotoxicity_score = c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY"),
  tam_m2_score = c("CD68", "CD163", "MRC1", "MSR1", "CSF1R", "C1QA", "C1QB", "C1QC"),
  treg_score = c("FOXP3", "IL2RA", "CTLA4", "CCR8", "IKZF2", "TIGIT"),
  caf_stroma_score = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRB"),
  mycaf_score = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"),
  icaf_score = c("IL6", "CXCL12", "CFD", "HAS1", "CXCL14", "PDGFRA"),
  endothelial_score = c("PECAM1", "VWF", "KDR", "CDH5", "ICAM1", "VCAM1", "SELE"),
  checkpoint_score = c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
)

parse_mean_number <- function(value) {
  if (is.na(value)) {
    return(NA_real_)
  }
  hits <- regmatches(as.character(value), gregexpr("-?\\d+(?:\\.\\d+)?", as.character(value), perl = TRUE))[[1]]
  if (length(hits) == 0 || hits[[1]] == "") {
    return(NA_real_)
  }
  mean(as.numeric(hits))
}

score_signature <- function(expr_mat, genes) {
  found <- intersect(genes, rownames(expr_mat))
  if (length(found) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_mat)), found = character()))
  }
  list(score = colMeans(zscore_rows_matrix(expr_mat[found, , drop = FALSE]), na.rm = TRUE), found = found)
}

mrna <- data.table::fread(file.path(RAW, "mRNA_RSEM_UQ_log2_Tumor.cct"))
data.table::setnames(mrna, 1, "gene")
sample_cols <- setdiff(names(mrna), "gene")
mrna[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
mrna <- mrna[, lapply(.SD, safe_median), by = gene, .SDcols = sample_cols]
expr_mat <- as.matrix(mrna[, ..sample_cols])
rownames(expr_mat) <- mrna$gene
storage.mode(expr_mat) <- "numeric"

out <- data.table::data.table(sample_id = sample_cols)
coverage <- list()
for (name in names(SIGNATURES)) {
  scored <- score_signature(expr_mat, SIGNATURES[[name]])
  out[[name]] <- scored$score
  coverage[[length(coverage) + 1]] <- data.frame(
    signature = name,
    n_genes = length(SIGNATURES[[name]]),
    n_found = length(scored$found),
    genes_found = paste(scored$found, collapse = ",")
  )
}

clinical <- data.table::fread(file.path(RAW, "clinical_table_140.tsv"))
clinical_cols <- c(
  "Neoplastic_cellularity", "Stromal_fraction", "Inflammation_fraction",
  "Acinar_fraction", "Islet_fraction", "Non_neoplastic_duct", "Fat_fraction"
)
for (col in intersect(clinical_cols, names(clinical))) {
  parsed <- vapply(clinical[[col]], parse_mean_number, numeric(1))
  lookup <- setNames(parsed, clinical$case_id)
  mean_name <- sprintf("%s_mean", tolower(col))
  z_name <- sprintf("%s_z", tolower(col))
  out[[mean_name]] <- unname(lookup[out$sample_id])
  out[[z_name]] <- zscore_vector(out[[mean_name]])
}

out$stromal_immune_composite <- (
  zscore_vector(out$caf_stroma_score) +
    zscore_vector(out$tam_m2_score) +
    zscore_vector(out$treg_score) +
    zscore_vector(out$stromal_fraction_mean) -
    zscore_vector(out$cd8_t_score) -
    zscore_vector(out$cytotoxicity_score)
) / 6

data.table::fwrite(out, file.path(PROCESSED, "immune_phenotypes.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(coverage, fill = TRUE), file.path(RESULTS, "immune_signature_gene_coverage.tsv"), sep = "\t")
cat("Wrote", file.path(PROCESSED, "immune_phenotypes.tsv"), "\n")
cat("Wrote", file.path(RESULTS, "immune_signature_gene_coverage.tsv"), "\n")
