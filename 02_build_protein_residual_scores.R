script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/")) else file.path(getwd(), "scripts")
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

RAW <- existing_dataset_dir("cptac", "CPTAC-PDAC")

parse_proteome <- function() {
  protein <- data.table::fread(file.path(RAW, "proteomics_gene_level_MD_abundance_tumor.cct"))
  data.table::setnames(protein, 1, "gene")
  sample_cols <- setdiff(names(protein), "gene")
  protein[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
  protein <- protein[, lapply(.SD, safe_median), by = gene, .SDcols = sample_cols]
  mat <- as.matrix(protein[, ..sample_cols])
  rownames(mat) <- protein$gene
  storage.mode(mat) <- "numeric"
  mat
}

add_scores <- function(score_frame, zmat, ann, prefix) {
  records <- list()
  add_record <- function(score_name, n_features) {
    records[[length(records) + 1]] <<- data.frame(score = score_name, n_features = n_features)
  }

  score_name <- sprintf("%s_score_all", prefix)
  n_name <- sprintf("%s_n_all", prefix)
  score_frame[[score_name]] <- colMeans(zmat, na.rm = TRUE)
  score_frame[[n_name]] <- colSums(!is.na(zmat))
  add_record(score_name, nrow(zmat))

  for (field in c("division", "category")) {
    values <- sort(unique(ann[[field]][!is.na(ann[[field]]) & trimws(ann[[field]]) != ""]))
    for (value in values) {
      idx <- which(ann[[field]] == value)
      subset <- zmat[idx, , drop = FALSE]
      name <- sprintf("%s_%s_%s", prefix, field, clean_name(value))
      score_col <- sprintf("%s_score", name)
      n_col <- sprintf("%s_n", name)
      score_frame[[score_col]] <- colMeans(subset, na.rm = TRUE)
      score_frame[[n_col]] <- colSums(!is.na(subset))
      add_record(score_col, nrow(subset))
    }
  }

  list(score_frame = score_frame, records = data.table::rbindlist(records, fill = TRUE))
}

build_layer <- function(layer, filename, feature_col, protein_mat) {
  layer_name <- layer
  mapping <- data.table::fread(file.path(PROCESSED, "matrisome_glyco_map.tsv"))
  mapping <- unique(mapping[layer == layer_name, .(feature_id, gene, division, category)], by = c("feature_id", "gene"))

  glyco <- data.table::fread(file.path(RAW, filename))
  data.table::setnames(glyco, feature_col, "feature_id")
  sample_cols <- setdiff(names(glyco), c("feature_id", "Gene"))
  glyco[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]

  merged <- merge(mapping, glyco, by = "feature_id", all = FALSE)
  sample_cols <- intersect(sample_cols, colnames(protein_mat))
  sample_cols <- intersect(sample_cols, names(merged))
  merged <- merged[gene %in% rownames(protein_mat)]
  values <- as.matrix(merged[, ..sample_cols])
  storage.mode(values) <- "numeric"
  merged <- merged[rowSums(!is.na(values)) >= 50]
  values <- as.matrix(merged[, ..sample_cols])
  storage.mode(values) <- "numeric"

  raw_z <- zscore_rows_matrix(values)
  residual <- matrix(NA_real_, nrow = nrow(merged), ncol = length(sample_cols), dimnames = list(NULL, sample_cols))
  for (i in seq_len(nrow(merged))) {
    gene <- merged$gene[[i]]
    y <- as.numeric(values[i, ])
    x <- as.numeric(protein_mat[gene, sample_cols])
    mask <- is.finite(y) & is.finite(x)
    if (sum(mask) < 50 || stats::sd(x[mask]) == 0) {
      next
    }
    fit <- stats::lm.fit(cbind(1, x[mask]), y[mask])
    alpha <- fit$coefficients[[1]]
    beta <- fit$coefficients[[2]]
    residual[i, ] <- y - (alpha + beta * x)
  }
  residual_z <- zscore_rows_matrix(residual)

  ann <- merged[, .(feature_id, gene, division, category)]
  score_frame <- data.frame(row.names = sample_cols)
  first <- add_scores(score_frame, raw_z, ann, sprintf("%s_raw", layer))
  second <- add_scores(first$score_frame, residual_z, ann, sprintf("%s_protein_residual", layer))
  records <- data.table::rbindlist(list(first$records, second$records), fill = TRUE)
  records[, layer := layer]
  list(scores = second$score_frame, ann = ann, records = records)
}

protein_mat <- parse_proteome()
sample_ids <- colnames(protein_mat)

matrisome <- data.table::fread(MATRISOME_FILE)
mat_genes <- sort(intersect(unique(matrisome$gene), rownames(protein_mat)))
protein_z <- zscore_rows_matrix(protein_mat[mat_genes, , drop = FALSE])

features <- data.frame(sample_id = sample_ids, row.names = sample_ids)
features$proteome_matrisome_score <- colMeans(protein_z, na.rm = TRUE)
features$proteome_matrisome_n <- colSums(!is.na(protein_z))
summary_records <- data.table::data.table(layer = "proteome", score = "proteome_matrisome_score", n_features = length(mat_genes))

layers <- list(
  list(layer = "nglyco_peptide", filename = "N-glycoproteomics_peptide_level_ratio_tumor.cct", feature_col = "Sequence"),
  list(layer = "nglyco_site", filename = "N-glycoproteomics_Site_level_ratio_tumor.cct", feature_col = "Modifications")
)

for (spec in layers) {
  built <- build_layer(spec$layer, spec$filename, spec$feature_col, protein_mat)
  features <- cbind(features, built$scores[rownames(features), , drop = FALSE])
  data.table::fwrite(built$ann, file.path(PROCESSED, sprintf("%s_matrisome_features_retained.tsv", spec$layer)), sep = "\t")
  summary_records <- data.table::rbindlist(list(summary_records, built$records), fill = TRUE)
}

data.table::fwrite(data.table::as.data.table(features), file.path(PROCESSED, "glyco_matrisome_features.tsv"), sep = "\t")
data.table::fwrite(summary_records, file.path(RESULTS, "glyco_matrisome_score_summary.tsv"), sep = "\t")
cat("Wrote", file.path(PROCESSED, "glyco_matrisome_features.tsv"), "\n")
cat("Wrote", file.path(RESULTS, "glyco_matrisome_score_summary.tsv"), "\n")
