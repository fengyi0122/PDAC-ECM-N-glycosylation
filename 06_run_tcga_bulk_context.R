script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

CANDIDATE_GENES <- read_gene_set("external_six_gene_module")
SIGNATURES <- list(
  cd8_t_score = c("CD8A", "CD8B", "TRAC", "CD3D", "CD3E"),
  cytotoxicity_score = c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY"),
  tam_m2_score = c("CD68", "CD163", "MRC1", "MSR1", "CSF1R", "C1QA", "C1QB", "C1QC"),
  treg_score = c("FOXP3", "IL2RA", "CTLA4", "CCR8", "IKZF2", "TIGIT"),
  caf_stroma_score = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRB"),
  mycaf_score = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"),
  mycaf_score_no_postn = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2"),
  icaf_score = c("IL6", "CXCL12", "CFD", "HAS1", "CXCL14", "PDGFRA"),
  endothelial_score = c("PECAM1", "VWF", "KDR", "CDH5", "ICAM1", "VCAM1", "SELE"),
  checkpoint_score = c("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT")
)

score_signature <- function(expr_mat, genes) {
  found <- intersect(genes, rownames(expr_mat))
  if (length(found) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_mat)), found = character()))
  }
  list(score = colMeans(zscore_rows_matrix(expr_mat[found, , drop = FALSE]), na.rm = TRUE), found = found)
}

fit_model <- function(data, outcome) {
  frame <- data.table::data.table(
    y = zscore_vector(data[[outcome]]),
    x = zscore_vector(data$candidate_ecm_glycoprotein_module)
  )
  frame <- frame[stats::complete.cases(frame)]
  result <- data.table::data.table(outcome = outcome, n = nrow(frame), beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, adj_r2 = NA_real_)
  if (nrow(frame) < 50 || stats::sd(frame$x) == 0) {
    return(result)
  }
  fit <- stats::lm(y ~ x, data = frame)
  coefficients <- summary(fit)$coefficients
  interval <- stats::confint(fit, "x", level = 0.95)
  result[, `:=`(
    beta = coefficients["x", "Estimate"],
    se = coefficients["x", "Std. Error"],
    ci_low = interval[[1]],
    ci_high = interval[[2]],
    p_value = coefficients["x", "Pr(>|t|)"],
    adj_r2 = summary(fit)$adj.r.squared
  )]
  result
}

expression <- read_xena_tcga_expression(unique(c(CANDIDATE_GENES, unlist(SIGNATURES))))
sample_columns <- setdiff(names(expression), "gene")
expression_matrix <- as.matrix(expression[, ..sample_columns])
rownames(expression_matrix) <- expression$gene
storage.mode(expression_matrix) <- "numeric"

scores <- data.table::data.table(sample_id = sample_columns)
candidate <- score_signature(expression_matrix, CANDIDATE_GENES)
scores$candidate_ecm_glycoprotein_module <- candidate$score
coverage <- list(data.table::data.table(signature = "candidate_ecm_glycoprotein_module", n_genes = length(CANDIDATE_GENES), n_found = length(candidate$found), genes_found = paste(candidate$found, collapse = ",")))

for (name in names(SIGNATURES)) {
  scored <- score_signature(expression_matrix, SIGNATURES[[name]])
  scores[[name]] <- scored$score
  coverage[[length(coverage) + 1]] <- data.table::data.table(signature = name, n_genes = length(SIGNATURES[[name]]), n_found = length(scored$found), genes_found = paste(scored$found, collapse = ","))
}

scores$stromal_immune_composite <- (
  zscore_vector(scores$caf_stroma_score) +
    zscore_vector(scores$tam_m2_score) +
    zscore_vector(scores$treg_score) +
    zscore_vector(scores$endothelial_score) -
    zscore_vector(scores$cd8_t_score) -
    zscore_vector(scores$cytotoxicity_score)
) / 6

model_outcomes <- c(names(SIGNATURES), "stromal_immune_composite")
models <- data.table::rbindlist(lapply(model_outcomes, function(outcome) fit_model(scores, outcome)), fill = TRUE)
models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::setorder(models, fdr, p_value)
data.table::fwrite(scores, file.path(PROCESSED, "tcga_paad_bulk_validation_scores.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(coverage, fill = TRUE), file.path(RESULTS, "tcga_paad_signature_coverage.tsv"), sep = "\t")
data.table::fwrite(models, file.path(RESULTS, "tcga_paad_bulk_validation_models.tsv"), sep = "\t")
cat("Wrote", file.path(RESULTS, "tcga_paad_bulk_validation_models.tsv"), "\n")
