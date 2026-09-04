script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

RAW <- existing_dataset_dir("geo_pdac", "GSE199102")
OUT <- file.path(RESULTS, "analysis_results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CANDIDATES <- read_gene_set("external_six_gene_module")
SIGNATURES <- list(
  candidate_module = CANDIDATES,
  caf_stroma_score = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "PDGFRB"),
  mycaf_score = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"),
  endothelial_score = c("PECAM1", "VWF", "KDR", "CDH5", "ICAM1", "VCAM1", "SELE"),
  immune_score = c("PTPRC", "CD3D", "CD3E", "CD8A", "CD68", "LST1", "NKG7")
)

normalize_id <- function(value) {
  gsub("-", ".", as.character(value), fixed = TRUE)
}

score_signature <- function(expression_matrix, genes) {
  found <- intersect(genes, rownames(expression_matrix))
  if (length(found) == 0) {
    return(list(score = rep(NA_real_, ncol(expression_matrix)), found = character()))
  }
  list(score = colMeans(zscore_rows_matrix(log1p(expression_matrix[found, , drop = FALSE])), na.rm = TRUE), found = found)
}

safe_mean <- function(x) {
  if (all(!is.finite(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

expression <- data.table::fread(file.path(RAW, "GSE199102_hPDAC_WTA_20210222T2101_Q3Norm_TargetCountMatrix.txt.gz"))
data.table::setnames(expression, 1, "gene")
sample_columns <- setdiff(names(expression), "gene")
expression[, (sample_columns) := lapply(.SD, as.numeric), .SDcols = sample_columns]
expression <- expression[, lapply(.SD, safe_median), by = gene, .SDcols = sample_columns]
expression_matrix <- as.matrix(expression[, ..sample_columns])
rownames(expression_matrix) <- expression$gene
colnames(expression_matrix) <- normalize_id(sample_columns)
storage.mode(expression_matrix) <- "numeric"

properties <- data.table::fread(file.path(RAW, "GSE199102_Broad_PDAC_WTA_AllSamples_SegmentProperties.txt.gz"))
properties[, sample_norm := normalize_id(Sample_ID)]
shared <- sort(intersect(colnames(expression_matrix), properties$sample_norm))
expression_matrix <- expression_matrix[, shared, drop = FALSE]
properties <- properties[match(shared, sample_norm)]

scores <- data.table::data.table(sample_id = shared)
coverage <- list()
for (name in names(SIGNATURES)) {
  scored <- score_signature(expression_matrix, SIGNATURES[[name]])
  scores[[name]] <- scored$score
  coverage[[length(coverage) + 1]] <- data.table::data.table(signature = name, n_genes = length(SIGNATURES[[name]]), n_found = length(scored$found), genes_found = paste(scored$found, collapse = ","))
}
for (column in intersect(c("Segment", "Patient", "TreatmentClass", "Status", "Not_malignant"), names(properties))) {
  scores[[column]] <- properties[[column]]
}
data.table::fwrite(scores, file.path(PROCESSED, "gse199102_geomx_segment_scores.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(coverage, fill = TRUE), file.path(RESULTS, "gse199102_signature_coverage.tsv"), sep = "\t")

score_columns <- names(SIGNATURES)
patient_segment <- scores[, lapply(.SD, safe_mean), by = .(Patient, Segment), .SDcols = score_columns]
candidate_wide <- data.table::dcast(patient_segment, Patient ~ Segment, value.var = "candidate_module")
set.seed(2026071204)
paired_records <- list()
for (comparison in list(c("CAF", "Epithelial"), c("CAF", "Immune"))) {
  group_a <- comparison[[1]]
  group_b <- comparison[[2]]
  if (!all(c(group_a, group_b) %in% names(candidate_wide))) {
    next
  }
  frame <- candidate_wide[is.finite(get(group_a)) & is.finite(get(group_b))]
  differences <- frame[[group_a]] - frame[[group_b]]
  bootstrap_means <- replicate(2000, mean(sample(differences, replace = TRUE)))
  paired_records[[length(paired_records) + 1]] <- data.table::data.table(
    comparison = paste0(group_a, "_vs_", group_b),
    n_patients = length(differences),
    mean_difference = mean(differences),
    median_difference = stats::median(differences),
    ci_low = stats::quantile(bootstrap_means, 0.025),
    ci_high = stats::quantile(bootstrap_means, 0.975),
    wilcoxon_p = stats::wilcox.test(differences, mu = 0, exact = FALSE)$p.value,
    positive_patients = sum(differences > 0)
  )
}
paired_effects <- data.table::rbindlist(paired_records, fill = TRUE)
data.table::fwrite(patient_segment, file.path(OUT, "geomx_patient_segment_means.tsv"), sep = "\t")
data.table::fwrite(paired_effects, file.path(OUT, "geomx_patient_level_paired_effects.tsv"), sep = "\t")

patient_correlation_records <- list()
for (score in setdiff(score_columns, "candidate_module")) {
  for (patient in unique(scores$Patient)) {
    frame <- scores[Patient == patient, .(candidate_module, score_value = get(score))]
    frame <- frame[stats::complete.cases(frame)]
    if (nrow(frame) < 4 || stats::sd(frame$candidate_module) == 0 || stats::sd(frame$score_value) == 0) {
      next
    }
    correlation <- suppressWarnings(stats::cor(frame$candidate_module, frame$score_value, method = "spearman"))
    patient_correlation_records[[length(patient_correlation_records) + 1]] <- data.table::data.table(Patient = patient, score = score, n_segments = nrow(frame), spearman_rho = correlation)
  }
}
patient_correlations <- data.table::rbindlist(patient_correlation_records, fill = TRUE)
set.seed(2026071205)
correlation_summary <- patient_correlations[, {
  bootstrap_medians <- replicate(2000, stats::median(sample(spearman_rho, replace = TRUE), na.rm = TRUE))
  positive_count <- sum(spearman_rho > 0, na.rm = TRUE)
  list(
    n_patients = .N,
    median_rho = stats::median(spearman_rho, na.rm = TRUE),
    ci_low = stats::quantile(bootstrap_medians, 0.025),
    ci_high = stats::quantile(bootstrap_medians, 0.975),
    positive_patients = positive_count,
    sign_test_p = stats::binom.test(positive_count, .N, p = 0.5, alternative = "greater")$p.value
  )
}, by = score]
correlation_summary[, fdr := stats::p.adjust(sign_test_p, method = "BH")]
data.table::fwrite(patient_correlations, file.path(OUT, "geomx_within_patient_correlations.tsv"), sep = "\t")
data.table::fwrite(correlation_summary, file.path(OUT, "geomx_patient_correlation_summary.tsv"), sep = "\t")
cat("Wrote", file.path(OUT, "geomx_patient_level_paired_effects.tsv"), "\n")
