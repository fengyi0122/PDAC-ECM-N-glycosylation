script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

GEOMX_RAW <- existing_dataset_dir("geo_pdac", "GSE199102")
KEY_OUTCOMES <- c("myCAF", "endothelial", "CAF_stroma", "stromal_immune_composite")
REFERENCE_CANDIDATES <- read_gene_set("external_six_gene_module")

fit_simple <- function(outcome, predictor) {
  frame <- data.table::data.table(y = zscore_vector(outcome), x = zscore_vector(predictor))
  frame <- frame[stats::complete.cases(frame)]
  result <- c(n = nrow(frame), beta = NA_real_, p_value = NA_real_)
  if (nrow(frame) < 20 || stats::sd(frame$x) == 0) {
    return(result)
  }
  fit <- stats::lm(y ~ x, data = frame)
  coefficients <- summary(fit)$coefficients
  c(n = nrow(frame), beta = coefficients["x", "Estimate"], p_value = coefficients["x", "Pr(>|t|)"])
}

feature_models <- data.table::fread(file.path(RESULTS, "analysis_results", "glycosite_feature_coherence_models.tsv"))
supported_features <- feature_models[outcome %in% KEY_OUTCOMES & beta > 0 & fdr < 0.10]
cptac_summary <- supported_features[, .(
  cptac_layers_positive = data.table::uniqueN(layer),
  cptac_outcomes_positive = data.table::uniqueN(outcome),
  cptac_positive_feature_outcome_associations = .N,
  cptac_best_fdr = min(fdr, na.rm = TRUE),
  cptac_max_beta = max(beta, na.rm = TRUE),
  cptac_supported_outcomes = paste(sort(unique(outcome)), collapse = ",")
), by = gene]
data.table::setorder(cptac_summary, cptac_best_fdr)
candidate_genes <- cptac_summary[cptac_layers_positive >= 2 & cptac_outcomes_positive >= 2, gene]
candidate_genes <- unique(c(candidate_genes, REFERENCE_CANDIDATES))

tcga_scores <- data.table::fread(file.path(PROCESSED, "tcga_paad_bulk_validation_scores.tsv"))
tcga_expression <- read_xena_tcga_expression(candidate_genes)
tcga_sample_columns <- setdiff(names(tcga_expression), "gene")

tcga_records <- list()
tcga_outcomes <- intersect(c("caf_stroma_score", "mycaf_score", "endothelial_score", "stromal_immune_composite"), names(tcga_scores))
for (row_index in seq_len(nrow(tcga_expression))) {
  gene <- tcga_expression$gene[[row_index]]
  expression <- as.numeric(tcga_expression[row_index, ..tcga_sample_columns])
  names(expression) <- tcga_sample_columns
  expression <- unname(expression[tcga_scores$sample_id])
  for (outcome in tcga_outcomes) {
    model <- fit_simple(tcga_scores[[outcome]], expression)
    tcga_records[[length(tcga_records) + 1]] <- data.table::data.table(gene = gene, outcome = outcome, n = model[["n"]], beta = model[["beta"]], p_value = model[["p_value"]])
  }
}
tcga_models <- data.table::rbindlist(tcga_records, fill = TRUE)
tcga_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(tcga_models, file.path(RESULTS, "c4_tcga_candidate_gene_models.tsv"), sep = "\t")
tcga_summary <- tcga_models[beta > 0, .(
  tcga_positive_outcomes = sum(fdr < 0.05, na.rm = TRUE),
  tcga_best_fdr = min(fdr, na.rm = TRUE),
  tcga_max_beta = max(beta, na.rm = TRUE),
  tcga_supported_outcomes = paste(sort(unique(outcome[fdr < 0.05])), collapse = ",")
), by = gene]

normalize_id <- function(value) {
  gsub("-", ".", as.character(value), fixed = TRUE)
}

geomx_expression <- data.table::fread(file.path(GEOMX_RAW, "GSE199102_hPDAC_WTA_20210222T2101_Q3Norm_TargetCountMatrix.txt.gz"))
data.table::setnames(geomx_expression, 1, "gene")
geomx_sample_columns <- setdiff(names(geomx_expression), "gene")
geomx_expression <- geomx_expression[gene %in% candidate_genes]
geomx_expression[, (geomx_sample_columns) := lapply(.SD, as.numeric), .SDcols = geomx_sample_columns]
geomx_expression <- geomx_expression[, lapply(.SD, safe_median), by = gene, .SDcols = geomx_sample_columns]
properties <- data.table::fread(file.path(GEOMX_RAW, "GSE199102_Broad_PDAC_WTA_AllSamples_SegmentProperties.txt.gz"))
properties[, sample_norm := normalize_id(Sample_ID)]
colnames(geomx_expression)[-1] <- normalize_id(colnames(geomx_expression)[-1])
shared <- sort(intersect(setdiff(names(geomx_expression), "gene"), properties$sample_norm))
properties <- properties[match(shared, sample_norm)]

geomx_records <- list()
for (row_index in seq_len(nrow(geomx_expression))) {
  gene <- geomx_expression$gene[[row_index]]
  expression <- zscore_vector(log1p(as.numeric(geomx_expression[row_index, ..shared])))
  segments <- data.table::data.table(Patient = properties$Patient, Segment = properties$Segment, expression = expression)[is.finite(expression)]
  patient_segment <- segments[, .(expression = mean(expression)), by = .(Patient, Segment)]
  patient_wide <- data.table::dcast(patient_segment, Patient ~ Segment, value.var = "expression")
  if (!all(c("CAF", "Epithelial") %in% names(patient_wide))) {
    next
  }
  paired <- patient_wide[is.finite(CAF) & is.finite(Epithelial)]
  differences <- paired$CAF - paired$Epithelial
  if (length(differences) < 10) {
    next
  }
  test <- stats::wilcox.test(differences, mu = 0, exact = FALSE)
  geomx_records[[length(geomx_records) + 1]] <- data.table::data.table(
    gene = gene,
    n_paired_patients = length(differences),
    geomx_caf_mean = mean(paired$CAF),
    geomx_epithelial_mean = mean(paired$Epithelial),
    geomx_caf_delta = mean(differences),
    geomx_positive_patients = sum(differences > 0),
    geomx_caf_p_value = test$p.value
  )
}
geomx_summary <- data.table::rbindlist(geomx_records, fill = TRUE)
geomx_summary[, geomx_caf_fdr := stats::p.adjust(geomx_caf_p_value, method = "BH")]
data.table::fwrite(geomx_summary, file.path(RESULTS, "c4_gse199102_candidate_gene_caf_enrichment.tsv"), sep = "\t")

core <- merge(cptac_summary, tcga_summary, by = "gene", all.x = TRUE)
core <- merge(core, geomx_summary, by = "gene", all.x = TRUE)
core[is.na(tcga_positive_outcomes), tcga_positive_outcomes := 0]
core[, core_status := data.table::fcase(
  cptac_layers_positive >= 2 & cptac_outcomes_positive >= 2 & tcga_positive_outcomes >= 2 & n_paired_patients >= 10 & geomx_caf_fdr < 0.05 & geomx_caf_delta > 0,
  "CORE",
  cptac_layers_positive >= 2 & (tcga_positive_outcomes >= 1 | (n_paired_patients >= 10 & geomx_caf_fdr < 0.10 & geomx_caf_delta > 0)),
  "SUPPORTING",
  default = "WEAK"
)]
core[, core_status_rank := match(core_status, c("CORE", "SUPPORTING", "WEAK"))]
data.table::setorder(core, core_status_rank, cptac_best_fdr)
core[, core_status_rank := NULL]
data.table::fwrite(core, file.path(RESULTS, "c4_candidate_core_set.tsv"), sep = "\t")
cat("Wrote", file.path(RESULTS, "c4_candidate_core_set.tsv"), "\n")
