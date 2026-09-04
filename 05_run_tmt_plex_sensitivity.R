script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table", "readxl"))
ensure_project_dirs()

RAW <- existing_dataset_dir("cptac", "CPTAC-PDAC")
mapping_file <- file.path(RAW, "S061_CPTAC_PDA_TMT11_Label_to_Sample_Mapping_File_JHU_r1_Feb2021.xlsx")
if (!file.exists(mapping_file)) {
  stop("The CPTAC TMT11 label-to-sample mapping file is unavailable.", call. = FALSE)
}

raw_mapping <- data.table::as.data.table(readxl::read_excel(mapping_file, sheet = "JHU_TMT11_Mapping"))
raw_mapping <- raw_mapping[grepl("Intact_Glycoproteome", raw_mapping[["Folder Name"]], fixed = TRUE)]
participant_columns <- grep("^TMT11-.* Participant ID$", names(raw_mapping), value = TRUE)

mapping <- data.table::rbindlist(lapply(participant_columns, function(participant_column) {
  channel <- sub("^TMT11-(.*) Participant ID$", "\\1", participant_column)
  specimen_column <- sub(" Participant ID$", " Specimen Label", participant_column)
  data.table::data.table(
    sample_id = trimws(as.character(raw_mapping[[participant_column]])),
    plex_dataset_name = as.character(raw_mapping[["Folder Name"]]),
    plex_number = suppressWarnings(as.integer(sub("^([0-9]+)CPTAC.*", "\\1", raw_mapping[["Folder Name"]]))),
    channel = channel,
    specimen_label = trimws(as.character(raw_mapping[[specimen_column]]))
  )
}), fill = TRUE)

mapping <- mapping[grepl("^C3[NL]-", sample_id)]
sample_map <- mapping[, .(
  n_channel_records = .N,
  n_distinct_plexes = data.table::uniqueN(plex_dataset_name),
  plex_dataset_name = paste(unique(plex_dataset_name), collapse = ";"),
  plex_number = paste(unique(plex_number), collapse = ";"),
  channels = paste(channel, collapse = ";"),
  specimen_labels = paste(specimen_label, collapse = ";")
), by = sample_id]
sample_map[, plex_id := ifelse(n_distinct_plexes == 1, sprintf("plex_%02d", as.integer(plex_number)), NA_character_)]
data.table::fwrite(sample_map, file.path(RESULTS, "pdc_pda_tmt11_sample_plex_map.tsv"), sep = "\t")

glyco <- data.table::fread(file.path(PROCESSED, "glyco_matrisome_features.tsv"))
phenotypes <- data.table::fread(file.path(PROCESSED, "immune_phenotypes.tsv"))
clinical <- data.table::fread(file.path(RAW, "clinical_table_140.tsv"))
clinical[, sample_id := case_id]
clinical[, histology_group := data.table::fifelse(histology_diagnosis == "PDAC", "PDAC", "Adenosquamous")]
analysis_data <- merge(glyco, phenotypes, by = "sample_id", all = FALSE)
analysis_data <- merge(analysis_data, clinical[, .(sample_id, histology_group)], by = "sample_id", all.x = TRUE)

expression <- data.table::fread(file.path(RAW, "mRNA_RSEM_UQ_log2_Tumor.cct"))
data.table::setnames(expression, 1, "gene")
sample_columns <- intersect(analysis_data$sample_id, names(expression))
expression[, (sample_columns) := lapply(.SD, as.numeric), .SDcols = sample_columns]
expression <- expression[, lapply(.SD, safe_median), by = gene, .SDcols = sample_columns]
expression_matrix <- as.matrix(expression[, ..sample_columns])
rownames(expression_matrix) <- expression$gene
storage.mode(expression_matrix) <- "numeric"

phenotype_definitions <- read_phenotype_gene_sets()
for (outcome_name in c("myCAF", "CAF_stroma", "endothelial")) {
  genes <- intersect(phenotype_definitions[[outcome_name]], rownames(expression_matrix))
  analysis_data[[outcome_name]] <- colMeans(zscore_rows_matrix(expression_matrix[genes, , drop = FALSE]), na.rm = TRUE)[match(analysis_data$sample_id, colnames(expression_matrix))]
}

for (layer in c("nglyco_peptide", "nglyco_site")) {
  cache <- readRDS(file.path(RESULTS, "analysis_results", sprintf("%s_analysis_layer_state.rds", layer)))
  score_table <- data.table::copy(cache$score_table[, .(sample_id, linear_score)])
  data.table::setnames(score_table, "linear_score", paste0(layer, "_primary_score"))
  analysis_data <- merge(analysis_data, score_table, by = "sample_id", all.x = TRUE)
}

analysis_data <- analysis_data[histology_group == "PDAC"]
analysis_data <- merge(analysis_data, sample_map[, .(sample_id, plex_id)], by = "sample_id", all.x = TRUE)

fit_model <- function(data, outcome, predictor, include_plex) {
  frame <- data.table::data.table(
    y = zscore_vector(data[[outcome]]),
    x = zscore_vector(data[[predictor]]),
    proteome = zscore_vector(data$proteome_matrisome_score),
    stroma = zscore_vector(data$stromal_fraction_z),
    purity = zscore_vector(data$neoplastic_cellularity_z),
    plex_id = factor(data$plex_id)
  )
  required <- c("y", "x", "proteome", "stroma", "purity", if (include_plex) "plex_id")
  frame <- frame[stats::complete.cases(frame[, ..required])]
  result <- data.table::data.table(n = nrow(frame), n_plex = if (include_plex) nlevels(frame$plex_id) else NA_integer_, beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, adj_r2 = NA_real_)
  if (nrow(frame) < 40 || (include_plex && nlevels(frame$plex_id) < 2)) {
    return(result)
  }
  formula <- if (include_plex) y ~ x + proteome + stroma + purity + plex_id else y ~ x + proteome + stroma + purity
  fit <- stats::lm(formula, data = frame)
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

predictors <- c(
  nglyco_peptide = "nglyco_peptide_primary_score",
  nglyco_site = "nglyco_site_primary_score"
)
outcomes <- c(myCAF = "myCAF", CAF_stroma = "CAF_stroma", endothelial = "endothelial")
records <- list()
for (layer in names(predictors)) {
  for (outcome_name in names(outcomes)) {
    for (model_name in c("base_adjusted", "plex_fixed_effect")) {
      model <- fit_model(analysis_data, outcomes[[outcome_name]], predictors[[layer]], model_name == "plex_fixed_effect")
      model[, `:=`(layer = layer, outcome = outcome_name, model = model_name)]
      records[[length(records) + 1]] <- model
    }
  }
}
models <- data.table::rbindlist(records, fill = TRUE)
models[, fdr := stats::p.adjust(p_value, method = "BH"), by = model]
data.table::setcolorder(models, c("layer", "outcome", "model", setdiff(names(models), c("layer", "outcome", "model"))))
data.table::fwrite(models, file.path(RESULTS, "cptac_batch_plex_sensitivity_models.tsv"), sep = "\t")
cat("Wrote", file.path(RESULTS, "cptac_batch_plex_sensitivity_models.tsv"), "\n")
