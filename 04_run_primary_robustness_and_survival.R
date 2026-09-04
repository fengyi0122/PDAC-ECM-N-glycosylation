script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/")) else file.path(getwd(), "scripts")
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table", "survival"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(survival))

set.seed(20260712)

OUT <- file.path(RESULTS, "analysis_results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

RAW <- existing_dataset_dir("cptac", "CPTAC-PDAC")

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- NA_real_
  sweep(sweep(mat, 1, means, "-"), 1, sds, "/")
}

safe_median <- function(x) {
  if (all(!is.finite(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

fit_effect <- function(data, outcome, predictor, covariates = c("proteome_matrisome_score", "stromal_fraction_z", "neoplastic_cellularity_z")) {
  cols <- unique(c(outcome, predictor, covariates))
  frame <- data[, ..cols]
  frame <- as.data.frame(frame)
  names(frame)[match(outcome, names(frame))] <- "outcome_value"
  names(frame)[match(predictor, names(frame))] <- "predictor_value"
  frame$outcome_value <- zscore(frame$outcome_value)
  frame$predictor_value <- zscore(frame$predictor_value)
  for (cov in covariates) {
    if (cov %in% names(frame) && is.numeric(frame[[cov]])) frame[[cov]] <- zscore(frame[[cov]])
  }
  frame <- frame[stats::complete.cases(frame), , drop = FALSE]
  if (nrow(frame) < 40) return(NULL)
  rhs <- c("predictor_value", covariates)
  rhs <- rhs[rhs %in% names(frame)]
  fit <- stats::lm(stats::as.formula(paste("outcome_value ~", paste(rhs, collapse = " + "))), data = frame)
  sm <- summary(fit)$coefficients
  if (!"predictor_value" %in% rownames(sm)) return(NULL)
  ci <- suppressMessages(stats::confint(fit, "predictor_value", level = 0.95))
  data.table::data.table(
    n = nrow(frame),
    beta = unname(sm["predictor_value", "Estimate"]),
    se = unname(sm["predictor_value", "Std. Error"]),
    ci_low = unname(ci[[1]]),
    ci_high = unname(ci[[2]]),
    p_value = unname(sm["predictor_value", "Pr(>|t|)"])
  )
}

parse_proteome <- function(sample_ids) {
  protein <- data.table::fread(file.path(RAW, "proteomics_gene_level_MD_abundance_tumor.cct"))
  data.table::setnames(protein, 1, "gene")
  sample_cols <- intersect(sample_ids, names(protein))
  protein[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
  protein <- protein[, lapply(.SD, safe_median), by = gene, .SDcols = sample_cols]
  mat <- as.matrix(protein[, ..sample_cols])
  rownames(mat) <- protein$gene
  storage.mode(mat) <- "numeric"
  mat
}

parse_mrna <- function(sample_ids) {
  expr <- data.table::fread(file.path(RAW, "mRNA_RSEM_UQ_log2_Tumor.cct"))
  data.table::setnames(expr, 1, "gene")
  sample_cols <- intersect(sample_ids, names(expr))
  expr[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
  expr <- expr[, lapply(.SD, safe_median), by = gene, .SDcols = sample_cols]
  mat <- as.matrix(expr[, ..sample_cols])
  rownames(mat) <- expr$gene
  storage.mode(mat) <- "numeric"
  mat
}

score_signature <- function(expr, genes) {
  found <- intersect(genes, rownames(expr))
  if (length(found) < 2) return(rep(NA_real_, ncol(expr)))
  colMeans(zscore_rows(expr[found, , drop = FALSE]), na.rm = TRUE)
}

extract_glycan_composition <- function(ids) {
  pattern <- "-N([0-9]+)H([0-9]+)F([0-9]+)S([0-9]+)G([0-9]+)$"
  hits <- regexec(pattern, ids)
  parts <- regmatches(ids, hits)
  get_part <- function(i) vapply(parts, function(x) if (length(x) == 6) as.numeric(x[[i]]) else NA_real_, numeric(1))
  data.table::data.table(N = get_part(2), H = get_part(3), F = get_part(4), S = get_part(5), G = get_part(6))
}

compute_ipw_scores <- function(residual_z, detect_prob, obs_n, feature_idx, sample_ids,
                               probability_lower = 0.05, probability_upper = 0.95,
                               weight_cap = 10, specification = "default") {
  scores <- rep(NA_real_, length(sample_ids))
  diagnostics <- vector("list", length(sample_ids))
  for (j in seq_along(sample_ids)) {
    idx <- feature_idx[is.finite(residual_z[feature_idx, j]) & is.finite(detect_prob[feature_idx, j])]
    if (length(idx) < 20) {
      diagnostics[[j]] <- data.table::data.table(
        sample_id = sample_ids[[j]], specification = specification,
        probability_lower = probability_lower, probability_upper = probability_upper,
        weight_cap = weight_cap, used_features = length(idx),
        effective_feature_count = NA_real_, mean_weight = NA_real_, median_weight = NA_real_,
        p95_weight = NA_real_, max_weight = NA_real_, probability_truncated_fraction = NA_real_,
        weight_capped_fraction = NA_real_
      )
      next
    }
    raw_probability <- detect_prob[idx, j]
    probability <- pmin(pmax(raw_probability, probability_lower), probability_upper)
    stabilizer <- obs_n[idx] / length(sample_ids)
    raw_weight <- stabilizer / probability
    weight <- pmin(raw_weight, weight_cap)
    scores[[j]] <- stats::weighted.mean(residual_z[idx, j], weight, na.rm = TRUE)
    diagnostics[[j]] <- data.table::data.table(
      sample_id = sample_ids[[j]], specification = specification,
      probability_lower = probability_lower, probability_upper = probability_upper,
      weight_cap = weight_cap, used_features = length(idx),
      effective_feature_count = sum(weight)^2 / sum(weight^2),
      mean_weight = mean(weight), median_weight = stats::median(weight),
      p95_weight = as.numeric(stats::quantile(weight, 0.95, names = FALSE)),
      max_weight = max(weight),
      probability_truncated_fraction = mean(raw_probability < probability_lower | raw_probability > probability_upper),
      weight_capped_fraction = mean(raw_weight > weight_cap)
    )
  }
  list(score = scores, diagnostics = data.table::rbindlist(diagnostics, fill = TRUE))
}

build_layer <- function(layer, filename, id_col, sample_ids, model_data, protein_mat, thresholds = c(40, 50, 70, 90, 110)) {
  ann <- data.table::fread(file.path(PROCESSED, sprintf("%s_matrisome_features_retained.tsv", layer)))
  ann <- unique(ann[category == "ECM Glycoproteins", .(feature_id, gene, division, category)], by = c("feature_id", "gene"))
  glyco <- data.table::fread(file.path(RAW, filename))
  data.table::setnames(glyco, id_col, "feature_id")
  if ("Gene" %in% names(glyco)) data.table::setnames(glyco, "Gene", "gene")
  sample_cols <- intersect(sample_ids, names(glyco))
  glyco[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
  merged <- merge(ann, glyco, by = c("feature_id", "gene"), all = FALSE)
  merged <- merged[gene %in% rownames(protein_mat)]
  values <- as.matrix(merged[, ..sample_cols])
  storage.mode(values) <- "numeric"
  protein_values <- t(vapply(merged$gene, function(g) protein_mat[g, sample_cols], numeric(length(sample_cols))))
  obs_n <- rowSums(is.finite(values) & is.finite(protein_values))
  keep <- obs_n >= min(thresholds)
  merged <- merged[keep]
  values <- values[keep, , drop = FALSE]
  protein_values <- protein_values[keep, , drop = FALSE]
  obs_n <- obs_n[keep]

  n_features <- nrow(values)
  residual_linear <- matrix(NA_real_, nrow = n_features, ncol = length(sample_cols), dimnames = list(NULL, sample_cols))
  residual_spline <- residual_linear
  detect_prob <- matrix(NA_real_, nrow = n_features, ncol = length(sample_cols), dimnames = list(NULL, sample_cols))
  detect_records <- vector("list", n_features)
  stroma <- model_data$stromal_fraction_z[match(sample_cols, model_data$sample_id)]
  purity <- model_data$neoplastic_cellularity_z[match(sample_cols, model_data$sample_id)]

  for (i in seq_len(n_features)) {
    y <- as.numeric(values[i, ])
    x <- as.numeric(protein_values[i, ])
    mask <- is.finite(y) & is.finite(x)
    residual_eligible <- sum(mask) >= min(thresholds) && is.finite(stats::sd(x[mask])) && stats::sd(x[mask]) > 0
    detect_record <- data.table::data.table(
      feature_id = merged$feature_id[[i]], gene = merged$gene[[i]], observed_n = obs_n[[i]],
      detection_rate = mean(is.finite(y)), residual_model_eligible = residual_eligible,
      detection_model_n = NA_integer_, detection_model_eligible = FALSE,
      detection_model_converged = FALSE, model_status = "residual_model_ineligible",
      protein_beta = NA_real_, stroma_beta = NA_real_, stroma_p = NA_real_,
      purity_beta = NA_real_, predicted_probability_complete_fraction = NA_real_
    )
    if (!residual_eligible) {
      detect_records[[i]] <- detect_record
      next
    }

    fit_linear <- stats::lm.fit(cbind(1, x[mask]), y[mask])
    pred_linear <- fit_linear$coefficients[[1]] + fit_linear$coefficients[[2]] * x
    residual_linear[i, mask] <- y[mask] - pred_linear[mask]

    spline_frame <- data.frame(y = y[mask], x = x[mask])
    fit_spline <- try(stats::lm(y ~ splines::ns(x, df = 3), data = spline_frame), silent = TRUE)
    if (!inherits(fit_spline, "try-error")) {
      pred_spline <- try(stats::predict(fit_spline, newdata = data.frame(x = x)), silent = TRUE)
      if (!inherits(pred_spline, "try-error")) residual_spline[i, mask] <- y[mask] - pred_spline[mask]
    }

    detected <- as.integer(is.finite(y))
    det_frame <- data.frame(detected = detected, protein = zscore(x), stroma = stroma, purity = purity)
    det_frame <- det_frame[stats::complete.cases(det_frame), , drop = FALSE]
    detection_eligible <- nrow(det_frame) >= 80 && length(unique(det_frame$detected)) == 2
    detect_record[, `:=`(
      detection_model_n = nrow(det_frame), detection_model_eligible = detection_eligible,
      model_status = if (detection_eligible) "fit_attempted" else "detection_model_ineligible"
    )]
    if (detection_eligible) {
      det_fit <- try(suppressWarnings(stats::glm(detected ~ protein + stroma + purity, family = stats::binomial(), data = det_frame)), silent = TRUE)
      if (!inherits(det_fit, "try-error") && isTRUE(det_fit$converged)) {
        coef_table <- summary(det_fit)$coefficients
        p_all <- try(stats::predict(det_fit, newdata = data.frame(protein = zscore(x), stroma = stroma, purity = purity), type = "response"), silent = TRUE)
        if (!inherits(p_all, "try-error")) detect_prob[i, ] <- as.numeric(p_all)
        detect_record[, `:=`(
          detection_model_converged = TRUE, model_status = "converged",
          protein_beta = if ("protein" %in% rownames(coef_table)) coef_table["protein", "Estimate"] else NA_real_,
          stroma_beta = if ("stroma" %in% rownames(coef_table)) coef_table["stroma", "Estimate"] else NA_real_,
          stroma_p = if ("stroma" %in% rownames(coef_table)) coef_table["stroma", "Pr(>|z|)"] else NA_real_,
          purity_beta = if ("purity" %in% rownames(coef_table)) coef_table["purity", "Estimate"] else NA_real_,
          predicted_probability_complete_fraction = if (!inherits(p_all, "try-error")) mean(is.finite(p_all)) else 0
        )]
      } else {
        detect_record[, model_status := "fit_failed_or_not_converged"]
      }
    }
    detect_records[[i]] <- detect_record
  }

  linear_z <- zscore_rows(residual_linear)
  spline_z <- zscore_rows(residual_spline)
  merged[, observed_n := obs_n]
  comp <- extract_glycan_composition(merged$feature_id)
  merged <- cbind(merged, comp)
  merged[, `:=`(
    fucosylated_composition = is.finite(F) & F > 0,
    sialylated_composition = (is.finite(S) & S > 0) | (is.finite(G) & G > 0),
    high_mannose_composition = is.finite(N) & N == 2 & is.finite(H) & H >= 5 & H <= 9 &
      is.finite(F) & F == 0 & is.finite(S) & S == 0 & is.finite(G) & G == 0
  )]

  default_idx <- which(obs_n >= 50)
  linear_score <- colMeans(linear_z[default_idx, , drop = FALSE], na.rm = TRUE)
  spline_score <- colMeans(spline_z[default_idx, , drop = FALSE], na.rm = TRUE)

  default_ipw <- compute_ipw_scores(
    linear_z, detect_prob, obs_n, default_idx, sample_cols,
    probability_lower = 0.05, probability_upper = 0.95, weight_cap = 10,
    specification = "p05_p95_cap10"
  )
  ipw_score <- default_ipw$score

  score_table <- data.table::data.table(
    sample_id = sample_cols,
    linear_score = linear_score,
    spline_score = spline_score,
    ipw_score = ipw_score,
    observed_feature_count = colSums(is.finite(values[default_idx, , drop = FALSE]))
  )

  class_definitions <- list(
    all_ecm = default_idx,
    fucosylated = default_idx[merged$fucosylated_composition[default_idx]],
    sialylated = default_idx[merged$sialylated_composition[default_idx]],
    high_mannose = default_idx[merged$high_mannose_composition[default_idx]],
    nonfucosylated = default_idx[!merged$fucosylated_composition[default_idx]],
    nonsialylated = default_idx[!merged$sialylated_composition[default_idx]]
  )
  class_scores <- data.table::data.table(sample_id = sample_cols)
  class_manifest <- list()
  for (class_name in names(class_definitions)) {
    idx <- class_definitions[[class_name]]
    if (length(idx) < 10) next
    class_scores[[class_name]] <- colMeans(linear_z[idx, , drop = FALSE], na.rm = TRUE)
    class_manifest[[length(class_manifest) + 1]] <- data.table::data.table(layer = layer, glycan_class = class_name, n_features = length(idx), n_genes = data.table::uniqueN(merged$gene[idx]))
  }

  threshold_scores <- list()
  for (threshold in thresholds) {
    idx <- which(obs_n >= threshold)
    if (length(idx) < 10) next
    threshold_scores[[length(threshold_scores) + 1]] <- data.table::data.table(
      sample_id = sample_cols, layer = layer, threshold = threshold,
      n_features = length(idx), score = colMeans(linear_z[idx, , drop = FALSE], na.rm = TRUE)
    )
  }

  detect_dt <- data.table::rbindlist(detect_records, fill = TRUE)
  if (nrow(detect_dt) > 0) detect_dt[, fdr := stats::p.adjust(stroma_p, method = "BH")]

  list(
    layer = layer, samples = sample_cols, annotation = merged, residual_linear_z = linear_z,
    score_table = score_table, class_scores = class_scores,
    class_manifest = data.table::rbindlist(class_manifest, fill = TRUE),
    threshold_scores = data.table::rbindlist(threshold_scores, fill = TRUE),
    detection_models = detect_dt, detection_probability = detect_prob,
    default_ipw_diagnostics = default_ipw$diagnostics
  )
}

clinical <- data.table::fread(file.path(RAW, "clinical_table_140.tsv"))
clinical[, sample_id := case_id]
clinical[, histology_group := data.table::fifelse(histology_diagnosis == "PDAC", "PDAC", "Adenosquamous")]
clinical[, follow_up_days_num := suppressWarnings(as.numeric(follow_up_days))]
clinical[, event := as.integer(tolower(vital_status) == "deceased")]
clinical[, age_num := suppressWarnings(as.numeric(age))]
clinical[, stage_group := data.table::fifelse(grepl("Stage I$|Stage IA|Stage IB|Stage IIA|Stage IIB", tumor_stage_pathological), "I-II",
                                               data.table::fifelse(grepl("Stage III|Stage IV", tumor_stage_pathological), "III-IV", NA_character_))]

phen <- data.table::fread(file.path(PROCESSED, "immune_phenotypes.tsv"))
glyco_scores <- data.table::fread(file.path(PROCESSED, "glyco_matrisome_features.tsv"))
model_data <- merge(phen, glyco_scores, by = "sample_id", all = FALSE)
model_data <- merge(model_data, clinical[, .(sample_id, histology_group, follow_up_days_num, event, age_num, sex, stage_group)], by = "sample_id", all.x = TRUE)
pdac_ids <- model_data[histology_group == "PDAC", sample_id]

expr <- parse_mrna(model_data$sample_id)
signature_definitions <- list(
  myCAF = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "POSTN"),
  myCAF_no_POSTN = c("ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2"),
  iCAF = c("IL6", "CXCL12", "CFD", "HAS1", "CXCL14", "PDGFRA"),
  apCAF = c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"),
  CD8_T_cell = c("CD3D", "CD3E", "CD8A", "CD8B"),
  cytotoxicity = c("NKG7", "GNLY", "GZMB", "PRF1"),
  TAM_M2 = c("CD68", "CD163", "MRC1", "MSR1"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
  CAF_stroma = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "ESAM")
)
for (name in names(signature_definitions)) {
  model_data[[name]] <- score_signature(expr, signature_definitions[[name]])[match(model_data$sample_id, colnames(expr))]
}

marker_table <- data.table::rbindlist(lapply(names(signature_definitions), function(name) {
  genes <- signature_definitions[[name]]
  data.table::data.table(signature = name, gene = genes, present_in_cptac_rna = genes %in% rownames(expr))
}))
data.table::fwrite(marker_table, file.path(OUT, "signature_gene_definitions.tsv"), sep = "\t")

protein_mat <- parse_proteome(pdac_ids)
pdac_model_data <- model_data[sample_id %in% pdac_ids]

layer_specs <- list(
  list(layer = "nglyco_peptide", filename = "N-glycoproteomics_peptide_level_ratio_tumor.cct", id_col = "Sequence"),
  list(layer = "nglyco_site", filename = "N-glycoproteomics_Site_level_ratio_tumor.cct", id_col = "Modifications")
)

use_cache <- identical(Sys.getenv("PDAC_USE_CACHE", unset = "0"), "1")
layers <- lapply(layer_specs, function(spec) {
  cache_file <- file.path(OUT, sprintf("%s_analysis_layer_state.rds", spec$layer))
  if (use_cache && file.exists(cache_file)) {
    message("Reading analysis layer state: ", spec$layer)
    return(readRDS(cache_file))
  }
  message("Building analysis layer state: ", spec$layer)
  built <- build_layer(spec$layer, spec$filename, spec$id_col, pdac_ids, pdac_model_data, protein_mat)
  saveRDS(built, cache_file)
  built
})
names(layers) <- vapply(layer_specs, `[[`, character(1), "layer")

for (layer_name in names(layers)) {
  score_cols <- copy(layers[[layer_name]]$score_table)
  data.table::setnames(score_cols, setdiff(names(score_cols), "sample_id"), paste0(layer_name, "_", setdiff(names(score_cols), "sample_id")))
  pdac_model_data <- merge(pdac_model_data, score_cols, by = "sample_id", all.x = TRUE)
}

figure2_state <- pdac_model_data[, .(
  sample_id,
  nglyco_peptide_residual = nglyco_peptide_linear_score,
  nglyco_site_residual = nglyco_site_linear_score,
  myCAF,
  CAF_stroma,
  iCAF,
  apCAF,
  endothelial,
  ecm_proteome = proteome_matrisome_score,
  stromal_fraction = stromal_fraction_z,
  neoplastic_cellularity = neoplastic_cellularity_z
)]
figure2_adjusted_cols <- c(
  "nglyco_peptide_residual", "nglyco_site_residual", "myCAF", "CAF_stroma",
  "iCAF", "apCAF", "endothelial"
)
figure2_covariate_cols <- c("ecm_proteome", "stromal_fraction", "neoplastic_cellularity")
for (state_column in figure2_adjusted_cols) {
  state_frame <- as.data.frame(figure2_state[, c(state_column, figure2_covariate_cols), with = FALSE])
  names(state_frame)[[1]] <- "state_value"
  complete <- stats::complete.cases(state_frame)
  adjusted_state <- rep(NA_real_, nrow(figure2_state))
  state_fit <- stats::lm(state_value ~ ecm_proteome + stromal_fraction + neoplastic_cellularity, data = state_frame[complete, , drop = FALSE])
  adjusted_state[complete] <- stats::residuals(state_fit)
  figure2_state[[state_column]] <- zscore(adjusted_state)
}
figure2_state[, (figure2_covariate_cols) := lapply(.SD, zscore), .SDcols = figure2_covariate_cols]
figure2_state[, residual_average := rowMeans(.SD, na.rm = TRUE), .SDcols = c("nglyco_peptide_residual", "nglyco_site_residual")]
data.table::setorder(figure2_state, residual_average)
figure2_state[, patient_order := seq_len(.N)]
data.table::fwrite(figure2_state, file.path(OUT, "figure2_patient_state_matrix.tsv"), sep = "\t")

outcomes <- c("myCAF", "myCAF_no_POSTN", "endothelial", "CAF_stroma", "stromal_immune_composite")
main_records <- list()
for (cohort in c("PDAC_135", "All_pancreatic_tumors_140")) {
  data_use <- if (cohort == "PDAC_135") pdac_model_data else model_data
  predictors <- if (cohort == "PDAC_135") {
    c(nglyco_peptide = "nglyco_peptide_linear_score", nglyco_site = "nglyco_site_linear_score")
  } else {
    c(
      nglyco_peptide = "nglyco_peptide_protein_residual_category_ecm_glycoproteins_score",
      nglyco_site = "nglyco_site_protein_residual_category_ecm_glycoproteins_score"
    )
  }
  for (layer_name in names(predictors)) {
    for (outcome in outcomes) {
      effect <- fit_effect(data_use, outcome, predictors[[layer_name]])
      if (is.null(effect)) next
      effect[, `:=`(cohort = cohort, layer = layer_name, outcome = outcome)]
      main_records[[length(main_records) + 1]] <- effect
    }
  }
}
main_models <- data.table::rbindlist(main_records, fill = TRUE)
main_models[, fdr := stats::p.adjust(p_value, method = "BH"), by = cohort]
data.table::fwrite(main_models, file.path(OUT, "primary_pdac_and_140_sensitivity_models.tsv"), sep = "\t")

threshold_records <- list()
method_records <- list()
detection_summary <- list()
class_records <- list()
feature_records <- list()
split_records <- list()
ipw_parameter_records <- list()
ipw_sample_diagnostic_records <- list()
detection_model_summaries <- list()

for (layer_name in names(layers)) {
  layer_obj <- layers[[layer_name]]
  threshold_long <- layer_obj$threshold_scores
  for (threshold_value in unique(threshold_long$threshold)) {
    tmp <- merge(pdac_model_data, threshold_long[threshold == threshold_value, .(sample_id, score, n_features)], by = "sample_id")
    for (outcome in outcomes) {
      effect <- fit_effect(tmp, outcome, "score")
      if (is.null(effect)) next
      effect[, `:=`(layer = layer_name, threshold = threshold_value, n_features = unique(tmp$n_features), outcome = outcome)]
      threshold_records[[length(threshold_records) + 1]] <- effect
    }
  }

  methods <- c(linear = paste0(layer_name, "_linear_score"), spline = paste0(layer_name, "_spline_score"), detection_IPW = paste0(layer_name, "_ipw_score"))
  for (method_name in names(methods)) {
    for (outcome in outcomes) {
      effect <- fit_effect(pdac_model_data, outcome, methods[[method_name]])
      if (is.null(effect)) next
      effect[, `:=`(layer = layer_name, method = method_name, outcome = outcome)]
      method_records[[length(method_records) + 1]] <- effect
    }
  }

  count_col <- paste0(layer_name, "_observed_feature_count")
  count_effect <- fit_effect(pdac_model_data, "stromal_fraction_z", count_col, covariates = c("proteome_matrisome_score", "neoplastic_cellularity_z"))
  if (!is.null(count_effect)) {
    count_effect[, `:=`(layer = layer_name, model = "observed_count_vs_stroma")]
    detection_summary[[length(detection_summary) + 1]] <- count_effect
  }
  detect_dt <- copy(layer_obj$detection_models)
  detect_dt[, layer := layer_name]
  data.table::fwrite(detect_dt, file.path(OUT, sprintf("%s_feature_detection_models.tsv", layer_name)), sep = "\t")

  converged_detection <- detect_dt[detection_model_converged == TRUE]
  detection_model_summaries[[length(detection_model_summaries) + 1]] <- data.table::data.table(
    layer = layer_name,
    features_in_ge40_robustness_universe = nrow(detect_dt),
    detection_model_eligible_features = sum(detect_dt$detection_model_eligible, na.rm = TRUE),
    converged_detection_models = nrow(converged_detection),
    median_detection_rate = stats::median(detect_dt$detection_rate, na.rm = TRUE),
    detection_rate_q1 = as.numeric(stats::quantile(detect_dt$detection_rate, 0.25, na.rm = TRUE, names = FALSE)),
    detection_rate_q3 = as.numeric(stats::quantile(detect_dt$detection_rate, 0.75, na.rm = TRUE, names = FALSE)),
    median_stroma_beta = stats::median(converged_detection$stroma_beta, na.rm = TRUE),
    stroma_beta_q1 = as.numeric(stats::quantile(converged_detection$stroma_beta, 0.25, na.rm = TRUE, names = FALSE)),
    stroma_beta_q3 = as.numeric(stats::quantile(converged_detection$stroma_beta, 0.75, na.rm = TRUE, names = FALSE)),
    positive_stroma_beta_features = sum(converged_detection$stroma_beta > 0, na.rm = TRUE),
    negative_stroma_beta_features = sum(converged_detection$stroma_beta < 0, na.rm = TRUE),
    stroma_fdr_below_0_05_features = sum(converged_detection$fdr < 0.05, na.rm = TRUE)
  )

  ipw_specifications <- data.table::data.table(
    specification = c("p05_p95_cap10", "p02_p98_cap10", "p05_p95_cap5", "p05_p95_cap20"),
    probability_lower = c(0.05, 0.02, 0.05, 0.05),
    probability_upper = c(0.95, 0.98, 0.95, 0.95),
    weight_cap = c(10, 10, 5, 20)
  )
  default_idx <- which(layer_obj$annotation$observed_n >= 50)
  for (specification_index in seq_len(nrow(ipw_specifications))) {
    specification_row <- ipw_specifications[specification_index]
    ipw_result <- compute_ipw_scores(
      layer_obj$residual_linear_z, layer_obj$detection_probability,
      layer_obj$annotation$observed_n, default_idx, layer_obj$samples,
      probability_lower = specification_row$probability_lower,
      probability_upper = specification_row$probability_upper,
      weight_cap = specification_row$weight_cap,
      specification = specification_row$specification
    )
    ipw_diagnostics <- copy(ipw_result$diagnostics)
    ipw_diagnostics[, layer := layer_name]
    ipw_sample_diagnostic_records[[length(ipw_sample_diagnostic_records) + 1]] <- ipw_diagnostics
    ipw_data <- merge(
      pdac_model_data,
      data.table::data.table(sample_id = layer_obj$samples, ipw_sensitivity_score = ipw_result$score),
      by = "sample_id", all.x = TRUE
    )
    for (outcome in c("myCAF", "CAF_stroma")) {
      effect <- fit_effect(ipw_data, outcome, "ipw_sensitivity_score")
      if (is.null(effect)) next
      effect[, `:=`(
        layer = layer_name, outcome = outcome,
        specification = specification_row$specification,
        probability_lower = specification_row$probability_lower,
        probability_upper = specification_row$probability_upper,
        weight_cap = specification_row$weight_cap
      )]
      ipw_parameter_records[[length(ipw_parameter_records) + 1]] <- effect
    }
  }

  class_scores <- copy(layer_obj$class_scores)
  for (class_name in setdiff(names(class_scores), "sample_id")) {
    tmp <- merge(pdac_model_data, class_scores[, .(sample_id, class_score = get(class_name))], by = "sample_id")
    for (outcome in outcomes) {
      effect <- fit_effect(tmp, outcome, "class_score")
      if (is.null(effect)) next
      effect[, `:=`(layer = layer_name, glycan_class = class_name, outcome = outcome)]
      class_records[[length(class_records) + 1]] <- effect
    }
  }

  residual_z <- layer_obj$residual_linear_z
  ann <- layer_obj$annotation
  default_idx <- which(ann$observed_n >= 50)
  for (outcome in outcomes) {
    y <- pdac_model_data[[outcome]][match(layer_obj$samples, pdac_model_data$sample_id)]
    cov <- pdac_model_data[match(layer_obj$samples, sample_id), .(proteome_matrisome_score, stromal_fraction_z, neoplastic_cellularity_z)]
    for (i in default_idx) {
      frame <- data.table::data.table(y = y, x = residual_z[i, ], pm = cov$proteome_matrisome_score, stroma = cov$stromal_fraction_z, purity = cov$neoplastic_cellularity_z)
      frame <- frame[stats::complete.cases(frame)]
      if (nrow(frame) < 40 || stats::sd(frame$x) == 0) next
      frame[, `:=`(y = zscore(y), x = zscore(x), pm = zscore(pm), stroma = zscore(stroma), purity = zscore(purity))]
      fit <- stats::lm(y ~ x + pm + stroma + purity, data = frame)
      sm <- summary(fit)$coefficients
      feature_records[[length(feature_records) + 1]] <- data.table::data.table(
        layer = layer_name, outcome = outcome, feature_id = ann$feature_id[[i]], gene = ann$gene[[i]],
        glycan_N = ann$N[[i]], glycan_H = ann$H[[i]], glycan_F = ann$F[[i]], glycan_S = ann$S[[i]], glycan_G = ann$G[[i]],
        n = nrow(frame), beta = sm["x", "Estimate"], p_value = sm["x", "Pr(>|t|)"]
      )
    }

    for (replicate in seq_len(100)) {
      shuffled <- sample(default_idx)
      half <- floor(length(shuffled) / 2)
      idx_a <- shuffled[seq_len(half)]
      idx_b <- shuffled[(half + 1):length(shuffled)]
      score_a <- colMeans(residual_z[idx_a, , drop = FALSE], na.rm = TRUE)
      score_b <- colMeans(residual_z[idx_b, , drop = FALSE], na.rm = TRUE)
      score_cor <- suppressWarnings(stats::cor(score_a, score_b, method = "spearman", use = "complete.obs"))
      for (half_name in c("A", "B")) {
        score <- if (half_name == "A") score_a else score_b
        tmp <- copy(pdac_model_data)
        tmp$split_score <- score[match(tmp$sample_id, layer_obj$samples)]
        effect <- fit_effect(tmp, outcome, "split_score")
        if (is.null(effect)) next
        effect[, `:=`(layer = layer_name, outcome = outcome, replicate = replicate, feature_half = half_name, split_score_rho = score_cor)]
        split_records[[length(split_records) + 1]] <- effect
      }
    }
  }
  data.table::fwrite(layer_obj$annotation, file.path(OUT, sprintf("%s_glycan_composition_manifest.tsv", layer_name)), sep = "\t")
  data.table::fwrite(layer_obj$class_manifest, file.path(OUT, sprintf("%s_glycan_class_manifest.tsv", layer_name)), sep = "\t")
}

threshold_models <- data.table::rbindlist(threshold_records, fill = TRUE)
threshold_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(threshold_models, file.path(OUT, "observation_threshold_sensitivity.tsv"), sep = "\t")

method_models <- data.table::rbindlist(method_records, fill = TRUE)
method_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(method_models, file.path(OUT, "residual_method_and_detection_ipw_sensitivity.tsv"), sep = "\t")

detection_summary_dt <- data.table::rbindlist(detection_summary, fill = TRUE)
detection_summary_dt[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(detection_summary_dt, file.path(OUT, "sample_detection_burden_models.tsv"), sep = "\t")

detection_model_summary_dt <- data.table::rbindlist(detection_model_summaries, fill = TRUE)
data.table::fwrite(detection_model_summary_dt, file.path(OUT, "detection_model_summary.tsv"), sep = "\t")

ipw_parameter_models <- data.table::rbindlist(ipw_parameter_records, fill = TRUE)
ipw_parameter_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(ipw_parameter_models, file.path(OUT, "detection_ipw_parameter_sensitivity.tsv"), sep = "\t")

ipw_sample_diagnostics <- data.table::rbindlist(ipw_sample_diagnostic_records, fill = TRUE)
data.table::fwrite(ipw_sample_diagnostics, file.path(OUT, "detection_ipw_sample_diagnostics.tsv"), sep = "\t")
ipw_weight_summary <- ipw_sample_diagnostics[, .(
  n_samples = .N,
  median_used_features = stats::median(used_features, na.rm = TRUE),
  median_effective_feature_count = stats::median(effective_feature_count, na.rm = TRUE),
  median_mean_weight = stats::median(mean_weight, na.rm = TRUE),
  median_p95_weight = stats::median(p95_weight, na.rm = TRUE),
  maximum_weight = max(max_weight, na.rm = TRUE),
  median_probability_truncated_fraction = stats::median(probability_truncated_fraction, na.rm = TRUE),
  median_weight_capped_fraction = stats::median(weight_capped_fraction, na.rm = TRUE)
), by = .(layer, specification, probability_lower, probability_upper, weight_cap)]
data.table::fwrite(ipw_weight_summary, file.path(OUT, "detection_ipw_weight_summary.tsv"), sep = "\t")

class_models <- data.table::rbindlist(class_records, fill = TRUE)
class_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(class_models, file.path(OUT, "glycan_composition_class_models.tsv"), sep = "\t")

feature_models <- data.table::rbindlist(feature_records, fill = TRUE)
feature_models[, fdr := stats::p.adjust(p_value, method = "BH"), by = .(layer, outcome)]
data.table::fwrite(feature_models, file.path(OUT, "glycosite_feature_coherence_models.tsv"), sep = "\t")

coherence_summary <- feature_models[, .(
  n_features = .N,
  positive_fraction = mean(beta > 0, na.rm = TRUE),
  median_beta = stats::median(beta, na.rm = TRUE),
  fdr_significant_fraction = mean(fdr < 0.05, na.rm = TRUE)
), by = .(layer, outcome)]
data.table::fwrite(coherence_summary, file.path(OUT, "glycosite_feature_coherence_summary.tsv"), sep = "\t")

gene_effects <- feature_models[, .(gene_beta = stats::median(beta, na.rm = TRUE), n_features = .N), by = .(layer, outcome, gene)]
gene_wide <- data.table::dcast(gene_effects, outcome + gene ~ layer, value.var = "gene_beta")
gene_concordance <- gene_wide[stats::complete.cases(gene_wide), .(
  n_genes = .N,
  spearman_rho = suppressWarnings(stats::cor(nglyco_peptide, nglyco_site, method = "spearman")),
  sign_concordance = mean(sign(nglyco_peptide) == sign(nglyco_site))
), by = outcome]
data.table::fwrite(gene_effects, file.path(OUT, "gene_level_cross_layer_effects.tsv"), sep = "\t")
data.table::fwrite(gene_concordance, file.path(OUT, "gene_level_cross_layer_concordance.tsv"), sep = "\t")

split_models <- data.table::rbindlist(split_records, fill = TRUE)
data.table::fwrite(split_models, file.path(OUT, "feature_split_half_models.tsv"), sep = "\t")
split_summary <- split_models[, .(
  median_split_score_rho = stats::median(split_score_rho, na.rm = TRUE),
  median_beta = stats::median(beta, na.rm = TRUE),
  beta_positive_fraction = mean(beta > 0, na.rm = TRUE),
  p_below_0_05_fraction = mean(p_value < 0.05, na.rm = TRUE)
), by = .(layer, outcome)]
data.table::fwrite(split_summary, file.path(OUT, "feature_split_half_summary.tsv"), sep = "\t")

caf_records <- list()
for (layer_name in names(layers)) {
  predictor <- paste0(layer_name, "_linear_score")
  for (caf_subtype in c("myCAF", "iCAF", "apCAF", "myCAF_no_POSTN")) {
    effect <- fit_effect(pdac_model_data, caf_subtype, predictor)
    if (is.null(effect)) next
    effect[, `:=`(layer = layer_name, caf_subtype = caf_subtype)]
    caf_records[[length(caf_records) + 1]] <- effect
  }
}
caf_models <- data.table::rbindlist(caf_records, fill = TRUE)
caf_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(caf_models, file.path(OUT, "caf_subtype_models.tsv"), sep = "\t")

postn_overlap_sensitivity <- main_models[
  cohort == "PDAC_135" & outcome %in% c("myCAF", "myCAF_no_POSTN")
]
postn_overlap_sensitivity[, sensitivity_fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(postn_overlap_sensitivity, file.path(OUT, "postn_overlap_sensitivity_models.tsv"), sep = "\t")

postn_predictor_records <- list()
postn_predictor_manifest <- list()
for (layer_name in names(layers)) {
  layer_obj <- layers[[layer_name]]
  full_idx <- which(layer_obj$annotation$observed_n >= 50)
  postn_idx <- full_idx[toupper(layer_obj$annotation$gene[full_idx]) == "POSTN"]
  postn_removed_idx <- setdiff(full_idx, postn_idx)
  predictor_scores <- list(
    full_residual_score = colMeans(layer_obj$residual_linear_z[full_idx, , drop = FALSE], na.rm = TRUE),
    postn_feature_removed_residual_score = colMeans(layer_obj$residual_linear_z[postn_removed_idx, , drop = FALSE], na.rm = TRUE)
  )
  postn_predictor_manifest[[length(postn_predictor_manifest) + 1]] <- data.table::data.table(
    layer = layer_name, full_feature_count = length(full_idx),
    postn_feature_count = length(postn_idx), postn_removed_feature_count = length(postn_removed_idx),
    postn_carrier_feature_fraction = length(postn_idx) / length(full_idx)
  )
  for (predictor_variant in names(predictor_scores)) {
    sensitivity_data <- merge(
      pdac_model_data,
      data.table::data.table(sample_id = layer_obj$samples, postn_sensitivity_score = predictor_scores[[predictor_variant]]),
      by = "sample_id", all.x = TRUE
    )
    for (outcome in c("myCAF", "myCAF_no_POSTN", "CAF_stroma")) {
      effect <- fit_effect(sensitivity_data, outcome, "postn_sensitivity_score")
      if (is.null(effect)) next
      effect[, `:=`(layer = layer_name, predictor_variant = predictor_variant, outcome = outcome)]
      postn_predictor_records[[length(postn_predictor_records) + 1]] <- effect
    }
  }
}
postn_predictor_models <- data.table::rbindlist(postn_predictor_records, fill = TRUE)
postn_predictor_models[, fdr := stats::p.adjust(p_value, method = "BH")]
data.table::fwrite(postn_predictor_models, file.path(OUT, "postn_feature_removed_predictor_models.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(postn_predictor_manifest, fill = TRUE), file.path(OUT, "postn_feature_removed_predictor_manifest.tsv"), sep = "\t")

survival_records <- list()
km_records <- list()
km_curve_records <- list()
km_risk_records <- list()
for (layer_name in names(layers)) {
  predictor <- paste0(layer_name, "_linear_score")
  count_col <- paste0(layer_name, "_observed_feature_count")
  surv_data <- copy(pdac_model_data)
  surv_data[, predictor_z := zscore(get(predictor))]
  surv_data[, age_z := zscore(age_num)]
  surv_data[, proteome_z := zscore(proteome_matrisome_score)]
  surv_data[, detection_z := zscore(get(count_col))]
  surv_data[, sex_factor := factor(sex)]
  surv_data[, stage_factor := factor(stage_group, levels = c("I-II", "III-IV"))]
  formulas <- list(
    unadjusted = survival::Surv(follow_up_days_num, event) ~ predictor_z,
    clinical_adjusted = survival::Surv(follow_up_days_num, event) ~ predictor_z + age_z + sex_factor + stage_factor,
    clinical_proteome_detection_adjusted = survival::Surv(follow_up_days_num, event) ~ predictor_z + age_z + sex_factor + stage_factor + proteome_z + detection_z
  )
  for (model_name in names(formulas)) {
    fit <- try(survival::coxph(formulas[[model_name]], data = surv_data, ties = "efron"), silent = TRUE)
    if (inherits(fit, "try-error")) next
    sm <- summary(fit)
    idx <- which(rownames(sm$coefficients) == "predictor_z")
    if (length(idx) != 1) next
    ph <- try(survival::cox.zph(fit), silent = TRUE)
    ph_global <- if (!inherits(ph, "try-error")) ph$table[nrow(ph$table), "p"] else NA_real_
    survival_records[[length(survival_records) + 1]] <- data.table::data.table(
      layer = layer_name, model = model_name, n = fit$n, events = fit$nevent,
      hazard_ratio = sm$conf.int[idx, "exp(coef)"], ci_low = sm$conf.int[idx, "lower .95"],
      ci_high = sm$conf.int[idx, "upper .95"], p_value = sm$coefficients[idx, "Pr(>|z|)"], ph_global_p = ph_global
    )
  }
  km_frame <- surv_data[is.finite(predictor_z) & is.finite(follow_up_days_num) & !is.na(event)]
  cutpoint <- stats::median(km_frame$predictor_z, na.rm = TRUE)
  km_frame[, group := factor(ifelse(predictor_z >= cutpoint, "High", "Low"), levels = c("Low", "High"))]
  logrank <- survival::survdiff(survival::Surv(follow_up_days_num, event) ~ group, data = km_frame)
  group_fit <- survival::coxph(survival::Surv(follow_up_days_num, event) ~ group, data = km_frame, ties = "efron")
  group_sm <- summary(group_fit)
  group_idx <- which(rownames(group_sm$coefficients) == "groupHigh")
  km_fit <- survival::survfit(survival::Surv(follow_up_days_num, event) ~ group, data = km_frame, conf.type = "log")
  km_sm <- summary(km_fit, censored = TRUE)
  km_curve <- data.table::data.table(
    layer = layer_name,
    group = sub("^group=", "", as.character(km_sm$strata)),
    time_days = km_sm$time,
    n_risk = km_sm$n.risk,
    n_event = km_sm$n.event,
    n_censor = km_sm$n.censor,
    survival = km_sm$surv,
    std_error = km_sm$std.err,
    ci_low = km_sm$lower,
    ci_high = km_sm$upper
  )
  km_initial <- km_frame[, .(n_risk = .N), by = group]
  km_initial[, `:=`(
    layer = layer_name,
    time_days = 0,
    n_event = 0,
    n_censor = 0,
    survival = 1,
    std_error = 0,
    ci_low = 1,
    ci_high = 1
  )]
  km_curve_records[[length(km_curve_records) + 1]] <- data.table::rbindlist(list(km_initial, km_curve), use.names = TRUE, fill = TRUE)
  risk_times <- seq(0, 4 * 365.25, by = 365.25)
  risk_sm <- summary(km_fit, times = risk_times, extend = TRUE)
  km_risk_records[[length(km_risk_records) + 1]] <- data.table::data.table(
    layer = layer_name,
    group = sub("^group=", "", as.character(risk_sm$strata)),
    time_days = risk_sm$time,
    time_years = risk_sm$time / 365.25,
    n_risk = risk_sm$n.risk
  )
  km_records[[length(km_records) + 1]] <- data.table::data.table(
    layer = layer_name, n = nrow(km_frame), events = sum(km_frame$event),
    high_n = sum(km_frame$group == "High"), low_n = sum(km_frame$group == "Low"),
    cutpoint_method = "median", cutpoint_z = cutpoint,
    grouped_hazard_ratio = group_sm$conf.int[group_idx, "exp(coef)"],
    grouped_ci_low = group_sm$conf.int[group_idx, "lower .95"],
    grouped_ci_high = group_sm$conf.int[group_idx, "upper .95"],
    grouped_cox_p = group_sm$coefficients[group_idx, "Pr(>|z|)"],
    logrank_p = stats::pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
  )
}
survival_models <- data.table::rbindlist(survival_records, fill = TRUE)
survival_models[, fdr := stats::p.adjust(p_value, method = "BH"), by = model]
km_summary <- data.table::rbindlist(km_records, fill = TRUE)
km_summary[, logrank_fdr := stats::p.adjust(logrank_p, method = "BH")]
data.table::fwrite(survival_models, file.path(OUT, "cptac_pdac_survival_models.tsv"), sep = "\t")
data.table::fwrite(km_summary, file.path(OUT, "cptac_pdac_km_summary.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(km_curve_records, fill = TRUE), file.path(OUT, "cptac_pdac_km_curve.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(km_risk_records, fill = TRUE), file.path(OUT, "cptac_pdac_km_risk_table.tsv"), sep = "\t")

decision_lines <- c(
  "PRIMARY ANALYSIS GATE",
  sprintf("PDAC-only primary cohort: n=%d; deaths=%d; survival records=%d", length(pdac_ids), clinical[histology_group == "PDAC", sum(event, na.rm = TRUE)], clinical[histology_group == "PDAC", sum(is.finite(follow_up_days_num))]),
  sprintf("Primary PDAC associations positive: %d/%d", main_models[cohort == "PDAC_135", sum(beta > 0, na.rm = TRUE)], main_models[cohort == "PDAC_135", .N]),
  sprintf("Primary PDAC associations FDR<0.05: %d/%d", main_models[cohort == "PDAC_135", sum(fdr < 0.05, na.rm = TRUE)], main_models[cohort == "PDAC_135", .N]),
  sprintf("Detection-IPW associations positive: %d/%d", method_models[method == "detection_IPW", sum(beta > 0, na.rm = TRUE)], method_models[method == "detection_IPW", .N]),
  sprintf("Detection-IPW associations FDR<0.05: %d/%d", method_models[method == "detection_IPW", sum(fdr < 0.05, na.rm = TRUE)], method_models[method == "detection_IPW", .N]),
  sprintf("Clinically adjusted survival associations p<0.05: %d/%d", survival_models[model == "clinical_adjusted", sum(p_value < 0.05, na.rm = TRUE)], survival_models[model == "clinical_adjusted", .N])
)
writeLines(decision_lines, file.path(OUT, "PRIMARY_ANALYSIS_GATE.txt"))
cat(paste(decision_lines, collapse = "\n"), "\n")
