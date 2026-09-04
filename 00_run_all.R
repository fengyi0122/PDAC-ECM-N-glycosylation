script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table", "survival", "readxl", "Matrix", "FNN", "CellChat", "digest", "SeuratObject"))

if (as.character(getRversion()) != EXPECTED_R_VERSION) {
  stop(sprintf("R version mismatch: found %s, expected %s.", as.character(getRversion()), EXPECTED_R_VERSION), call. = FALSE)
}

manifest <- data.table::fread(file.path(script_dir, "input_files.tsv"))
required_manifest_columns <- c("kind", "dataset", "relative_path", "source_url", "bytes", "sha256")
if (!all(required_manifest_columns %in% names(manifest))) {
  stop("The input manifest is missing required identity fields.", call. = FALSE)
}
resolve_input <- function(kind, dataset, relative_path) {
  if (kind == "bundled") {
    return(file.path(script_dir, relative_path))
  }
  file.path(existing_dataset_dir(kind, dataset), relative_path)
}
manifest[, resolved_path := mapply(resolve_input, kind, dataset, relative_path, USE.NAMES = FALSE)]
missing_inputs <- manifest[!file.exists(resolved_path)]
if (nrow(missing_inputs) > 0) {
  stop(sprintf("Missing required inputs: %s", paste(missing_inputs$relative_path, collapse = ", ")), call. = FALSE)
}
input_sizes <- file.info(manifest$resolved_path)$size
size_mismatches <- manifest[input_sizes != as.numeric(bytes)]
if (nrow(size_mismatches) > 0) {
  stop(sprintf("Input size mismatch: %s", paste(size_mismatches$relative_path, collapse = ", ")), call. = FALSE)
}
input_hashes <- vapply(manifest$resolved_path, digest::digest, character(1), file = TRUE, algo = "sha256", serialize = FALSE)
hash_mismatches <- manifest[tolower(input_hashes) != tolower(sha256)]
if (nrow(hash_mismatches) > 0) {
  stop(sprintf("Input SHA-256 mismatch: %s", paste(hash_mismatches$relative_path, collapse = ", ")), call. = FALSE)
}

if (length(read_gene_set("external_six_gene_module")) != 6 || length(read_gene_set("visium_frozen_24_gene_module")) != 24) {
  stop("Frozen candidate gene-set definitions are incomplete.", call. = FALSE)
}
if (!all(c("myCAF", "CAF_stroma", "endothelial") %in% names(read_phenotype_gene_sets()))) {
  stop("Frozen primary phenotype definitions are incomplete.", call. = FALSE)
}

Sys.setenv(
  PDAC_PROJECT_ROOT = ROOT,
  PDAC_DATABASE_ROOT = DATABASE_ROOT,
  PDAC_USE_CACHE = "0"
)

python_executable <- Sys.getenv("PDAC_PYTHON", unset = "python")
python_status <- system2(python_executable, file.path(script_dir, "00_validate_python_environment.py"))
if (!identical(python_status, 0L)) {
  stop("Python environment validation failed.", call. = FALSE)
}

analysis_steps <- c(
  "01_prepare_matrisome_mapping.R",
  "02_build_protein_residual_scores.R",
  "03_build_transcriptomic_phenotypes.R",
  "04_run_primary_robustness_and_survival.R",
  "05_run_tmt_plex_sensitivity.R",
  "06_run_tcga_bulk_context.R",
  "07_run_scrna_cell_context.R",
  "08_run_gse202051_patient_context.py",
  "09_run_geomx_patient_context.R",
  "10_run_candidate_integration.R",
  "11_run_visium_spatial_context.R",
  "12_run_enzyme_and_spatial_neighborhood.R",
  "13_run_score_contribution_decomposition.R",
  "14_run_candidate_ecm_receptor_context.R"
)

for (step in analysis_steps) {
  if (grepl("\\.R$", step)) {
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  } else {
    executable <- python_executable
  }
  status <- system2(executable, file.path(script_dir, step))
  if (!identical(status, 0L)) {
    stop(sprintf("Analysis step failed: %s", step), call. = FALSE)
  }
}

required_outputs <- c(
  file.path(RESULTS, "analysis_results", "nglyco_peptide_analysis_layer_state.rds"),
  file.path(RESULTS, "analysis_results", "nglyco_site_analysis_layer_state.rds"),
  file.path(RESULTS, "analysis_results", "primary_pdac_and_140_sensitivity_models.tsv"),
  file.path(RESULTS, "analysis_results", "detection_model_summary.tsv"),
  file.path(RESULTS, "analysis_results", "detection_ipw_parameter_sensitivity.tsv"),
  file.path(RESULTS, "analysis_results", "postn_feature_removed_predictor_models.tsv"),
  file.path(RESULTS, "analysis_results", "cptac_pdac_survival_models.tsv"),
  file.path(RESULTS, "analysis_results", "cptac_pdac_km_summary.tsv"),
  file.path(RESULTS, "analysis_results", "gse202051_signature_correlation_summary.tsv"),
  file.path(RESULTS, "analysis_results", "geomx_patient_level_paired_effects.tsv"),
  file.path(RESULTS, "analysis_results", "glycan_class_enzyme_expression_models.tsv"),
  file.path(RESULTS, "analysis_results", "gse272362_neighborhood_summary.tsv"),
  file.path(RESULTS, "analysis_results", "residual_score_module_contributions.tsv"),
  file.path(RESULTS, "analysis_results", "receptor_ligand", "gse154778_candidate_ecm_interaction_summary.tsv"),
  file.path(RESULTS, "analysis_results", "receptor_ligand", "gse154778_candidate_ecm_enrichment_by_patient.tsv"),
  file.path(RESULTS, "analysis_results", "receptor_ligand", "gse154778_postn_interaction_summary.tsv"),
  file.path(RESULTS, "cptac_batch_plex_sensitivity_models.tsv"),
  file.path(RESULTS, "tcga_paad_bulk_validation_models.tsv"),
  file.path(RESULTS, "gse154778_candidate_module_by_inferred_cell_type.tsv"),
  file.path(RESULTS, "c4_candidate_core_set.tsv"),
  file.path(RESULTS, "gse272362_spatial_module_summary.tsv")
)
missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs) > 0) {
  stop(sprintf("Required analysis outputs were not generated: %s", paste(missing_outputs, collapse = ", ")), call. = FALSE)
}
empty_outputs <- required_outputs[file.info(required_outputs)$size <= 0]
if (length(empty_outputs) > 0) {
  stop(sprintf("Required analysis outputs are empty: %s", paste(empty_outputs, collapse = ", ")), call. = FALSE)
}
