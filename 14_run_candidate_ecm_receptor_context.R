script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table", "Matrix", "CellChat"))
ensure_project_dirs()

set.seed(20260721)
RAW <- existing_dataset_dir("geo_pdac", "GSE154778")
OUT <- file.path(RESULTS, "analysis_results", "receptor_ligand")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
use_cache <- identical(Sys.getenv("PDAC_USE_CACHE", unset = "0"), "1")

sample_map <- data.table::data.table(patient = sprintf("P%02d", 1:10))
labels <- data.table::fread(file.path(PROCESSED, "gse154778_marker_inferred_cell_scores.tsv"))
labels <- labels[patient %in% sample_map$patient]

frozen_24 <- read_gene_set("visium_frozen_24_gene_module")
frozen_6 <- read_gene_set("external_six_gene_module")

suppressPackageStartupMessages(library(CellChat))
data(CellChatDB.human)
db_use <- CellChatDB.human
db_interactions <- data.table::as.data.table(db_use$interaction)
db_interactions <- db_interactions[annotation == "ECM-Receptor" | pathway_name == "PERIOSTIN"]
db_interactions <- unique(db_interactions, by = "interaction_name")
db_interactions[, interaction_name := as.character(interaction_name)]
db_use$interaction <- as.data.frame(db_interactions)
rownames(db_use$interaction) <- db_use$interaction$interaction_name
candidate_db <- db_interactions[ligand %in% frozen_24]
candidate_db[, module_class := data.table::fcase(
  ligand == "POSTN", "POSTN branch",
  ligand %in% frozen_6, "Six-gene module",
  default = "Distributed 24-gene program"
)]

required_genes <- unique(CellChat:::extractGeneSubsetFromPair(db_interactions, complex_input = db_use$complex, geneInfo = db_use$geneInfo))
gene_manifest <- file.path(OUT, "cellchat_required_genes.tsv")
data.table::fwrite(data.table::data.table(gene = required_genes), gene_manifest, sep = "\t")
selected_matrix_file <- file.path(PROCESSED, "gse154778_cellchat_selected_counts.tsv.gz")
library_file <- file.path(PROCESSED, "gse154778_cell_library_sizes.tsv")
extraction_audit <- file.path(OUT, "gse154778_cellchat_extraction_audit.tsv")
python <- Sys.getenv("PDAC_PYTHON", unset = "python")
extractor <- file.path(script_dir, "14a_extract_gse154778_cellchat_matrix.py")
if (!use_cache || !file.exists(selected_matrix_file) || !file.exists(library_file) || !file.exists(extraction_audit)) {
  status_extract <- system2(python, c(
    extractor,
    "--input", file.path(RAW, "GSE154778_dgeMtx.csv.gz"),
    "--genes", gene_manifest,
    "--output-matrix", selected_matrix_file,
    "--output-library", library_file,
    "--output-audit", extraction_audit
  ))
  if (!identical(status_extract, 0L)) stop("GSE154778 CellChat matrix extraction failed", call. = FALSE)
}
selected_counts <- data.table::fread(selected_matrix_file)
library_sizes <- data.table::fread(library_file)
selected_genes <- selected_counts$gene

read_sample <- function(patient_id) {
  cell_ids <- grep(sprintf("^%s:", patient_id), names(selected_counts), value = TRUE)
  patient_labels <- labels[patient == patient_id][match(cell_ids, cell_id)]
  if (length(cell_ids) != nrow(patient_labels) || anyNA(patient_labels$cell_id)) stop(sprintf("Cell metadata mismatch for %s", patient_id), call. = FALSE)
  counts <- as.matrix(selected_counts[, ..cell_ids])
  rownames(counts) <- selected_genes
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  lib <- library_sizes[match(cell_ids, cell_id), library_size]
  if (any(!is.finite(lib) | lib <= 0)) stop(sprintf("Invalid library sizes for %s", patient_id), call. = FALSE)
  list(counts = counts, library_size = lib, meta = patient_labels[, .(cell_id, cell_type = inferred_cell_type)])
}

communication_records <- list()
cell_count_records <- list()
patient_status <- list()
patient_cache <- file.path(OUT, "patient_cellchat")
dir.create(patient_cache, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(sample_map))) {
  patient_id <- sample_map$patient[[i]]
  sample_data <- read_sample(patient_id)
  counts_by_type <- sample_data$meta[, .(n_cells = .N), by = cell_type]
  counts_by_type[, patient := patient_id]
  cell_count_records[[length(cell_count_records) + 1]] <- counts_by_type
  caf_n <- counts_by_type[cell_type == "CAF", n_cells]
  if (!length(caf_n)) caf_n <- 0
  receivers <- counts_by_type[cell_type != "CAF" & n_cells >= 10, cell_type]
  if (caf_n < 20 || !length(receivers)) {
    patient_status[[length(patient_status) + 1]] <- data.table::data.table(patient = patient_id, caf_cells = caf_n, eligible_receivers = length(receivers), status = "excluded")
    next
  }
  keep_types <- c("CAF", receivers)
  keep_cells <- sample_data$meta[cell_type %in% keep_types, cell_id]
  cache_file <- file.path(patient_cache, sprintf("%s_communications.tsv", patient_id))
  if (use_cache && file.exists(cache_file)) {
    comm <- data.table::fread(cache_file)
    if (nrow(comm)) communication_records[[length(communication_records) + 1]] <- comm
    patient_status[[length(patient_status) + 1]] <- data.table::data.table(patient = patient_id, caf_cells = caf_n, eligible_receivers = length(receivers), status = "analyzed")
    next
  }
  counts <- sample_data$counts[, keep_cells, drop = FALSE]
  lib <- sample_data$library_size[match(keep_cells, sample_data$meta$cell_id)]
  norm <- Matrix::t(Matrix::t(counts) / lib * 10000)
  norm@x <- log1p(norm@x)
  meta <- as.data.frame(sample_data$meta[match(keep_cells, cell_id), .(cell_type)])
  meta$samples <- factor(patient_id)
  rownames(meta) <- keep_cells
  cellchat <- createCellChat(object = norm, meta = meta, group.by = "cell_type")
  cellchat@DB <- db_use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0, thresh.fc = 0, thresh.p = 1)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = FALSE, distance.use = FALSE, nboot = 100, seed.use = 20260721 + i)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  comm <- data.table::as.data.table(subsetCommunication(cellchat, sources.use = "CAF", targets.use = receivers, thresh = 1))
  if (nrow(comm)) {
    comm[, interaction_name := as.character(interaction_name)]
    comm[, patient := patient_id]
    communication_records[[length(communication_records) + 1]] <- comm
  }
  data.table::fwrite(comm, cache_file, sep = "\t")
  patient_status[[length(patient_status) + 1]] <- data.table::data.table(patient = patient_id, caf_cells = caf_n, eligible_receivers = length(receivers), status = "analyzed")
  rm(cellchat, norm, counts, sample_data)
  gc()
}

cell_counts <- data.table::rbindlist(cell_count_records, fill = TRUE)
status <- data.table::rbindlist(patient_status, fill = TRUE)
communication <- data.table::rbindlist(communication_records, fill = TRUE)
communication[, interaction_name := as.character(interaction_name)]
numeric_interaction <- grepl("^[0-9]+$", communication$interaction_name)
if (any(numeric_interaction)) {
  interaction_index <- as.integer(communication$interaction_name[numeric_interaction])
  if (any(interaction_index < 1 | interaction_index > nrow(db_interactions))) stop("Invalid CellChat interaction index", call. = FALSE)
  communication[numeric_interaction, interaction_name := db_interactions$interaction_name[interaction_index]]
}
eligible_contexts <- cell_counts[n_cells >= 10 & cell_type != "CAF", .(target = cell_type), by = patient]
eligible_contexts <- merge(eligible_contexts, status[status == "analyzed", .(patient)], by = "patient")

eligible_contexts[, join_key := 1L]
db_grid <- db_interactions[, .(interaction_name, ligand, receptor, pathway_name, annotation)]
db_grid[, join_key := 1L]
background_grid <- merge(eligible_contexts, db_grid, by = "join_key", allow.cartesian = TRUE)
background_grid[, join_key := NULL]
observed <- communication[, .(patient, target, interaction_name, probability = prob, permutation_p = pval)]
background_grid <- merge(background_grid, observed, by = c("patient", "target", "interaction_name"), all.x = TRUE)
background_grid[!is.finite(probability), probability := 0]
background_grid[!is.finite(permutation_p), permutation_p := 1]
background_grid[, within_context_percentile := data.table::fifelse(probability > 0, data.table::frank(probability, ties.method = "average") / .N, 0), by = .(patient, target)]
background_grid[, significant := permutation_p < 0.05]
background_grid[, candidate_ligand := ligand %in% frozen_24]
background_grid[, six_gene_ligand := ligand %in% frozen_6]
background_grid[, postn_branch := ligand == "POSTN"]

candidate_long <- background_grid[candidate_ligand == TRUE]
candidate_long <- merge(candidate_long, candidate_db[, .(interaction_name, module_class)], by = "interaction_name", all.x = TRUE)
candidate_summary <- candidate_long[, .(
  eligible_patients = .N,
  positive_patients = sum(probability > 0),
  significant_patients = sum(significant),
  conservation_fraction = mean(significant),
  median_probability = stats::median(probability),
  median_within_context_percentile = stats::median(within_context_percentile),
  top_quartile_fraction = mean(within_context_percentile >= 0.75)
), by = .(interaction_name, ligand, receptor, pathway_name, module_class, target)]
candidate_summary[, recurrent_top_ranked := eligible_patients >= 3 & conservation_fraction >= 0.5 & median_within_context_percentile >= 0.75]
data.table::setorder(candidate_summary, -recurrent_top_ranked, -conservation_fraction, -median_within_context_percentile)

postn_long <- candidate_long[postn_branch == TRUE]
postn_summary <- candidate_summary[ligand == "POSTN"]

enrichment_by_context <- background_grid[, {
  candidate_values <- within_context_percentile[candidate_ligand]
  background_values <- within_context_percentile[!candidate_ligand]
  test <- tryCatch(stats::wilcox.test(candidate_values, background_values, alternative = "greater", exact = FALSE), error = function(e) NULL)
  list(
    candidate_median_percentile = stats::median(candidate_values),
    background_median_percentile = stats::median(background_values),
    percentile_difference = stats::median(candidate_values) - stats::median(background_values),
    enrichment_p = if (is.null(test)) NA_real_ else test$p.value
  )
}, by = .(patient, target)]
enrichment_by_context[, enrichment_fdr := stats::p.adjust(enrichment_p, method = "BH")]

enrichment_summary <- enrichment_by_context[, .(
  eligible_patients = .N,
  positive_patients = sum(percentile_difference > 0),
  median_percentile_difference = stats::median(percentile_difference),
  fdr_significant_patients = sum(enrichment_fdr < 0.05)
), by = target]

recurrent <- candidate_summary[recurrent_top_ranked == TRUE]
gate_pass <- status[status == "analyzed", .N] >= 8 &&
  data.table::uniqueN(candidate_db$ligand) >= 5 &&
  nrow(recurrent) >= 3 &&
  nrow(postn_summary) >= 2

data.table::fwrite(status, file.path(OUT, "gse154778_cellchat_patient_status.tsv"), sep = "\t")
data.table::fwrite(cell_counts, file.path(OUT, "gse154778_cellchat_cell_counts.tsv"), sep = "\t")
data.table::fwrite(db_interactions, file.path(OUT, "cellchat_ecm_receptor_database.tsv"), sep = "\t")
data.table::fwrite(candidate_db, file.path(OUT, "cptac_candidate_ligand_receptor_database.tsv"), sep = "\t")
data.table::fwrite(communication, file.path(OUT, "gse154778_cellchat_observed_communications.tsv"), sep = "\t")
data.table::fwrite(candidate_long, file.path(OUT, "gse154778_candidate_ecm_interactions_by_patient.tsv"), sep = "\t")
data.table::fwrite(candidate_summary, file.path(OUT, "gse154778_candidate_ecm_interaction_summary.tsv"), sep = "\t")
data.table::fwrite(postn_long, file.path(OUT, "gse154778_postn_interactions_by_patient.tsv"), sep = "\t")
data.table::fwrite(postn_summary, file.path(OUT, "gse154778_postn_interaction_summary.tsv"), sep = "\t")
data.table::fwrite(enrichment_by_context, file.path(OUT, "gse154778_candidate_ecm_enrichment_by_patient.tsv"), sep = "\t")
data.table::fwrite(enrichment_summary, file.path(OUT, "gse154778_candidate_ecm_enrichment_summary.tsv"), sep = "\t")

writeLines(c(
  "CANDIDATE ECM LIGAND-RECEPTOR ANALYSIS GATE",
  sprintf("Primary tumors analyzed: %d", status[status == "analyzed", .N]),
  sprintf("CellChat ECM-receptor background interactions: %d", nrow(db_interactions)),
  sprintf("Frozen 24-gene candidate ligands represented: %d", data.table::uniqueN(candidate_db$ligand)),
  sprintf("Candidate ligand-receptor interactions represented: %d", nrow(candidate_db)),
  sprintf("Recurrent top-ranked interaction-recipient contexts: %d", nrow(recurrent)),
  sprintf("POSTN interaction-recipient contexts evaluated: %d", nrow(postn_summary)),
  sprintf("Decision: %s", if (gate_pass) "PASS" else "STOP"),
  "Interpretation ceiling: RNA-inferred receptor context for CPTAC-supported candidate ECM ligands, not glycoform-specific binding, protein secretion, or functional communication"
), file.path(OUT, "CANDIDATE_ECM_CELLCHAT_GATE.txt"))
cat(readLines(file.path(OUT, "CANDIDATE_ECM_CELLCHAT_GATE.txt")), sep = "\n")
