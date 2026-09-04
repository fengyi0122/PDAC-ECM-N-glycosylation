script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/")) else file.path(getwd(), "scripts")
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table", "FNN"))
suppressPackageStartupMessages(library(data.table))

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

safe_median <- function(x) {
  if (all(!is.finite(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

fit_adjusted <- function(data, outcome, predictor) {
  frame <- data[, .(
    y = zscore(get(outcome)),
    x = zscore(get(predictor)),
    proteome = zscore(proteome_matrisome_score),
    stroma = zscore(stromal_fraction_z),
    purity = zscore(neoplastic_cellularity_z)
  )]
  frame <- frame[stats::complete.cases(frame)]
  if (nrow(frame) < 40) return(NULL)
  fit <- stats::lm(y ~ x + proteome + stroma + purity, data = frame)
  sm <- summary(fit)$coefficients
  ci <- stats::confint(fit, "x")
  data.table::data.table(
    n = nrow(frame), beta = sm["x", "Estimate"], se = sm["x", "Std. Error"],
    ci_low = ci[[1]], ci_high = ci[[2]], p_value = sm["x", "Pr(>|t|)"]
  )
}

clinical <- data.table::fread(file.path(RAW, "clinical_table_140.tsv"))
clinical <- clinical[histology_diagnosis == "PDAC", .(sample_id = case_id)]
phen <- data.table::fread(file.path(PROCESSED, "immune_phenotypes.tsv"))
glyco <- data.table::fread(file.path(PROCESSED, "glyco_matrisome_features.tsv"))
model_data <- merge(merge(phen, glyco, by = "sample_id"), clinical, by = "sample_id")

expr <- data.table::fread(file.path(RAW, "mRNA_RSEM_UQ_log2_Tumor.cct"))
data.table::setnames(expr, 1, "gene")
sample_cols <- intersect(model_data$sample_id, names(expr))
expr[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]
expr <- expr[, lapply(.SD, safe_median), by = gene, .SDcols = sample_cols]
expr_mat <- as.matrix(expr[, ..sample_cols])
rownames(expr_mat) <- expr$gene
storage.mode(expr_mat) <- "numeric"

enzyme_groups <- list(
  fucosylation = c("FUT3", "FUT8", "FUT11"),
  sialylation = c("ST3GAL1", "ST3GAL3", "ST3GAL4", "ST6GAL1"),
  mannose_processing = c("MAN1A1", "MAN1A2", "MAN2A1", "MAN2A2", "MGAT1", "MGAT2"),
  branching = c("MGAT3", "MGAT5", "B4GALT1", "B4GALT4")
)
enzyme_manifest <- data.table::rbindlist(lapply(names(enzyme_groups), function(group) {
  data.table::data.table(enzyme_group = group, gene = enzyme_groups[[group]], present_in_cptac_rna = enzyme_groups[[group]] %in% rownames(expr_mat))
}))
data.table::fwrite(enzyme_manifest, file.path(OUT, "glycosylation_enzyme_manifest.tsv"), sep = "\t")

for (gene in enzyme_manifest[present_in_cptac_rna == TRUE, gene]) {
  values <- as.numeric(expr_mat[gene, model_data$sample_id])
  model_data[[paste0("enzyme_", gene)]] <- values
}

class_records <- list()
enzyme_records <- list()
for (layer_name in c("nglyco_peptide", "nglyco_site")) {
  cache_file <- file.path(OUT, sprintf("%s_analysis_layer_state.rds", layer_name))
  if (!file.exists(cache_file)) stop("Missing analysis layer state: ", cache_file)
  layer <- readRDS(cache_file)
  class_scores <- copy(layer$class_scores)
  class_scores[, layer := layer_name]
  class_records[[length(class_records) + 1]] <- class_scores
  layer_data <- merge(model_data, class_scores, by = "sample_id", all.x = TRUE)

  class_to_group <- list(
    fucosylated = "fucosylation",
    sialylated = "sialylation",
    high_mannose = "mannose_processing"
  )
  for (glycan_class in names(class_to_group)) {
    group <- class_to_group[[glycan_class]]
    genes <- enzyme_manifest[enzyme_group == group & present_in_cptac_rna == TRUE, gene]
    for (gene in genes) {
      outcome <- paste0("enzyme_", gene)
      effect <- fit_adjusted(layer_data, outcome, glycan_class)
      if (is.null(effect)) next
      effect[, `:=`(layer = layer_name, glycan_class = glycan_class, enzyme_group = group, enzyme = gene)]
      enzyme_records[[length(enzyme_records) + 1]] <- effect
    }
  }
}
class_score_table <- data.table::rbindlist(class_records, fill = TRUE)
data.table::fwrite(class_score_table, file.path(OUT, "glycan_composition_sample_scores.tsv"), sep = "\t")
enzyme_models <- data.table::rbindlist(enzyme_records, fill = TRUE)
enzyme_models[, fdr := stats::p.adjust(p_value, method = "BH"), by = .(layer, glycan_class)]
data.table::fwrite(enzyme_models, file.path(OUT, "glycan_class_enzyme_expression_models.tsv"), sep = "\t")

spatial_scores <- data.table::fread(file.path(PROCESSED, "gse272362_spatial_module_scores.tsv"))
spatial_rds <- file.path(existing_dataset_dir("geo_pdac", "GSE272362"), "Zenodo", "PDAC_Updated.rds")
if (!file.exists(spatial_rds)) stop("Missing GSE272362 spatial RDS: ", spatial_rds)
message("Reading spatial object for CD8-NK neighborhood analysis")
obj <- readRDS(spatial_rds)
rctd <- as.matrix(methods::slot(obj@assays[["rctd_fullfinal"]], "data"))
if (!"CD8-NK cells" %in% rownames(rctd)) stop("CD8-NK cells are unavailable in rctd_fullfinal")
cd8_nk <- data.table::data.table(spot = colnames(rctd), rctd_cd8_nk = as.numeric(rctd["CD8-NK cells", ]))
spatial_scores <- merge(spatial_scores, cd8_nk, by = "spot", all.x = TRUE)
data.table::fwrite(spatial_scores[, .(spot, image, Origin, ecm_glyco_module, caf_module, rctd_cd8_nk, ecm_caf_colocalized, imagecol, imagerow)],
                   file.path(OUT, "gse272362_spatial_neighborhood_scores.tsv"), sep = "\t")

section_records <- list()
primary <- spatial_scores[Origin == "Pancreas" & is.finite(imagecol) & is.finite(imagerow) & is.finite(rctd_cd8_nk)]
for (section_name in unique(primary$image)) {
  sec <- primary[image == section_name]
  if (nrow(sec) < 100 || sum(sec$ecm_caf_colocalized, na.rm = TRUE) < 10) next
  coords <- as.matrix(sec[, .(imagecol, imagerow)])
  knn <- FNN::get.knnx(coords, coords, k = 7)
  neighbor_index <- knn$nn.index[, -1, drop = FALSE]
  sec[, neighbor_cd8_nk := apply(neighbor_index, 1, function(idx) mean(rctd_cd8_nk[idx], na.rm = TRUE))]
  cd8_threshold <- stats::quantile(sec$rctd_cd8_nk, 0.75, na.rm = TRUE)
  cd8_high <- which(sec$rctd_cd8_nk >= cd8_threshold)
  if (length(cd8_high) < 5) next
  nearest_cd8 <- FNN::get.knnx(coords[cd8_high, , drop = FALSE], coords, k = 1)$nn.dist[, 1]
  sec[, nearest_cd8_nk_distance := nearest_cd8]
  coloc <- sec$ecm_caf_colocalized
  section_records[[length(section_records) + 1]] <- data.table::data.table(
    image = section_name,
    n_spots = nrow(sec),
    n_colocalized = sum(coloc, na.rm = TRUE),
    neighbor_cd8_delta = mean(sec$neighbor_cd8_nk[coloc], na.rm = TRUE) - mean(sec$neighbor_cd8_nk[!coloc], na.rm = TRUE),
    nearest_cd8_distance_delta = stats::median(sec$nearest_cd8_nk_distance[coloc], na.rm = TRUE) - stats::median(sec$nearest_cd8_nk_distance[!coloc], na.rm = TRUE)
  )
}
section_neighborhood <- data.table::rbindlist(section_records, fill = TRUE)
data.table::fwrite(section_neighborhood, file.path(OUT, "gse272362_section_neighborhood_effects.tsv"), sep = "\t")

neighborhood_summary <- data.table::data.table(
  metric = c("neighbor_cd8_delta", "nearest_cd8_distance_delta"),
  n_sections = nrow(section_neighborhood),
  median_effect = c(stats::median(section_neighborhood$neighbor_cd8_delta, na.rm = TRUE), stats::median(section_neighborhood$nearest_cd8_distance_delta, na.rm = TRUE)),
  positive_sections = c(sum(section_neighborhood$neighbor_cd8_delta > 0, na.rm = TRUE), sum(section_neighborhood$nearest_cd8_distance_delta > 0, na.rm = TRUE)),
  wilcoxon_p = c(
    stats::wilcox.test(section_neighborhood$neighbor_cd8_delta, mu = 0, exact = FALSE)$p.value,
    stats::wilcox.test(section_neighborhood$nearest_cd8_distance_delta, mu = 0, exact = FALSE)$p.value
  )
)
neighborhood_summary[, fdr := stats::p.adjust(wilcoxon_p, method = "BH")]
data.table::fwrite(neighborhood_summary, file.path(OUT, "gse272362_neighborhood_summary.tsv"), sep = "\t")

cat("GLYCAN ENZYME MODELS\n")
print(enzyme_models[order(fdr)][1:min(.N, 10)])
cat("\nSPATIAL NEIGHBORHOOD SUMMARY\n")
print(neighborhood_summary)
