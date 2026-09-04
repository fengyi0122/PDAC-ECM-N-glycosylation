script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/")) else file.path(getwd(), "scripts")
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table", "Matrix", "SeuratObject"))
ensure_project_dirs()

GSE272362 <- existing_dataset_dir("geo_pdac", "GSE272362")
RDS <- file.path(GSE272362, "Zenodo", "PDAC_Updated.rds")
if (!file.exists(RDS)) {
  stop(sprintf("Missing GSE272362 RDS: %s", RDS), call. = FALSE)
}

out_scores <- file.path(PROCESSED, "gse272362_spatial_module_scores.tsv")
out_summary <- file.path(RESULTS, "gse272362_spatial_module_summary.tsv")
out_gene_coverage <- file.path(RESULTS, "gse272362_ecm_glyco_module_gene_coverage.tsv")
out_correlations <- file.path(RESULTS, "figure4_v3_2_gse272362_module_correlations.tsv")

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  t(apply(mat, 1, zscore))
}

score_rows <- function(mat, features, fallback = NULL) {
  found <- intersect(features, rownames(mat))
  if (length(found) >= 2) {
    return(list(score = Matrix::colMeans(zscore_rows(mat[found, , drop = FALSE]), na.rm = TRUE), found = found, source = "gene_set"))
  }
  if (!is.null(fallback) && fallback %in% rownames(mat)) {
    return(list(score = zscore(as.numeric(mat[fallback, ])), found = fallback, source = "fallback"))
  }
  stop(sprintf("Could not score features. Found only %d requested genes and fallback `%s` is unavailable.", length(found), fallback), call. = FALSE)
}

message("Reading ", RDS)
obj <- readRDS(RDS)

core_genes <- read_gene_set("visium_frozen_24_gene_module")

spatial <- methods::slot(obj@assays[["Spatial"]], "data")
ecm_scored <- score_rows(spatial, core_genes)
gene_coverage <- data.table::data.table(
  module = "CPTAC_CORE_ECM_glyco_module",
  requested_genes = length(core_genes),
  found_genes = length(ecm_scored$found),
  score_source = ecm_scored$source,
  genes_found = paste(ecm_scored$found, collapse = ",")
)
data.table::fwrite(gene_coverage, out_gene_coverage, sep = "\t")

fges <- methods::slot(obj@assays[["fges"]], "data")
if (is.null(dim(fges))) {
  fges <- obj@assays[["fges"]]@data
}
fges <- as.matrix(fges)

get_fges <- function(feature) {
  if (!feature %in% rownames(fges)) {
    stop(sprintf("Missing fges feature `%s`.", feature), call. = FALSE)
  }
  zscore(as.numeric(fges[feature, ]))
}

immune_features <- intersect(c("Effector-cells", "Tcells", "NK-cells", "Bcells", "Myeloid-cells-traffic"), rownames(fges))
if (length(immune_features) == 0) {
  immune_features <- intersect(c("Effector-cells", "Tcells"), rownames(fges))
}
immune_mat <- fges[immune_features, , drop = FALSE]
immune_score <- Matrix::colMeans(zscore_rows(immune_mat), na.rm = TRUE)

rctd <- as.matrix(methods::slot(obj@assays[["rctd_fullfinal"]], "data"))
rctd_rows <- rownames(rctd)
get_rctd <- function(feature) {
  if (!feature %in% rctd_rows) return(rep(NA_real_, ncol(rctd)))
  as.numeric(rctd[feature, ])
}

scores <- data.table::data.table(
  spot = colnames(spatial),
  ecm_glyco_module = zscore(ecm_scored$score),
  caf_module = get_fges("Cancer-associated-fibroblasts"),
  immune_module = zscore(immune_score),
  tumor_epithelial_module = get_fges("Tumor-proliferation-rate"),
  fges_ecm = get_fges("ECM"),
  rctd_mycaf = get_rctd("myCAF"),
  rctd_icaf = get_rctd("iCAF"),
  rctd_tumor_epithelial = get_rctd("Tumor Epithelial cells"),
  rctd_immune = rowSums(rctd[intersect(c("B cells", "CD4+ cells", "CD8-NK cells", "DCs", "C1Q-TAM", "FCN1-TAM", "SPP1-TAM", "Proliferative T cells"), rctd_rows), , drop = FALSE])
)
scores[, rctd_caf := rowSums(cbind(rctd_mycaf, rctd_icaf), na.rm = TRUE)]

meta <- data.table::as.data.table(obj@meta.data, keep.rownames = "spot")
keep_meta <- intersect(
  c("spot", "Sample_ID2", "SlideName", "Origin", "Treatment", "spot_class", "first_type", "second_type", "first_class", "second_class", "patient", "orig.ident"),
  names(meta)
)
meta <- meta[, ..keep_meta]

coords <- data.table::rbindlist(lapply(names(obj@images), function(image_name) {
  co <- data.table::as.data.table(methods::slot(obj@images[[image_name]], "coordinates"), keep.rownames = "spot")
  co[, image := image_name]
  co
}), fill = TRUE)

dt <- merge(scores, meta, by = "spot", all.x = TRUE)
dt <- merge(dt, coords, by = "spot", all.x = TRUE)
dt <- dt[!is.na(imagecol) & !is.na(imagerow)]
dt[, Sample_ID2 := trimws(Sample_ID2)]
dt[, `:=`(
  ecm_high = FALSE,
  caf_high = FALSE,
  immune_high = FALSE
)]
dt[, ecm_high := ecm_glyco_module >= stats::quantile(ecm_glyco_module, 0.75, na.rm = TRUE), by = image]
dt[, caf_high := caf_module >= stats::quantile(caf_module, 0.75, na.rm = TRUE), by = image]
dt[, immune_high := immune_module >= stats::quantile(immune_module, 0.75, na.rm = TRUE), by = image]
dt[, ecm_caf_colocalized := ecm_high & caf_high]

summary <- dt[, .(
  n_spots = .N,
  origin = paste(sort(unique(stats::na.omit(Origin))), collapse = ";"),
  median_ecm = stats::median(ecm_glyco_module, na.rm = TRUE),
  median_caf = stats::median(caf_module, na.rm = TRUE),
  median_immune = stats::median(immune_module, na.rm = TRUE),
  coloc_spots = sum(ecm_caf_colocalized, na.rm = TRUE),
  coloc_fraction = mean(ecm_caf_colocalized, na.rm = TRUE),
  first_type_mode = names(sort(table(first_type), decreasing = TRUE))[1]
), by = .(image, Sample_ID2)]
data.table::setorder(summary, -coloc_spots, -coloc_fraction)

primary <- dt[Origin == "Pancreas"]
correlations <- data.table::rbindlist(lapply(split(primary, primary$image), function(section) {
  data.table::data.table(
    image = section$image[[1]],
    Sample_ID2 = section$Sample_ID2[[1]],
    module = c("CAF module", "Immune module", "Tumor epithelial"),
    spearman_rho = c(
      suppressWarnings(stats::cor(section$ecm_glyco_module, section$caf_module, method = "spearman", use = "complete.obs")),
      suppressWarnings(stats::cor(section$ecm_glyco_module, section$immune_module, method = "spearman", use = "complete.obs")),
      suppressWarnings(stats::cor(section$ecm_glyco_module, section$tumor_epithelial_module, method = "spearman", use = "complete.obs"))
    )
  )
}), fill = TRUE)

data.table::fwrite(dt, out_scores, sep = "\t")
data.table::fwrite(summary, out_summary, sep = "\t")
data.table::fwrite(correlations, out_correlations, sep = "\t")

cat("Wrote", out_scores, "\n")
cat("Wrote", out_summary, "\n")
cat("Wrote", out_gene_coverage, "\n")
cat("Wrote", out_correlations, "\n")
