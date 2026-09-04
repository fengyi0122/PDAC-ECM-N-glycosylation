script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

RAW <- existing_dataset_dir("cptac", "CPTAC-PDAC")
matrisome <- unique(data.table::fread(MATRISOME_FILE)[, .(gene, division, category)])

specifications <- list(
  list(layer = "nglyco_peptide", filename = "N-glycoproteomics_peptide_level_ratio_tumor.cct", feature_col = "Sequence"),
  list(layer = "nglyco_site", filename = "N-glycoproteomics_Site_level_ratio_tumor.cct", feature_col = "Modifications")
)

records <- lapply(specifications, function(spec) {
  features <- data.table::fread(file.path(RAW, spec$filename), select = c(spec$feature_col, "Gene"))
  data.table::setnames(features, c(spec$feature_col, "Gene"), c("feature_id", "gene"))
  features <- unique(features[!is.na(feature_id) & !is.na(gene), .(feature_id, gene)])
  mapped <- merge(features, matrisome, by = "gene", all = FALSE)
  mapped[, layer := spec$layer]
  data.table::setcolorder(mapped, c("layer", "feature_id", "gene", "division", "category"))
  mapped
})

mapping <- data.table::rbindlist(records, fill = TRUE)
data.table::setorder(mapping, layer, division, category, gene, feature_id)
data.table::fwrite(mapping, file.path(PROCESSED, "matrisome_glyco_map.tsv"), sep = "\t")
cat("Wrote", file.path(PROCESSED, "matrisome_glyco_map.tsv"), "\n")
