script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_file, winslash = "/"))
source(file.path(script_dir, "00_paths_and_functions.R"))
ensure_packages(c("data.table"))
ensure_project_dirs()

set.seed(20260717)
out_dir <- file.path(RESULTS, "analysis_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_programs <- c(
  "CAF/matricellular" = "score_decomposition_caf_matricellular",
  "Basement membrane/laminin" = "score_decomposition_basement_laminin",
  "Elastic/fibulin" = "score_decomposition_elastic_fibulin",
  "Vascular ECM" = "score_decomposition_vascular_ecm"
)
module_map <- unlist(lapply(names(module_programs), function(label) {
  stats::setNames(rep(label, length(read_gene_set(module_programs[[label]]))), read_gene_set(module_programs[[label]]))
}))

covariance_share <- function(component, total) {
  mask <- is.finite(component) & is.finite(total)
  if (sum(mask) < 3 || stats::var(total[mask]) == 0) return(NA_real_)
  stats::cov(component[mask], total[mask]) / stats::var(total[mask])
}

bootstrap_share <- function(component, total, replicates = 2000L) {
  n <- length(total)
  values <- replicate(replicates, {
    idx <- sample.int(n, n, replace = TRUE)
    covariance_share(component[idx], total[idx])
  })
  stats::quantile(values, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

decompose_layer <- function(layer) {
  cache <- readRDS(file.path(out_dir, sprintf("%s_analysis_layer_state.rds", layer)))
  annotation <- data.table::as.data.table(cache$annotation)
  keep <- which(annotation$observed_n >= 50)
  annotation <- annotation[keep]
  residual <- cache$residual_linear_z[keep, , drop = FALSE]
  denominator <- colSums(is.finite(residual))
  total <- colSums(residual, na.rm = TRUE) / denominator
  annotation[, module := unname(module_map[gene])]
  annotation[is.na(module), module := "Other ECM glycoproteins"]

  module_rows <- lapply(unique(annotation$module), function(label) {
    idx <- which(annotation$module == label)
    component <- colSums(residual[idx, , drop = FALSE], na.rm = TRUE) / denominator
    interval <- bootstrap_share(component, total)
    data.table::data.table(
      layer = layer,
      component = label,
      contribution_share = covariance_share(component, total),
      ci_low = interval[[1]],
      ci_high = interval[[2]],
      retained_features = length(idx),
      retained_genes = data.table::uniqueN(annotation$gene[idx])
    )
  })

  gene_rows <- lapply(unique(annotation$gene), function(gene_name) {
    idx <- which(annotation$gene == gene_name)
    component <- colSums(residual[idx, , drop = FALSE], na.rm = TRUE) / denominator
    interval <- bootstrap_share(component, total)
    data.table::data.table(
      layer = layer,
      gene = gene_name,
      module = annotation$module[idx][[1]],
      contribution_share = covariance_share(component, total),
      ci_low = interval[[1]],
      ci_high = interval[[2]],
      retained_features = length(idx),
      measured_values = sum(is.finite(residual[idx, , drop = FALSE]))
    )
  })

  list(
    modules = data.table::rbindlist(module_rows),
    genes = data.table::rbindlist(gene_rows),
    check = data.table::data.table(
      layer = layer,
      tumors = length(total),
      retained_features = nrow(annotation),
      residual_score_variance = stats::var(total),
      module_share_sum = sum(vapply(module_rows, function(x) x$contribution_share, numeric(1))),
      gene_share_sum = sum(vapply(gene_rows, function(x) x$contribution_share, numeric(1)))
    )
  )
}

layers <- lapply(c("nglyco_peptide", "nglyco_site"), decompose_layer)
module_results <- data.table::rbindlist(lapply(layers, `[[`, "modules"))
gene_results <- data.table::rbindlist(lapply(layers, `[[`, "genes"))
checks <- data.table::rbindlist(lapply(layers, `[[`, "check"))

data.table::setorder(module_results, layer, -contribution_share)
data.table::setorder(gene_results, layer, -contribution_share)
data.table::fwrite(module_results, file.path(out_dir, "residual_score_module_contributions.tsv"), sep = "\t")
data.table::fwrite(gene_results, file.path(out_dir, "residual_score_gene_contributions.tsv"), sep = "\t")
data.table::fwrite(checks, file.path(out_dir, "residual_score_contribution_checks.tsv"), sep = "\t")

cat("Wrote descriptive residual-score contribution tables to", out_dir, "\n")
