main_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    stop("Run analysis scripts with Rscript.", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

script_dir <- dirname(main_script_path())
project_root_env <- Sys.getenv("PDAC_PROJECT_ROOT", unset = "")
project_candidates <- c(
  project_root_env,
  file.path(script_dir, "..", ".."),
  file.path(script_dir, "..")
)
project_candidates <- unique(project_candidates[nzchar(project_candidates)])
project_matches <- project_candidates[dir.exists(project_candidates)]
if (length(project_matches) == 0) {
  stop("Set PDAC_PROJECT_ROOT to an existing project directory.", call. = FALSE)
}
ROOT <- normalizePath(project_matches[[1]], winslash = "/", mustWork = TRUE)
matrisome_candidates <- c(
  file.path(script_dir, "resources", "matrisome_experimental_current.tsv"),
  file.path(ROOT, "data", "external", "matrisome", "matrisome_experimental_current.tsv")
)
matrisome_matches <- matrisome_candidates[file.exists(matrisome_candidates)]
if (length(matrisome_matches) == 0) {
  stop("The versioned Matrisome annotation is unavailable.", call. = FALSE)
}
MATRISOME_FILE <- normalizePath(matrisome_matches[[1]], winslash = "/", mustWork = TRUE)

database_root_env <- Sys.getenv("PDAC_DATABASE_ROOT", unset = "")
database_candidates <- c(
  database_root_env,
  file.path(ROOT, "Database"),
  file.path(ROOT, "..", "Database"),
  file.path(ROOT, "..", "..", "Database")
)
database_candidates <- unique(database_candidates[nzchar(database_candidates)])
database_matches <- database_candidates[dir.exists(database_candidates)]
if (length(database_matches) == 0) {
  stop("Set PDAC_DATABASE_ROOT to the directory containing GEO, CPTAC, and TCGA-PAAD.", call. = FALSE)
}
DATABASE_ROOT <- normalizePath(database_matches[[1]], winslash = "/", mustWork = TRUE)
GEO_PDAC_ROOT <- file.path(DATABASE_ROOT, "GEO", "PDAC")
CPTAC_ROOT <- file.path(DATABASE_ROOT, "CPTAC")
TCGA_PAAD_ROOT <- file.path(DATABASE_ROOT, "TCGA-PAAD")
RAW_LOCAL <- file.path(ROOT, "data", "raw")
PROCESSED <- file.path(ROOT, "data", "processed")
RESULTS <- file.path(ROOT, "results", "tables")
FIGURES <- file.path(ROOT, "results", "figures")

dataset_dir <- function(kind, dataset, create = FALSE) {
  path <- switch(
    kind,
    "geo_pdac" = file.path(GEO_PDAC_ROOT, dataset),
    "cptac" = file.path(CPTAC_ROOT, dataset),
    "tcga_paad" = if (identical(dataset, "UCSC_Xena")) TCGA_PAAD_ROOT else file.path(TCGA_PAAD_ROOT, dataset),
    "local_raw" = file.path(RAW_LOCAL, dataset),
    stop(sprintf("Unsupported dataset kind: %s", kind), call. = FALSE)
  )
  if (create) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

existing_dataset_dir <- function(kind, dataset) {
  path <- dataset_dir(kind, dataset)
  if (!dir.exists(path)) {
    stop(sprintf("Required dataset %s was not found in the configured database directory.", dataset), call. = FALSE)
  }
  path
}

ensure_project_dirs <- function() {
  dir.create(PROCESSED, recursive = TRUE, showWarnings = FALSE)
  dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
  dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
}

EXPECTED_R_VERSION <- "4.5.3"
EXPECTED_R_PACKAGE_VERSIONS <- c(
  "data.table" = "1.18.2.1",
  "survival" = "3.8.6",
  "readxl" = "1.4.5",
  "Matrix" = "1.7.5",
  "FNN" = "1.1.4.1",
  "CellChat" = "2.2.0.9001",
  "digest" = "0.6.39",
  "SeuratObject" = "5.3.0"
)

ensure_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required R packages: %s. Restore the frozen environment before running the analysis.", paste(missing, collapse = ", ")), call. = FALSE)
  }
  versioned <- intersect(packages, names(EXPECTED_R_PACKAGE_VERSIONS))
  mismatched <- versioned[vapply(versioned, function(package) {
    as.character(utils::packageVersion(package)) != EXPECTED_R_PACKAGE_VERSIONS[[package]]
  }, logical(1))]
  if (length(mismatched) > 0) {
    details <- vapply(mismatched, function(package) {
      sprintf("%s=%s expected=%s", package, as.character(utils::packageVersion(package)), EXPECTED_R_PACKAGE_VERSIONS[[package]])
    }, character(1))
    stop(sprintf("R package version mismatch: %s", paste(details, collapse = "; ")), call. = FALSE)
  }
  invisible(TRUE)
}

read_gene_set <- function(program) {
  path <- file.path(script_dir, "candidate_gene_sets.tsv")
  sets <- data.table::fread(path)
  selected_program <- program
  genes <- sets[program == selected_program][order(order_index), gene]
  if (length(genes) == 0) {
    stop(sprintf("Gene set not found: %s", program), call. = FALSE)
  }
  genes
}

read_phenotype_gene_sets <- function() {
  path <- file.path(script_dir, "phenotype_gene_sets.tsv")
  sets <- data.table::fread(path)
  signatures <- unique(sets$signature)
  stats::setNames(lapply(signatures, function(name) {
    selected_name <- name
    sets[signature == selected_name][order(order_index), gene]
  }), signatures)
}

zscore_vector <- function(x, zero_if_constant = FALSE) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    fill <- if (zero_if_constant) 0 else NA_real_
    return(rep(fill, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

zscore_rows_matrix <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, stats::sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- NA_real_
  sweep(sweep(mat, 1, means, "-"), 1, sds, "/")
}

clean_name <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  value <- gsub("_+", "_", value)
  gsub("^_|_$", "", value)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

read_xena_tcga_expression <- function(genes = NULL) {
  raw_dir <- existing_dataset_dir("tcga_paad", "UCSC_Xena")
  expression <- data.table::fread(file.path(raw_dir, "TCGA-PAAD.star_tpm.tsv"))
  mapping <- data.table::fread(
    file.path(raw_dir, "gencode.v36.annotation.gtf.gene.probemap"),
    select = c("id", "gene")
  )
  data.table::setnames(expression, 1, "id")
  if (!is.null(genes)) {
    mapping <- mapping[gene %in% genes]
  }
  expression <- expression[mapping, on = "id", nomatch = 0]
  sample_columns <- setdiff(names(expression), c("id", "gene"))
  sample_columns <- sample_columns[substr(sample_columns, 14, 15) == "01"]
  expression <- expression[!is.na(gene) & nzchar(gene), c("gene", sample_columns), with = FALSE]
  expression[, (sample_columns) := lapply(.SD, as.numeric), .SDcols = sample_columns]
  expression[, lapply(.SD, safe_median), by = gene, .SDcols = sample_columns]
}
