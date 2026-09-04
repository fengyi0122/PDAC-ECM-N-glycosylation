PDAC GLYCOSYLATED ECM ANALYSIS CODE

This directory contains the complete core R and Python code for the PDAC glycosylated ECM analyses. Presentation-layer formatting is outside the computational scope.

CONFIGURATION

Set PDAC_PROJECT_ROOT to an existing project directory where processed data and result tables will be written.
Set PDAC_DATABASE_ROOT to the database directory containing GEO/PDAC, CPTAC, and TCGA-PAAD.

The required database layout is:
PDAC_DATABASE_ROOT/CPTAC/CPTAC-PDAC
PDAC_DATABASE_ROOT/TCGA-PAAD
PDAC_DATABASE_ROOT/GEO/PDAC/GSE154778
PDAC_DATABASE_ROOT/GEO/PDAC/GSE199102
PDAC_DATABASE_ROOT/GEO/PDAC/GSE202051
PDAC_DATABASE_ROOT/GEO/PDAC/GSE272362

ENVIRONMENT

Use R 4.5.3 and Python 3.12.10. Install renv, restore the R environment from renv.lock, create an isolated Python environment, and install requirements.txt before running the analysis. SeuratObject 5.3.0 is required to read the GSE272362 Seurat object. CellChat 2.2.0.9001 is used for patient-stratified candidate ECM receptor-context inference. Exact software versions are listed in environment_versions.tsv.

R environment setup:
Rscript -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::restore()"

Python environment setup:
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

Set PDAC_PYTHON to the Python executable inside the isolated environment when it is not available as python.

INPUTS

Required dataset accessions, filenames, source URLs, byte sizes, and SHA-256 values are listed in input_files.tsv. Public datasets are not redistributed with the code. The runner verifies every input before starting the analysis.
The UCSC Xena TCGA-PAAD STAR-TPM download is distributed as a gzip archive and should be decompressed to the filename listed in input_files.tsv. The associated GENCODE v36 probe map is used directly for Ensembl-to-gene-symbol mapping.

GENE SETS

candidate_gene_sets.tsv freezes the six-gene external RNA module, the 24-gene Visium spatial module, and the ECM module definitions used for descriptive residual-score decomposition. phenotype_gene_sets.tsv freezes the CPTAC primary phenotype definitions shared by the primary and TMT11 plex models. dataset_signature_gene_sets.tsv records the complete dataset-specific marker definitions used in CPTAC, TCGA-PAAD, GSE199102, and GSE154778 analyses. fdr_analysis_families.tsv records each Benjamini-Hochberg correction family, its correction unit, the models sharing the correction, the number of tests, and the source result file. The Visium module contains 24 genes, of which 23 are detected in GSE272362. The versioned Matrisome annotation used by the analysis is bundled under resources.

EXECUTION

Run Rscript 00_run_all.R from this directory. The runner validates input identity and all required R and Python versions, rebuilds analysis intermediates, and executes the listed R and Python scripts in dependency order. The GSE154778 helper streams the dense expression matrix to recover CellChat-required genes and full-library sizes before the R interaction analysis. Individual scripts can also be run separately after their upstream outputs have been generated.

OUTPUTS

Processed scores are written to data/processed. Statistical result tables are written to results/tables and results/tables/analysis_results. The analysis outputs include the 135-PDAC primary and all-140-tumor sensitivity models, feature-level detection diagnostics, four detection-IPW specifications, predictor-side and outcome-side POSTN-removal sensitivity, continuous and median-split survival analyses, GSE202051 patient-level analyses, patient-level GeoMx analyses, MAN2A2 and spatial-neighborhood analyses, descriptive residual-score decomposition, and patient-stratified CellChat candidate ECM receptor-context and candidate-set enrichment analyses. In c4_candidate_core_set.tsv, cptac_positive_feature_outcome_associations counts positive model rows with beta > 0 and BH FDR < 0.10 across the assessed CPTAC glycofeature layer-outcome combinations; it is not a count of unique glycofeatures. Detection diagnostics use the at-least-40-observation robustness feature universe, whereas the primary residual scores use the at-least-50-observation feature set. The runner verifies that all key outputs exist and are nonempty before successful completion.

CHECKSUMS

SHA256SUMS.txt records the SHA-256 checksum of every file in this release except the checksum file itself.

LICENSE

The original R and Python code and project-authored documentation and configuration files are released under the MIT License. See LICENSE. The bundled Matrisome annotation is third-party-derived material and is not relicensed under the MIT License; see THIRD_PARTY_NOTICES.md.

CITATION

If you use or adapt this code, or if this code contributes materially to analyses or results reported in a publication, please cite both the versioned software release and the associated article. Citation metadata for the software release are provided in CITATION.cff and on the Zenodo record.

Modified or derivative implementations that retain or substantially follow the analytical framework of this repository should clearly acknowledge their derivation from this codebase and cite the corresponding software version. Use the published version of record for the associated article citation when it becomes available.

PUBLIC RELEASE

GitHub repository: https://github.com/fengyi0122/PDAC-ECM-N-glycosylation
Version: v1.0.0
GitHub release: https://github.com/fengyi0122/PDAC-ECM-N-glycosylation/releases/tag/v1.0.0
Zenodo DOI: https://doi.org/10.5281/zenodo.21536074

Version v1.0.0 contains the core R and Python analysis code associated with the article "Protein-residual ECM N-glycosylation delineates a composition-resolved myCAF-enriched stromal state in pancreatic ductal adenocarcinoma." The versioned GitHub release is archived in Zenodo under the DOI listed above.
