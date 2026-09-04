from pathlib import Path
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


CANDIDATE_GENES = ["POSTN", "MXRA5", "PXDN", "MFAP4", "LTBP2", "TNC"]
MIN_NUCLEI = 20
MIN_CORRELATION_NUCLEI = 50
BOOTSTRAPS = 2000
SEED = 20260719
SIGNATURES = [
    ("MALIGNANT CELLS", "Malignant cells", "Epithelial"),
    ("DUCTAL", "Ductal", "Epithelial"),
    ("myCAF", "myCAF", "CAF"),
    ("Tuveson_mCAF", "Tuveson mCAF", "CAF"),
    ("iCAF", "iCAF", "CAF"),
    ("Tuveson_iCAF", "Tuveson iCAF", "CAF"),
    ("apCAF", "apCAF", "CAF"),
    ("CAF", "CAF", "CAF"),
    ("PanCAF", "PanCAF", "CAF"),
    ("FIBROBLASTS", "Fibroblasts", "CAF"),
    ("Pan_Immune", "Pan-immune", "Immune"),
    ("IMMUNE", "Immune", "Immune"),
    ("M2", "M2", "Immune"),
    ("Macrophage", "Macrophage", "Immune"),
    ("CD8_Tcells", "CD8 T cells", "Immune"),
    ("CD4_regulatory", "Regulatory CD4 T cells", "Immune"),
    ("ENDOTHELIAL", "Endothelial", "Endothelial"),
]


def first_column(columns, candidates):
    lower = {str(column).lower(): column for column in columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def normalize_context(value):
    label = str(value).strip().upper()
    if "FIBRO" in label or label in {"CAF", "CAFS"}:
        return "CAF"
    if "ENDOTHEL" in label:
        return "Endothelial"
    if "IMMUNE" in label:
        return "Immune"
    if any(token in label for token in ["MALIGNANT", "DUCTAL", "ACINAR", "EPITHEL"]):
        return "Epithelial"
    return np.nan


def normalize_status(value, patient):
    label = str(value).strip().lower() if pd.notna(value) else ""
    patient_label = str(patient).strip().lower()
    if any(token in label for token in ["naive", "untreated", "treatment-naive"]):
        return "Untreated"
    if label in {"crt", "crtl", "crtx", "crtn", "rt", "gart"} or any(token in label for token in ["treated", "therapy"]):
        return "Treated"
    if any(token in label for token in ["normal", "control", "non-malignant"]):
        return "Non-tumor"
    if any(token in patient_label for token in ["normal", "control"]):
        return "Non-tumor"
    if "_u_" in patient_label or patient_label.startswith("u_") or patient_label.startswith("pdac_u"):
        return "Untreated"
    if "_t_" in patient_label or patient_label.startswith("t_") or patient_label.startswith("pdac_t"):
        return "Treated"
    if patient_label.startswith("u") and patient_label[1:].isdigit():
        return "Untreated"
    if patient_label.startswith("t") and patient_label[1:].isdigit():
        return "Treated"
    return "Unknown"


def bootstrap_ci(values, rng, statistic=np.mean):
    values = np.asarray(values, dtype=float)
    estimates = np.empty(BOOTSTRAPS, dtype=float)
    for index in range(BOOTSTRAPS):
        estimates[index] = statistic(rng.choice(values, size=len(values), replace=True))
    return np.quantile(estimates, [0.025, 0.975])


def main():
    project_root = Path(os.environ.get("PDAC_PROJECT_ROOT", Path(__file__).resolve().parents[2]))
    database_root_value = os.environ.get("PDAC_DATABASE_ROOT")
    if not database_root_value:
        raise RuntimeError("Set PDAC_DATABASE_ROOT to the directory containing GEO, CPTAC, and TCGA-PAAD")
    database_root = Path(database_root_value)
    input_path = database_root / "GEO" / "PDAC" / "GSE202051" / "GSE202051_totaldata-final-toshare.h5ad"
    output_dir = project_root / "results" / "tables" / "analysis_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_path, backed="r")
    obs = adata.obs.copy()
    patient_column = first_column(obs.columns, ["pid", "patient", "patient_id", "sample", "sample_id"])
    context_column = first_column(obs.columns, ["broad_celltypes", "cells_labels", "cell_type", "celltype", "broad_celltype"])
    status_column = first_column(obs.columns, ["treatment_status", "new_treatment", "status", "treatment", "condition", "therapy"])
    if patient_column is None or context_column is None:
        raise RuntimeError(f"Required patient or cell-context columns are absent: {list(obs.columns)}")

    patient = obs[patient_column].astype(str)
    context = obs[context_column].map(normalize_context)
    if status_column is None:
        status_values = pd.Series(np.nan, index=obs.index)
    else:
        status_values = obs[status_column]
    status = pd.Series([normalize_status(value, pid) for value, pid in zip(status_values, patient)], index=obs.index)
    keep = context.notna() & status.isin(["Untreated", "Treated"])

    var_names = pd.Index(adata.var_names.astype(str))
    found_genes = [gene for gene in CANDIDATE_GENES if gene in var_names]
    if len(found_genes) != len(CANDIDATE_GENES):
        raise RuntimeError(f"Candidate genes missing from GSE202051: {sorted(set(CANDIDATE_GENES) - set(found_genes))}")
    gene_indices = [var_names.get_loc(gene) for gene in found_genes]
    selected_indices = np.flatnonzero(keep.to_numpy())
    if len(selected_indices) == 0:
        raise RuntimeError("No PDAC nuclei remained after cell-context and treatment-status filtering")
    expression = adata[:, gene_indices].X
    if sparse.issparse(expression):
        expression = expression.toarray()
    expression = np.asarray(expression, dtype=float)
    expression = expression[selected_indices, :]
    transform = "dataset_log1p_normalized"
    if np.nanmax(expression) > 30:
        expression = np.log1p(expression)
        transform = "log1p_applied_to_dataset_normalized_counts"

    module = np.nanmean(expression, axis=1)
    module_sd = np.nanstd(module, ddof=1)
    if not np.isfinite(module_sd) or module_sd == 0:
        raise RuntimeError("Candidate module has zero or undefined variance")
    module_z = (module - np.nanmean(module)) / module_sd
    selected_obs = pd.DataFrame(
        {
            "patient": patient.iloc[selected_indices].to_numpy(),
            "treatment_status": status.iloc[selected_indices].to_numpy(),
            "cell_context": context.iloc[selected_indices].to_numpy(),
            "candidate_module_z": module_z,
            "candidate_detected": np.any(expression > 0, axis=1).astype(float),
        }
    )
    missing_signatures = [name for name, _, _ in SIGNATURES if name not in obs.columns]
    if missing_signatures:
        raise RuntimeError(f"Author-provided GSE202051 signatures are absent: {missing_signatures}")
    for signature, _, _ in SIGNATURES:
        selected_obs[signature] = pd.to_numeric(obs.iloc[selected_indices][signature], errors="coerce").to_numpy()

    patient_context = (
        selected_obs.groupby(["patient", "treatment_status", "cell_context"], observed=True)
        .agg(
            n_nuclei=("candidate_module_z", "size"),
            candidate_module_z_mean=("candidate_module_z", "mean"),
            candidate_module_z_median=("candidate_module_z", "median"),
            detected_fraction=("candidate_detected", "mean"),
        )
        .reset_index()
    )
    patient_context["eligible"] = patient_context["n_nuclei"] >= MIN_NUCLEI
    patient_context.to_csv(output_dir / "gse202051_patient_cell_context.tsv", sep="\t", index=False)

    eligible = patient_context.loc[patient_context["eligible"]].copy()
    rng = np.random.default_rng(SEED)
    summary_rows = []
    for cell_context, frame in eligible.groupby("cell_context", observed=True):
        low, high = bootstrap_ci(frame["candidate_module_z_mean"].to_numpy(), rng)
        summary_rows.append(
            {
                "cell_context": cell_context,
                "n_patients": frame["patient"].nunique(),
                "n_nuclei": int(frame["n_nuclei"].sum()),
                "mean_patient_score": frame["candidate_module_z_mean"].mean(),
                "median_patient_score": frame["candidate_module_z_mean"].median(),
                "ci_low": low,
                "ci_high": high,
                "mean_detected_fraction": frame["detected_fraction"].mean(),
            }
        )
    cell_summary = pd.DataFrame(summary_rows).sort_values("mean_patient_score", ascending=False)
    cell_summary.to_csv(output_dir / "gse202051_cell_context_summary.tsv", sep="\t", index=False)

    wide = eligible.pivot_table(
        index=["patient", "treatment_status"],
        columns="cell_context",
        values="candidate_module_z_mean",
        aggfunc="mean",
    ).reset_index()
    comparison_targets = {"Epithelial": "CAF_vs_Epithelial", "Immune": "CAF_vs_Immune", "Endothelial": "CAF_vs_Endothelial"}
    contrast_rows = []
    for target, comparison in comparison_targets.items():
        if "CAF" not in wide.columns or target not in wide.columns:
            continue
        frame = wide[["patient", "treatment_status", "CAF", target]].dropna().copy()
        frame["difference"] = frame["CAF"] - frame[target]
        if len(frame) < 3:
            continue
        low, high = bootstrap_ci(frame["difference"].to_numpy(), rng)
        test = wilcoxon(frame["difference"], alternative="two-sided", zero_method="wilcox", method="auto")
        untreated = frame.loc[frame["treatment_status"] == "Untreated", "difference"]
        treated = frame.loc[frame["treatment_status"] == "Treated", "difference"]
        contrast_rows.append(
            {
                "comparison": comparison,
                "n_patients": len(frame),
                "mean_difference": frame["difference"].mean(),
                "median_difference": frame["difference"].median(),
                "ci_low": low,
                "ci_high": high,
                "wilcoxon_p": test.pvalue,
                "positive_patients": int((frame["difference"] > 0).sum()),
                "positive_fraction": float((frame["difference"] > 0).mean()),
                "untreated_n": len(untreated),
                "untreated_mean": untreated.mean() if len(untreated) else np.nan,
                "treated_n": len(treated),
                "treated_mean": treated.mean() if len(treated) else np.nan,
            }
        )
    contrasts = pd.DataFrame(contrast_rows)
    if contrasts.empty:
        raise RuntimeError("No eligible patient-paired GSE202051 contrasts were available")
    contrasts["fdr"] = multipletests(contrasts["wilcoxon_p"], method="fdr_bh")[1]
    contrasts.to_csv(output_dir / "gse202051_patient_paired_contrasts.tsv", sep="\t", index=False)

    signature_rows = []
    for (patient_id, treatment_status), frame in selected_obs.groupby(["patient", "treatment_status"], observed=True):
        candidate_values = frame["candidate_module_z"].to_numpy(dtype=float)
        for signature, display_label, family in SIGNATURES:
            signature_values = frame[signature].to_numpy(dtype=float)
            complete = np.isfinite(candidate_values) & np.isfinite(signature_values)
            n_complete = int(complete.sum())
            rho = np.nan
            if n_complete >= MIN_CORRELATION_NUCLEI:
                candidate_complete = candidate_values[complete]
                signature_complete = signature_values[complete]
                if np.nanstd(candidate_complete) > 0 and np.nanstd(signature_complete) > 0:
                    rho = spearmanr(candidate_complete, signature_complete).statistic
            signature_rows.append(
                {
                    "patient": patient_id,
                    "treatment_status": treatment_status,
                    "signature": signature,
                    "display_label": display_label,
                    "signature_family": family,
                    "n_nuclei": n_complete,
                    "spearman_rho": rho,
                }
            )
    patient_signatures = pd.DataFrame(signature_rows)
    patient_signatures.to_csv(output_dir / "gse202051_patient_signature_correlations.tsv", sep="\t", index=False)

    signature_summary_rows = []
    for signature, display_label, family in SIGNATURES:
        frame = patient_signatures.loc[
            (patient_signatures["signature"] == signature) & patient_signatures["spearman_rho"].notna()
        ].copy()
        values = frame["spearman_rho"].to_numpy(dtype=float)
        if len(values) < 3:
            continue
        low, high = bootstrap_ci(values, rng, statistic=np.median)
        if np.allclose(values, 0):
            p_value = 1.0
        else:
            p_value = wilcoxon(values, alternative="two-sided", zero_method="wilcox", method="auto").pvalue
        untreated = frame.loc[frame["treatment_status"] == "Untreated", "spearman_rho"]
        treated = frame.loc[frame["treatment_status"] == "Treated", "spearman_rho"]
        signature_summary_rows.append(
            {
                "signature": signature,
                "display_label": display_label,
                "signature_family": family,
                "n_patients": len(frame),
                "median_rho": np.median(values),
                "mean_rho": np.mean(values),
                "ci_low": low,
                "ci_high": high,
                "wilcoxon_p": p_value,
                "positive_patients": int((values > 0).sum()),
                "positive_fraction": float((values > 0).mean()),
                "untreated_n": len(untreated),
                "untreated_median": untreated.median() if len(untreated) else np.nan,
                "treated_n": len(treated),
                "treated_median": treated.median() if len(treated) else np.nan,
            }
        )
    signature_summary = pd.DataFrame(signature_summary_rows)
    if signature_summary.empty:
        raise RuntimeError("No eligible patient-stratified GSE202051 signature correlations were available")
    signature_summary["fdr"] = multipletests(signature_summary["wilcoxon_p"], method="fdr_bh")[1]
    signature_summary.to_csv(output_dir / "gse202051_signature_correlation_summary.tsv", sep="\t", index=False)

    gate_rows = []
    for comparison in ["CAF_vs_Epithelial", "CAF_vs_Immune"]:
        row = contrasts.loc[contrasts["comparison"] == comparison]
        if row.empty:
            gate_rows.append({"criterion": f"{comparison}_available", "passed": False, "value": "absent", "threshold": "required"})
            continue
        row = row.iloc[0]
        checks = [
            ("n_patients", row["n_patients"] >= 8, row["n_patients"], ">=8"),
            ("mean_difference", row["mean_difference"] > 0, row["mean_difference"], ">0"),
            ("positive_fraction", row["positive_fraction"] >= 0.75, row["positive_fraction"], ">=0.75"),
            ("fdr", row["fdr"] < 0.05, row["fdr"], "<0.05"),
            ("treatment_direction", row["untreated_mean"] > 0 and row["treated_mean"] > 0, f"{row['untreated_mean']:.6g};{row['treated_mean']:.6g}", "both >0"),
        ]
        for name, passed, value, threshold in checks:
            gate_rows.append({"criterion": f"{comparison}_{name}", "passed": bool(passed), "value": value, "threshold": threshold})
    gate = pd.DataFrame(gate_rows)
    overall_pass = bool(gate["passed"].all())
    gate = pd.concat(
        [gate, pd.DataFrame([{"criterion": "overall_reanalysis_gate", "passed": overall_pass, "value": overall_pass, "threshold": "all core criteria pass"}])],
        ignore_index=True,
    )
    gate.to_csv(output_dir / "gse202051_reanalysis_gate.tsv", sep="\t", index=False)

    audit = pd.DataFrame(
        [
            {
                "dataset": "GSE202051",
                "input_file": input_path.name,
                "n_total_nuclei": adata.n_obs,
                "n_selected_pdac_nuclei": len(selected_obs),
                "n_selected_patients": selected_obs["patient"].nunique(),
                "patient_column": patient_column,
                "cell_context_column": context_column,
                "status_column": status_column if status_column is not None else "inferred_from_patient",
                "expression_transform": transform,
                "candidate_genes": ",".join(found_genes),
                "minimum_nuclei_per_patient_context": MIN_NUCLEI,
                "minimum_nuclei_per_patient_signature_correlation": MIN_CORRELATION_NUCLEI,
                "author_provided_signatures": len(SIGNATURES),
                "minimum_patients_per_signature": int(signature_summary["n_patients"].min()),
                "reanalysis_gate_passed": overall_pass,
            }
        ]
    )
    audit.to_csv(output_dir / "gse202051_full_object_audit.tsv", sep="\t", index=False)
    adata.file.close()
    print(f"GSE202051 reanalysis gate: {'PASS' if overall_pass else 'FAIL'}")


if __name__ == "__main__":
    main()
