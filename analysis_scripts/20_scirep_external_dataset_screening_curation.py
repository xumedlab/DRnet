from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


ROWS = [
    {
        "dataset": "GSE160306",
        "organism": "Homo sapiens",
        "specimen_context": "postmortem retina, macula and retinal periphery",
        "assay_platform": "bulk RNA-seq, Illumina HiSeq 4000",
        "sample_summary": "79 GEO samples; macula-only discovery subset used 39 samples",
        "candidate_gene_overlap": "7/7 in discovery matrix",
        "same_compartment_macula_retina": "yes_for_discovery_subset",
        "stage_label_compatibility": "yes",
        "status": "used_as_discovery_dataset",
        "reason_or_use": "Primary macula-restricted stage-aware discovery cohort; not an external validation dataset.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160306",
    },
    {
        "dataset": "GSE102485",
        "organism": "Homo sapiens",
        "specimen_context": "PDR neovascular proliferative membrane, non-diabetic membrane controls, normal retina",
        "assay_platform": "bulk RNA-seq, Illumina NextSeq 500",
        "sample_summary": "30 GEO samples; directional comparison used 20 diabetic membrane, 5 non-diabetic membrane, and 3 normal retina samples",
        "candidate_gene_overlap": "7/7 checked",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "partial_PDR_only",
        "status": "used_as_cross_compartment_directional_support",
        "reason_or_use": "Full seven-gene coverage and human PDR relevance, but specimen type differs from macular retina; not formal same-compartment validation.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102485",
    },
    {
        "dataset": "GSE94019",
        "organism": "Homo sapiens",
        "specimen_context": "CD31+ endothelial cells from fibrovascular membranes and control retina",
        "assay_platform": "bulk RNA-seq, Illumina HiSeq 2000",
        "sample_summary": "13 samples: 9 fibrovascular membrane endothelial-cell samples and 4 retinal endothelial-cell controls",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "PDR_membrane_endothelial_only",
        "status": "screened_out_for_primary_expansion",
        "reason_or_use": "Human PDR-relevant endothelial-cell dataset, but cell-enriched membrane context is not comparable with bulk macular retina.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94019",
    },
    {
        "dataset": "GSE179568",
        "organism": "Homo sapiens",
        "specimen_context": "retinal neovascularization membranes, macular pucker, and macular hole controls",
        "assay_platform": "bulk RNA-seq, Illumina HiSeq 1000",
        "sample_summary": "24 samples: 7 PDR retinal neovascularization membranes, 10 macular pucker, and 7 macular hole samples",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "PDR_membrane_only",
        "status": "screened_out_for_primary_expansion",
        "reason_or_use": "Useful orthogonal membrane biology resource, but controls and specimen type do not match macular retinal tissue.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179568",
    },
    {
        "dataset": "GSE60436",
        "organism": "Homo sapiens",
        "specimen_context": "fibrovascular membranes from PDR and commercial retina controls",
        "assay_platform": "Illumina HumanWG-6 v3.0 expression beadchip",
        "sample_summary": "9 samples: 6 fibrovascular membranes and 3 retina controls",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "PDR_membrane_only",
        "status": "screened_out_for_primary_expansion",
        "reason_or_use": "Small array dataset in a fibrovascular membrane context; not suitable for same-compartment bulk macular RNA-seq replication.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60436",
    },
    {
        "dataset": "GSE53257",
        "organism": "Homo sapiens",
        "specimen_context": "DR, diabetes, and normal samples on custom mitoscriptome array",
        "assay_platform": "custom Human 8x15k mitoscriptome microarray",
        "sample_summary": "16 samples",
        "candidate_gene_overlap": "1/7 checked",
        "same_compartment_macula_retina": "unclear",
        "stage_label_compatibility": "limited",
        "status": "screened_out_for_gene_overlap",
        "reason_or_use": "Custom mitoscriptome platform has insufficient overlap with the current seven-gene candidate set.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53257",
    },
    {
        "dataset": "GSE221521",
        "organism": "Homo sapiens",
        "specimen_context": "peripheral blood leukocytes",
        "assay_platform": "bulk RNA-seq, Illumina NovaSeq 6000",
        "sample_summary": "193 samples: 50 healthy controls, 74 diabetic without DR, and 69 DR",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "blood_DR_vs_DM_available",
        "status": "screened_out_for_tissue_context",
        "reason_or_use": "Large human cohort but blood leukocytes are not comparable with macular retinal tissue; suitable only for separate circulating-biomarker exploration.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521",
    },
    {
        "dataset": "GSE185011",
        "organism": "Homo sapiens",
        "specimen_context": "peripheral blood mononuclear cells",
        "assay_platform": "bulk RNA-seq, Illumina NovaSeq 6000",
        "sample_summary": "25 samples: 5 each in HC, T2DM, DR, DPN, and DN groups",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "blood_DR_group_only",
        "status": "screened_out_for_tissue_context",
        "reason_or_use": "PBMC dataset with small DR subgroup; useful only for separate blood-compartment exploration, not retinal replication.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185011",
    },
    {
        "dataset": "GSE189005",
        "organism": "Homo sapiens",
        "specimen_context": "whole blood cells",
        "assay_platform": "Affymetrix Human Clariom D array",
        "sample_summary": "45 samples including control, T2DM, T2DR, and T2DN subgroups",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "blood_T2DR_subgroups_available",
        "status": "screened_out_for_tissue_context",
        "reason_or_use": "Whole-blood array dataset is not comparable with bulk macular retina and includes complication/time subgroups rather than retinal stage labels.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189005",
    },
    {
        "dataset": "GSE191210",
        "organism": "Homo sapiens",
        "specimen_context": "vitreous samples from PDR and idiopathic macular hole controls",
        "assay_platform": "Affymetrix Human Clariom D array and non-coding RNA profiling",
        "sample_summary": "9 samples: 3 controls, 3 anti-VEGF pretreated PDR, and 3 untreated PDR",
        "candidate_gene_overlap": "not_used_for_quantitative_extension",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "PDR_vitreous_only",
        "status": "screened_out_for_tissue_context",
        "reason_or_use": "Vitreous fluid array dataset is biologically adjacent but not retinal tissue; small sample size and treatment strata limit interpretation.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191210",
    },
    {
        "dataset": "GSE236333",
        "organism": "Rattus norvegicus",
        "specimen_context": "rat whole retina in STZ-diabetes and lixisenatide treatment model",
        "assay_platform": "bulk RNA-seq, Illumina HiSeq 4000",
        "sample_summary": "18 rat retina samples",
        "candidate_gene_overlap": "not_applicable_species_difference",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "animal_intervention_model",
        "status": "screened_out_for_species_and_design",
        "reason_or_use": "Animal intervention dataset can inform future mechanistic hypotheses but cannot increase the human macular discovery sample size.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236333",
    },
    {
        "dataset": "GSE178121",
        "organism": "Mus musculus",
        "specimen_context": "mouse retinal single-cell suspensions in STZ-diabetes model",
        "assay_platform": "single-cell RNA-seq, 10x Chromium",
        "sample_summary": "2 pooled samples: control retina and STZ diabetic retina",
        "candidate_gene_overlap": "not_applicable_species_and_scRNA",
        "same_compartment_macula_retina": "no",
        "stage_label_compatibility": "animal_model_only",
        "status": "screened_out_for_species_and_data_type",
        "reason_or_use": "Useful for future cell-type attribution, but mouse scRNA-seq cannot be merged with human bulk macular retina discovery data.",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178121",
    },
]


def write_screening_table() -> None:
    df = pd.DataFrame(ROWS)

    output_paths = [
        ROOT / "supplementary_tables" / "Supplementary_Table_S6_external_dataset_screening.csv",
        ROOT / "analysis_data" / "results_tables" / "external_comparison_dataset_screening.csv",
        ROOT / "analysis_data" / "results_tables" / "external_validation_dataset_screening.csv",
    ]
    for path in output_paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)

    workbook = ROOT / "Scientific_Reports_supplementary_tables.xlsx"
    if workbook.exists():
        with pd.ExcelWriter(workbook, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            df.to_excel(writer, sheet_name="S6_dataset_screen", index=False)


if __name__ == "__main__":
    write_screening_table()
