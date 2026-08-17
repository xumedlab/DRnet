# Locked independent validation protocol for P2RX4

Protocol frozen on 2026-08-13 before downloading or inspecting gene-level values from GSE276892.

## Fixed target and rationale

- The only confirmatory gene is `P2RX4`, selected as the rank-1 candidate in the GSE160306 discovery analysis.
- No other gene will be promoted on the basis of GSE276892.
- GSE276892 did not participate in candidate selection, model development, threshold choice, or cell-type localization.

## Validation cohort

- Dataset: GSE276892.
- Assay: bulk RNA sequencing of fluorescence-sorted human vitreous hyalocytes.
- Disease group: eight independent patients with proliferative diabetic retinopathy (PDR).
- Control group: nine independent patients with macular hole or macular pucker, including seven control profiles reused from GSE147657 and two newly deposited controls.
- Primary analysis unit: patient. Technical lanes, if present, must be summed before normalization and must never be treated as independent observations.

## Fixed primary hypothesis and analysis

- Directional hypothesis: `P2RX4` expression is higher in PDR hyalocytes than in control hyalocytes.
- Primary effect: log2 fold change from a negative-binomial count model with disease group as the sole design term.
- Primary test: one-sided Wald test for a positive disease coefficient.
- Success criterion: positive log2 fold change and one-sided P < 0.05.
- Because exactly one gene was fixed before data inspection, no across-gene multiplicity correction applies to this confirmatory test.

## Fixed sensitivity analyses

- Two-sided Wald P value and 95% confidence interval.
- Patient-level normalized-expression difference with a two-sided Mann-Whitney test and Cliff's delta.
- Leave-one-patient-out direction stability.
- Separate comparisons against macular-hole and macular-pucker controls, reported as descriptive sensitivities because each subgroup is small.
- If age and sex are available for all patients, an age/sex-adjusted negative-binomial sensitivity model will be added; otherwise missing covariates will be reported and not imputed.

## Failure and interpretation rules

- A non-positive effect or one-sided P >= 0.05 is a failed independent validation and must be reported as such.
- A positive but non-significant effect is directional support only, not validation.
- A significant result in another candidate cannot replace a failed P2RX4 test.
- This compartment-specific result may support disease-state regulation in ocular myeloid cells, but it cannot by itself establish whole-retina abundance, microglial exclusivity, causality, diagnostic performance, or therapeutic efficacy.

## Provenance

- GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276892
- Primary article: https://pubmed.ncbi.nlm.nih.gov/39543723/
- Discovery cohort: GSE160306.
