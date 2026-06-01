# Suggested Bioinformatics Sections For The Manuscript

## Methods 2.1-2.6

Transcriptomic data were obtained from GEO under accession number GSE160306. Only macular samples were retained for formal bioinformatic analysis. The final cohort included 10 healthy controls, 10 diabetic samples without retinopathy progression, 9 NPDR samples, and 10 NPDR/PDR + DME samples. Raw count matrices were used for differential expression analysis, whereas log2(CPM+1) values recomputed from the same raw counts were used for PCA, ssGSEA, LASSO modeling, and visualization.

Differential expression for the primary comparison (healthy control vs NPDR/PDR + DME) and the supportive pairwise comparisons was performed using DESeq2 implemented through the PyDESeq2 framework in Python. Genes with adjusted P < 0.05 and |log2 fold change| >= 0.5 were considered significant. A severity-trend analysis across the four disease stages was additionally performed using Spearman correlation, and positively correlated genes with FDR < 0.1 were retained as progression-related candidates.

Inflammation-related activity was quantified by ssGSEA using the HALLMARK_INFLAMMATORY_RESPONSE gene set from MSigDB. Inflammation-associated candidate genes were defined as the intersection between significant primary DEGs and the Hallmark inflammatory response gene set. LASSO logistic regression with repeated stratified cross-validation was then applied to the candidate genes in the primary comparison cohort, and model performance was summarized with out-of-fold predicted probabilities and ROC analysis.

For functional interpretation, each selected gene was used as an anchor for Spearman-ranked preranked GSEA against the Hallmark collection. Relative immune activity was estimated with ssGSEA using 28 curated immune cell signatures, and group differences as well as gene-immune correlations were evaluated with nonparametric statistics and multiple-testing correction.

## Results 3.1-3.5

After restricting the analysis to macular samples, 39 samples were included in the final transcriptomic cohort, consisting of 10 healthy controls, 10 diabetic samples, 9 NPDR samples, and 10 advanced DR samples (NPDR/PDR + DME). DESeq2 identified 468 significant genes in the primary contrast of healthy control versus advanced DR. Intersecting these DEGs with the Hallmark inflammatory response signature yielded 11 inflammation-associated core genes (MSR1, TIMP1, OPRK1, LYN, FZD5, TLR3, NDP, NMI, CLEC5A, CYBB, CMKLR1). Among them, 8 genes also showed a positive severity trend across the four clinical groups (MSR1, TIMP1, OPRK1, LYN, FZD5, TLR3, NMI, CMKLR1).

LASSO logistic regression retained 7 genes in the final signature (MSR1, NMI, FZD5, TIMP1, CMKLR1, LYN, TLR3). The combined signature achieved an out-of-fold AUC of 0.870, indicating moderate-to-good discrimination within the discovery cohort while still requiring independent external validation. The inflammation ssGSEA score showed a primary-comparison P value of 0.0963, and the four-group severity trend yielded a Spearman rho of 0.340.

Gene-centered preranked GSEA highlighted recurrent Hallmark programs related to INFLAMMATORY RESPONSE, EPITHELIAL MESENCHYMAL TRANSITION, INTERFERON ALPHA RESPONSE, INTERFERON GAMMA RESPONSE, COAGULATION. Immune ssGSEA further suggested that the most perturbed immune signatures included Immature_B_cell, Macrophage, Natural_killer_cell, Memory_B_cell, Regulatory_T_cell, supporting an association between the selected inflammatory genes and the retinal immune microenvironment.

## Do Not Auto-Fill

The current repository does not contain the raw data for clinical validation, ELISA/RT-qPCR assays, or cell experiments. Therefore, Methods 2.7 onward and Results 3.6 onward should be completed only after the corresponding wet-lab data are provided.
