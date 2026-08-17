# DRnet: prioritizing P2RX4 in human diabetic retinopathy

DRnet is the analysis repository for the Research Article **“Donor-aware cross-cohort transcriptomics prioritizes P2RX4 in human diabetic retinopathy.”** It asks a focused question: which inflammation-related gene is most consistently prioritized by donor-aware analysis of human diabetic macula, and how well does that candidate hold up in separate ocular transcriptomic datasets?

The answer from this analysis is **P2RX4**. In the discovery cohort, P2RX4 ranked first for both the overall diabetic-retinopathy association and the DME-conditioned association. The follow-up datasets are deliberately treated as stress tests, not as proof of replication: the raw-count GSE276892 reconstruction was directionally compatible but imprecise, while the stronger GSE179568 membrane signal could not be separated cleanly from age and tissue-composition differences.

In practical terms, this repository supports **candidate prioritization**. It does not claim that P2RX4 is a diagnostic marker, that its RNA abundance measures receptor activity, or that the observed associations are causal.

## What is included

- `analysis_scripts/` contains the donor-aware discovery analysis, external P2RX4 checks, raw-count negative-binomial analysis, normal-retina localization, and the single command-line entry point.
- `project_inputs/` contains the frozen discovery inputs used by the analysis.
- `analysis_data/` contains the processed public-cohort inputs, accession manifests, author-derived cell-type mapping, and the checksum for the reconstructed-count archive.
- `analysis_results/` contains the compact tables and JSON summaries behind the reported results. Large all-gene intermediate tables are regenerated rather than stored in Git.
- `remote_linux_reconstruction/` contains the STAR/featureCounts workflow used to rebuild GSE276892 counts from 31 FASTQ lanes.
- `tests/` checks the statistical helpers, cohort mappings, reconstruction workflow, and release-independent analysis code.

Manuscript source, journal templates, rendered figures, submission PDFs, and duplicate archive manifests are intentionally kept out of this code repository.

## Quick start

The tested environment uses Python 3.10. From a fresh clone:

```bash
git clone https://github.com/xumedlab/DRnet.git
cd DRnet
python -m venv .venv
```

Activate the environment and install the pinned dependencies:

```bash
# Linux or macOS
source .venv/bin/activate

# Windows PowerShell
# .\.venv\Scripts\Activate.ps1

python -m pip install -r requirements.txt
```

The raw-count analysis uses a 444,319,883-byte reconstruction archive that is distributed separately on the matching [GitHub Releases](https://github.com/xumedlab/DRnet/releases) page. Place it at:

```text
analysis_data/independent_validation/remote_results/DRnet_GSE276892_remote_results.tar.gz
```

Expected SHA-256:

```text
CF6784D6852D30A7A2A67FDFC0C7FB93E0A3ABBF3207CE775F860709795D65A6
```

Run a short smoke test first:

```bash
python analysis_scripts/run_final_research_article_pipeline.py --quick
```

Run the formal analysis used for the article:

```bash
python analysis_scripts/run_final_research_article_pipeline.py
```

The formal run uses random seed `20260813`, 2,000 donor bootstraps, 4,999 studentized wild-bootstrap-t draws, and 10,000 external-cohort bootstrap draws. The `--quick` output is only a software check and must not be interpreted as the reported analysis.

## Test the code

```bash
python -m pytest -q
python -m ruff check analysis_scripts tests remote_linux_reconstruction/scripts
```

## Data and analysis flow

| Stage | Main input | Main output |
|---|---|---|
| Donor-aware discovery | GSE160306-derived frozen macular matrix and donor manifest | total and DME-conditioned rankings, bootstrap and leave-one-donor-out stability |
| Processed-cohort checks | GSE276892, GSE179568, and GSE94019 deposited data | P2RX4 effect estimates, robust sensitivities, QC associations, and overlap audit |
| Raw-count reconstruction | GSE276892 plus reused GSE147657 FASTQ lanes | STAR/featureCounts counts and PyDESeq2 negative-binomial models |
| Normal-retina context | six GSE130636 donor-region expression matrices plus the Voigt author mapping | donor-aggregated cell-class and regional localization summaries |

The discovery inputs and the smaller processed public-cohort files are stored in Git so that the main analysis is inspectable. FASTQ files, reference genomes, STAR indices, BAM files, and the 444 MB reconstructed-count archive are excluded because they are large and independently downloadable or reproducible from the included manifests and Linux workflow.

## Reading the results responsibly

- The primary discovery unit is the donor, not the individual tissue specimen.
- No primary estimand survived correction across all 158 screened genes.
- P2RX4 was stable to leave-one-donor-out analysis, but its bootstrap top-five frequency was 61.65%; this is prioritization under uncertainty, not a fixed multigene signature.
- GSE276892 includes seven reused controls and provides limited source overlap. Its raw-count disease estimate was imprecise (`log2 fold change = 0.295`, two-sided `P = 0.615`).
- GSE179568 has no age common support between the main disease and control groups. Adjusted coefficients therefore depend on extrapolation.
- Separate GEO accessions do not prove that patient identities do not overlap.

These limitations are part of the result, not post hoc exclusions.

## Rebuilding the raw counts

The exact Linux workflow is documented in `remote_linux_reconstruction/README.md`. It uses STAR 2.7.8a, featureCounts 2.0.1, and GENCODE GRCh38.p13 release 42 in two reference configurations. Technical lanes are counted separately and then summed to 17 biological samples before statistical modelling.

## Questions

For reproducibility questions or suspected code defects, please open a GitHub issue and include the command used, operating system, Python version, and full error message.
