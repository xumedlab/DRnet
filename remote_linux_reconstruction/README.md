# Rebuilding GSE276892 counts on Linux

This workflow reconstructs gene-level counts for the GSE276892 validation dataset from the original FASTQ files. It uses the official Linux builds of STAR 2.7.8a and Subread/featureCounts 2.0.1. The formal analysis was run on Linux because the supported STAR release does not provide a native Windows binary comparable to the authors' environment.

Two workflows are run:

1. `all_regions_source_reconstruction` uses the GENCODE GRCh38.p13 release 42 all-regions FASTA and GTF. This is the closest available reconstruction of the source description.
2. `primary_assembly_sensitivity` uses the release 42 primary-assembly FASTA and GTF as a standard-reference sensitivity analysis.

Both workflows run unstranded exon-level featureCounts aggregation by `gene_id`. The 21 GSE147657 FASTQ files are three technical lanes from each of seven reused controls. Lanes are processed separately and summed by donor, producing 17 biological samples: eight PDR and nine controls.

## Server requirements

- Linux x86_64
- at least 32 GiB RAM; 64 GiB is recommended
- 16–32 CPU threads recommended
- at least 180 GiB free when inputs and work files share a filesystem with `KEEP_BAM=0`
- `bash`, `python3`, `tar`, `gzip`, and GNU coreutils
- no root access or Conda installation is required

Use a local SSD or NVMe filesystem for the STAR indices and temporary alignments. High-latency network storage is suitable for fixed inputs but not for index construction.

## Required input layout

Keep these relative paths under one reconstruction directory:

```text
DRnet_GSE276892_reconstruction/
├── remote_linux_reconstruction/
│   ├── README.md
│   ├── preflight.sh
│   ├── run_reconstruction.sh
│   ├── collect_results.sh
│   ├── config/
│   │   ├── sample_sheet.tsv
│   │   ├── input_files.tsv
│   │   └── provenance_sources.tsv
│   └── scripts/
│       ├── verify_remote_inputs.py
│       └── aggregate_featurecounts.py
├── analysis_data/independent_validation/raw_reads/
│   ├── GSE276892/                 # 10 FASTQ files
│   └── GSE147657/                 # 21 FASTQ lanes
└── tools/
    ├── downloads/
    │   ├── STAR-2.7.8a.tar.gz
    │   └── subread-2.0.1-Linux-x86_64.tar.gz
    └── reference/
        ├── gencode_v42_all_regions/
        │   ├── GRCh38.p13.genome.fa.gz
        │   └── gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz
        └── gencode_v42_primary_assembly/
            ├── GRCh38.primary_assembly.genome.fa.gz
            └── gencode.v42.primary_assembly.annotation.gtf.gz
```

The included manifests record the expected file names, sizes, MD5 values, and SHA-256 values. Do not start alignment until the preflight check passes.

## Run the reconstruction

```bash
cd /path/to/DRnet_GSE276892_reconstruction
chmod +x remote_linux_reconstruction/*.sh

env \
  WORK_ROOT=/local/ssd/DRnet_GSE276892_work \
  MIN_FREE_GB=100 \
  bash remote_linux_reconstruction/preflight.sh
```

A successful check ends with `Preflight PASS`. Start the long run only after that message:

```bash
nohup env \
  WORK_ROOT=/local/ssd/DRnet_GSE276892_work \
  MIN_FREE_GB=100 \
  THREADS=16 \
  KEEP_BAM=0 \
  CLEAN_GENERATED_REFERENCE_AND_INDEX=1 \
  RUN_PRIMARY_ASSEMBLY_SENSITIVITY=1 \
  bash remote_linux_reconstruction/run_reconstruction.sh \
  > remote_linux_reconstruction/nohup_reconstruction.log 2>&1 &

echo $! > remote_linux_reconstruction/reconstruction.pid
```

Follow progress with:

```bash
tail -n 80 -f remote_linux_reconstruction/nohup_reconstruction.log
```

The scripts resume from completed indices, lane counts, and STAR logs. Do not start a second instance while the first process is running.

## Collect the results

After a successful run, the workflow creates:

```text
remote_work/DRnet_GSE276892_remote_results.tar.gz
remote_work/DRnet_GSE276892_remote_results.tar.gz.sha256
```

The archive contains both 17-sample gene-count matrices, P2RX4 counts, featureCounts summaries, STAR QC logs, sample-level QC, environment records, tool versions, input verification, and output checksums. STAR indices, decompressed references, and BAM files do not need to be retained when `KEEP_BAM=0`.

Place the downloaded archive and checksum at:

```text
analysis_data/independent_validation/remote_results/
```

The local Python pipeline then runs the P2RX4 negative-binomial models, source and covariate sensitivities, two-new-control direction check, QC correlations, leave-one-out analysis, and reference sensitivity.

## Fixed analysis choices

- `SJDB_OVERHANG=100`
- featureCounts `-t exon -g gene_id -s 0`
- `KEEP_BAM=0` after validated counts and summaries are written
- `CLEAN_GENERATED_REFERENCE_AND_INDEX=1` after each validated workflow
- no ComBat adjustment, because disease status and data source have insufficient design overlap

The FASTQ-to-count reconstruction was added after inspection of the processed matrix and is reported as a protocol deviation, not as a retrospectively confirmatory analysis.

## Upstream resources

- [STAR](https://github.com/alexdobin/STAR)
- [Subread and featureCounts](https://subread.sourceforge.net/)
- [GENCODE human release 42](https://www.gencodegenes.org/human/release_42.html)
- [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/)
