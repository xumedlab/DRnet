#!/usr/bin/env python3
"""Aggregate lane-level featureCounts outputs to biological samples and QC."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_gtf_gene_map(path: Path) -> dict[str, tuple[str, str]]:
    gene_map: dict[str, tuple[str, str]] = {}
    attribute_pattern = re.compile(r'(\S+)\s+"([^"]+)"')
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(attribute_pattern.findall(fields[8]))
            gene_id = attrs.get("gene_id")
            if gene_id:
                gene_map[gene_id] = (
                    attrs.get("gene_name", ""),
                    attrs.get("gene_type", attrs.get("gene_biotype", "")),
                )
    return gene_map


def parse_featurecounts(path: Path) -> tuple[list[dict[str, str]], dict[str, int]]:
    with path.open("r", encoding="utf-8") as handle:
        rows = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(rows, delimiter="\t")
    if reader.fieldnames is None or len(reader.fieldnames) < 7:
        raise ValueError(f"Malformed featureCounts file: {path}")
    count_column = reader.fieldnames[-1]
    annotation_fields = reader.fieldnames[:6]
    annotation: list[dict[str, str]] = []
    counts: dict[str, int] = {}
    for row in reader:
        gene_id = row["Geneid"]
        annotation.append({field: row[field] for field in annotation_fields})
        counts[gene_id] = int(row[count_column])
    return annotation, counts


def parse_star_log(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = line.split("|", 1)
            output[key.strip()] = value.strip().rstrip("%")
    return output


def parse_featurecounts_summary(path: Path) -> dict[str, int]:
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"Empty featureCounts summary: {path}")
    count_column = [name for name in rows[0] if name != "Status"]
    if len(count_column) != 1:
        raise ValueError(f"Unexpected featureCounts summary columns: {path}")
    return {row["Status"]: int(row[count_column[0]]) for row in rows}


def aggregate(
    sample_sheet_path: Path,
    workflow_root: Path,
    gtf_path: Path,
    output_dir: Path,
) -> None:
    sample_rows = read_tsv(sample_sheet_path)
    if len(sample_rows) != 31:
        raise ValueError(f"Expected 31 run rows, found {len(sample_rows)}")
    gene_map = parse_gtf_gene_map(gtf_path)
    sample_order: list[str] = []
    sample_metadata: dict[str, dict[str, str]] = {}
    sample_counts: dict[str, defaultdict[str, int]] = {}
    sample_input_reads: defaultdict[str, int] = defaultdict(int)
    sample_unique_reads: defaultdict[str, int] = defaultdict(int)
    sample_assigned_reads: defaultdict[str, int] = defaultdict(int)
    lane_qc: list[dict[str, object]] = []
    reference_annotation: list[dict[str, str]] | None = None
    reference_gene_ids: list[str] | None = None

    for row in sample_rows:
        run = row["run_accession"]
        sample = row["sample_id"]
        if sample not in sample_order:
            sample_order.append(sample)
            sample_metadata[sample] = row
            sample_counts[sample] = defaultdict(int)
        fc_path = workflow_root / "lane_counts" / f"{run}.featureCounts.txt"
        summary_path = Path(f"{fc_path}.summary")
        star_log_path = workflow_root / "align" / run / "Log.final.out"
        annotation, counts = parse_featurecounts(fc_path)
        gene_ids = [item["Geneid"] for item in annotation]
        if reference_annotation is None:
            reference_annotation = annotation
            reference_gene_ids = gene_ids
        elif gene_ids != reference_gene_ids:
            raise ValueError(f"Gene order differs in {fc_path}")
        for gene_id, value in counts.items():
            sample_counts[sample][gene_id] += value

        star = parse_star_log(star_log_path)
        fc_summary = parse_featurecounts_summary(summary_path)
        input_reads = int(star["Number of input reads"])
        unique_reads = int(star["Uniquely mapped reads number"])
        assigned_reads = fc_summary.get("Assigned", 0)
        sample_input_reads[sample] += input_reads
        sample_unique_reads[sample] += unique_reads
        sample_assigned_reads[sample] += assigned_reads
        lane_qc.append(
            {
                "run_accession": run,
                "sample_id": sample,
                "input_reads": input_reads,
                "uniquely_mapped_reads": unique_reads,
                "uniquely_mapped_percent": star.get("Uniquely mapped reads %", ""),
                "multimapped_percent": star.get("% of reads mapped to multiple loci", ""),
                "too_short_percent": star.get("% of reads unmapped: too short", ""),
                "featurecounts_assigned_reads": assigned_reads,
                "featurecounts_assigned_percent": (
                    100.0 * assigned_reads / input_reads if input_reads else 0.0
                ),
            }
        )

    if reference_annotation is None or reference_gene_ids is None:
        raise ValueError("No featureCounts rows were read")
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "sample_gene_counts.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = [
            "Geneid",
            "gene_name",
            "gene_type",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
            *sample_order,
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for annotation in reference_annotation:
            gene_id = annotation["Geneid"]
            gene_name, gene_type = gene_map.get(gene_id, ("", ""))
            writer.writerow(
                {
                    **annotation,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    **{sample: sample_counts[sample][gene_id] for sample in sample_order},
                }
            )

    lane_fields = list(lane_qc[0])
    with (output_dir / "lane_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=lane_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(lane_qc)

    sample_qc: list[dict[str, object]] = []
    for sample in sample_order:
        counts = sample_counts[sample]
        total_counts = sum(counts.values())
        input_reads = sample_input_reads[sample]
        metadata = sample_metadata[sample]
        sample_qc.append(
            {
                "sample_id": sample,
                "diagnosis": metadata["diagnosis"],
                "source_dataset": metadata["source_dataset"],
                "is_reused_control": metadata["is_reused_control"],
                "number_of_lanes": metadata["lanes_per_sample"],
                "input_reads": input_reads,
                "uniquely_mapped_reads": sample_unique_reads[sample],
                "uniquely_mapped_percent": (
                    100.0 * sample_unique_reads[sample] / input_reads if input_reads else 0.0
                ),
                "featurecounts_assigned_reads": sample_assigned_reads[sample],
                "featurecounts_assigned_percent": (
                    100.0 * sample_assigned_reads[sample] / input_reads if input_reads else 0.0
                ),
                "total_gene_counts": total_counts,
                "detected_genes": sum(value > 0 for value in counts.values()),
            }
        )
    with (output_dir / "sample_qc.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample_qc[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_qc)

    p2rx4_rows = []
    for annotation in reference_annotation:
        gene_id = annotation["Geneid"]
        gene_name, gene_type = gene_map.get(gene_id, ("", ""))
        if gene_name == "P2RX4":
            p2rx4_rows.append(
                {
                    "Geneid": gene_id,
                    "gene_name": gene_name,
                    "gene_type": gene_type,
                    **{sample: sample_counts[sample][gene_id] for sample in sample_order},
                }
            )
    if len(p2rx4_rows) != 1:
        raise ValueError(f"Expected one P2RX4 gene row, found {len(p2rx4_rows)}")
    with (output_dir / "P2RX4_counts.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["Geneid", "gene_name", "gene_type", *sample_order]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(p2rx4_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-sheet", type=Path, required=True)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    aggregate(
        args.sample_sheet.resolve(),
        args.workflow_root.resolve(),
        args.gtf.resolve(),
        args.output_dir.resolve(),
    )


if __name__ == "__main__":
    main()
