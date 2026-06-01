import csv
import runpy
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline_utils import ensure_dirs, log_message

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']


def read_manifest(path):
    with open(path, encoding='utf-8') as f:
        return list(csv.DictReader(f))


def read_expression_table(path):
    df = pd.read_csv(path, sep=None, engine='python')
    first_col = df.columns[0]
    if first_col != 'ensemblID':
        df = df.rename(columns={first_col: 'ensemblID'})
    return df


def write_tsv(df, path):
    df.to_csv(path, sep='\t', index=False)


def main():
    ensure_dirs(PROC_DIR, RESULT_DIR / 'logs')

    mac = read_manifest(PROC_DIR / 'manifest_macula_4groups.csv')
    primary = read_manifest(PROC_DIR / 'manifest_primary_binary.csv')
    mac_ids = [row['sample_id'] for row in mac]
    primary_ids = [row['sample_id'] for row in primary]

    counts_raw = read_expression_table(RAW_DIR / 'GSE160306_human_retina_DR_totalRNA_counts.txt')
    cpm_raw = read_expression_table(RAW_DIR / 'GSE160306_human_retina_DR_totalRNA_normalized_cpm.txt')

    missing_counts = [sample for sample in mac_ids if sample not in counts_raw.columns]
    if missing_counts:
        raise ValueError(f'Raw counts missing required samples: {missing_counts}')

    missing_cpm = [sample for sample in mac_ids if sample not in cpm_raw.columns]

    counts_mac = counts_raw[['ensemblID', *mac_ids]].copy()
    keep_mask = (counts_mac[mac_ids] >= cfg['MIN_COUNT']).sum(axis=1) >= cfg['MIN_SAMPLES']
    counts_mac = counts_mac.loc[keep_mask].copy()
    counts_primary = counts_mac[['ensemblID', *primary_ids]].copy()

    library_sizes = counts_raw[mac_ids].sum(axis=0)
    cpm_values = counts_mac[mac_ids].div(library_sizes, axis=1) * 1_000_000.0
    log2cpm_mac = pd.concat(
        [
            counts_mac[['ensemblID']].copy(),
            pd.DataFrame(
                np.log2(cpm_values.astype(float) + 1.0),
                columns=mac_ids,
                index=counts_mac.index,
            ),
        ],
        axis=1,
    )

    write_tsv(counts_mac, PROC_DIR / 'counts_macula_4groups.tsv')
    write_tsv(counts_primary, PROC_DIR / 'counts_primary_binary.tsv')
    write_tsv(log2cpm_mac, PROC_DIR / 'log2cpm_macula_4groups.tsv')

    for name, rows in [('pheno_macula_4groups.csv', mac), ('pheno_primary_binary.csv', primary)]:
        with open(PROC_DIR / name, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    with open(PROC_DIR / 'gene_filtering_summary.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                'genes_before_filter',
                'genes_after_filter',
                'macula_samples',
                'primary_samples',
                'normalized_cpm_missing_samples',
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                'genes_before_filter': int(len(counts_raw)),
                'genes_after_filter': int(len(counts_mac)),
                'macula_samples': len(mac_ids),
                'primary_samples': len(primary_ids),
                'normalized_cpm_missing_samples': ';'.join(missing_cpm) if missing_cpm else 'none',
            }
        )

    log_lines = [
        '# Expression preparation summary',
        '',
        f'- macula samples retained: {len(mac_ids)}',
        f'- primary comparison samples retained: {len(primary_ids)}',
        f'- genes before filter: {len(counts_raw)}',
        f'- genes after filter: {len(counts_mac)}',
        '- log2CPM source: recomputed from raw counts to keep count/log2CPM sample universe consistent',
        f"- normalized CPM file missing macula samples: {', '.join(missing_cpm) if missing_cpm else 'none'}",
    ]
    (RESULT_DIR / 'logs' / '02_prepare_expression_details.md').write_text(
        '\n'.join(log_lines) + '\n',
        encoding='utf-8',
    )
    log_message(
        '02_prepare_expression',
        f'genes_after_filter={len(counts_mac)} macula_samples={len(mac_ids)} missing_normalized_cpm={len(missing_cpm)}',
    )


if __name__ == '__main__':
    main()
