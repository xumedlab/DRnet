#!/usr/bin/env python3
"""
离线构建 Ensembl->Gene Symbol 映射文件。

输出: data_raw/ensembl_to_symbol.csv
字段: ensembl_id,gene_symbol

优先级:
1) --mapping-file (CSV/TSV，需包含 ensembl/symbol 列)
2) --gtf (Ensembl GTF，解析 gene_id + gene_name)
3) 二者同时提供时合并（mapping-file 优先覆盖）
"""

import argparse
import csv
from pathlib import Path


def detect_delim(path: Path) -> str:
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        head = f.read(2048)
    return '\t' if head.count('\t') >= head.count(',') else ','


def norm_ensembl(x: str) -> str:
    return (x or '').split('.')[0].strip()


def read_mapping_table(path: Path):
    delim = detect_delim(path)
    out = {}
    with open(path, 'r', encoding='utf-8', errors='ignore', newline='') as f:
        r = csv.DictReader(f, delimiter=delim)
        if not r.fieldnames:
            return out
        fields = {c.lower().strip(): c for c in r.fieldnames}
        ekeys = [k for k in fields if 'ensembl' in k and ('id' in k or 'gene' in k)]
        skeys = [k for k in fields if ('symbol' in k) or ('hgnc' in k) or (k == 'gene') or ('gene_name' in k)]
        if not ekeys or not skeys:
            raise ValueError(f'{path} 缺少可识别的 ensembl/symbol 列。')
        ecol, scol = fields[ekeys[0]], fields[skeys[0]]
        for row in r:
            eid = norm_ensembl(row.get(ecol, ''))
            sym = (row.get(scol, '') or '').strip().upper()
            if eid.startswith('ENSG') and sym:
                out[eid] = sym
    return out


def parse_gtf_attributes(attr_text: str):
    attrs = {}
    # key "value";
    for part in attr_text.strip().split(';'):
        part = part.strip()
        if not part or ' ' not in part:
            continue
        k, v = part.split(' ', 1)
        attrs[k.strip()] = v.strip().strip('"')
    return attrs


def read_gtf(path: Path):
    out = {}
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature != 'gene':
                continue
            attrs = parse_gtf_attributes(cols[8])
            eid = norm_ensembl(attrs.get('gene_id', ''))
            gname = (attrs.get('gene_name', '') or '').strip().upper()
            if eid.startswith('ENSG') and gname:
                out[eid] = gname
    return out


def write_mapping(path: Path, mapping: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['ensembl_id', 'gene_symbol'])
        w.writeheader()
        for eid, sym in sorted(mapping.items()):
            w.writerow({'ensembl_id': eid, 'gene_symbol': sym})


def main():
    ap = argparse.ArgumentParser(description='离线构建 ensembl_to_symbol.csv')
    ap.add_argument('--mapping-file', type=Path, help='已有映射表 CSV/TSV')
    ap.add_argument('--gtf', type=Path, help='Ensembl GTF 文件路径')
    ap.add_argument('--out', type=Path, default=Path('data_raw/ensembl_to_symbol.csv'))
    args = ap.parse_args()

    mapping = {}
    if args.gtf:
        if not args.gtf.exists():
            raise FileNotFoundError(f'GTF 不存在: {args.gtf}')
        mapping.update(read_gtf(args.gtf))
    if args.mapping_file:
        if not args.mapping_file.exists():
            raise FileNotFoundError(f'mapping-file 不存在: {args.mapping_file}')
        # mapping-file 优先级更高，后覆盖
        mapping.update(read_mapping_table(args.mapping_file))

    if not mapping:
        raise RuntimeError('未构建到任何映射。请提供 --gtf 或 --mapping-file。')

    write_mapping(args.out, mapping)
    print(f'Wrote {len(mapping)} mappings -> {args.out}')


if __name__ == '__main__':
    main()
