import csv, math, json
from pathlib import Path
from urllib import request, parse

GROUPS = ["healthy control", "diabetic", "NPDR", "NPDR/PDR + DME"]
SEVERITY_MAP = {"healthy control": 0, "diabetic": 1, "NPDR": 2, "NPDR/PDR + DME": 3}


def ensure_dirs(*dirs):
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)


def log_message(step, msg, logs_dir="results/logs"):
    ensure_dirs(logs_dir)
    with open(Path(logs_dir) / f"{step}.log", "a", encoding="utf-8") as f:
        f.write(msg.rstrip() + "\n")


def read_gmt(path):
    sets = {}
    with open(path, encoding="utf-8", errors="ignore") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"Invalid GMT line {ln} in {path}")
            sets[parts[0]] = [g.strip() for g in parts[2:] if g.strip()]
    return sets


def write_gmt(path, sets):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        for k, v in sets.items():
            w.writerow([k, "na", *v])


def parse_series_matrix(path):
    rows = {}
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            rec = next(csv.reader([line], delimiter="\t"))
            key = rec[0].replace("!Sample_", "")
            vals = [x.strip().strip('"') for x in rec[1:]]
            rows.setdefault(key, []).append(vals)
    titles = rows["title"][0]
    data = [{"sample_id": s} for s in titles]
    for k, vv in rows.items():
        if k == "title":
            continue
        if k == "characteristics_ch1":
            for vals in vv:
                for i, v in enumerate(vals):
                    if ": " in v:
                        kk, vv2 = v.split(": ", 1)
                        data[i][kk.strip()] = vv2.strip()
            continue
        vals = vv[0]
        for i, v in enumerate(vals):
            data[i][k] = v
    return data


def read_table(path):
    with open(path, encoding="utf-8", newline="") as f:
        r = csv.reader(f)
        header = next(r)
        rows = [row for row in r if row]
    return header, rows


def write_csv(path, rows, fieldnames):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def bh_adjust(pvals):
    n = len(pvals)
    if n == 0:
        return []
    idx = sorted(range(n), key=lambda i: pvals[i])
    out = [1.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(idx), 1):
        k = n - rank + 1
        val = min(prev, pvals[i] * n / k)
        out[i] = val
        prev = val
    return out


def mann_whitney_u(x, y):
    vals = [(v, 0) for v in x] + [(v, 1) for v in y]
    vals.sort(key=lambda z: z[0])
    ranks = [0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j < len(vals) and vals[j][0] == vals[i][0]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            ranks[k] = r
        i = j
    rx = sum(r for r, (_, g) in zip(ranks, vals) if g == 0)
    nx, ny = len(x), len(y)
    u = rx - nx * (nx + 1) / 2
    mu = nx * ny / 2
    sigma = math.sqrt(nx * ny * (nx + ny + 1) / 12) if nx * ny > 0 else 1
    z = (u - mu) / (sigma + 1e-9)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u, p


def spearman(x, y):
    def rank(v):
        s = sorted((val, i) for i, val in enumerate(v))
        r = [0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j < len(v) and s[j][0] == s[i][0]:
                j += 1
            rr = (i + j + 1) / 2
            for k in range(i, j):
                r[s[k][1]] = rr
            i = j
        return r

    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    rho = num / (dx * dy + 1e-12)
    n = len(x)
    t = abs(rho) * math.sqrt((n - 2) / (1 - rho * rho + 1e-12)) if n > 2 else 0
    p = 2 * (1 - 0.5 * (1 + math.erf(t / math.sqrt(2))))
    return rho, p


def ssgsea_score(expr_dict, geneset):
    genes = [g for g in geneset if g in expr_dict]
    if not genes:
        return 0.0
    return sum(expr_dict[g] for g in genes) / len(genes)


def normalize_ensembl_id(gid):
    return str(gid).split('.')[0].strip()


def _read_mapping_file(path):
    rows = []
    with open(path, encoding='utf-8', errors='ignore', newline='') as f:
        sample = f.read(2048)
        f.seek(0)
        delim = '\t' if sample.count('\t') >= sample.count(',') else ','
        r = csv.DictReader(f, delimiter=delim)
        for row in r:
            rows.append({(k or '').lower().strip(): (v or '').strip() for k, v in row.items()})
    if not rows:
        return {}
    ensembl_keys = [k for k in rows[0].keys() if 'ensembl' in k and ('id' in k or 'gene' in k)]
    symbol_keys = [k for k in rows[0].keys() if ('symbol' in k) or (k == 'gene') or ('hgnc' in k)]
    if not ensembl_keys or not symbol_keys:
        return {}
    ek, sk = ensembl_keys[0], symbol_keys[0]
    out = {}
    for row in rows:
        e, s = normalize_ensembl_id(row.get(ek, '')), row.get(sk, '').upper()
        if e.startswith('ENSG') and s and s != 'NAN':
            out[e] = s
    return out


def _query_mygene_batch(ensembl_ids):
    out = {}
    if not ensembl_ids:
        return out
    # MyGene batch endpoint
    chunk = 400
    for i in range(0, len(ensembl_ids), chunk):
        ids = ensembl_ids[i:i + chunk]
        payload = parse.urlencode({
            'ids': ','.join(ids),
            'fields': 'symbol',
            'species': 'human'
        }).encode('utf-8')
        req = request.Request('https://mygene.info/v3/gene', data=payload, method='POST')
        req.add_header('Content-Type', 'application/x-www-form-urlencoded')
        try:
            with request.urlopen(req, timeout=40) as resp:
                arr = json.loads(resp.read().decode('utf-8'))
            if isinstance(arr, list):
                for rec in arr:
                    q = normalize_ensembl_id(rec.get('query', ''))
                    s = str(rec.get('symbol', '')).upper()
                    if q.startswith('ENSG') and s and s != 'NAN':
                        out[q] = s
        except Exception:
            return out
    return out


def _query_ensembl_xrefs(ensembl_ids):
    # slower fallback but robust if mygene is blocked
    out = {}
    for eid in ensembl_ids[:2000]:  # avoid unbounded runtime
        url = f'https://rest.ensembl.org/xrefs/id/{eid}?external_db=HGNC;content-type=application/json'
        try:
            req = request.Request(url, headers={'Accept': 'application/json'})
            with request.urlopen(req, timeout=20) as resp:
                arr = json.loads(resp.read().decode('utf-8'))
            if isinstance(arr, list) and arr:
                label = str(arr[0].get('display_id', '')).upper()
                if label:
                    out[eid] = label
        except Exception:
            continue
    return out




def _query_biomart(ensembl_ids):
    out = {}
    if not ensembl_ids:
        return out
    chunk = 400
    for i in range(0, len(ensembl_ids), chunk):
        ids = ensembl_ids[i:i+chunk]
        values = ','.join(ids)
        xml = f"""<?xml version='1.0' encoding='UTF-8'?>
<!DOCTYPE Query>
<Query virtualSchemaName='default' formatter='TSV' header='0' uniqueRows='1' count='' datasetConfigVersion='0.6'>
  <Dataset name='hsapiens_gene_ensembl' interface='default'>
    <Filter name='ensembl_gene_id' value='{values}'/>
    <Attribute name='ensembl_gene_id'/>
    <Attribute name='hgnc_symbol'/>
  </Dataset>
</Query>"""
        try:
            q = parse.urlencode({'query': xml})
            url = f'https://www.ensembl.org/biomart/martservice?{q}'
            with request.urlopen(url, timeout=60) as resp:
                txt = resp.read().decode('utf-8', errors='ignore')
            for line in txt.splitlines():
                parts = line.strip().split('	')
                if len(parts) >= 2:
                    e = normalize_ensembl_id(parts[0])
                    g = parts[1].strip().upper()
                    if e.startswith('ENSG') and g:
                        out[e] = g
        except Exception:
            return out
    return out


def load_ensembl_symbol_mapping(raw_dir='data_raw', proc_dir='data_processed', ensembl_ids=None):
    raw_dir = Path(raw_dir)
    proc_dir = Path(proc_dir)
    ensure_dirs(proc_dir)
    cache = proc_dir / 'ensembl_to_symbol_mapping.csv'

    mapping = {}
    if cache.exists():
        mapping = _read_mapping_file(cache)

    if not mapping:
        for p in [
            raw_dir / 'ensembl_to_symbol.csv',
            raw_dir / 'ensembl_to_symbol.tsv',
            raw_dir / 'gene_annotation.csv',
            raw_dir / 'gene_annotation.tsv',
        ]:
            if p.exists():
                mapping = _read_mapping_file(p)
                if mapping:
                    break

    if ensembl_ids:
        need = sorted({normalize_ensembl_id(x) for x in ensembl_ids if normalize_ensembl_id(x).startswith('ENSG')})
        missing = [x for x in need if x not in mapping]
        if missing:
            online = _query_mygene_batch(missing)
            mapping.update(online)
            missing2 = [x for x in missing if x not in mapping]
            if missing2:
                mapping.update(_query_biomart(missing2))
                missing3 = [x for x in missing2 if x not in mapping]
                if missing3:
                    mapping.update(_query_ensembl_xrefs(missing3))

    if mapping:
        with open(cache, 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=['ensembl_id', 'gene_symbol'])
            w.writeheader()
            for e, s in sorted(mapping.items()):
                w.writerow({'ensembl_id': e, 'gene_symbol': s})
    return mapping


def matrix_ensembl_to_symbol(gene_values, mapping):
    symbol_values = {}
    dedup_log = []
    for gid, vals in gene_values.items():
        eid = normalize_ensembl_id(gid)
        sym = mapping.get(eid)
        if not sym:
            continue
        if sym in symbol_values:
            old_mean = sum(symbol_values[sym]) / len(symbol_values[sym])
            new_mean = sum(vals) / len(vals)
            if new_mean > old_mean:
                symbol_values[sym] = vals
                dedup_log.append({'gene_symbol': sym, 'kept_ensembl': eid, 'reason': 'higher_mean_expression'})
        else:
            symbol_values[sym] = vals
    return symbol_values, dedup_log
