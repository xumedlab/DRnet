import csv
from pathlib import Path
import shutil
from pipeline_utils import ensure_dirs, read_gmt, write_gmt, log_message
import runpy
cfg=runpy.run_path('00_config.py')
RAW_DIR,PROC_DIR,RESULT_DIR=cfg['RAW_DIR'],cfg['PROC_DIR'],cfg['RESULT_DIR']

def clean_name(n):
    n=n.strip().replace('memeory','memory').replace('  ',' ')
    return '_'.join(n.split())

def main():
    ensure_dirs(RAW_DIR, PROC_DIR, RESULT_DIR/'logs', RESULT_DIR/'tables')
    # bootstrap raw files if they are still in repo root
    for fn in ['hallmark_human_gene_symbols.gmt','hallmark_inflammatory_response.gmt','immune_28_signatures_tisidb.gmt']:
        src=Path(fn); dst=RAW_DIR/fn
        if src.exists() and not dst.exists(): shutil.copyfile(src, dst)

    hall=read_gmt(RAW_DIR/'hallmark_human_gene_symbols.gmt')
    infl=read_gmt(RAW_DIR/'hallmark_inflammatory_response.gmt')
    immune=read_gmt(RAW_DIR/'immune_28_signatures_tisidb.gmt')
    key='HALLMARK_INFLAMMATORY_RESPONSE'
    include=key in hall
    same=set(hall.get(key,[]))==set(infl.get(key,[]))

    cleaned={clean_name(k):sorted(set(v)) for k,v in immune.items()}
    out_gmt=PROC_DIR/'immune_28_signatures_charoentong2017_human_cleaned.gmt'
    write_gmt(out_gmt, cleaned)

    summary=[
        {'source':'hallmark_human_gene_symbols.gmt','set_count':len(hall)},
        {'source':'hallmark_inflammatory_response.gmt','set_count':len(infl)},
        {'source':'immune_28_signatures_tisidb.gmt','set_count':len(immune)},
        {'source':'immune_cleaned','set_count':len(cleaned)},
    ]
    with open(RESULT_DIR/'tables'/'gene_set_summary.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['source','set_count']);w.writeheader();w.writerows(summary)
    md=RESULT_DIR/'logs'/'gene_set_source_and_cleaning_log.md'
    md.write_text(f"# Gene set validation\n\n- Hallmark includes {key}: {include}\n- Standalone inflammatory set equals hallmark member: {same}\n- Immune set count: {len(immune)}\n- Cleaned immune set count: {len(cleaned)}\n",encoding='utf-8')
    log_message('00a_validate_and_clean_gene_sets', f'cleaned immune sets: {len(cleaned)}')

if __name__=='__main__': main()
