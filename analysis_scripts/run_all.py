import argparse
import datetime as dt
import json
import subprocess
import sys
import time

import runpy

from pipeline_utils import ensure_dirs

cfg = runpy.run_path('00_config.py')
RESULT_DIR = cfg['RESULT_DIR']
LOG_DIR = RESULT_DIR / 'logs'
STEPS = [
    '00a_validate_and_clean_gene_sets.py',
    '01_parse_manifest.py',
    '02_prepare_expression.py',
    '03_qc_and_pca.py',
    '04_differential_expression.py',
    '05_inflammation_scoring.py',
    '06_candidate_selection.py',
    '07_lasso_signature.py',
    '08_enrichment_analysis.py',
    '09_immune_infiltration.py',
    '10_make_figures.py',
    '11_export_master_summary.py',
]


def ts_now():
    return dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def emit(msg, log_file):
    line = f"[{ts_now()}] {msg}"
    print(line, flush=True)
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def load_state(statef):
    if not statef.exists():
        return {'completed': []}

    try:
        state = json.loads(statef.read_text(encoding='utf-8'))
    except json.JSONDecodeError:
        return {'completed': []}

    completed = []
    seen = set()
    for step in state.get('completed', []):
        if step in STEPS and step not in seen:
            completed.append(step)
            seen.add(step)
    return {'completed': completed}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    args = ap.parse_args()

    ensure_dirs(LOG_DIR)
    statef = LOG_DIR / 'run_state.json'
    run_log = LOG_DIR / 'run_all.log'

    if args.resume:
        state = load_state(statef)
    else:
        state = {'completed': []}
        statef.write_text(json.dumps(state, indent=2), encoding='utf-8')

    emit(f'Pipeline start. total_steps={len(STEPS)} resume={args.resume}', run_log)

    total_start = time.perf_counter()
    for i, step in enumerate(STEPS, 1):
        if args.resume and step in state['completed']:
            emit(f'SKIP step={i}/{len(STEPS)} script={step} reason=already_completed', run_log)
            continue

        emit(f'START step={i}/{len(STEPS)} script={step}', run_log)
        step_start = time.perf_counter()
        try:
            subprocess.run([sys.executable, step], check=True)
        except subprocess.CalledProcessError as e:
            elapsed = time.perf_counter() - step_start
            emit(
                f'FAIL step={i}/{len(STEPS)} script={step} '
                f'elapsed_sec={elapsed:.2f} returncode={e.returncode}',
                run_log,
            )
            raise

        elapsed = time.perf_counter() - step_start
        if step not in state['completed']:
            state['completed'].append(step)
        statef.write_text(json.dumps(state, indent=2), encoding='utf-8')
        emit(f'DONE step={i}/{len(STEPS)} script={step} elapsed_sec={elapsed:.2f}', run_log)

    total_elapsed = time.perf_counter() - total_start
    emit(f'Pipeline finished. elapsed_sec={total_elapsed:.2f}', run_log)


if __name__ == '__main__':
    main()
