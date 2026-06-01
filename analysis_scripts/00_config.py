from pathlib import Path

RAW_DIR = Path('data_raw')
PROC_DIR = Path('data_processed')
RESULT_DIR = Path('results')
TABLE_DIR = RESULT_DIR / 'tables'
LOG_DIR = RESULT_DIR / 'logs'
FIG_DIR = RESULT_DIR / 'figures'

GROUPS = ["healthy control", "diabetic", "NPDR", "NPDR/PDR + DME"]
PRIMARY_CTRL = "healthy control"
PRIMARY_CASE = "NPDR/PDR + DME"
SEVERITY_MAP = {"healthy control":0, "diabetic":1, "NPDR":2, "NPDR/PDR + DME":3}

PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 0.5
MIN_COUNT = 10
MIN_SAMPLES = 5
RANDOM_SEED = 202501

FIG_DPI = 300
FIG_FORMATS = ("png", "pdf")
