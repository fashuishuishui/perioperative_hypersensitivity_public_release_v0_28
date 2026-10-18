from pathlib import Path
import json,os
ROOT=Path(__file__).resolve().parents[1]
CONFIG=Path(os.environ.get('PV_INPUT_CONFIG',ROOT/'config/paths.json'))
INPUTS=json.loads(CONFIG.read_text(encoding='utf-8'))
def path(key):
 p=Path(INPUTS[key]);return p if p.is_absolute() else (CONFIG.parent/p).resolve()
RAW=path('raw_sqlite');LATEST=path('latest_sqlite');LICENSE=path('meddra_pt_smq');REFERENCE=path('meddra_llt')
for name in ['work','aggregate_data','metadata']:(ROOT/name).mkdir(exist_ok=True)
