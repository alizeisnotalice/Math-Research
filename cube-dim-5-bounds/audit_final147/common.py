"""Qualified138 helpers with frozen data configuration;no H/log on import."""
import importlib.util,hashlib,json
from pathlib import Path
from fractions import Fraction as F
p=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit123-sol/cube-dim-5-bounds/audit_residual138/common.py')
assert hashlib.sha256(p.read_bytes()).hexdigest()=='917b7f647859f5d91ff3d97bfa72160aab02099746a3103579e508c7ba8a264a'
spec=importlib.util.spec_from_file_location('readonly138_general_helpers',p)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
for name in ('sha','digest','write','nodup','read','lines','check_locks','Stream','node','box','tree','geometry','assert_pids_absent'):
    globals()[name]=getattr(base,name)
HERE=Path(__file__).resolve().parent
CONFIG=json.loads((HERE/'config.json').read_text())
assert CONFIG['audit_id']==HERE.name and CONFIG['workers']==15 and CONFIG['epsilon']=='1/128' and CONFIG['target']=='11/100000'
DESK=base.DESK
SCRATCH=DESK/'cube-finite-dim-scratch'/HERE.name
PROD=Path(CONFIG['production_dir'])
AUTHOR=Path(CONFIG['production_code_dir'])
SOURCE_RESULT=Path(CONFIG['source_result_manifest'])
FAN=base.FAN
FAN_SHA=base.FAN_SHA
KERNEL_SHA=base.KERNEL_SHA
MODEL=base.MODEL
TARGET=base.TARGET
