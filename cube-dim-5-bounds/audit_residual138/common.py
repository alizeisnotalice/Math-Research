"""Read-only qualified123 general helpers;138 paths only;no H/log on import."""
import importlib.util
from pathlib import Path
from fractions import Fraction as F
import hashlib
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-audit123-sol/cube-dim-5-bounds/audit_g123_sol')
p=OLD/'common.py'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='59d53a9f318490b8782156a60728a3e7079194737eaa809bcf699c392bd2e4f6'
spec=importlib.util.spec_from_file_location('readonly123_general_helpers',p)
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
for name in ('sha','digest','write','nodup','read','lines','check_locks','Stream','node','box','tree','geometry'):
    globals()[name]=getattr(base,name)
HERE=Path(__file__).resolve().parent
DESK=base.DESK
SCRATCH=DESK/'cube-finite-dim-scratch/audit_residual138'
PROD=DESK/'cube-finite-dim-scratch/upper_parallel137/residual_preflight'
AUTHOR=DESK/'Math-Research-cube5-upper/cube-dim-5-bounds/upper_parallel137'
FAN=base.FAN
FAN_SHA=base.FAN_SHA
KERNEL_SHA=base.KERNEL_SHA
MODEL=base.MODEL
TARGET=base.TARGET

def assert_pids_absent(pids):
    import subprocess
    if not pids:return
    proc=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid=,ppid=,command='],capture_output=True,text=True)
    assert not proc.stdout.strip(),('production/audit PID present;inspect before PF',proc.stdout)
