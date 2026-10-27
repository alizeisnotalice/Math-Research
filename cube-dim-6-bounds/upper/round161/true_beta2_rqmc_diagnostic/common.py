"""C161: standard-library-only I/O and pinning; importing performs no science."""
from pathlib import Path
import hashlib,json,os
P=Path(__file__).resolve().parent
SOURCE_FILES=('common.py','prepare.py','rqmc_core.py','run_prepared.py','PROOF.md','OPERATOR.md','source_pins.json','plan.json')
def need(v,message):
 if not v:raise ArithmeticError(message)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fsync_dir(p):
 fd=os.open(p,os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)
def write_new(path,value):
 path=Path(path)
 with path.open('x') as f:
  json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 fsync_dir(path.parent)
def check_sources(pins):
 for label,row in pins['sources'].items():
  path=Path(row['path']);need(path.stat().st_size==row['bytes']and sha(path)==row['sha256'],'SOURCE_CHANGED:'+label)
def source_freeze():return {name:sha(P/name)for name in SOURCE_FILES}
def check_freeze(frozen):
 need(set(frozen)==set(SOURCE_FILES),'EXACT_SOURCE_FREEZE_LIST')
 for name,value in frozen.items():need(sha(P/name)==value,'FROZEN_SOURCE:'+name)
