"""Opt-in future replay: restore exact frozen gate inputs in an exclusive directory.
Not executed during final delivery; no source check is patched or weakened.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,sys
p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if not a.verify:p.error('Use --verify for a new mathematical audit execution.')
P=Path(__file__).resolve().parents[1];reg=P/'round020/final_bounds_gate_registration.json';r=json.loads(reg.read_bytes());assert not a.output.exists();a.output.mkdir(parents=True)
for rel,h in r['sources'].items():
 source=P/'round020/final_bounds_ledger_at_gate.json' if rel=='results/bounds_ledger.json' else P/rel
 assert hashlib.sha256(source.read_bytes()).hexdigest()==h,rel
 dest=a.output/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
for rel in ['round020/final_bounds_gate.py','round020/final_bounds_gate_registration.json']:
 dest=a.output/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(P/rel,dest)
assert hashlib.sha256((a.output/'round020/final_bounds_gate.py').read_bytes()).hexdigest()==r['code_sha256']
raise SystemExit(subprocess.run([sys.executable,str(a.output/'round020/final_bounds_gate.py')],cwd=a.output).returncode)
