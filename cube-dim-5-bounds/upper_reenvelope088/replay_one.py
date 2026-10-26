from pathlib import Path
import argparse,gzip,json,hashlib
import engine as E
H=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--rank',type=int,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists();a.output.mkdir(parents=True);E.S=a.output.resolve();reg=json.loads((H/'registration.json').read_text());fp=hashlib.sha256((E.sha(H/'engine.py')+E.sha(H/'registration.json')).encode()).hexdigest();row=json.loads(gzip.decompress((H/'nodes.json.gz').read_bytes()))[a.rank];E.init(fp,reg['time_endpoints']);r=E.worker(row);assert r['status']=='complete';print(json.dumps(r,indent=2))
