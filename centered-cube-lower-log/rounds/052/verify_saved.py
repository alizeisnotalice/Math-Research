"""Round052 isolated reproduction; keep outputs away from saved payloads."""
from pathlib import Path
import tempfile,shutil,subprocess,os,sys,json,hashlib,time,math,argparse
HERE=Path(__file__).resolve().parent

def compare(a,b,path=''):
 if isinstance(a,dict):
  assert set(a)==set(b),(path,'keys')
  for k in a:
   if k in ('seconds','elapsed_seconds'):continue
   compare(a[k],b[k],path+'/'+k)
 elif isinstance(a,list):
  assert len(a)==len(b),(path,'length')
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
 elif isinstance(a,(int,float)) and not isinstance(a,bool):assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-12),(path,a,b)
 else:assert a==b,(path,a,b)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();t0=time.perf_counter()
 p=json.loads((HERE/'provenance.json').read_text())
 for rel,h in p['payload_sha256'].items():assert hashlib.sha256((HERE/rel).read_bytes()).hexdigest()==h,rel
 for rel,h in p['reference_sha256'].items():assert hashlib.sha256((HERE.parent/rel).read_bytes()).hexdigest()==h,rel
 scripts=['flat_budget/probe.py','flat_budget/check_order.py','amplitude/probe.py','tailgate/probe.py','tailgate/verify_exact.py']
 results=['flat_budget/RESULTS.json','flat_budget/ORDER_CHECKS.json','amplitude/RESULTS.json','tailgate/RESULTS.json','tailgate/EXACT_CHECKS.json']
 with tempfile.TemporaryDirectory(prefix='cube052-') as td:
  dst=Path(td)/'052';shutil.copytree(HERE,dst);env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
  for rel in scripts:
   r=subprocess.run([sys.executable,str(dst/rel)],cwd=dst,env=env,capture_output=True,text=True)
   assert r.returncode==0,(rel,r.stderr[-4000:])
  for rel in results:compare(json.loads((HERE/rel).read_text()),json.loads((dst/rel).read_text()),rel)
 out={'status':'PASS','scripts_run':len(scripts),'results_compared':len(results),'payload_hashes_checked':len(p['payload_sha256']),'reference_hashes_checked':len(p['reference_sha256']),'elapsed_seconds':time.perf_counter()-t0,'python':sys.version.split()[0]}
 if args.output:args.output.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
