"""Hash-checked isolated reproduction of round051; no writes to saved payloads."""
from pathlib import Path
import json,hashlib,tempfile,shutil,subprocess,os,sys,math,time,argparse
HERE=Path(__file__).resolve().parent

def compare(a,b,path=''):
 if isinstance(a,dict):
  assert set(a)==set(b),(path,'keys')
  for k in a:
   if k in ('seconds','elapsed_seconds','inputs_dir'):continue
   compare(a[k],b[k],path+'/'+k)
 elif isinstance(a,list):
  assert len(a)==len(b),(path,'length')
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
 elif isinstance(a,(int,float)) and not isinstance(a,bool):
  assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-12),(path,a,b)
 else:assert a==b,(path,a,b)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();start=time.perf_counter()
 p=json.loads((HERE/'provenance.json').read_text())
 for rel,h in p['payload_sha256'].items():assert hashlib.sha256((HERE/rel).read_bytes()).hexdigest()==h,rel
 for rel,h in p['reference_sha256'].items():assert hashlib.sha256((HERE.parent/rel).read_bytes()).hexdigest()==h,rel
 with tempfile.TemporaryDirectory(prefix='cube051-') as td:
  rounds=Path(td);dst=rounds/'051';shutil.copytree(HERE,dst)
  for rel in ('046/round046_oracle/source_oracle.py','046/round046_gf/MANIFEST.json'):
   target=rounds/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(HERE.parent/rel,target)
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',CUBE_ORACLE_DIR=str(rounds/'046/round046_oracle'),CUBE_INPUT_DIR=str(dst/'channel/inputs'))
  steps=[('countgate/probe.py',[]),('countgate/cap_probe.py',[]),('countgate/check_exact.py',[]),('softgate/probe.py',[]),('channel/channel.py',[]),('channel/channel.py',['--high-pilot','--manifest',str(rounds/'046/round046_gf/MANIFEST.json')])]
  for rel,opts in steps:
   r=subprocess.run([sys.executable,str(dst/rel),*opts],cwd=dst,env=env,capture_output=True,text=True)
   assert r.returncode==0,(rel,r.stderr[-4000:])
  resultfiles=['countgate/RESULTS.json','countgate/CAP_RESULTS.json','countgate/CHECKS.json','softgate/RESULTS.json','channel/RESULTS.json','channel/HIGH_PILOT.json']
  for rel in resultfiles:compare(json.loads((HERE/rel).read_text()),json.loads((dst/rel).read_text()),rel)
 result={'status':'PASS','scripts_run':len(steps),'result_summaries_compared':len(resultfiles),'payload_hashes_checked':len(p['payload_sha256']),'reference_hashes_checked':len(p['reference_sha256']),'elapsed_seconds':time.perf_counter()-start,'python':sys.version.split()[0]}
 if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
