from pathlib import Path
from fractions import Fraction as F
import json,hashlib,math,os,shutil,subprocess,sys,tempfile,argparse,time
import numpy as np
H=Path(__file__).resolve().parent
SKIP={'seconds','elapsed_seconds','profile_npz_sha256'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cmp(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP,(set(a),set(b))
  for k in set(a)-SKIP:cmp(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):cmp(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-10),(a,b)
 else:assert a==b,(a,b)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists();st=time.time();prov=json.loads((H/'provenance.json').read_text());runs=[]
 for n,hh in prov['payload_SHA256'].items():assert sha(H/n)==hh,n
 for n,hh in prov['reference_SHA256'].items():assert sha(H.parent/n)==hh,n
 with tempfile.TemporaryDirectory(prefix='cube48-')as td:
  root=Path(td);env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
  for d in ['round048_banks','round048_process','round048_roughness']:shutil.copytree(H/d,root/d)
  old=root/'round046_gf';old.mkdir();shutil.copyfile(H.parent/'046/round046_gf/MANIFEST.json',old/'MANIFEST.json')
  b=root/'round048_banks';shutil.copyfile(H.parent/'046/round046_gf/gf_check.c',b/'gf_check.c')
  assert sha(b/'gf_check.c')==json.loads((b/'FROZEN_PLAN.json').read_text())['source_C_sha256']
  subprocess.run(['/usr/bin/clang','-O3','-std=c11',str(b/'gf_check.c'),'-o',str(b/'gf_check')],check=True)
  def run(d,s,out,compare=True,extra=()):
   p=subprocess.run([sys.executable,str(root/d/s),*map(str,extra)],cwd=root/d,env=env,capture_output=True,text=True,timeout=180)
   if p.returncode:raise RuntimeError(s+'\n'+p.stderr[-6000:])
   data=json.loads((root/d/out).read_text())
   if compare:cmp(json.loads((H/d/out).read_text()),data)
   runs.append(d+'/'+s);return data
  m=run('round048_banks','build_banks.py','MANIFEST.json');r=run('round048_banks','ranks.py','RANKS.json');co=run('round048_banks','covariance.py','COVARIANCE.json')
  npverify=root/'round048_roughness/old_verify_gf.py';shutil.copyfile(H.parent/'046/round046_independent/verify_gf.py',npverify)
  independently=[]
  for bank in m['banks']:
   source=root/f"bank{bank['bank']}.json";source.write_text(json.dumps(dict(m,basis=bank['basis'])))
   dest=root/'round048_roughness'/f"bank{bank['bank']}_verified.json"
   data=run('round048_roughness','old_verify_gf.py',dest.name,False,('--manifest',source,'--output',dest));assert data['status']=='PASS';independently.append(bank['bank'])
  outward=[]
  for row in co['scales']:
   if row['integer_scale']:assert row['correlation']is None;continue
   lo,hi=map(F,row['Pearson_correlation_rational_bounds']);dl,dh=map(F,row['Pearson_correlation_outward_decimal_9_places']);assert dl<=lo<=hi<=dh;outward.append(row['L'])
  g=run('round048_process','run_process.py','RESULTS.json',False)
  sm=run('round048_process','check_and_package.py','SUMMARY.json')
  assert hashlib.sha256(json.dumps(g['rows'],sort_keys=True,separators=(',',':')).encode()).hexdigest()==sm['canonical_numeric_rows_sha256']
  with np.load(H/'round048_process/PROFILE_MIN.npz')as x,np.load(root/'round048_process/PROFILE_MIN.npz')as y:
   assert x.files==y.files
   for k in x.files:assert np.array_equal(x[k],y[k]),k
  rough=run('round048_roughness','probe.py','RESULTS.json');density=run('round048_roughness','density_check.py','DENSITY_CHECK.json')
  out=dict(status='PASS',runs=runs,independent_NumPy_full_GF_banks=independently,bank_count=len(m['banks']),rational_outward_intervals_checked=outward,exact_small_covariance=co['small_exact_test']['passed'],cross_rank_histogram=r['cross_endpoint_pair_rank_histogram'],process_configs=len(g['rows']),all_process_profile_arrays_exact=True,process_density_checks=sm['independent_density_checks'],root_density_cases=density['full_covariance_cases'],roughness_cases=len(rough['rows']),new_independent_queries=0,actual_HL_process_transfer_proved=False,new_Cn_log_gain_proved=False,seconds=time.time()-st)
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
