from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,os,shutil,subprocess,sys,tempfile,time
H=Path(__file__).resolve().parent
SKIP={'runtime_python','seconds','elapsed_seconds','numpy'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cmp(a,b):
 if isinstance(a,dict):
  assert set(a)-SKIP==set(b)-SKIP
  for k in set(a)-SKIP:cmp(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):cmp(x,y)
 elif isinstance(a,float):assert math.isclose(a,b,rel_tol=3e-8,abs_tol=3e-9),(a,b)
 else:assert a==b,(a,b)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();assert not args.output.exists();st=time.time();prov=json.loads((H/'provenance.json').read_text());runs=[]
 for n,hh in prov['payload_SHA256'].items():assert sha(H/n)==hh,n
 for n,hh in prov['reference_SHA256'].items():assert sha(H.parent/n)==hh,n
 with tempfile.TemporaryDirectory(prefix='cube49-')as td:
  root=Path(td)
  for d in ['round049_affinity','round049_minorant','round049_multibank_cf']:shutil.copytree(H/d,root/d)
  old=root/'round048_banks';old.mkdir();shutil.copyfile(H.parent/'048/round048_banks/COVARIANCE.json',old/'COVARIANCE.json')
  def run(d,s,out):
   p=subprocess.run([sys.executable,str(root/d/s)],cwd=root/d,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=180)
   if p.returncode:raise RuntimeError(s+'\n'+p.stderr[-6000:])
   a=json.loads((root/d/out).read_text());cmp(json.loads((H/d/out).read_text()),a);runs.append(d+'/'+s);return a
  affinity=run('round049_affinity','probe.py','RESULTS.json');finite=run('round049_affinity','finite_check.py','FINITE.json');minor=run('round049_minorant','probe.py','RESULTS.json');multi=run('round049_multibank_cf','probe.py','RESULTS.json')
  grid=json.loads((H/'round049_affinity/GRID.json').read_text());total=F(0)
  for r in affinity['rows']:
   E=F(r['E_exact']);m=E.numerator//E.denominator;assert m==r['floor_E'] and E>0
   prob=F(12,2**m);assert prob==F(r['overlap_probability_rational_upper']);total+=prob
  assert total==F(grid['sum_delta_exact'])
  for r in minor['uniform_Gaussian_nonpositivity_certificates']:
   assert F(r['uniform_Gaussian_bridge_rational_upper'])<0 and r['all_s_gaussian_reference_strictly_negative_certified']
  out=dict(status='PASS',runs=runs,true_source_scales=len(affinity['rows']),true_source_probability_caps_exact=True,finite_probability_cases=finite['cases'],grid_error_exact=str(total),Gaussian_calibration_all_s_negative_widths=[r['width']for r in minor['uniform_Gaussian_nonpositivity_certificates']],multibank_result_keys=list(multi),new_independent_queries=0,actual_HL_new_log_gain_proved=False,seconds=time.time()-st)
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
