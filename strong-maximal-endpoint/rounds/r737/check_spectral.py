import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
import json,time,re,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import mpmath as mp
HERE=Path(__file__).resolve().parent

def one(h):
 start=time.monotonic();term=1;total=0
 for j in range(h//2):
  total+=term;term=term*(h-j)*3//(j+1)
 numerator=2*total+term;denominator=2*4**h
 vals=[]
 for prec in [80,120]:
  with mp.workdps(prec):
   v=mp.log(denominator)-mp.log(numerator);cap=h*mp.log(2/mp.sqrt(3));gap=v-cap-mp.log(h)/2
   vals.append([mp.nstr(x,70) for x in [v,cap,gap]])
 assert vals[0]==vals[1]
 return dict(h=h,dimension=3*h,NP_numerator_bits=numerator.bit_length(),NP_denominator_bits=denominator.bit_length(),V=vals[1][0],fixed_capacity=vals[1][1],centered_gap=vals[1][2],precision_pair=[80,120],saved_digits=70,elapsed=time.monotonic()-start)
if __name__=='__main__':
 start=time.monotonic();wl=json.loads((HERE/'check_formula_wolfram.json').read_text());hs=[r['h'] for r in wl['spectral_diagnostics']]
 with ProcessPoolExecutor(max_workers=12) as p:res=list(p.map(one,hs))
 with mp.workdps(100):
  errors=[]
  for a,b in zip(res,wl['spectral_diagnostics']):
   for key in ['V','fixed_capacity','centered_gap']:
    val=b[key].split('`')[0].replace('*^','e');errors.append(abs(mp.mpf(a[key])-mp.mpf(val)))
  error=max(errors);assert error<mp.mpf('1e-60');err=mp.nstr(error,12)
 out=dict(status='complete',workers=12,method='exact integer NP tail then independent 80/120 digit logarithms',cross_software_max_absolute_error=err,results=res,elapsed=time.monotonic()-start)
 (HERE/'check_spectral.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','cross_software_max_absolute_error','elapsed']}))
