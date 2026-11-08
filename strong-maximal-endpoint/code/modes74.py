import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from columns71 import ratio_sign

def run(j):
 n=j['N'];c=j['c'];s=j['S'];last=min(n,math.ceil(s/c)-1)
 def lr(k):
  g=s-c*k
  return -math.inf if g<=c else math.log(n-k)-math.log(k+1)+(n-1)*math.log1p(-c/g)
 lo,hi=0,last
 while lo<hi:
  m=(lo+hi)//2
  if lr(m)>0:lo=m+1
  else:hi=m
 mode=lo;ss=F.from_float(s);cc=F.from_float(c);assert ss-cc*last>0 and (last==n or ss-cc*(last+1)<=0)
 a=ratio_sign(n,ss,cc,mode,'nonpositive') if mode<last else dict(certified=True,endpoint=True);b=ratio_sign(n,ss,cc,mode-1,'nonnegative') if mode else dict(certified=True,endpoint=True);assert a['certified'] and b['certified']
 return dict(j,mode=mode,right=a,left=b)
if __name__=='__main__':
 inp=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True);d=json.load(open(inp));unique={(r['N'],c,node['S']) for r in d['rows'] for c in [r['c_actual'],r['c_reference']] for node in r['quadrature_nodes']};jobs=[dict(N=n,c=c,S=s) for n,c,s in sorted(unique)];start=time.monotonic()
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,jobs,chunksize=4))
 (out/'modes74.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,workers=12,candidate_SHA256=__import__('hashlib').sha256(inp.read_bytes()).hexdigest(),rows=rows,scope='exact binary-rational c,S posterior modes only; no floating integral certificate'),indent=2)+'\n');(out/'modes_summary74.json').write_text(json.dumps(dict(count=len(rows),all_certified=True,maximum_terms=max(max(r['right'].get('terms',0),r['left'].get('terms',0)) for r in rows),candidate_SHA256=__import__('hashlib').sha256(inp.read_bytes()).hexdigest(),scope='exact binary-rational mode brackets; full per-node output reproducible with modes74.py'),indent=2)+'\n');print('PASS',len(rows),'exact posterior mode certificates')
