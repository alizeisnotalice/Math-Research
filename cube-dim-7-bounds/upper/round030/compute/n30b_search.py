#!/usr/bin/env python3
"""Finite fixed-t N30B exploration: 54 signed evaluations max per dimension."""
import os,sys,json,time,math,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
from scipy.integrate import quad
root=Path(__file__).resolve().parents[3]
rows=json.loads((root/'revision_v3/numeric/table.json').read_text())

def run(n):
 start=time.monotonic();cpu=time.process_time();tq=Q(next(x['t'] for x in rows if x['n']==n));t=float(tq);y=math.log(t);h=(1+1/t)/2;A=(1+1/t)**n
 B=1-h**n+n*quad(lambda z:1-.5*(1+1/(2*t)-.5*math.exp(-z))**(n-1),0,y,epsabs=1e-11)[0]
 calls=0
 def g(R):
  nonlocal calls
  calls+=1
  if calls>64: raise ValueError('call cap')
  s=R/t;D=s**n-((s+1/t)/2)**n
  P=n*quad(lambda z:1-.5*(1+1/(2*t)-.5*math.exp(-z))**(n-1),math.log(R),y,epsabs=1e-11)[0]
  return P-D-math.sqrt(A*(B-D))
 low=1.;high=t;gl=g(low);gh=g(high)
 if not gl>0 or not gh<0:return {'n':n,'pid':os.getpid(),'status':'no endpoint sign bracket','left':gl,'right':gh,'calls':calls}
 for _ in range(52):
  mid=(low+high)/2;gm=g(mid)
  if gm>0:low=mid
  else:high=mid
 # 1% retreat; fixed dyadic denominator, exact candidate condition checked separately.
 R=Q(round((1+.99*(low-1))*2**20),2**20);s=R/tq
 assert 1<R<tq and Q(1,tq)<s<1
 return {'n':n,'t':str(tq),'root_bracket_float':[low,high],'R_candidate':str(R),'sigma_candidate':str(s),'candidate_float_g':g(float(R)),'calls':calls,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'floating bracket exploration only, no uniqueness/global claim; separate exact certification'}
if __name__=='__main__':
 start=time.monotonic()
 with concurrent.futures.ProcessPoolExecutor(max_workers=7) as ex: r=list(ex.map(run,(4,5,6,7,8,16,64)))
 out={'batch':'N30B','stage':'exploratory search','workers':7,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'wall_seconds':time.monotonic()-start,'results':r};Path(__file__).with_name('n30b_search.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(r,indent=2))
