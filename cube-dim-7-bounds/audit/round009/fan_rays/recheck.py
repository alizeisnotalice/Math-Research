#!/usr/bin/env python3
"""Independent exact R009 cofactor audit; stdout only; frozen sources read only."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
from itertools import product, combinations, islice, accumulate
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor
from collections import Counter,defaultdict
from pathlib import Path
import hashlib, json, sys, time, math, resource, argparse
from flint import fmpz_mat
import flint
HERE=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,s):
 if not x: raise ValueError(s)
def build():
 body=[(1,)+x for x in product((-1,1),repeat=6)]
 prefix=[tuple(accumulate(v)) for v in body]
 nn=[v for v in prefix if min(v)>=0]; mixed=[v for v in prefix if min(v)<0]
 need(len(body)==64 and len(nn)==35 and len(mixed)==29 and len(set(prefix))==64,'prefix classification')
 def vals(s): return tuple(1 if c=='+' else -1 for c in s)
 words=[vals('+--')+x for x in product((-1,1),repeat=4)]
 words += [vals(p)+x for p in ('++---','+-+--') for x in product((-1,1),repeat=2)]
 words += [vals(p) for p in ('+++----','++-+---','++--+--','+-++---','+-+-+--')]
 ordered=[tuple(accumulate(v)) for v in words]
 need(len(ordered)==29 and set(ordered)==set(mixed),'first-hit walk reconstruction')
 walls=[tuple(int(i==j) for j in range(7)) for i in range(7)]
 return body,prefix,walls+ordered
N=None
def init(n):
 global N
 N=n
def work(job):
 start,batch=job; t=time.process_time(); c=Counter(); rays={}; strata=[0]*7
 for ids in batch:
  a=[N[i] for i in ids]; co=[]
  for j in range(7):
   co.append((-1 if j%2 else 1)*int(fmpz_mat([[r[k] for k in range(7) if k!=j] for r in a]).det()))
  strata[sum(i<7 for i in ids)]+=1
  if not any(co): c['rank_below_6']+=1; continue
  if min(co)<0<max(co): c['rank6_mixed']+=1; continue
  if max(co)<=0: co=[-x for x in co]
  if min(co)<0 or max(co)<=0: raise ArithmeticError('orientation')
  g=math.gcd(*co); w=tuple(x//g for x in co)
  if any(sum(r[j]*w[j] for j in range(7)) for r in a): raise ArithmeticError('cofactor kernel')
  c['rank6_feasible']+=1
  if w not in rays: rays[w]=ids
 return {'start':start,'end':start+len(batch),'counts':dict(c),'strata':strata,'rays':[(w,ids) for w,ids in rays.items()],'pid':os.getpid(),'cpu':time.process_time()-t,'rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
def main():
 p=argparse.ArgumentParser(); p.add_argument('--candidate',type=Path,required=True); p.add_argument('--theory',type=Path,required=True); a=p.parse_args()
 lock=json.loads((HERE/'input_lock.json').read_text()); cand=a.candidate.resolve(); theory=a.theory.resolve()
 for name,h in lock['candidate'].items(): need(sha(cand/name)==h,'candidate hash '+name)
 for name,h in lock['theory'].items(): need(sha(theory/name)==h,'theory hash '+name)
 body,prefix,normals=build(); inp=json.loads((cand/'inputs.json').read_text()); th=json.loads((theory/'fan_normals.json').read_text())
 need([list(x) for x in body]==inp['body_vectors'] and [list(x) for x in prefix]==inp['prefix_normals'],'body input mismatch')
 need([list(x) for x in normals]==inp['arrangement_normals']==th['facet_normals']+th['internal_prefix_normals'],'ordered normal mismatch')
 total=math.comb(36,6); need(total==1947792,'subset total')
 starttime=time.perf_counter(); it=combinations(range(36),6)
 def jobs():
  start=0
  while True:
   batch=list(islice(it,2048))
   if not batch:break
   yield start,batch; start+=len(batch)
 results=[]
 with ProcessPoolExecutor(max_workers=15,initializer=init,initargs=(normals,)) as pool:
  for out in pool.map(work,jobs(),chunksize=1,buffersize=30):results.append(out)
 elapsed=time.perf_counter()-starttime
 counts=Counter(); rays={}; cpus=defaultdict(float); rss={}; strata=[0]*7; end=0
 for r in results:
  need(r['start']==end and r['end']>end,'chunk coverage'); end=r['end']; counts.update(r['counts']); cpus[str(r['pid'])]+=r['cpu']; rss[str(r['pid'])]=max(rss.get(str(r['pid']),0),r['rss'])
  strata=[x+y for x,y in zip(strata,r['strata'])]
  for w,ids in r['rays']:
   w=tuple(w); ids=tuple(ids)
   if w not in rays or ids<rays[w]:rays[w]=ids
 need(end==total and sum(counts.values())==total,'incomplete enumeration')
 need(strata==[math.comb(7,i)*math.comb(29,6-i) for i in range(7)],'strata mismatch')
 cat=[]
 for idx,(w,ids) in enumerate(sorted(rays.items())):
  x=[sum(w[i:]) for i in range(7)]; den=sum(x); active=[i for i,n in enumerate(normals) if sum(a*b for a,b in zip(n,w))==0]
  signs=[(t>0)-(t<0) for n in prefix for t in [sum(a*b for a,b in zip(n,w))]]
  minor=None
  for omit in range(7):
   cols=[i for i in range(7) if i!=omit]; d=int(fmpz_mat([[normals[j][k] for k in cols] for j in ids]).det())
   if d:minor={'columns':cols,'determinant':d};break
  need(minor is not None and set(ids)<=set(active),'witness')
  cat.append(dict(id=idx,w_primitive=list(w),xi_primitive=x,normalization_denominator=den,xi_normalized=[str(Fraction(i,den)) for i in x],active_normals=active,normal_incidence=[int(i in active) for i in range(36)],body_dot_signs=signs,witness_rows=list(ids),witness_minor=minor))
 claim=json.loads((cand/'rays.json').read_text()); rc=json.loads((cand/'receipt.json').read_text())
 need(cat==claim,'full ray data mismatch')
 for k in ('rank_below_6','rank6_mixed','rank6_feasible'):need(counts[k]==rc[k],'category mismatch '+k)
 need(len(cat)==rc['unique_rays'] and strata==rc['wall_strata'],'receipt mismatch')
 for name,h in lock['candidate'].items(): need(sha(cand/name)==h,'candidate changed '+name)
 for name,h in lock['theory'].items(): need(sha(theory/name)==h,'theory changed '+name)
 return {'status':'PASS','scope':'all 1947792 six-subsets independently computed using seven FLINT integer determinants','combination_count':end,'category_counts':dict(counts),'wall_strata':strata,'unique_rays':len(cat),'complete_rays':cat,'worker_cpu_seconds':dict(sorted(cpus.items())),'worker_peak_rss_bytes':rss,'wall_seconds':elapsed,'chunks':len(results),'workers_requested':15,'flint_version':flint.__version__,'executable':sys.executable,'source_hashes_unchanged':True}
if __name__=='__main__':
 try: print(json.dumps(main(),indent=2));sys.exit(0)
 except BaseException as e: print(json.dumps({'status':'FAIL','error':repr(e)}),flush=True);sys.exit(1)
