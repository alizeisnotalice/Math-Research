from pathlib import Path
import numpy as np,json,itertools,time
P=Path(__file__).resolve().parent;results=[];start=time.perf_counter()
def rank(vs):
 bas={}
 for x in vs:
  v=int(x)
  while v:
   bit=v.bit_length()-1
   if bit in bas:v^=bas[bit]
   else:bas[bit]=v;break
 return len(bas)
for d,poly in [(10,0x409),(12,0x1053)]:
 K=(1<<d)-1;ex=[];z=1
 for _ in range(K):
  ex.append(z);z<<=1
  if z&(1<<d):z^=poly
 assert z==1 and len(set(ex))==K
 ex=np.array(ex);log=np.full(K+1,-1);log[ex]=np.arange(K)
 def sub(l):
  a=int(ex[l]);b=int(ex[(l+log[6])%K]);return [a,b,a^b]
 # Count every defective third lag for each pair-safe (0,b), without sampling.
 defect_counts=[]
 for second in range(1,K):
  u=sub(0);v=sub(second)
  if rank(u[:2]+v[:2])!=4:continue
  span=set(u+v+[a^b for a in u for b in v])
  candidates3={(int(log[z])-int(log[w]))%K for z in span for w in [1,6,7]}
  badthird=[]
  for third in candidates3:
   ww=sub(third)
   if third not in [0,second] and rank(u[:2]+ww[:2])==4 and rank(v[:2]+ww[:2])==4:
    assert rank(u[:2]+v[:2]+ww[:2])<6
    badthird.append(third)
  defect_counts.append(len(badthird))
 print('defective thirds',d,min(defect_counts),max(defect_counts),float(np.mean(defect_counts)),flush=True)
 defect_summary=dict(pair_count=len(defect_counts),min=min(defect_counts),max=max(defect_counts),mean=float(np.mean(defect_counts)),histogram={str(k):defect_counts.count(k) for k in sorted(set(defect_counts))})
 for spanname,candidates in [('window',np.arange(K//9+1)),('full',np.arange(K))]:
  trials=[]
  for trial in range(8):
   order=candidates if trial==0 else np.random.default_rng(np.random.SeedSequence([2026102901,d,trial])).permutation(candidates)
   chosen=[];spaces=[];bad=set()
   for k in order:
    u=sub(int(k))
    if any(v in bad for v in u):continue
    for v in spaces:
     bad.update(u);bad.update(v);bad.update(a^b for a in u for b in v)
    bad.update(u);spaces.append(u);chosen.append(int(k))
   for ij in itertools.combinations(range(len(chosen)),2):assert rank([v for i in ij for v in spaces[i][:2]])==4
   for ij in itertools.combinations(range(len(chosen)),3):assert rank([v for i in ij for v in spaces[i][:2]])==6
   trials.append(dict(trial=trial,size=len(chosen),lags=chosen))
  result=dict(d=d,K=K,defective_thirds=defect_summary,domain=spanname,candidates=len(candidates),trials=trials,best=max(t['size'] for t in trials));results.append(result)
  print(d,spanname,len(candidates),[t['size'] for t in trials],flush=True)
(P/'results.json').write_text(json.dumps(dict(seconds=time.perf_counter()-start,seed=2026102901,results=results),indent=2))
