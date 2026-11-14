"""Exact integer grid-witness certificate for R=4, breakpoints 0,1,2,4.
Specified rectangles only, not complete halo. No floating-point classification.
"""
import json,math,time,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent

def mul(poly,a):
 out=[a*poly[0]]
 for j in range(1,len(poly)):out.append(a*poly[j]+2*poly[j-1])
 out.append(2*poly[-1]);return out

def one(n):
 started=time.monotonic();difference=[0]*(n+2);base=[1];classified=0
 for c in range(n+1):
  poly=base[:]
  for b in range(n-c+1):
   a=n-b-c;den=3**b*5**c;tail=den;last=a
   for j in range(1,len(poly)):
    tail-=poly[j-1]
    if 2*tail>den:last=a+j
    else:break
   weight=math.comb(n,c)*math.comb(n-c,b)*2**c
   difference[0]+=weight;difference[last+1]-=weight;classified+=1
   poly=mul(poly,1)
  base=mul(base,3)
 vols=[];v=0
 for k in range(n+1):v+=difference[k];vols.append(v)
 source=[0]*(n+1);running=0
 for k in range(n,-1,-1):running+=math.comb(n,k)*3**(n-k);source[k]=running
 best=max(range(1,n+1),key=lambda k:F(vols[k],source[k]));ratio=F(vols[best],source[best])
 return dict(n=n,best_k=best,witness_positive_orthant_volume=str(vols[best]),source_positive_orthant_volume=str(source[best]),ratio=str(ratio),numerical_log_ratio=math.log(ratio.numerator)-math.log(ratio.denominator),classified_cells=classified,elapsed=time.monotonic()-started)
if __name__=='__main__':
 start=time.monotonic()
 results=[]
 with ProcessPoolExecutor(max_workers=12) as pool:
  futures=[pool.submit(one,n) for n in [8,16,32,64,128,256,320,384,448,512,640,768]]
  for future in as_completed(futures):
   rec=future.result();results.append(rec);(HERE/f'exact_grid_n{rec["n"]}.json').write_text(json.dumps(rec)+'\n');print(json.dumps({k:rec[k] for k in ['n','best_k','numerical_log_ratio','elapsed']}),flush=True)
 results.sort(key=lambda r:r['n'])
 out=dict(status='complete',integer_exact_classification=True,workers=12,grid=[0,1,2,4],R=4,threshold='1/2 strict',scope='Certified union of count-grid cells with explicit rectangles; not full halo nor general endpoint coefficient',results=results,elapsed=time.monotonic()-start)
 (HERE/'exact_grid.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='complete',elapsed=out['elapsed'],summary=[{k:r[k] for k in ['n','best_k','numerical_log_ratio','classified_cells','elapsed']} for r in results])))
