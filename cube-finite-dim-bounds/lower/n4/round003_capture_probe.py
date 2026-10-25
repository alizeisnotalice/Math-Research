"""Refresh unrestricted captured-prefix witnesses after joint parameter exploration."""
from pathlib import Path
from itertools import product
import json,time
import numpy as np
from scipy.stats import qmc
from round003_stair_search import full_boxes,volume
P=Path(__file__).resolve().parent

def probe(z,power=12):
    st=time.monotonic();a=np.array(z[:4]);f,h,w=z[4:];q=f*a[1];L=a.copy();L[0]*=2
    b,_=full_boxes(z);base=volume(b)[0]/2
    centers=qmc.Sobol(4,scramble=True,seed=2026102504).random_base2(power)*L
    covered=np.zeros(len(centers),bool)
    for box in b:covered|=np.all(centers>=box[::2],axis=1)&np.all(centers<=box[1::2],axis=1)
    ids=np.array(list(product(range(-3,5),repeat=4)));parity=ids[:,0]%2;sign=1-2*parity
    points=ids*a;points[:,0]=2*a[0]*(ids[:,0]//2)+parity*2*a[0]*h;points[:,1]+=sign*q/2
    weights=1+sign*w;added={}
    for c in centers[~covered]:
        dist=np.max(abs(points-c),axis=1);idx=np.argsort(dist)[:256];mass=np.cumsum(weights[idx]);score=mass/(2*dist[idx])**4;best=np.argmax(score)
        if score[best]<=1+1e-10:continue
        ss=idx[:best+1];key=tuple(sorted(map(int,ss)))
        if key in added:continue
        r=mass[best]**.25/2;lo=np.maximum(0,np.max(points[ss],axis=0)-r);hi=np.minimum(L,np.min(points[ss],axis=0)+r)
        if np.all(hi>lo):added[key]=np.stack([lo,hi],axis=1).ravel()
    union=b if not added else np.concatenate([b,np.array(list(added.values()))])
    enriched=float(volume(union)[0]/2)
    return dict(z=list(map(float,z)),baseline_ratio=float(base),enriched_ratio=enriched,sobol_centers=2**power,seed=2026102504,uncovered_centers=int(sum(~covered)),added_subsets=len(added),subsets=[list(k) for k in added],seconds=time.monotonic()-st)

def main():
    data=json.loads((P/'round003_search_results.json').read_text());rows=data['rows']
    choices=[max(rows,key=lambda r:r['ratio']),max((r for r in rows if r['z'][4]>=.03),key=lambda r:r['ratio']),max((r for r in rows if r['z'][4]>=.03 and (abs(r['z'][5]-.5)>=.03 or abs(r['z'][6])>=.03)),key=lambda r:r['ratio'])]
    out=[]
    for x in choices:
        r=probe(x['z']);out.append(r);print({k:v for k,v in r.items() if k!='subsets'},flush=True)
    (P/'round003_probe_results.json').write_text(json.dumps(dict(status='numerical finite witness refresh only',selection='best overall, best f>=.03, best f>=.03 with |h-.5|>=.03 or |w|>=.03',rows=out),indent=2)+'\n')
if __name__=='__main__':main()
