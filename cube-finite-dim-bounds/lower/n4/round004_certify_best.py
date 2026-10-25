"""Exact outer cover of the best selected finite satellite witness union only."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json,time,hashlib
import numpy as np
from scipy.stats import qmc
from round004_satellite_search import A as AF,PRIMARY,capture
from certify_scan import A,scale,radius,floor,ceil,volume
P=Path(__file__).resolve().parent;H=10**12

def main():
    st=time.monotonic();results=json.loads((P/'round004_results.json').read_text())
    row=max((r for r in results['rows'] if r['epsilon']>0),key=lambda r:r['ratio'])
    assert row['a']==AF.tolist() and row['centers']==4096
    eps=F(str(row['epsilon']));phi=list(map(F,row['phi']))
    extra=capture(AF,np.array(row['phi']),float(eps),qmc.Sobol(4,scramble=True,seed=20261025405).random_base2(12))
    extra|=capture(AF,np.array(row['phi']),float(eps),qmc.Sobol(4,scramble=True,seed=20261025404).random_base2(9))
    keys=set(PRIMARY)|extra
    added=sorted(extra-set(PRIMARY));(P/'round004_best_keys.json').write_text(json.dumps([list(k) for k in added],separators=(',',':'))+'\n')
    pos={(orb,j,k):scale(A[j],F(k)+(phi[j] if orb else 0)) for orb,j,k in product(range(2),range(4),range(-8,9))}
    boxes=set();extrema=[]
    for key in keys:
        ll=[key[:4],key[4:8]];hh=[key[8:12],key[12:]]
        cnt=[max(0,np.prod([max(0,hh[o][j]-ll[o][j]+1) for j in range(4)]).item()) for o in range(2)]
        rr=radius(F(cnt[0])+eps*int(cnt[1]));box=[];valid=True
        active=[o for o in range(2) if cnt[o]>0]
        for j in range(4):
            lower=max(F(0),max(pos[o,j,hh[o][j]][0] for o in active)-rr[1])
            upper=min(A[j][1],min(pos[o,j,ll[o][j]][1] for o in active)+rr[1])
            if upper<=lower:valid=False;break
            box.extend([floor(lower*H),ceil(upper*H)])
        if valid:boxes.add(tuple(box))
        for o in active:extrema.extend(ll[o]+hh[o])
    assert min(extrema)>=-8 and max(extrema)<=8
    # Integer outer boxes, exact joint union. No floating-point containment pruning.
    b=np.array(sorted(boxes),dtype=np.int64);lo=b[:,::2];hi=b[:,1::2];keep=[]
    for i in range(len(b)):
        included=np.all(lo<=lo[i],axis=1)&np.all(hi>=hi[i],axis=1);included[i]=False
        if not included.any():keep.append(tuple(map(int,b[i])))
    upper=F(volume(tuple(keep)),H**4)/(1+eps)
    old=json.loads((P.parents[1]/'baseline/original_certificates.json').read_text())[1];oldlo=F(old['lower_bound_interval'][0]);assert upper<oldlo
    out=dict(scope='Exact upper enclosure of selected finite witness union only; NOT an upper bound for the periodic input weak ratio or the family.',a='exact inherited R4 formulas',epsilon=str(eps),phase_exact=[str(x) for x in phi],key_count=len(keys),additional_keys=len(added),outer_boxes=len(keep),rounding_denominator=H,finite_union_ratio_upper=str(upper),old_R4_lower=str(oldlo),strictly_below_old=True,index_min=int(min(extrema)),index_max=int(max(extrema)),transfer_margin=8,transfer_mass='(1+epsilon)*(2N+17)^4',elapsed=time.monotonic()-st,input_sha256={name:hashlib.sha256((P/name).read_bytes()).hexdigest() for name in ['round004_results.json','round004_best_keys.json','round004_preregister.md']})
    (P/'round004_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
