"""Stdlib exact dual verification for one finite lattice model only."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations_with_replacement, permutations
import json

ROOT=Path(__file__).resolve().parents[1]

def verify():
    rec=json.loads((ROOT/'results/round005_lp_dual.json').read_text())
    q=rec['q'];assert q==4 and rec['closed_obstacle']
    points=list(combinations_with_replacement(range(3*q),3))
    lookup={p:i for i,p in enumerate(points)}
    c=[F(len(set(permutations(p)))*2**sum(v>0 for v in p),q**3) for p in points]
    z=c[:];ys=list(map(F,rec['y']));assert len(ys)==len(points)-1
    for p,y in zip(points[1:],ys):
        assert y<=0
        z[lookup[p]]-=6*y
        for j in range(3):
            for d in (-1,1):
                v=list(p);v[j]+=d;v=tuple(sorted(abs(x) for x in v))
                if v[-1]<3*q:z[lookup[v]]+=y
    assert min(z)>=0
    lower=sum((v/F(8) for p,v in zip(points,z) if p[-1]<=q),F(0))
    assert lower==F(rec['dual_bound']) and lower>8
    return {'status':'PASS exact finite LP dual','bound':str(lower),'minimum_slack':str(min(z)),'continuous_bound':False}

if __name__=='__main__':print(json.dumps(verify(),indent=2))
