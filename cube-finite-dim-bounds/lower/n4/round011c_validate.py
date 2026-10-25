"""Exact affine-extrema proof for fixed L/hull paths, not a union derivative."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
from certify_scan import A as IA,radius
from round011_embedding import split_block
from round011_chess_engine import BASE
P=Path(__file__).resolve().parent

def main():
    source=json.loads((P/'round011_gain_certificate.json').read_text());key=source['new_key'];pts=[tuple(k) for k in product(range(2),range(3),range(2),range(3))];L=[k for k in pts if (k[0]%2,k[1]%2)!=(1,1)];assert len(L)==30 and len(pts)==36
    groups=[[],[],[],[]]
    for k in L:groups[2*(k[0]%2)+k[1]%2].append((k[0]//2,k[1]//2,k[2],k[3]))
    ll=[];hh=[]
    for group in groups:
        ll+=list(map(min,zip(*group))) if group else [0]*4;hh+=list(map(max,zip(*group))) if group else [-1]*4
    assert ll+hh==key
    proof=[]
    for path in range(3):
        H=[F(1,50),F(1,50) if path==0 else F(0),F(1,100) if path==0 else F(0),F(0)];spans=[]
        for name,atoms in [('L',L),('hull',pts)]:
            axes=[]
            for j in range(4):
                linear=[(F(k[j]),H[j]*(1 if (k[0]+k[1])%2==0 else -1)) for k in atoms]
                minform=min(linear,key=lambda q:(q[0],q[1]));maxform=max(linear,key=lambda q:(q[0],q[1]))
                # Affine dominance at both endpoints is an exact all-t proof.
                assert all(minform[0]+t*minform[1]<=q[0]+t*q[1]<=maxform[0]+t*maxform[1] for q in linear for t in [0,1])
                axes.append(dict(min=list(map(str,minform)),max=list(map(str,maxform)),span=[str(maxform[k]-minform[k]) for k in range(2)]))
            spans.append(dict(name=name,axes=axes))
        proof.append(dict(path=path,affine_extrema=spans))
    assert len([k for k in pts if (k[0]+k[1])%2==0])==18
    assert [len(g) for g in groups]==[12,6,12,0]
    hulls=[]
    for shift in product(range(-1,2),repeat=4):
        lo=tuple(shift[j]*(2 if j<2 else 1) for j in range(4));hi=tuple(lo[j]+[1,2,1,2][j] for j in range(4));assert split_block(lo,hi) in BASE;hulls.append(lo)
    # All t and all3 paths: period factors<=1, spans bounded by A-path extrema.
    maxspan=[F(26,25),F(51,25),F(51,50),F(2)];rmin=radius(F(147,5))[0];margin=min(2*rmin-IA[j][1]*maxspan[j] for j in range(4));assert margin>0
    out=dict(status='exact fixed-key affine extrema and containing-hull identities only',Lkey_matches_original_integer_points=True,L_count=30,hull_count=36,L_counts=[12,6,12,0],L_delta=-6,hull_delta=0,all_t_interval=['0','1'],affine_extrema_proof=proof,all81_Q_translated_hulls_in_BASE=True,positive_receiver_width_lower=str(margin),gate='F=2P1*h1+rL-r36; F<=0 implies L receiver subset hull receiver; F>0 only permits escape, other old boxes can shield',derivative_scope='No assertion on full old-union derivative or global path sign')
    (P/'round011c_validation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status']);print('81 hulls,24 affine extrema systems passed; positive width lower',float(margin))
if __name__=='__main__':main()
