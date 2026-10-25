"""vg11 exact32-key/count/mass agreement and actual four-moving-coset geometry."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import sys,json,time,hashlib,numpy as np
from scipy.stats import qmc
from round011_embedding import blocks,split_block,COSETS,SIGNS
from round011_chess_engine import decode,capture,boxes
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();vg=P.parent/'variation_geometry'
    if not(vg/'vg11_chess.py').exists():vg=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-variation/cube-finite-dim-bounds/lower/variation_geometry')
    sys.path.insert(0,str(vg));import vg11_chess as chess
    checked=0;classes={}
    for lo,hi,ell in blocks():
        ours=split_block(lo,hi);other=chess.split_block(lo,hi);assert ours==other['key'] and other['N']==int(np.prod(ell));classes.setdefault(other['counts'],ours);checked+=1
    masschecks=0;delta_signs=set()
    for counts,key in classes.items():
        n=sum(counts);delta=sum(s*k for s,k in zip(SIGNS,counts));delta_signs.add((delta>0)-(delta<0))
        for e in [F(-1,20),F(1,20),F(1,10)]:assert chess.finite_mass(key,e)==n+e*delta;masschecks+=1
        d=chess.radius_descriptor(key,F(1,10),F(1,2));assert d['complete_supercell_mass']==4 and d['complete_mass_slope']==0 and d['radius_velocity_factor']==F(delta,40*n)
    h=[F(1,100),F(1,50),F(1,100),F(0)];x=np.r_[np.zeros(4),.05,list(map(float,h))];periods,Q,ph,w=decode(x);exact=chess.coset_phases(h);assert np.allclose(ph,np.array(exact,float),atol=1e-15)
    anchored=chess.coset_phases(h,anchor=True);u=chess.super_shift_direction(h)
    for i,((p,q),s) in enumerate(zip(COSETS,SIGNS)):
        base=(F(p,2),F(q,2),F(0),F(0));assert anchored[i]==tuple(base[j]+(s-1)*u[j] for j in range(4))
    keys=capture(Q,ph,w,qmc.Sobol(4,scramble=True,seed=20261025113).random_base2(7));gchecks=0
    for key in sorted(keys)[:96]:
        ll=np.array(key[:16]).reshape(4,4);hh=np.array(key[16:]).reshape(4,4);points=[];mass=F(0)
        for orbit,((p,q),s) in enumerate(zip(COSETS,SIGNS)):
            inds=list(product(*(range(int(ll[orbit,j]),int(hh[orbit,j])+1) for j in range(4))))
            for index in inds:points.append(np.array([2*index[0]+p,2*index[1]+q,index[2],index[3]])*periods+s*periods*np.array(list(map(float,h))))
            mass+=len(inds)*(1+s*F(1,20))
        assert mass==chess.finite_mass(key,F(1,20));points=np.array(points);r=float(mass)**.25/2;lo=np.maximum(0,points.max(axis=0)-r);hi=np.minimum(Q,points.min(axis=0)+r);bb,_=boxes(Q,ph,w,[key])
        if np.all(hi>lo+1e-14):assert len(bb)==1 and np.allclose(bb[0],np.stack([lo,hi],axis=1).ravel(),atol=1e-13,rtol=1e-13)
        else:assert not len(bb)
        gchecks+=1
    out=dict(original_blocks_checked=checked,distinct_count_classes=len(classes),exact_mass_checks=masschecks,signed_delta_cases=sorted(delta_signs),all_four_cosets_move_correctly=True,anchoring_negative_cosets_minus2u_passed=True,actual_capture_keys=len(keys),explicit_physical_point_keys=gchecks,vg11_source_sha256=hashlib.sha256((vg/'vg11_chess.py').read_bytes()).hexdigest(),all_passed=True,seconds=time.monotonic()-st)
    (P/'round011_validation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
