"""Independent vg8 parity/mass interface and physical displacement checks."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import sys,json,time,hashlib,numpy as np
from scipy.stats import qmc
from round008_embedding import original_blocks,split_block
from round008_superperiod import A,decode,capture,boxes
P=Path(__file__).resolve().parent

def main():
    st=time.monotonic();vg=P.parent/'variation_geometry'
    if not(vg/'vg8_parity.py').exists():vg=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-variation/cube-finite-dim-bounds/lower/variation_geometry')
    sys.path.insert(0,str(vg));import vg8_parity as parity
    checked=0;classes={}
    for lo,hi,ell in original_blocks():
        ours=split_block(lo,hi);other=parity.split_block(lo,hi);assert ours==other['key'] and other['N']==int(np.prod(ell));checked+=1;classes.setdefault((other['N_even'],other['N_odd']),ours)
    masschecks=0;signed=set()
    for (ne,no),key in classes.items():
        for e in [F(1,20),F(-1,20),F(1,10)]:assert parity.finite_mass(key,e)==ne+no+e*(ne-no);masschecks+=1
        d=parity.radius_descriptor(key,F(1,10),F(1,2));assert d['complete_supercell_mass']==2 and d['complete_mass_slope']==0;assert d['radius_velocity_factor']==F(ne-no,40*(ne+no));signed.add((ne>no)-(ne<no))
    x=np.r_[np.zeros(4),.05,[.02,.02,.02,0.]];periods,Q,phases,weights=decode(x);keys=capture(Q,phases,weights,qmc.Sobol(4,scramble=True,seed=20261025810).random_base2(7));geometrychecks=0
    for key in sorted(keys)[:128]:
        ll=np.array(key[:8]).reshape(2,4);hh=np.array(key[8:]).reshape(2,4);points=[];mass=F(0)
        for orb in [0,1]:
            inds=list(product(*(range(int(ll[orb,j]),int(hh[orb,j])+1) for j in range(4))))
            for index in inds:
                oldindex=np.array([2*index[0]+orb,*index[1:]]);pos=oldindex*periods+(periods*np.array(x[5:]) if orb else 0);points.append(pos)
            mass+=len(inds)*(1+F(1,20) if orb==0 else 1-F(1,20))
        assert mass==parity.finite_mass(key,F(1,20));points=np.array(points);r=float(mass)**.25/2;lo=np.maximum(0,points.max(axis=0)-r);hi=np.minimum(Q,points.min(axis=0)+r);bb,_=boxes(Q,phases,weights,[key])
        if np.all(hi>lo+1e-14):assert len(bb)==1 and np.allclose(bb[0],np.stack([lo,hi],axis=1).ravel(),atol=1e-13,rtol=1e-13)
        else:assert not len(bb)
        geometrychecks+=1
    out=dict(original_blocks_checked=checked,distinct_parity_count_classes=len(classes),exact_mass_checks=masschecks,signed_delta_cases=sorted(signed),actual_capture_keys=len(keys),explicit_physical_point_checks=geometrychecks,vg8_source_sha256=hashlib.sha256((vg/'vg8_parity.py').read_bytes()).hexdigest(),all_passed=True,seconds=time.monotonic()-st)
    (P/'round008_validation_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
