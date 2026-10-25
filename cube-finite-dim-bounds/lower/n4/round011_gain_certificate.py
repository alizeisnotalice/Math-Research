"""One existing chess point: exact new-witness cube outside all old embedded boxes."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import json,time,numpy as np
from scipy.stats import qmc
from round011_chess_engine import BASE,decode,capture,boxes
from round011_embedding import COSETS,SIGNS
from certify_scan import A as INTERVAL_A,scale,radius
P=Path(__file__).resolve().parent

def find_uncovered(target,old):
    lo,hi=target[0],target[1];cuts=sorted({lo,hi}|{float(v) for b in old for v in b[:2] if lo<v<hi})
    for a,b in zip(cuts,cuts[1:]):
        if b<=a+1e-14:continue
        active=[z[2:] for z in old if z[0]<=a and b<=z[1]]
        if not active:return [(a+b)/2]+[(target[j]+target[j+1])/2 for j in range(2,len(target),2)]
        if len(target)>2:
            answer=find_uncovered(target[2:],active)
            if answer is not None:return [(a+b)/2]+answer
    return None

def coefficients(keys,h100):
    data=np.array(list(keys),int);ll=data[:,:16].reshape(-1,4,4);hh=data[:,16:].reshape(-1,4,4);counts=np.prod(np.maximum(0,hh-ll+1),axis=2);offset=np.array([[p,q,0,0] for p,q in COSETS]);multipliers=np.array([2,2,1,1]);sh=np.array(SIGNS)[:,None]*h100
    low=(ll*multipliers+offset)*100+sh;high=(hh*multipliers+offset)*100+sh;low[counts==0]=10**9;high[counts==0]=-10**9
    return low.min(axis=1),high.max(axis=1),counts

@lru_cache(None)
def position(j,c):return scale(INTERVAL_A[j],F(int(c),100))

def main():
    st=time.monotonic();reg=json.loads((P/'round011_preregister.json').read_text());spec=reg['directions'][6];eps=F(spec['epsilon']);h=list(map(F,spec['h']));h100=np.array([int(q*100) for q in h]);x=np.r_[np.zeros(4),float(eps),list(map(float,h))];periods,Q,ph,w=decode(x);extra=capture(Q,ph,w,qmc.Sobol(4,scramble=True,seed=20261025111).random_base2(8));oldb,_=boxes(Q,ph,w,BASE);newb,newkeys=boxes(Q,ph,w,BASE|extra);chosen=None
    for b,key in zip(newb,newkeys):
        if key in BASE:continue
        point=find_uncovered(tuple(b),[tuple(z) for z in oldb])
        if point is not None:chosen=(key,point);break
    assert chosen is not None
    key,pfloat=chosen;point=list(map(lambda z:F(str(z)),pfloat));cmn,cmx,counts=coefficients([key],h100);n=int(sum(counts[0]));delta=int(counts[0]@np.array(SIGNS));m=F(n)+eps*delta;rlo,rhi=radius(m);margins=[]
    for j in range(4):
        lower_upper=position(j,int(cmx[0,j]))[1]-rlo;upper_lower=position(j,int(cmn[0,j]))[0]+rlo;margins.extend([point[j]-lower_upper,upper_lower-point[j],point[j],INTERVAL_A[j][0]*(2 if j<2 else 1)-point[j]])
    assert min(margins)>0
    oldkeys=sorted(BASE);oldmin,oldmax,oldcounts=coefficients(oldkeys,h100);oldn=oldcounts.sum(axis=1);olddelta=oldcounts@np.array(SIGNS);oldmass=oldn+float(eps)*olddelta;rf=oldmass**.25/2;lofloat=oldmax/100*periods-rf[:,None];hifloat=oldmin/100*periods+rf[:,None];pv=np.array(pfloat);axis=np.argmax(np.maximum(lofloat-pv,pv-hifloat),axis=1);minimum=None;hist=[0]*4
    for i in range(len(oldkeys)):
        mass=F(int(oldn[i]))+eps*int(olddelta[i]);_,rh=radius(mass);ordered=[int(axis[i])]+[j for j in range(4) if j!=axis[i]];gap=None
        for j in ordered:
            left=position(j,int(oldmax[i,j]))[0]-rh-point[j];right=point[j]-(position(j,int(oldmin[i,j]))[1]+rh)
            if max(left,right)>0:gap=max(left,right);hist[j]+=1;break
        assert gap is not None
        minimum=gap if minimum is None else min(minimum,gap)
    margin=min([minimum,*margins,rlo]);limit=3*margin/4;D=10**12;eta=F((limit*D).numerator//(limit*D).denominator,D);assert 0<eta<margin
    residual=margin-eta;contraction=F(((residual/2)*D).numerator//((residual/2)*D).denominator,D);assert 0<contraction<residual and contraction<rlo
    assert rlo**4<=m/16<rhi**4
    gain=(2*eta)**4/4
    out=dict(status='strict positive selected-witness increment at existing case6; NOT a new R4/input/family bound',case=6,epsilon=str(eps),h=list(map(str,h)),periods='original exact R4 formulas',new_key=list(map(int,key)),new_counts=list(map(int,counts[0])),new_mass=str(m),new_radius_interval=[str(rlo),str(rhi)],new_radius_fourth_power=str(m/16),rational_radius_contraction=str(contraction),rational_center=list(map(str,point)),cube_infinity_radius=str(eta),old_embedded_boxes_excluded=len(oldkeys),old_exclusion_axis_histogram=hist,old_minimum_exclusion_margin=str(minimum),new_inside_and_domain_minimum_margin=str(min(margins)),strict_finite_normalized_gain_lower=str(gain),strict_gain_float=float(gain),radius_root_enclosure='Fraction fourth-root brackets denominator2*10^16',proof='All points within eta of center remain strictly in new receiver box and strictly outside every old raw receiver box; mass denominator4. Contract common radius by the positive rational amount, smaller than residual coordinate margin, to obtain strict threshold1.',seconds=time.monotonic()-st)
    (P/'round011_gain_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if 'margin' not in k},indent=2))
if __name__=='__main__':main()
