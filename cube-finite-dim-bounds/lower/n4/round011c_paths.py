"""Finite analytic chess paths: exact single-hull gate and full finite unions."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import json,time,numpy as np
from scipy.stats import qmc
from certify_scan import A as IA,radius
from round011_chess_engine import BASE,decode,capture,boxes,A
from round007_union_fast import volume
P=Path(__file__).resolve().parent;OLD=1.9017426111009554

def path(which,t):
    y=np.zeros(9);y[4]=float(t*(F(1,10) if which==0 else F(-1,5)));y[5]=float(t/50)
    if which==0:y[6]=float(t/50);y[7]=float(t/100)
    if which==2:y[0]=-np.log1p(float(t/25))
    return y

def gate(which,t):
    """F=2P1h1+rL-r36; all endpoints enclosed with Fraction."""
    alpha=F(1,10) if which==0 else F(-1,5);a,b=IA[0]
    if which==2:a/=1+t/25;b/=1+t/25
    r=radius(30-6*alpha*t);r36=radius(F(36))
    return 2*a*t/50+r[0]-r36[1],2*b*t/50+r[1]-r36[0]

def threshold(which):
    assert gate(which,F(0))[1]<0 and gate(which,F(1))[0]>0
    lo,hi=F(0),F(1)
    while hi-lo>F(1,2**40):
        mid=(lo+hi)/2;l,u=gate(which,mid)
        if u<0:lo=mid
        elif l>0:hi=mid
        else:raise AssertionError('root brackets too wide for requested isolation')
    # Uniform positive derivative: A has negative radius velocity; B/C radius increases.
    derivative=2*IA[0][0]/50
    if which==0:derivative-=F(3,5)/(64*radius(F(147,5))[0]**3)
    if which==2:derivative*=F(25,26)**2
    assert derivative>0
    return dict(theta_interval=list(map(str,(lo,hi))),theta_float=(float(lo)+float(hi))/2,F_at_lower=list(map(str,gate(which,lo))),F_at_upper=list(map(str,gate(which,hi))),monotone_derivative_lower=str(derivative),rL_always_below_r36=True,first_grid_m=int(hi*20)+1)

def translated_L(key):
    ll=np.array(key[:16],int).reshape(4,4);hh=np.array(key[16:],int).reshape(4,4);active=np.all(hh>=ll,axis=1);out=set()
    for shift in product(range(-1,2),repeat=4):
        l=ll.copy();h=hh.copy();l[active]+=shift;h[active]+=shift;out.add(tuple(map(int,np.r_[l.ravel(),h.ravel()])))
    return out

def ratio(qq,ph,ww,keys):
    b,selected=boxes(qq,ph,ww,keys);return volume(b)[0]/4,len(selected)

def main():
    st=time.monotonic();reg=json.loads((P/'round011c_preregister.json').read_text());source=json.loads((P/'round011_gain_certificate.json').read_text());lk=translated_L(source['new_key']);assert len(lk)==81
    thresholds=[threshold(k) for k in range(3)];rows=[];candidate=None;centers=qmc.Sobol(4,scramble=True,seed=20261025113).random_base2(8)
    for k in range(3):
        print('GATE',k,thresholds[k],flush=True)
        for m in range(1,21):
            t=F(m,20);x=path(k,t);pp,qq,ph,ww=decode(x);assert abs(sum(ww)-4)<1e-14
            old,nb=ratio(qq,ph,ww,BASE);withL,nl=ratio(qq,ph,ww,BASE|lk);extra=capture(qq,ph,ww,centers);full,nf=ratio(qq,ph,ww,BASE|lk|extra)
            row=dict(path=k,m=m,t=str(t),old_ratio=old,L_ratio=withL,full_ratio=full,old_union_signed_loss=OLD-old,L_gain=withL-old,event_gain=full-withL,total_signed_gain=full-OLD,hull_gate_sign=('positive' if gate(k,t)[0]>0 else 'negative' if gate(k,t)[1]<0 else 'unresolved'),selected_boxes=[nb,nl,nf],captured_keys=len(extra))
            rows.append(row)
            if m in [1,thresholds[k]['first_grid_m'],20]:print('POINT',row,flush=True)
            if full>OLD+1e-11:candidate=dict(row=row,x=x.tolist());break
        if candidate is not None:break
    rlo,rhi=radius(F(30));slo,shi=radius(F(36));out=dict(status='fixed analytic paths only; no input/family upper bound',model=reg['model'],L_counts=[12,6,12,0],L_mass='30-6*epsilon',hull_mass=36,baseline_r36_minus_r30_interval=[str(slo-rhi),str(shi-rlo)],thresholds=thresholds,rows=rows,candidate=candidate,stop='potential improvement requires exact priority' if candidate else '60 registered candidates complete; no new bound; no further optimization',elapsed=time.monotonic()-st)
    (P/'round011c_results.json').write_text(json.dumps(out,indent=2)+'\n');print('END',out['stop'],'elapsed',out['elapsed'],flush=True)
if __name__=='__main__':main()
