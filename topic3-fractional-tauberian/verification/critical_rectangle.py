#!/usr/bin/env python3
"""Exact scalar regressions and high-precision strict-phase checks.

The all-input peak theorem is proved analytically by endpoint reduction;
these finite samples check its formulas and explicit attaining crosses.
"""
from fractions import Fraction as F
from pathlib import Path
import random,json,platform


def squared_peak(w,d):
    if d>=1-w*w:return F(1)
    if d>=1-w:return (1-d)/(2*(1-d)-w*w)
    raise ValueError('use rational root parametrization in the first branch')


def exact_checks():
    rng=random.Random(20261105);counts=[0,0];min_slack=None;attainers=[]
    for branch in range(2):
        for k in range(2048):
            w=F(rng.randrange(1,64),64)
            if branch==0:
                ell0=w+(1-w)*F(rng.randrange(65),64);d=ell0*(ell0-w);p2=ell0*ell0/(1+d);target=ell0
            else:
                d=1-w+w*(1-w)*F(rng.randrange(65),64);p2=squared_peak(w,d);target=(1-d)/w
            W=1/w
            if k%3==0:ell=w*F(rng.randrange(1,65),64)
            elif k%3==1:ell=w+(W-w)*F(rng.randrange(65),64)
            else:ell=W*(1+F(rng.randrange(65),64))
            amax=min(ell*ell,w*ell,F(1));a=amax*F(rng.randrange(65),64)
            blo=max(F(0),d-1+a);bhi=min(d,ell*ell-a);assert blo<=bhi
            b=blo+(bhi-blo)*F(rng.randrange(65),64)
            value2=(a+b)**2/(ell*ell*(1-d+2*b));slack=p2-value2;assert slack>=0
            min_slack=slack if min_slack is None else min(min_slack,slack)
            # A concrete cross: add its whole vertical caps, trim equal outer x strips.
            a0=w*target;added=target*target-a0;removed=d-added
            assert 0<=removed<=1-a0 and added>=0
            horizontal_width=W-removed/w;mass=w*horizontal_width+target*target-w*target
            assert horizontal_width>=target>=w and mass==1-removed+added
            assert target*target/mass==p2 and removed+added==d
            counts[branch]+=1
            if k<3:attainers.append({'branch':branch+1,'w':str(w),'delta':str(d),'square_side':str(target),'trimmed_width':str(horizontal_width),'squared_peak':str(p2)})
    return {'samples_by_branch':counts,'minimum_exact_upper_slack':str(min_slack),'exact_attainment_checks':sum(counts),'attainer_examples':attainers}


def precision_checks():
    import mpmath as mp
    mp.mp.dps=90;rows=[]
    for A0 in ['0','0.25','1','4']:
        A=mp.mpf(A0)
        for n in [5,7,9,11]:
            t=mp.mpf(2)**(-n);alpha=1-t;d=A*t*t
            ws=alpha*mp.sqrt(1+d)-d/(alpha*mp.sqrt(1+d));pc=-mp.log(ws)/t
            values=[]
            for sign in [-1,0,1]:
                p=pc+sign*t**3;w=mp.exp(-p*t);ell=(w+mp.sqrt(w*w+4*d))/2
                gap=ell/mp.sqrt(1+d)-alpha
                if sign:assert gap*sign<0
                else:assert abs(gap)<mp.mpf('1e-85')
                assert d<=1-w
                values.append(mp.nstr(gap,40))
            cut=(-mp.log(alpha)-mp.log(1-d)/2)/t
            if A:assert pc>cut
            else:assert abs(pc-cut)<mp.mpf('1e-85')
            rows.append({'A':A0,'t':mp.nstr(t,30),'critical_p':mp.nstr(pc,45),'cut_p':mp.nstr(cut,45),'peak_minus_alpha_below_at_above':values})
    return {'mpmath':mp.__version__,'precision_digits':90,'parameter_rows':rows,'scope':'strict phase signs and equality at high precision are regressions; exact boundary is derived analytically'}


def main():
    import sys
    out={'python':platform.python_version(),'seed':20261105,'exact_arithmetic':'Fraction','exact':exact_checks()}
    if '--exact-only' not in sys.argv:out['high_precision']=precision_checks()
    target=Path(__file__).resolve().parents[1]/'results/critical_rectangle_python.json'
    if '--exact-only' not in sys.argv:target.write_text(json.dumps(out,indent=2)+'\n')
    print({'exact':out['exact'],'high_precision_rows':len(out.get('high_precision',{}).get('parameter_rows',[]))})
if __name__=='__main__':main()
