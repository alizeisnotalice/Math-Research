#!/usr/bin/env python3
"""Independent CDF polynomial interpolation and 224-bit integer-root audit."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import hashlib,json,os,time
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def value(p,x):return sum((F(c)*x**k for k,c in enumerate(p)),F(0))
def integ(p,l,r):return sum((F(c)*(r**(k+1)-l**(k+1))/F(k+1) for k,c in enumerate(p)),F(0))
def sqrt224(x):
    scale=1<<224;n=isqrt(x.numerator*scale*scale//x.denominator)
    lo,hi=F(n,scale),F(n+1,scale)
    assert lo*lo<=x<hi*hi
    return lo,hi
start=time.monotonic();print('verify111_pid',os.getpid(),flush=True)
r=read(HERE/'results.json');reg=read(HERE/'registration.json');cat=read(HERE/'ratio_catalog.json')
assert r['catalog_sha256']==sha(HERE/'ratio_catalog.json')
for path,h in reg['sources'].items():assert sha(path)==h,path
for row in r['new_ratio_points']:
    x=F(row['x']);a,b=sqrt224(1+x*x/16);c,d=sqrt224(1+x*x/256)
    lower=(a+x/4)**2*(c+x/16)**8;upper=(b+x/4)**2*(d+x/16)**8
    assert F(row['lower'])<=lower<upper<=F(row['upper'])
old103=read(next(p for p in reg['sources'] if p.endswith('upper_adaptive103/ratio_catalog.json')))
old091=read(next(p for p in reg['sources'] if p.endswith('upper_shift091/results.json')))
old091={x['x']:x['R_new'] for x in old091['ratios']}
for x,row in cat.items():
    if row['origin']=='new111':continue
    source={'lower':'1','upper':'1'} if x=='0' else old103[x] if row['origin']=='saved103' else old091[x]
    assert row['lower']==source['lower'] and row['upper']==source['upper']
plans=read(HERE/'geometry_plans.json');countpieces=0;countcells=0
for plan,case,saved in zip(plans,reg['finite_cases'],r['cases']):
    spans=[tuple(map(F,z)) for z in case['offsets']];pieces=plan['profile']['pieces'];a=F(plan['a'])
    knots=sorted({F(0)}|{abs(z) for ab in spans for z in ab})
    assert [(F(p['lo']),F(p['hi'])) for p in pieces]==list(zip(knots,knots[1:]))
    J0=F(0);previousW=None
    for p in pieces:
        l,u=F(p['lo']),F(p['hi']);cdf=list(map(F,p['cdf']));W=list(map(F,p['W']))
        assert len(cdf)<=5
        # Degree<=4, five distinct rational values certify the entire CDF
        # polynomial on this interval because every coordinate factor is affine.
        for j in range(1,6):
            x=l+(u-l)*F(j,6);actual=F(1)
            for lo,hi in spans:
                actual*=F(int(x>=abs(lo))) if lo==hi else max(F(0),min(hi,x)-max(lo,-x))/(hi-lo)
            assert value(cdf,x)==actual
        derivative=[F(k)*W[k] for k in range(1,len(W))]
        expected=cdf[:];expected[0]-=1
        assert derivative==expected
        if previousW is not None:assert previousW==value(W,l)
        previousW=value(W,u);assert value(W,l)>=0 and previousW>=0
        J0+=integ(W,l,u);countpieces+=1
    assert previousW in (None,F(0)) and J0==F(plan['profile']['J0'])
    weights=F(0);previous_right=F(0);low=high=F(0)
    for cell in plan['cells']:
        l,u,w=map(F,(cell['left'],cell['right'],cell['weight']))
        assert previous_right<=l<u and w>0;previous_right=u
        candidates=[p for p in pieces if F(p['lo'])<=l<u<=F(p['hi'])];assert len(candidates)==1
        assert w==integ(candidates[0]['W'],l,u)
        assert F(cell['x_left'])==l/a and F(cell['x_right'])==u/a
        low+=w*F(cat[cell['x_left']]['lower']);high+=w*F(cat[cell['x_right']]['upper'])
        weights+=w;countcells+=1
    assert weights==J0
    h=knots[-1];out=saved['result']
    if h:
        R=F(cat[str(h/a)]['upper']);old=J0*R
        assert F(out['J_lower'])==max(low,J0)
        assert F(out['J_upper'])==min(high,old)
        assert F(out['J_upper'])<old
        mu=max(abs((lo+hi)/2) for lo,hi in spans)
        assert F(out['factor_upper'])==min(R,1+mu/a+F(25,24)*F(out['J_upper'])/a**2)
    else:assert out['factor_upper']=='1'
v={'status':'passed','pid':os.getpid(),'wall_seconds':time.monotonic()-start,
   'new_x_checked_by_224bit_grid':len(r['new_ratio_points']),'old_ratio_values_recomputed':0,
   'CDF_polynomial_pieces':countpieces,'positive_weight_cells':countcells,
   'new_H_jet_log_calls':0,'results_sha256':sha(HERE/'results.json'),'code_sha256':sha(__file__)}
(HERE/'verification.json').write_text(json.dumps(v,sort_keys=True,indent=2)+'\n');print(json.dumps(v,sort_keys=True))
