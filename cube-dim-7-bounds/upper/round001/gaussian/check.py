#!/usr/bin/env python3
"""Deterministic symbolic/Arb checks, no new shape draws."""
import hashlib,json,math
from pathlib import Path
from fractions import Fraction
import sympy as s
from flint import arb,ctx
ROOT=Path(__file__).resolve().parent

def main():
    ctx.prec=256
    pi=arb.pi(); q=arb(7)/6
    C=1/(2*pi**(arb(7)/2)); H=16*pi**3/15
    D=2*(arb(65)/2).gamma()/arb(32).gamma()
    mg=arb(1)/4*((arb(65)/128).gamma()/(arb(1)/2).gamma())**64
    K=1/(arb(105)*arb(2)**(arb(25)/6)*pi**(arb(13)/12))
    iu=K*mg**q
    gneed=7*K**(1/q)*mg/(arb(128)**(1/q))
    n=s.Integer(7); cq=s.gamma(3)/(4*s.pi**(n/2)); hq=2*s.pi**(n/2)/s.gamma(n/2)
    symbolic_ratio=s.simplify(cq**s.Rational(7,6)*hq/896*(105*2**s.Rational(25,6)*s.pi**s.Rational(13,12)))
    assert symbolic_ratio==1
    V=s.Matrix([[1]+[1 if (b>>k)&1 else -1 for k in range(6)] for b in range(64)]).T
    assert V*V.T==64*s.eye(7)
    anchors=[]
    for k in range(7):
        a=s.Matrix([1]+[1]*k+[0]*(6-k))
        v=s.Matrix([1]*7); w=s.Matrix([1]+[1]*k+[-1]*(6-k))
        assert a==(v+w)/2
        indices=[i for i in range(64) if V[:,i]==v or V[:,i]==w]
        remain=V[:,[i for i in range(64) if i not in indices]]
        assert remain.rank()==7
        anchors.append(dict(k=k,body_span_rank=s.Matrix.hstack(v,w).rank(),remainder_rank=7,potential='infinity by positive convolution lower bound'))
    generic=[]
    for j in range(8):
        x=[Fraction(1)]+[Fraction(j+k,j+7) for k in range(1,7)]
        assert len({abs(a) for a in x})==7
        generic.append(dict(j=j,coordinates=[str(a) for a in x],distinct_absolute_values=7,outside_S3=True))
    # Fixed exact-shape checks are not new random samples.
    assert s.det(sum((64*V[:,i]*V[:,i].T for i in range(64)),s.zeros(7)))==s.Integer(4096)**7
    intervals=dict(C=str(C),H=str(H),D=str(D),M_lower=str(1/(32*pi.sqrt())),M_Gamma=str(mg),M_simple_upper=str(D/64),I_coefficient=str(K),I_Gamma_upper=str(iu),g_required_for_U128_if_using_M_Gamma=str(gneed),coefficient_ratio_symbolic=str(symbolic_ratio))
    out=dict(status='pass',arb_bits=ctx.prec,certificate_scope='rigorous analytic M and complete I upper enclosure only; no global g and no C7 bound',intervals=intervals,anchors=anchors,generic_points=generic)
    if (ROOT/'results.json').exists():
        import numpy as np
        r=json.loads((ROOT/'results.json').read_text()); a=np.load(ROOT/'observables.npz')
        assert r['status']=='completed' and r['shape_count']==8192
        assert a['Z'].shape==(8192,) and a['T'].shape==(8192,8)
        assert np.all(np.isfinite(a['Z'])) and np.all(np.isfinite(a['T']))
        assert np.all((a['Z']>0)&(a['Z']<=1/64)) and np.all(a['T']>0)
        assert len(r['blocks'])==8 and all(b['samples']==1024 for b in r['blocks'])
        assert r['diagnostics']['clipped']==r['diagnostics']['dropped']==r['diagnostics']['resampled']==0
        assert abs(float(a['Z'].mean())-r['Z']['mean'])<1e-18
        out['actual_results_check']='pass: 8192 retained rows, eight blocks, no clipping/dropping/resampling'
    (ROOT/'check-results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
