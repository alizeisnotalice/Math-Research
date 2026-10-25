#!/usr/bin/env python3
"""Arb certificates on an auxiliary METHOD BUDGET, not on the actual I."""
import json
import os
import time
import sympy as sp
from flint import arb, ctx
from rank2_angle import ROOT, TRIPLES, geometry, digest, dump

START = time.monotonic()
ctx.prec=256

def a_rational(q):
    q=sp.Rational(q)
    return arb(int(q.p))/int(q.q)

def main():
    reg=ROOT/'certificate_registration.json'
    assert json.loads(reg.read_text())['status']=='registered_before_computation'
    beta=arb(3)/4
    b=arb(3)/20
    rho=arb(247)/256
    s=b*rho
    r=2*s
    theta=[arb(3)/16,arb(3)/8,arb(7)/16]
    moment=lambda t:(1+t/beta).gamma()/(1+t).gamma()
    rows=[]
    lowers=[arb(1),arb(1)]
    for triple,weight in zip(TRIPLES,theta):
        W,gram=geometry(triple)
        C=arb(0)
        h_real_numerator=arb(0)
        h_imag_exact=sp.Integer(0)
        for i in range(16):
            x,y=W.row(i)
            norm=sp.simplify(x*x+y*y)
            if norm==0:
                continue
            scale=a_rational(norm/16)**beta
            C+=scale
            h_real_numerator+=scale*a_rational(sp.simplify((x*x-y*y)/norm))
            h_imag_exact+=sp.simplify((norm/16)**sp.Rational(3,4)*2*x*y/norm)
        assert sp.simplify(h_imag_exact)==0
        h_squared=(h_real_numerator/C)**2
        mu1=((1-beta)+(r+beta)*h_squared)/(r+1)
        prefactor=arb(2)**r*moment(r)*C**(-r/beta)
        partials=[prefactor,prefactor*(1+s*mu1)]
        for k in (0,1):
            lowers[k]*=partials[k]**(weight/rho)
        rows.append({'triple':triple,'gram_determinant':gram,'C':str(C),
                    'h_squared':str(h_squared),'mu1':str(mu1),
                    'negative_moment_partial_sums':[str(v) for v in partials]})
    factor=2*(arb(7)/4).gamma()*moment(s)**(3/rho)
    c=3/(arb(2)**(arb(15)/2)*arb.pi()**(arb(5)/2))
    js=[factor*x for x in lowers]
    is_=[c**(arb(10)/7)*arb.pi()**2/60*j**(arb(10)/7) for j in js]
    threshold=arb(118)/10_000_000
    checks=[bool(x>threshold) for x in is_]
    assert all(checks)
    result={'status':'certified_lower_bounds_on_default_rank2_Holder_budget_only',
            'script':digest(__file__),'registration':digest(reg),'source_geometry':digest(ROOT/'rank2_angle.py'),
            'parameters':{'b':'3/20','rho':'247/256','s':'741/5120','r':'741/2560'},
            'representatives':rows,'budget_J_floors':[str(x) for x in js],
            'budget_I_floors':[str(x) for x in is_],
            'floor_strictly_above_118e_minus7':checks,
            'not_a_lower_bound_on_actual_J_or_I':True,'new_C5_upper_bound':False,
            'runtime':{'seconds':time.monotonic()-START,'pid':os.getpid(),'precision_bits':ctx.prec}}
    dump('budget_floor.json',result)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    main()
