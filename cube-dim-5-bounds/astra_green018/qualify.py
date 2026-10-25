"""Finite DD4 checks and one new evaluation against a saved old full-G ball."""
from pathlib import Path
from fractions import Fraction as F
import json,time,os,sys
import sympy as sp
from flint import arb,ctx
from kernel import GreenKernel,dd4,ar,sha,TABLE
P=Path(__file__).resolve().parent
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-5/cube-dim-5-bounds/numerics/green_generic.json')
def exact_dd(nodes,values):
    values=list(values)
    for order in range(1,len(nodes)):
        values=[(values[i+1]-values[i])/(nodes[i+order]-nodes[i]) for i in range(len(nodes)-order)]
    return values[0]
def main():
    began=time.monotonic();ctx.prec=256
    nodes=list(map(F,[-3,-1,2,4,7]));t=F(2,3)
    lhs=exact_dd(nodes,[t/(t*t+b*b) for b in nodes]);re,im=F(1),F(0)
    for b in nodes:
        u,v=t/(t*t+b*b),b/(t*t+b*b);re,im=re*u-im*v,re*v+im*u
    assert lhs==re
    # Five coincident nonzero nodes exercise the highest Hermite derivative.
    for node in (-2,2):assert dd4((node,)*5).contains(ar(F(-1,4*node**4)))
    # An independent rational confluent interpolation solve checks a mixed pattern.
    specifications=[(-3,0),(-3,1),(1,0),(2,0),(2,1)]
    matrix=sp.Matrix([[sp.diff(sp.Symbol('z')**j,sp.Symbol('z'),order).subs(sp.Symbol('z'),node) for j in range(5)] for node,order in specifications])
    values=[arb(3).log(),ar(F(-1,3)),arb(0),arb(2).log(),ar(F(1,2))]
    inv=matrix.inv();direct=sum((ar(F(str(inv[4,j])))*values[j] for j in range(5)),arb(0))
    assert (dd4((-3,-3,1,2,2))-direct).contains(0)
    engine=GreenKernel();rejected=False
    try:engine.projections((1,1,1,1,0))
    except ValueError:rejected=True
    assert rejected
    old=json.loads(OLD.read_text());saved=old['points'][0]
    point=tuple(map(F,saved['point']));assert point==(F(1),F(936,1009),F(673,1009),F(68,1009),F(917,1009))
    result=engine.evaluate(point);newball=arb(result['G_ball']);oldball=arb(saved['G'])
    assert newball.overlaps(oldball) and (newball-oldball).contains(0)
    budget=json.loads((P.parent/'astra_upper014'/'budget.json').read_text())['cases'][0]
    assert F(result['G_upper_dyadic'])<F(budget['g_required_for_cost_32_using_exact_expression']['rational_lower'])
    out={'status':'new3504_kernel_qualification_pass','DD4_resolvent_derivative_sign_identity_exact':True,
         'pure_fivefold_nonzero_Hermite_checks':2,'mixed_Hermite_independent_interpolation_check':True,'zero_projection_guard_pass':True,
         'new_value':result,'saved_old_G_ball':saved['G'],'difference_from_saved_old':str(newball-oldball),
         'old26880_baseline_executed':False,'old_point_below_current014_required_g':True,
         'new_global_g':False,'new_C5_bound':False,
         'source_hashes':{str(x):sha(x) for x in [P/'registration.md',P/'kernel.py',Path(__file__),TABLE,OLD]},
         'pid':os.getpid(),'seconds':time.monotonic()-began,'python':sys.executable}
    (P/'qualification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},indent=2))
if __name__=='__main__':main()
