"""Audit013 regression and repaired raw-branch guards; fixed cases only."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
import sympy as sy
from vg6_bounds import Endpoint as E,guarded_extreme

P=Path(__file__).resolve().parent


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def main():
    raw=[E(0,0),E(F(-1,64),1)]
    failed=False
    try:
        guarded_extreme(raw,F(1,16))
    except ValueError:
        failed=True
    assert failed
    bad_epsilon=F(1,32)
    bad_value=max(F(0),F(-1,64)+bad_epsilon)
    assert bad_value==F(1,64)>0
    repaired,guards=guarded_extreme(raw,F(1,256))
    assert repaired==E(0,0,0)
    assert guards==[dict(gap=F(1,64),movement_upper=F(1,256),strict_slack=F(3,256))]
    for t in [F(0),F(1,512),F(1,256)]:
        assert max(F(0),F(-1,64)+t)==0
    # Alternative full-cap repair: B=16 is exact, with equality at epsilon=1/32.
    x=sy.symbols('x')
    square=sy.expand(16*(x-sy.Rational(1,32))**2)
    assert sy.expand(16*x*x-(x-sy.Rational(1,64))-square)==0
    assert bad_value==16*bad_epsilon**2
    tied=[E(0,1,1),E(0,-1,2),E(-1,10,1)]
    tied_max,tied_guards=guarded_extreme(tied,F(1,64),'max')
    tied_min,min_guards=guarded_extreme([E(-e.base,-e.velocity,e.remainder) for e in tied],F(1,64),'min')
    assert tied_max==E(0,1,2) and tied_min==E(0,-1,2)
    assert tied_guards==min_guards
    for t in [F(0),F(1,128),F(1,64)]:
        values=[e.base+e.velocity*t+e.remainder*t*t for e in tied]
        assert abs(max(values)-t)<=2*t*t
        assert abs(min(-v for v in values)+t)<=2*t*t
    out=dict(status='audit013 raw-clipping regression and repaired guards PASS',
             prior_defect='1-Lipschitz clipping alone does not preserve a fixed baseline B through later branch changes',
             bad_cap=F(1,16),bad_epsilon=bad_epsilon,bad_value=bad_value,
             original_claimed_B=0,large_cap_rejected=failed,repaired_cap=F(1,256),
             repaired_endpoint=dict(base=repaired.base,velocity=repaired.velocity,B=repaired.remainder),
             guards=guards,full_cap_alternative_B=16,nonnegative_square=str(square),
             tied_max=dict(base=tied_max.base,velocity=tied_max.velocity,B=tied_max.remainder),
             tied_min=dict(base=tied_min.base,velocity=tied_min.velocity,B=tied_min.remainder),
             tied_guards=tied_guards,
             source_sha256={name:hashlib.sha256((P/name).read_bytes()).hexdigest()
                            for name in ['vg6_bounds.py','vg6_clipping_check.py']},
             scope='new fixed guard regression only; old 512/MATLAB/phase searches not rerun')
    (P/'vg6_clipping_results.json').write_text(json.dumps(encode(out),indent=2)+'\n')
    print(json.dumps(encode(out),indent=2))


if __name__=='__main__':main()
