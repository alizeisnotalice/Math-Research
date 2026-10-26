#!/usr/bin/env python3
"""New symbolic/certificate-interface toys only. No real H or old grid run."""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import sympy as s

OUT = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    assert not (OUT/'symbolic_results.json').exists()
    registration = {'registered_utc': datetime.now(timezone.utc).isoformat(), 'pid': os.getpid(),
        'code_sha256': sha(__file__), 'finite_new_identities': 3, 'finite_artificial_endpoint_toys': 1,
        'new_H_gradient_kernel_log_exp_calls': 0, 'real_endpoint_or_grid_calls': 0}
    (OUT/'symbolic_registration.json').write_text(json.dumps(registration,sort_keys=True,indent=2)+'\n')
    e0,e1,e2,m0,m1,m2 = s.symbols('e0 e1 e2 m0 m1 m2',positive=True)
    h0,h1,h2 = m0+m1+m2,m1+m2,m2
    step_bound = h0/e0**2-(e0**-2-e1**-2)*h1-(e1**-2-e2**-2)*h2
    assert s.simplify(step_bound-(m0/e0**2+m1/e1**2+m2/e2**2)) == 0
    L,U = s.symbols('L U',positive=True)
    q128 = (1-s.Rational(3,4)*L/U)/(16*s.Rational(1,128)**2)
    assert s.simplify(q128-(1024-768*L/U)) == 0
    norm2,base = s.symbols('norm2 base',nonnegative=True)
    assert s.simplify((base-q128*norm2)-(base-1024*norm2)-768*L/U*norm2) == 0
    # Artificial same-anchor interval: H128/g >=8/7, U128=140micro,
    # L64=105micro. log(8/7)>=2/15>1/8 by the elementary log lower bound.
    lower_mass = F(21,200000)
    upper_mass = F(7,50000)
    delta = F(1,8)
    q = F(1024)-F(768)*lower_mass/upper_mass
    assert q == 448
    radius2 = 4*F(1,128)**2
    old_margin = delta-1024*radius2
    new_margin = delta-q*radius2
    assert old_margin == -F(1,8) and new_margin == F(1,64)
    assert F(2,15) > delta and F(11,100000)*F(8,7) < upper_mass
    result = {'status': 'PASS_new_symbolic_identities_and_artificial_rational_toy',
        'pid': os.getpid(), 'registration_sha256': sha(OUT/'symbolic_registration.json'),
        'cutoff_step_telescoping': True, 'epsilon128_q_identity': True, 'margin_improvement_identity': True,
        'artificial_toy': {'L64':str(lower_mass),'U128':str(upper_mass),'saved_Delta':str(delta),
            'loggradient':'zero','box_offset':'[-1/128,1/128]^4','q':str(q),
            'old_margin':str(old_margin),'new_margin':str(new_margin),'actual_H_claim':False},
        'new_H_gradient_kernel_log_exp_calls':0,'real_endpoint_or_grid_calls':0}
    (OUT/'symbolic_results.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
