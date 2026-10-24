#!/usr/bin/env python3
"""Reproduce the single analytic corner-cap constant; no box search or LP."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse
import importlib.util
import json
import sys
import sympy as sp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
CORE_HASH = 'ed7a2460b2b71333072093843aa5818b822d11d07bbf7d6ba70b841d1c0485e2'
EXPECTED = Q('2399675517088728818701176616244058049538048286599298645832413152820217373259687257513350102964429/44514041409384860089196449672854782875219825655221939086914062500000000000000000000000000000000000')


def core_path(explicit):
    options = [explicit] if explicit else [ROOT.parent/'round008/certify.py', ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments/round008/certify.py']
    for path in options:
        if path.is_file() and sha256(path.read_bytes()).hexdigest() == CORE_HASH:
            return path
    raise FileNotFoundError('Pass --core <reviewed round008/certify.py>.')


def verify(source):
    spec = importlib.util.spec_from_file_location('round017_cap_repro_log', source)
    c = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = c
    spec.loader.exec_module(c)
    R=Q(15,8)
    E={2:Q(69,100),6:Q(19,50),7:Q(11,10),8:Q(3,4)}
    eps=Q(2,199); tau=1-eps; bmax=Q(299,20000); rlo=Q(99,100)
    log2=c.log_point(2)
    H={2:8-8*log2,6:Q(1334,5)-384*log2,7:Q(3108,5)-896*log2,8:Q(149104,105)-2048*log2}
    omega={2:4*eps*(c.log_point(199).hi-1)}
    M={2:Q(4)}
    L={}
    for p in (6,7,8):
        m=p-2; n=p-1
        K=1/(tau*(p-1))+sum((eps**k/Q(m-k) for k in range(m)),Q(0))
        K+=eps**m*c.log_point(1/eps).hi
        L[p]=max(Q(0),2*p*(p-1)*K-4*p)
        omega[p]=eps*L[p]
        M[p]=2*p*(sum((Q(2**k+1,n-k) for k in range(n)),Q(0))+Q(1,2*n))
    U={p:H[p].hi+omega[p]+bmax**(p-1)*M[p] for p in E}
    F=sum((rlo**j/((j+2)*R**(j+1)) for j in range(29)),Q(0))
    tail=2*(2*F+1/R)/(R-1)
    margin=tail-sum((E[p]*U[p] for p in E),Q(0))
    if not (margin > 0 and margin == EXPECTED):
        raise ValueError('fixed cap constant did not exactly reproduce')
    t=sp.symbols('t',positive=True)
    def jp(m):
        return sum((1+t)**k/sp.Integer(m-k) for k in range(m))+(1+t)**m*(sp.log(t)-sp.log(1+t))
    def jm(m):
        return sum((1-t)**k/sp.Integer(m-k) for k in range(m))+(1-t)**m*(sp.log(t)-sp.log(1-t))
    for p in (6,7,8):
        h=2*p*(jp(p-1)+jm(p-1))
        rhs=2*p*(p-1)*(jp(p-2)-jm(p-2))+4*p/t
        if sp.cancel(sp.expand(sp.diff(h,t)-rhs)) != 0:
            raise ValueError('derivative identity failed')
    return {'status':'one_fixed_support_corner_cap_strictly_positive',
            'execution':'portable_verifier', 'candidate_R':'15/8', 'cap_r_q':['99/100','1'],
            'epsilon_t':str(eps),'b_max':str(bmax),
            'derivative_L_upper_display':{str(p):float(x) for p,x in L.items()},
            'global_h_upper_M':{str(p):str(x) for p,x in M.items()},
            'modulus_upper_display':{str(p):float(x) for p,x in omega.items()},
            'Apsi_upper_display':{str(p):float(x) for p,x in U.items()},
            'tail_lower_display':float(tail),
            'strict_cap_generator_margin_fraction':str(margin),
            'strict_cap_generator_margin_display':float(margin),
            'new_LPs':0,'new_generator_points':0,'new_adaptive_runs':0,
            'new_symbolic_derivative_identities':3,'core_sha256':CORE_HASH,
            'proof_sha256':sha256((ROOT/'boundary_cap.md').read_bytes()).hexdigest(),
            'verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope':'One fixed analytic corner cap only; not a global certificate or new general upper bound.'}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--core',type=Path)
    args=parser.parse_args()
    result=verify(core_path(args.core))
    (ROOT/'boundary_cap_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
