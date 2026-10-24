#!/usr/bin/env python3
"""Exact auxiliary checks only; no LP, sampling, quadrature or certification.

Prove a common small-radius bound using the reviewed Fraction log core, and
check the homogeneity identities underlying the unscaled two-slice formula.
"""
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
REVIEWED_CORE_SHA256 = 'ed7a2460b2b71333072093843aa5818b822d11d07bbf7d6ba70b841d1c0485e2'


def find_core(explicit):
    if explicit is not None:
        return explicit.resolve()
    for path in [ROOT.parent/'round008/certify.py',
                 ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments/round008/certify.py']:
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError('Pass --core <reviewed round008/certify.py>.')


def verify(core_path):
    if sha256(core_path.read_bytes()).hexdigest() != REVIEWED_CORE_SHA256:
        raise ValueError('reviewed Fraction/log kernel hash mismatch')
    spec = importlib.util.spec_from_file_location('round016_auxiliary_fraction_log', core_path)
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    harmonic7 = sum((Q(1, k) for k in range(1, 8)), Q(0))
    delta = Q(1, 100)
    bound = (2*(1+delta)**6*harmonic7+core.log_point(Q(1, 101))
             +(1-delta)**7*core.log_point(Q(1, 99)))
    if not bound.hi < 0:
        raise ValueError('common small-radius upper bound is not negative')
    # For n=1..7, k=0..n, this proves the two polynomial homogeneity
    # identities; k=n also proves all coefficients multiplying logarithms.
    a, x = sp.symbols('a x', positive=True)
    homogeneity_checks = 0
    for n in range(1, 8):
        for k in range(n+1):
            for sign in (1, -1):
                residual = sp.expand(a**n*(1+sign*x/a)**k-a**(n-k)*(a+sign*x)**k)
                if residual != 0:
                    raise ValueError('unscaled homogeneity identity failed')
                homogeneity_checks += 1
    u, v = sp.symbols('u v', positive=True)
    aa, bb, d, s = 1-v, 1-u, u-v, 1-u-v
    substitutions = [aa+u-(1+d), aa-u-s, bb+v-(1-d), bb-v-s]
    if any(sp.expand(value) != 0 for value in substitutions):
        raise ValueError('unscaled coordinate substitutions failed')
    return {
        'status': 'exact_small_radius_and_unscaled_identity_checks_passed',
        'new_LPs': 0, 'new_generator_points': 0, 'new_quadrature_points': 0,
        'small_radius': {
            'powers': list(range(2, 9)), 'slice_range': '0<t<=1/100',
            'spatial_range': '0<r<=1/100 in the open quadrants',
            'harmonic7': str(harmonic7),
            'formula': 'B=2*(101/100)^6*Harm7+log(1/101)+(99/100)^7*log(1/99)',
            'B_interval': bound.pair(), 'B_upper_display': float(bound.hi),
            'strictly_negative': True,
            'conclusion': 'h_p(t)/(2p)<=B<0 for p=2..8; positive slice scaling preserves the sign',
            'analytic_assumptions': 'All log terms are negative; polynomial terms use upper bounds and log coefficients use lower bounds. t=u/(1-v)<=r and its symmetric counterpart.'
        },
        'unscaled_identity': {
            'homogeneity_identities_exactly_zero': homogeneity_checks,
            'coordinate_substitutions_exactly_zero': 4,
            'domain': 'u>0,v>0,u+v<1; n=p-1>=1',
            'log_step': 'In this positive domain log((u/a)/(1+u/a))=log(u/(a+u)) and the minus/symmetric identities hold exactly.',
            'scope': 'Verifies the algebraic substitution underlying continuous_interface.md; endpoint inclusion still requires the stated one-sided interval rules.'
        },
        'core_sha256': REVIEWED_CORE_SHA256,
        'interface_sha256': sha256((ROOT/'continuous_interface.md').read_bytes()).hexdigest(),
        'verifier_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Auxiliary exact lemmas only; no candidate-wide certificate or new general upper bound.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--core', type=Path)
    args = parser.parse_args()
    result = verify(find_core(args.core))
    (ROOT/'auxiliary_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
