"""Minimal exact transcription identities; no search or prior verifier reuse."""
import hashlib
import json
import os
from pathlib import Path
import sympy as S

HERE = Path(__file__).resolve().parent
q, h, eps = S.symbols('q h eps', positive=True)
s = S.Rational(3, 2)
R = s*q-h
a0 = (1-eps)*q
g = R-a0
a = 1-g/(2*R)
ell = (1+a0/R)/2
ca, E, T, mass, Mk = S.symbols('ca E T mass Mk', positive=True)
A, B = S.symbols('A B', positive=True)
x, L, p, arate = S.symbols('x L p arate', positive=True)
tau = S.sqrt(A*B)
sigma, chi = S.Rational(25,256), S.Rational(63,1280)
qmax = S.Rational(11569,12800)
eta = (1-qmax)/(2*qmax)

checks = {
    'relative_guard_identity': S.simplify(g-R/3-(eps*q-2*h/3)) == 0,
    'finite_h_mean_identity': S.simplify(a-S.Rational(5,6)-h/(3*R)+eps*q/(2*R)) == 0,
    'negative_pin_derivative': S.simplify(S.diff(ell,q)+(1-eps)*h/(2*R**2)) == 0,
    'strong_log_transport': 3*S.Rational(3,2)*S.Rational(4,3) == 6,
    'strong_J1_transport': S.Rational(1,2)*S.Rational(3,2)*S.Rational(4,3) == 1,
    'strong_J2_transport': S.Rational(1,10)*S.Rational(3,2)*S.Rational(4,3) == S.Rational(1,5),
    'distinct_source_cap_rate': (eta-chi)*(1-sigma) == S.Rational(13627383,3790929920),
    'distinct_tau_rate': sigma/2-chi*(1-sigma) == S.Rational(1447,327680),
    'finite_a_energy_remainder': S.simplify((ca**2*E-T*ca*mass)/(ca*Mk)-(E-T*mass)/Mk-(ca-1)*E/Mk) == 0,
    'complete_root_identity': S.simplify((1+A/tau)*(tau+B)-(S.sqrt(A)+S.sqrt(B))**2) == 0,
    'Cantelli_monotonic_derivative': S.simplify(S.diff((arate*x-p)**2/x,x)-(arate**2-p**2/x**2)) == 0,
    'complete_price_slack': S.Rational(291,1000)-S.Rational(257,2500)-S.Rational(47,250) == S.Rational(1,5000),
}
assert all(checks.values()), checks
inputs = {
    'R073_proof': Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round073/cantelli_frontier/PROOF.md'),
    'R073_review': HERE.parents[1]/'round073/cantelli_frontier/review.md',
    'relative_J': Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md'),
    'J1': Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round066/first_profile/LEMMA.md'),
    'J2': Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile/PROOF.md'),
}
receipt = {
    'stage': 'pre-freeze transcription review only',
    'pid': os.getpid(), 'sympy': S.__version__, 'checks': checks,
    'sources': {key: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for key,path in inputs.items()},
    'scope': 'Twelve exact symbolic/rational adaptation identities only. No old certificate, finite-gate batch, parameter search, MATLAB, or manuscript acceptance.',
}
(HERE/'adaptation_identity_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
print(json.dumps({'pid': receipt['pid'], 'sympy': receipt['sympy'], 'checks': checks}, indent=2))
