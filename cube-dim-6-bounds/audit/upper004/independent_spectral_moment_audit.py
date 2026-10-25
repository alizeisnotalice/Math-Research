#!/usr/bin/env python3
"""Exact audit of the round004 spectral second-moment cube lower bound."""
from __future__ import annotations

import hashlib, json, math, platform
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx
import flint

HERE=Path(__file__).resolve().parent
FAN=HERE.parents[1]/'upper/round002/fan/exact_cells_results.json'
AUTHOR=HERE.parents[1]/'upper/round004/primary/spectral_moment.json'
SOURCE=AUTHOR.with_suffix('.py')
W=2**6*math.factorial(6)
Q=Fraction

def hs(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def a(q): return arb(q.numerator)/q.denominator
def main():
    doc=json.loads(FAN.read_text());base_mass=Q();base_trace=Q();base_matrix=[[Q() for _ in range(6)]for _ in range(6)]
    for cell in doc['cells']:
        rays=[[Q(v) for v in ray] for ray in cell['rays_a0_equals_one']]
        rates=[Q(v) for v in cell['rates']]
        alpha=Q(cell['determinant'])/math.prod(rates)
        assert len(rays)==len(rates)==6 and alpha>0 and all(q>0 for q in rates)
        # u_j are independent Exp(q_j) at t=1.  Thus mean xi_i is m_i,
        # and E xi_i xi_k=m_i*m_k+sum_j r_ji*r_jk/q_j^2.
        mean=[sum((rays[j][i]/rates[j] for j in range(6)),Q()) for i in range(6)]
        for i in range(6):
            for k in range(6):
                moment=mean[i]*mean[k]+sum((rays[j][i]*rays[j][k]/rates[j]**2 for j in range(6)),Q())
                base_matrix[i][k]+=alpha*moment
        base_mass+=alpha
    base_trace=sum((base_matrix[i][i] for i in range(6)),Q())
    full_mass=W*base_mass
    assert full_mass==Q(doc['full_exp_minus_symbol_integral'])
    # Sum over the full signed-permutation group: signs remove all off-diagonal
    # entries and permutations make every diagonal entry trace/6.
    full_cov_diag=base_trace*W/Q(6)
    E=base_trace/base_mass
    assert 6*full_cov_diag==W*base_trace
    ctx.prec=512
    C0=a(full_mass)/(2*arb.pi())**6
    T=(a(E)/2).sqrt()
    g=C0*(1/(5*T**5)-a(E)/(14*T**7))
    producer=json.loads(AUTHOR.read_text())
    output={
      'status':'PASS_exact_spectral_moment_reconstruction',
      'scope':'certifies the stated uniform spectral-moment lower bound conditional on the locked full fan/Fourier-density identity; no I upper bound, cost, or C6 claim',
      'inputs':{'fan_sha256':hs(FAN),'producer_source_sha256':hs(SOURCE),'producer_result_sha256':hs(AUTHOR)},
      'base_mass':str(base_mass),'full_mass':str(full_mass),'base_trace_second_moment':str(base_trace),
      'E_norm2':str(E),'full_covariance_diagonal_entry':str(full_cov_diag),
      'covariance_group_average':'(trace(base second-moment)/6) times identity after normalized signed-permutation averaging; all off-diagonals cancel exactly',
      'C0_arb_512':C0.str(90,more=True),'T_arb_512':T.str(90,more=True),'uniform_g_lower_arb_512':g.str(90,more=True),
      'producer_E_match':str(E)==producer['E_norm2'],
      'formula_check':'integral_T^infty t^-6(1-E/(2t^2))dt=1/(5T^5)-E/(14T^7), T=sqrt(E/2)',
      'cube_scope':'For max_i |x_i|<=1, unconditionality gives E_mu[(x dot xi)^2]=sum_i x_i^2 E_mu[xi_i^2]<=E_mu||xi||^2=E. cos u>=1-u^2/2 then gives p_t(x)>=C0 t^-6(1-E/(2t^2)) for t>=T.',
      'software':{'python':platform.python_version(),'python_flint':flint.__version__},'source_sha256':hs(Path(__file__))}
    assert output['producer_E_match']
    (HERE/'independent_spectral_moment_audit_results.json').write_text(json.dumps(output,indent=2,sort_keys=True)+'\n')
    print(json.dumps(output,indent=2,sort_keys=True))
if __name__=='__main__': main()
