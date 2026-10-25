#!/usr/bin/env python3
"""Exact/MP audit of the alpha=3/2 stable short-time face interface."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
from itertools import permutations
from math import prod
import mpmath as mp

mp.mp.dps=80
alpha=mp.mpf(3)/2
# Levy density c |z|^(-1-alpha) for characteristic exponent |xi|^alpha.
c=mp.gamma(alpha+1)*mp.sin(mp.pi*alpha/2)/mp.pi
# Big jumps above x plus Chebyshev for the compensated small-jump process.
tail_constant=4*c/(alpha*(2-alpha))
# Symmetric unimodality: P(|X|>=a*x)>=2(1-a)x*s_t(x), optimized at a=alpha/(alpha+1).
K=2*c*(alpha+1)**(alpha+1)/(alpha**(alpha+1)*(2-alpha))
A=mp.gamma(1/alpha)/(mp.pi*alpha)  # s_1(0)
# Conservative quotient certificate: the four projected singleton directions
# -2 e_i have determinant 16.  For each choice of the one large variable,
# use this quotient density; summing the 16 choices cancels that determinant.
def det_int(rows):
    """Tiny exact determinant audit; rows/columns have the same determinant."""
    n=len(rows)
    return sum(
        (-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
        * prod(rows[i][p[i]] for i in range(n))
        for p in permutations(range(n))
    )

# After a signed coordinate reflection any selected body direction is v.
# The following four listed body directions are a complement.  Subtracting v
# from them and expanding along coordinate one identifies the exact quotient
# determinant with the displayed four-by-four determinant.
v=((1,1,1,1,1),)
complement=((1,-1,-1,-1,1),(1,-1,-1,1,-1),(1,-1,1,-1,-1),(1,1,-1,-1,-1))
full_det=det_int(v+complement)
quotient_rows=tuple(tuple(w[j]-v[0][j] for j in range(1,5)) for w in complement)
quotient_det=det_int(quotient_rows)
assert full_det == quotient_det == -48

# The exact 48 determinant gives a sharper face certificate.  The quotient
# singleton certificate (determinant 16) is retained below as an audit only.
Cface=(mp.mpf(1)/3)*K*16**(1+alpha)*A**4
Cquotient=3*Cface
q=mp.mpf(10)/7
short_coefficient=16**(-1/q)*(mp.mpf(5)/6)*Cface**(1-1/q)*(16*K)**(1/q)
epsilons=[mp.mpf(1)/4096,mp.mpf(1)/16384,mp.mpf(1)/65536]
def fmt(x):return mp.nstr(x,40)
payload={
 "alpha":"3/2","one_dimensional_characteristic_function":"exp(-t*|xi|^(3/2))",
 "levy_density_constant":"c=Gamma(5/2)*sin(3*pi/4)/pi=3*sqrt(2)/(8*sqrt(pi))","c_decimal":fmt(c),
 "tail_decomposition":"P(|Z_t|>x) <= 4*c*t/[alpha*(2-alpha)]*x^(-alpha)=2*sqrt(2/pi)*t*x^(-3/2)","tail_constant_decimal":fmt(tail_constant),
 "density_tail_constant":"K=2*c*(alpha+1)^(alpha+1)/[alpha^(alpha+1)*(2-alpha)]=(25/6)*sqrt(10/(3*pi))", "K_decimal":fmt(K),
 "density_peak":"A=s_1(0)=2*Gamma(2/3)/(3*pi)","A_decimal":fmt(A),
 "explicit_complement":{"v":v[0],"other_directions":complement,"full_5_by_5_det":str(int(full_det)),"quotient_rows":quotient_rows,"quotient_4_by_4_det":str(int(quotient_det))},
 "face_pointwise":"p_t(1,y)<=Cf*t^(1-4/alpha)=Cf*t^(-5/3)","Cf_formula":"(1/3)*K*16^(5/2)*A^4","Cf_decimal":fmt(Cface),
 "singleton_quotient_crosscheck":"The four quotient singleton directions have determinant 16 and give 3*Cf, hence a weaker independent certificate.","Cquotient_decimal":fmt(Cquotient),
 "marginal":"integral_R4 p_t(1,y)dy=s_(16t)(1)<=16*K*t", "marginal_constant_decimal":fmt(16*K),
 "q":"10/7","interpolation":"||p_t(1,.)||_q <= Cf^(3/10)*(16K)^(7/10)*t^(1/5)",
 "positive_face_shorttime":"||int_0^eps p_t(1,.)dt||_L^q([0,infty)^4) <= B*eps^(6/5)","B_formula":"16^(-7/10)*(5/6)*Cf^(3/10)*(16K)^(7/10)","B_decimal":fmt(short_coefficient),
 "finite_epsilon_error_budget":[{"epsilon":"1/"+str(int(1/e)),"bound_B_epsilon_6_over_5":fmt(short_coefficient*e**(mp.mpf(6)/5))} for e in epsilons],
 "scope":"All displayed short-time inequalities follow from the stated independent 16-direction stable model and the global one-dimensional K bound. They do not prove the alpha=3/2 rank-stratum assertions, finite-source transfer, g/I, or a C5 improvement."
}
raw=json.dumps(payload,sort_keys=True,indent=2).encode()+b'\n';Path(__file__).with_name('round009_alpha_shorttime_results.json').write_bytes(raw);print(raw.decode(),end='');print('sha256',hashlib.sha256(raw).hexdigest())
