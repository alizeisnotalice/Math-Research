"""Exact admissibility of two selected rational H4 parameters; NOT a cost bound."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json

rows=[]
for n,A,B,lo,hi,target,power in (
    (3,Q('0.39660079267158965'),Q(3,1000),Q(31547,10000),Q(31549,10000),Q(625,2),5),
    (4,Q('0.1289749174365545'),Q(1,1000),Q(35495,10000),Q(35497,10000),Q(2000),6)):
    assert lo**power<target<hi**power
    assert lo*lo>n
    amin_upper=Q(n+2,4)*lo**(2-n)
    outer_margin=A-amin_upper
    # sqrt(3)<7/4; fourth-dimensional corner expression is rational.
    cube_lower=Q(4,7)+Q(12,5)*B-A if n==3 else Q(1,4)+4*B-A
    assert outer_margin>0 and cube_lower>0
    core=Q(11,5) if n==3 else Q(5,2)
    core_margin=core**(2-n)-A-B*Q(n-1,n+2)*core**4
    assert core<lo and core*core>n and core_margin>0
    pi_lower=Q(157,50)  # standard strict pi>3.14
    radial_core=core**2/2-A*core**n/n
    inner_mass_lower=(4*pi_lower if n==3 else 2*pi_lower**2)*radial_core
    reciprocal_root3_upper=Q(577351,1000000)
    assert 3*reciprocal_root3_upper**2>1
    cube_upper=reciprocal_root3_upper+Q(12,5)*B-A if n==3 else Q(1,4)+4*B-A
    cost_lower=inner_mass_lower/(2**n*cube_upper)
    assert cost_lower>(8 if n==3 else 16)
    rows.append({'n':n,'A_exact':str(A),'B_exact':str(B),
                 'R_polynomial':f'R^{power}={target}',
                 'R_interval':[str(lo),str(hi)],
                 'outer_negative_margin_lower_exact':str(outer_margin),
                 'outer_negative_margin_lower_decimal':float(outer_margin),
                 'cube_positive_margin_lower_exact':str(cube_lower),
                 'cube_positive_margin_lower_decimal':float(cube_lower),
                 'positive_core_radius':str(core),'core_positive_margin_exact':str(core_margin),
                 'mass_lower_from_core_exact':str(inner_mass_lower),
                 'cube_min_upper_exact':str(cube_upper),
                 'specific_kernel_cost_lower_exact':str(cost_lower),
                 'specific_kernel_cost_lower_decimal':float(cost_lower),
                 'scope':'exact admissibility and selected-kernel cost > target; NOT a C_n lower bound or whole-family obstacle'})
print(json.dumps({'rows':rows,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
