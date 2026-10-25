"""Exact selected-kernel failure certificates; not C_n or whole-family bounds.

Uses complete H6 sphere ranges proved/enumerated by cb7: n3 max=32/693,
n4 max=1/16; H4 max=(n-1)/(n+2), both harmonics have sphere mean zero.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json

rows=[]
for n,A,B4,B6,R in (
    (3,Q('.4089479272201547'),Q('.005004477192161995'),Q('.0009386738283819751'),Q(16,5)),
    (4,Q('.18157763153667977'),Q('.004'),Q('.0005'),Q(18,5))):
    core=Q(21,10);p4=Q(n-1,n+2);p6=Q(32,693) if n==3 else Q(1,16)
    assert B4>0 and B6>0 and core<R and core*core>n
    core_margin=core**(2-n)-A-B4*p4*core**4-B6*p6*core**6
    assert core_margin>0
    # The radial bound is decreasing, so the entire core ball is positive.
    h4corner=Q(n)-Q(3,n+2)*n*n
    h6corner=Q(n)-Q(15,n+8)*n*n+Q(30,(n+8)*(n+4))*n**3
    inv_upper=Q(577351,1000000) if n==3 else Q(1,4)
    if n==3:assert 3*inv_upper**2>1
    corner_upper=inv_upper-A-B4*h4corner-B6*h6corner
    pi_lower=Q(157,50)
    mass_lower=(4*pi_lower if n==3 else 2*pi_lower**2)*(core*core/2-A*core**n/n)
    cost_lower=mass_lower/(2**n*corner_upper)
    assert corner_upper>0 and cost_lower>(8 if n==3 else 16)
    rows.append({'n':n,'A_exact':str(A),'B4_exact':str(B4),'B6_exact':str(B6),'R_exact':str(R),
                 'positive_core_radius':str(core),'core_margin_lower_exact':str(core_margin),
                 'core_margin_lower_display':float(core_margin),'mass_lower_exact':str(mass_lower),
                 'cube_min_upper_exact':str(corner_upper),'specific_kernel_cost_lower_exact':str(cost_lower),
                 'specific_kernel_cost_lower_display':float(cost_lower),
                 'scope':'rejects only this normalized kernel; no C_n lower bound and no whole-family obstacle'})
print(json.dumps({'rows':rows,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
