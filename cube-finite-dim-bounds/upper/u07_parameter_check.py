"""Exact selected free-R guards, cube monotonicity and cost rejection.
Uses cb7's complete sphere enumeration and exact H6 upper ranges.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json

def sphere_min(n,b4,b6,r):
    d=b6*r**6;e=b4*r**4-Q(15,n+8)*d;c=-Q(3,n+2)*b4*r**4+Q(30,(n+8)*(n+4))*d
    values=[]
    for m in range(1,n+1):
        values.append(d/Q(m*m)+e/Q(m)+c)
        if not d:continue
        T=-2*e/(3*d)
        for p in range(1,m):
            q=m-p
            if p==q:continue # any stationary continuum has the uniform representative value
            u=(1-q*T)/(p-q);v=T-u
            if u>0 and v>0:values.append(d*(p*u**3+q*v**3)+e*(p*u*u+q*v*v)+c)
    return min(values),len(values)

rows=[]
for n,A,b4,b6,R,core in (
    (3,Q('.41061707799980696'),Q('.0052535036230922875'),Q('.0009893885508363638'),Q('3.1643099997673794'),Q(21,10)),
    (4,Q('.14344247730484339'),Q('.002778738033100787'),Q('.0003790027706453'),Q('3.3162574578337103'),Q(23,10))):
    assert min(A,b4,b6)>0 and R*R>n
    sphere,nvalues=sphere_min(n,b4,b6,R);eta=A-R**(2-n)+sphere;assert eta>0
    # p*n^(-n/2) >=2/21 (n3, sqrt3<7/4), or exactly1/16 (n4).
    radial_lower=Q(2,21) if n==3 else Q(1,16)
    gradient_margin=radial_lower-2*Q(3,n+2)*b4*n-3*Q(15,n+8)*b6*n
    assert gradient_margin>0
    p4=Q(n-1,n+2);p6=Q(32,693) if n==3 else Q(1,16)
    core_margin=core**(2-n)-A-b4*p4*core**4-b6*p6*core**6
    assert core<R and core*core>n and core_margin>0
    h4corner=Q(n)-Q(3,n+2)*n*n
    h6corner=Q(n)-Q(15,n+8)*n*n+Q(30,(n+8)*(n+4))*n**3
    inv_upper=Q(577351,1000000) if n==3 else Q(1,4)
    corner_upper=inv_upper-A-b4*h4corner-b6*h6corner;assert corner_upper>0
    pi_lower=Q(157,50)
    mass_lower=(4*pi_lower if n==3 else 2*pi_lower**2)*(core*core/2-A*core**n/n)
    lower=mass_lower/(2**n*corner_upper);assert lower>(8 if n==3 else 16)
    rows.append({'n':n,'A_exact':str(A),'B4_exact':str(b4),'B6_exact':str(b6),'R_exact':str(R),
                 'sphere_candidate_count':nvalues,'outer_negative_margin_exact':str(eta),'outer_negative_margin_display':float(eta),
                 'cube_coordinate_derivative_margin_exact':str(gradient_margin),'cube_coordinate_derivative_margin_display':float(gradient_margin),
                 'core_radius':str(core),'core_margin_exact':str(core_margin),'core_margin_display':float(core_margin),
                 'mass_lower_exact':str(mass_lower),'cube_min_upper_exact':str(corner_upper),
                 'selected_kernel_cost_lower_exact':str(lower),'selected_kernel_cost_lower_display':float(lower)})
print(json.dumps({'status':'exact arithmetic selected-kernel certificate, not C_n or whole-family bound','rows':rows,
                  'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
