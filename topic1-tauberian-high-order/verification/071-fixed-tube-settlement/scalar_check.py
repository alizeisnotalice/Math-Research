"""R071 narrow independent Fraction check. No geometry is encoded."""
from fractions import Fraction as F
import json

a=F(20,81)
gradient_first=-a+2*F(10,81)
assert gradient_first==0
# Along r=1/3+a*v, Phi's v^2 coefficient as polynomial in n.
peak_v2=[-a*a+F(40,81)*a-F(10,243)-F(13,486),F(2,81),-F(1,81)]
closed=[-F(91,13122),F(324,13122),-F(162,13122)]
qnorm_v2=[-F(10,243)-F(13,486)+F(400,6561),F(2,81),-F(1,81)]
assert peak_v2==closed==qnorm_v2
# Exact Hessian action on the symmetric and antisymmetric directions.
H=[[F(-2,3),F(-1,3)],[F(-1,3),F(-2,3)]]
for vec,lam in [([1,1],F(-1)),([1,-1],-F(1,3))]:
    assert [sum(H[i][j]*vec[j] for j in range(2)) for i in range(2)]==[lam*x for x in vec]
high=[]
for q in range(3,25):
    beta=F(q*(2*q-1),q-1)
    long_gap=beta+2-F(19,2)
    boundary_gap=q+F(q,q-1)+5-F(19,2)
    assert long_gap==F((q-3)*(4*q-5),2*(q-1))>=0
    assert boundary_gap==F((q-3)*(2*q-3),2*(q-1))>=0
    high.append({'q':q,'long_exponent':str(beta+2),
                 'short_exponent':q+6,'boundary_exponent':str(q+F(q,q-1)+5)})
reserves={'M':8-1-1-1,'raw':8-1-1,'slot':8-1}
assert reserves=={'M':5,'raw':6,'slot':7}
result={'arithmetic':'Python standard-library Fraction exact',
        'moving_gradient_v_residual':str(gradient_first),
        'peak_v2_coefficients_in_n':list(map(str,peak_v2)),
        'peak_minus_qnorm_v2_coefficients':['0','0','0'],
        'base_hessian_eigenvalues':['-1','-1/3'],
        'remaining_reserve_units':reserves,
        'shell_total_v_exponent':'(2-n)+(n-1)+8 = 9',
        'normal_width_exponent':str(6+F(7,2)),
        'high_support_samples':high,
        'N9_error_exponents':['9','9','10','little-o(9)','10'],
        'max_all_m_exponent':str(1+F(9,2)),
        'sigma_error_exponent':6,
        'geometry_certified':False,'all_assertions_passed':True}
print(json.dumps(result,indent=2))
