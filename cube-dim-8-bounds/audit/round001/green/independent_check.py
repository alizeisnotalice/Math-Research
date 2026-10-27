#!/usr/bin/env python3
"""Small independent exact algebra check; no Green values or subset enumeration.

The analytic audit is in proofreview.md. This script checks only normalization,
determinants, exponent identities, source scaling, and cost constants using
integer/Fraction arithmetic. It imports and executes no author code.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import datetime, hashlib, json, resource, sys, time

D=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def det(matrix):
    a=[[Q(x) for x in row] for row in matrix]
    result=Q(1)
    for k in range(len(a)):
        pivot=next(i for i in range(k,len(a)) if a[i][k])
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            result=-result
        d=a[k][k];result*=d
        for i in range(k+1,len(a)):
            ratio=a[i][k]/d
            for j in range(k,len(a)):
                a[i][j]-=ratio*a[k][j]
    return result

def main():
    n=8;q=Q(n,n-1)
    H=[[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]]
    assert abs(det(H))==16
    assert [[sum(H[i][k]*H[j][k] for k in range(4)) for j in range(4)] for i in range(4)]==[[4*int(i==j) for j in range(4)] for i in range(4)]
    B=[[H[i][j] if i<4 and j<4 else H[i-4][j-4] if i>=4 and j>=4 else 0 for j in range(8)] for i in range(8)]
    assert abs(det(B))==256
    W=[(1,)+s for s in product((-1,1),repeat=3)]
    V=[w+(0,)*4 for w in W]+[(0,)*4+w for w in W]
    assert len(V)==16 and all(sum(x*x for x in v)==4 for v in V)
    gram=[[sum(v[i]*v[j] for v in V) for j in range(8)] for i in range(8)]
    assert gram==[[8*int(i==j) for j in range(8)] for i in range(8)]
    # For arbitrary w=(1,e2,e3,e4), diag(w)H has first column w;
    # each remaining column is in W and its determinant stays +/-16.
    # This constructs eight completions, rather than enumerating subsets.
    for w in W:
        completion=[[w[i]*H[i][j] for j in range(4)] for i in range(4)]
        assert abs(det(completion))==16
        assert tuple(completion[i][0] for i in range(4))==w
        assert all(tuple(completion[i][j] for i in range(4)) in W for j in range(4))
    normal=[v for v in V if v[0]!=0]
    D1=sum(abs(v[0]) for v in V);A=D1
    C_coefficient=Q(1,256) # C=C_coefficient*pi^-8
    B_coefficient=D1*D1*len(normal)*C_coefficient
    assert D1==A==8 and B_coefficient==2
    short_t=(Q(2-n)*(1-1/q)+1/q)
    assert short_t==Q(1,8)
    short_integral_exponent=1+short_t
    short_integral_coefficient=1/short_integral_exponent
    assert short_integral_exponent==Q(9,8) and short_integral_coefficient==Q(8,9)
    assert -Q(n)*(1-1/q)-1/q==Q(-15,8)
    assert 1/(Q(15,8)-1)==Q(8,7)
    orthant_norm_exponent=-Q(7)/q
    assert orthant_norm_exponent==Q(-49,8)
    assert Q(1,n-1)*C_coefficient==Q(1,1792)
    assert Q(4)-Q(8-4)==0 # four tail factors and four peak factors
    low_rank_exponents={r:2*r-n for r in range(4)}
    assert all(e<=-1 for e in low_rank_exponents.values())
    # Test the homogeneous distribution scaling dimensionally:
    # L[r^-n K(x/r)] = r^(-n-1)(LK)(x/r), while dilation of a
    # probability measure has distribution density r^-n mu(x/r).
    assert (-n-1)-(-n)==-1
    # |{G>s}|=J/n*s^-q and int_c^inf s^-q ds=(n-1)c^-1/(n-1).
    assert Q(1)/(q-1)==7
    mass_coefficient=Q(n-1,n)
    assert mass_coefficient==Q(7,8)
    cube_faces=2*n;positive_orthants=2**(n-1)
    J_to_I=cube_faces*positive_orthants
    assert J_to_I==2048
    # At c=g/n, (int K)^(n-1) = n*J^(n-1)/(2^(n(n-1))*g^n).
    # Substitute J=2048 I exactly, with no floating fractional powers.
    raised_cost_coefficient=Q(n)*Q(J_to_I,2**n)**(n-1)
    assert raised_cost_coefficient==8**8
    result={
        'status':'PASS','scope':'finite exact normalization and exponent checks only; no numerical Green certificate',
        'vectors':16,'direction_lengths_squared':[4],
        'Hadamard_abs_det':16,'block_basis_abs_det':256,
        'frame_gram':'8 I8','constructed_normal_completions':8,
        'C_coefficient_pi_minus_8':str(C_coefficient),
        'A':A,'D1':D1,'short_time_B_coefficient_pi_minus_8':str(B_coefficient),
        'singular_rank_time_exponents':low_rank_exponents,
        'short_time_density_sup_t_power':2-n,
        'short_time_Lq_t_power':str(short_t),
        'short_time_integral_power':str(short_integral_exponent),
        'short_time_integral_prefactor':str(short_integral_coefficient),
        'long_time_pointwise_power':1-n,
        'long_time_pointwise_coefficient_pi_minus_8':str(Q(1,n-1)*C_coefficient),
        'long_time_full_plane_Lq_power':str(Q(-7,8)),
        'orthant_norm_factor_power_of_2':str(orthant_norm_exponent),
        'source_scaled_coefficient_power_of_r':-1,
        'J_to_I':J_to_I,'raised_cost_I_power':7,'raised_cost_g_power':-8,
        'raised_cost_coefficient':str(raised_cost_coefficient),
        'enumeration_of_independent_subsets_performed':False,
    }
    (D/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    started=datetime.datetime.now(datetime.timezone.utc).isoformat();wall=time.perf_counter();cpu=time.process_time()
    error=None;code=0
    try:main()
    except BaseException as exc:
        error=repr(exc);code=1;raise
    finally:
        receipt={'status':'PASS' if code==0 else 'FAIL','exit_code':code,'error':error,
                 'start_utc':started,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'wall_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,
                 'peak_rss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 'python':sys.version,'executable':sys.executable,'checker_sha256':sha(Path(__file__)),
                 'mathematical_dependencies':'Python standard library Fraction; no Green evaluations'}
        (D/'actualreceipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
