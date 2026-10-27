#!/usr/bin/env python3
"""Exact algebraic audit of one fixed probe and the analytic cluster family.

No pair enumeration, optimizer, author import, or new candidate search.
Geometric inequalities are certified from endpoint bounds and sign patterns.
"""
from fractions import Fraction as Q
from pathlib import Path
import datetime, hashlib, json, resource, sys, time

D=Path(__file__).resolve().parent
A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube8-sol-theory/cube-dim-8-bounds/upper/round003/coupling')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    plan=json.loads((A/'probe_plan.json').read_text())
    saved=json.loads((A/'verification.json').read_text())
    assert plan['candidate_count']==saved['probe_count']==1
    assert saved['probe_plan_sha256']==sha(A/'probe_plan.json')
    n=8;d=2**n-1;a=Q(8);t=Q(2);b=a*(1+1/t)
    M=a**(-n);tau=Q(1,2**30);S=d*M;D2=d*M*M
    # Any two distinct vertices of {-7,+7}^8 differ in some coordinate
    # by exactly 14. The sup-distance is 14, regardless of how many
    # signs differ. Thus all pairwise assertions follow without a loop.
    vertex_abs=Q(7);distance=2*vertex_abs
    assert vertex_abs<a and distance>b and distance>2
    assert a+vertex_abs<16 and 0<=tau
    assert b==12 and tau>0 and a>2*t and 1<t
    E=Q(0);calibrated=S*S/(2*d)-D2/2
    stronger=S*S/(2*(d-1))-D2/2
    assert calibrated==E==0 and stronger==Q(255,142989288169013248)>0
    assert D2==M*S
    assert Q(saved['weight_M'])==M and Q(saved['unpaid_mass_S'])==S
    assert Q(saved['sum_squared_weights'])==D2 and Q(saved['exact_energy'])==E
    assert Q(saved['tau'])==tau and Q(saved['false_stronger_rhs_D254'])==stronger
    assert saved['family_size']==2**n and saved['unpaid_size']==d
    assert saved['pair_checks']==(2**n)*(2**n-1)//2==32640
    # One algebraic family, not additional geometric probes. For t=2,L=8:
    # min abs=6m+1, max abs=8m-1, b=12m. For m>=1 the margins
    # across clusters and from the anchor are 2 and 2m respectively.
    # Verify polynomial coefficients of these affine identities exactly.
    L=Q(8);q=1+1/t
    assert L>2*t and L*(2-q)>=4
    cross_slope=2*(L-2)-q*L
    cross_intercept=Q(2)
    anchor_margin_slope=(2*L-2)-q*L
    assert cross_slope==0 and cross_intercept==2 and anchor_margin_slope==2
    # Same-class diameter=2m-2; the b-star margin=10m+2>0.
    assert q*L-2==10
    # For arbitrary N=m^8, each class has N weights M=(8m)^-8.
    # Algebraic identities after dividing by M^2 use only integer N,d:
    # d*N*(N-1)/2 = (d*N)^2/(2d) - d*N/2.
    # Compare coefficients of the polynomial in N.
    assert [Q(-d,2),Q(d,2)]==[Q(-d,2),Q(d*d,2*d)]
    fixed_positive_mass=Q(d,8**8)
    # tau=(N-1)M is positive for m>=2, and S=d(tau+M) exactly.
    assert d==255 and Q(2*d+1)==511 and d+1==256
    rough_cap_comparison=(Q(81,16)+3)**2
    assert rough_cap_comparison==Q(16641,256)<256<511
    result={
        'status':'PASS','scope':'selector-model geometry/energy obstruction only; no new C8 bound or C8 lower bound',
        'fixed_registered_probe':{'n':8,'a':str(a),'t':str(t),'b':str(b),'tau':str(tau),
          'minimum_distinct_center_distance':str(distance),'star_margin_at_zero':str(a-vertex_abs),
          'anchor_incoming_weight':'0','all_selector_tests':'0',
          'all_small_cube_increments':'entire open cube',
          'M':str(M),'S':str(S),'D2':str(D2),'E':'0','calibrated_rhs':str(calibrated),
          'false_stronger_rhs_d254':str(stronger)},
        'analytic_family_endpoint_coefficients':{'a':'8m','min_abs':'6m+1','max_abs':'8m-1',
          'cross_cluster_distance_lower':'12m+2','cross_capture_margin':'2',
          'anchor_distance_lower':'14m','anchor_capture_margin':'2m',
          'within_cluster_diameter':'2m-2','within_capture_margin':'10m+2',
          'cluster_size':'m^8','positive_fixed_mass':str(fixed_positive_mass),
          'energy_identity':'E=S^2/510-M*S/2',
          'exact_eligibility_threshold_for_m_ge_2':'tau=(m^8-1)*(8m)^-8',
          'last_index_linear_equality':'S=255*(tau+M)'},
        'rough_existing_cap_family_comparison':str(rough_cap_comparison),
        'coarse_scalar_infima_energy_then_last_index':[511,256],
        'pair_enumeration_performed':False,'new_geometric_candidate_search_performed':False,
        'source_plan_sha256':sha(A/'probe_plan.json'),
        'source_verification_sha256':sha(A/'verification.json'),
    }
    (D/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    start=datetime.datetime.now(datetime.timezone.utc).isoformat();wall=time.perf_counter();cpu=time.process_time()
    code=0;error=None
    try:main()
    except BaseException as exc:code=1;error=repr(exc);raise
    finally:
        receipt={'status':'PASS' if code==0 else 'FAIL','exit_code':code,'error':error,
          'start_utc':start,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'wall_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,
          'peak_rss_native_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          'python':sys.version,'executable':sys.executable,'checker_sha256':sha(Path(__file__)),
          'mathematical_dependencies':'Python standard-library Fraction only'}
        (D/'actualreceipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
