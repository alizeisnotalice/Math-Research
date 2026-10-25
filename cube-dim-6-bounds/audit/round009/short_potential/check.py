#!/usr/bin/env python3
"""F009: exact old-data reprocessing; eight times, zero kernel evaluations."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
import sympy as s
from flint import arb,ctx,fmpz
ctx.prec=384
START=time.monotonic()
HERE=Path(__file__).resolve().parent
MAIN=HERE.parents[2]
GP=MAIN/'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json'
HP=MAIN/'upper/round007/g_time_moments/results.json'
g=json.loads(GP.read_text());h=json.loads(HP.read_text())
assert g['status']=='PASS_EXACT_FMPQ_LOCAL_REFERENCE'
assert g['point_numerators']==h['center_numerators']
assert g['denominator']==h['denominator']
assert g['point_numerators'][0]==g['denominator']
assert len(set(g['point_numerators']))==6 and all(a>0 for a in g['point_numerators'])
# Any point in a span of two body directions has <=2 distinct absolute coordinate values.
# Therefore this saved point is off S2, where K_T is continuous.
G=arb(fmpz(g['Q_over_pi5_numerator']))/arb(fmpz(g['Q_over_pi5_denominator']))/arb.pi()**5
assert G.overlaps(arb(g['G_arb_384']))
def ar(q):
    q=Q(q);return arb(q.numerator)/q.denominator
def sr(a):return a.str(85)
expected=[Q(1,1000)*2**i for i in range(8)]
assert [Q(a['tau']) for a in h['time_nodes']]==expected
rows=[];possible=[];excluded=[]
for node in h['time_nodes']:
    Tq=Q(node['tau']);T=ar(Tq)
    H=arb(node['H_arb384']);K=G-H
    assert K>0 and H>0 and H<G
    threshold=T/4096
    point_cost=T/(64*K)
    if K<threshold:
        label='EXCLUDED_FOR_STRICT_C6_BELOW_64_BY_THIS_CONTINUOUS_POINT'
        excluded.append(str(Tq));assert point_cost>64
    elif K>threshold:
        label='PASSES_ONLY_THIS_POINT_NECESSARY_TEST_GLOBAL_MINIMUM_UNCERTIFIED'
        possible.append(str(Tq));assert point_cost<64
        assert K>T/4000
    else:raise AssertionError('Threshold comparison unresolved')
    rows.append({'T':str(Tq),'H_saved':sr(H),'K_point':sr(K),'K_over_T':sr(K/T),
      'required_global_minimum_threshold_exact':str(Tq/4096),
      'point_surplus_over_threshold':sr(K-threshold),
      'point_cost_T_over_64K':sr(point_cost),
      'optional_global_target_for_cost_62_point_5':str(Tq/4000),
      'true_model_cost_is_at_least_point_cost':True,'diagnostic':label})
assert possible==['4/125','8/125'] and len(excluded)==6

# Exact multiplier identities: no invocation of a kernel.
a,r,T,t=s.symbols('a r T t',positive=True)
khat=(1-s.exp(-T*a))/a
assert s.simplify(s.integrate(-a*s.exp(-t*a),(t,0,T))-(s.exp(-T*a)-1))==0
assert s.simplify(-a*khat.subs(a,r*a)-(s.exp(-r*T*a)-1)/r)==0
assert s.limit(khat,a,0)==T
assert s.simplify(khat.subs(a,r*a)-(1-s.exp(-r*T*a))/(r*a))==0
assert s.simplify(T/(64*(T/4096))-64)==0
assert s.simplify(T/(64*(T/4000))-s.Rational(125,2))==0
paths=[GP,HP,GP.with_name('exact_local_reference.py'),MAIN/'upper/theory/001_n6_body_frame_and_singular_strata.md']
out={'status':'PASS_CONDITIONAL_SHORT_POTENTIAL_TRANSFER_AND_EIGHT_OLD_TIME_DIAGNOSTICS',
 'new_C6_bound':False,'kernel_evaluations':0,'new_time_nodes':0,'old_times_reprocessed':8,'MC_samples':0,'matlab_calls':0,
 'same_exact_point':g['point_numerators'],'point_denominator':g['denominator'],
 'point_off_S2_by_six_distinct_positive_coordinates':True,
 'G_reconstructed_from_exact_coefficient':sr(G),'G_overlaps_saved_interval':True,
 'multiplier_identities':'L K_T=p_T-delta_0; dilation endpoint r*T; L(K_T)_r=(p_(rT)-delta_0)/r',
 'excluded_times':excluded,'times_not_excluded_by_this_point_only':possible,
 'whole_cube_essential_minimum_certified':False,'rows':rows,
 'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'elapsed_seconds':time.monotonic()-START,'exit_code':0}
(HERE/'check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'G':sr(G),'rows':[{'T':a['T'],'K_point':a['K_point'],
  'point_cost':a['point_cost_T_over_64K'],'diagnostic':a['diagnostic']} for a in rows],
  'excluded_times':excluded,'remaining_times':possible,'elapsed_seconds':out['elapsed_seconds']},indent=2))
