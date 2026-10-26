"""Primary binding of complete103/112 integrals and their exact cost comparisons.

Analytic review is recorded separately. This finite acceptance does not evaluate
H, R, log, or an old numerical baseline and does not assume g110 is proved.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json,gzip,os,time
H=Path(__file__).resolve().parent
P=H.parent


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fourth(q):
    D=1<<320
    n=isqrt(isqrt(q.numerator*D**4//q.denominator))
    assert n**4*q.denominator<=q.numerator*D**4<(n+1)**4*q.denominator
    return F(n,D),F(n+1,D)


def check_budget(b,U,E,g):
    B,I,C=map(F,[b['B_upper'],b['I_upper'],b['cost_upper']])
    assert B>E>0 and I>0 and C>0
    assert (B-E)**5>=U**4 and I**4>=B**5 and C**4>=(5*B/g)**5
    return B,I,C


def main():
    t=time.monotonic();rp=H/'primary_integrals_registration.json';r=json.loads(rp.read_text())
    for p,h in r['locks'].items():assert sha(P/p)==h,p
    a=json.loads((P/'audit_adaptive113/results.json').read_text())
    assert a['status']=='pass' and a['new_children']==28427 and a['new_endpoint_inputs']==2075171
    assert a['all_new_range_mean_bins_exact'] and a['rational_sum_fallbacks']==0
    led=a['full_ledger'];assert led['final_nodes']==97972 and led['physical_frontier']==94580
    assert led['volume']=='1' and led['all_parent_prefixes']==2048
    g=F(21,200000);E=F(json.loads((P/'upper_adaptive103/registration.json').read_text())['E076_upper'])
    b103=json.loads(gzip.decompress((P/'upper_adaptive103/g105_prefix_budgets.json.gz').read_bytes()))[-1]
    assert all(b103[k]==led['budget_g105'][k] for k in ('I_upper','cost_upper','strict32'))
    U103=F(led['budget_g100']['UH'])
    B103,I103,C103=check_budget(b103,U103,E,g)
    cap103=F(51807888,10**12);cap095=F(54702517,10**12)
    assert I103<=cap103<cap095 and F(led['I_simple'])==cap103
    floor=json.loads((P/'audit_adaptive113/floor_results.json').read_text())
    assert floor['strict_family_exclusion'] and floor['not_actual_I_or_model_C5_lower_bound']
    assert F(floor['independent_I_expression_lower'])>F(floor['T_upper'])
    assert F(floor['family_cost_lower'])>32 and F(floor['family_S_upper'])<0
    w=json.loads((P/'audit_path121/whole_results.json').read_text())
    assert w['status']=='pass' and w['nodes']==97972 and w['positive_bands']==7151956
    assert w['blocks']==25 and w['mean_exact_bins']==w['root256_same_final_bin']==97972
    assert w['mean_improved']==w['self_improved']==97971
    assert w['all_full_offset_and_old_input_links'] and w['new_H_log_Arb_calls']==0
    last=w['final'];assert last['volume']=='1' and last['physical_frontier']==94580
    assert last['all_lower_and_topology_unchanged'] and last['all_unprocessed_self_paid_old']
    b=json.loads(gzip.decompress((P/'upper_path112/numeric_hierarchy.json.gz').read_bytes()))['budget']
    assert all(b[k]==last['final_budget'][k] for k in ('UH','I_upper','cost_upper','strict32'))
    U=F(b['UH']);B,I,C=check_budget(b,U,E,g)
    cap=F(47619360,10**12);assert I<=cap<cap103 and F(last['I_simple'])==cap
    assert not b['strict32'] and F(last['S_upper'])<0
    IL=F(last['LH']);assert 0<IL<U<U103
    gu=json.loads((P/'audit_zero041/result.json').read_text())
    assert gu['point']==['1','1','1','2/3','1/16']
    gU=F(gu['pure_Fraction_G_interval'][1]);assert g<gU
    model_low=fourth(3125*IL**4/gU**5)[0]
    cost=fourth(3125*cap**4/g**5)
    assert model_low<32<cost[0]
    threshold=[32*x for x in fourth((g/5)**5)]
    margin=[threshold[0]-cap,threshold[1]-cap];assert margin[1]<0
    hypothetical_g=F(11,100000);cond=fourth(3125*cap**4/hypothetical_g**5)
    assert cond[1]<32 and 32**4*hypothetical_g**5-3125*cap**4>0
    out={'status':'primary_accepts_complete_alpha1_I103_then_I112',
         'qualified_g':str(g),'I103_cap':str(cap103),'best_complete_I_cap':str(cap),
         'exact_I095_cap_minus_new':str(cap095-cap),
         'exact_I103_cap_minus_new':str(cap103-cap),
         'exact_I092_cap_minus_new':str(F(6010267,10**11)-cap),
         'exact_I010_cap_minus_new':str(F(71555899,10**12)-cap),
         'complete_H512_integral_upper':str(U),'complete_G_integral_lower':str(IL),
         'shorttime_E_norm':str(E),'B_upper':str(B),
         'no_extra_tube':True,'full_face_all_positive_time':True,
         'general_L1_all_common_positive_scales':'same029 retained; no model/time mixing',
         'certified_simple_cap_cost_interval':list(map(str,cost)),
         'certified_simple_cap_cost_display':list(map(float,cost)),
         'fixed_model_cost_interval':[str(model_low),str(cost[1])],
         'fixed_model_excluded':False,'model_cost_not_a_C5_lower_bound':True,
         'S_certificate_interval_for_simple_cap':list(map(str,margin)),
         'old_current_mean_variance_only_family':'strictly excluded under113 fixedM_U/E/g/tree;111 changesM_U and lies outside that family',
         'conditional_g110_pair':{'g_not_yet_qualified':True,'g':str(hypothetical_g),
            'cost_interval':list(map(str,cond)),'display_only':list(map(float,cond)),
            'exact_positive_gate_margin':str(32**4*hypothetical_g**5-3125*cap**4),
            'not_a_C5_bound':True},
         'formal_C5_lower_changed':False,'formal_C5_upper':'32',
         'formal_C5_upper_changed':False,'dual_breakthrough_time':None,'post_dual_rounds':0,
         'primary_review':'111 report/pathmean/secant;112 report/production_v2;113 report/floor/changed ledger;121 report/scalar/full_v2 read; full positive path, convex chord, source bindings, roots, means, full volume and029 reviewed',
         'independence_limits':'original H Arb/fan and075/085/091 shared;121 independently rederived new scalar weights/mean/chord/DFS and source identity; no H rerun',
         'latest_model_policy':'no new/follow-up Astra dispatch; archived evidence only',
         'new_H_R_log_or_old_baseline_calls':0,'pid':os.getpid(),'seconds':time.monotonic()-t,
         'source_sha256':sha(__file__),'registration_sha256':sha(rp),'input_locks':r['locks']}
    p=H/'primary_integrals_acceptance.json';assert not p.exists();p.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','best_complete_I_cap','certified_simple_cap_cost_display','formal_C5_upper_changed','pid','seconds']}))


if __name__=='__main__':main()
