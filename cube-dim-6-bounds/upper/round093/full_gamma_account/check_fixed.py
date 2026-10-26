#!/usr/bin/env python3
"""A093 exact rational full-G upper-Hessian ledger; no density or native code."""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import os
import resource
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parent
sys.set_int_max_str_digits(100000)
resource.setrlimit(resource.RLIMIT_CPU,(27,27))
START,CPU=time.monotonic(),time.process_time()
OUT={"schema":"A093-fixed-results-v1","pid":os.getpid(),"status":"RUNNING",
     "configurations_started":0,"configurations_completed":0,"configurations":[],
     "actual_work":{"new_exact_arithmetic_configurations":0,"Hessian_component_polynomials":0,
                    "nonzero_component_coefficients":0,"time_candidates_evaluated":1,
                    "time_scan_points":0,"body_four_subsets_enumerated":0,"body_matrices_enumerated":0,
                    "old085_checks_rerun":0,"full6D_density_evaluations":0,"native_calls":0,
                    "actual_G_vertices":0,"actual_F_or_gradient_fields":0,"old_baseline_reruns":0,
                    "matlab_calls":0,"recursive_agents":0,"automatic_retries":0}}
COMP={}


def save():
    OUT["wall_seconds"]=time.monotonic()-START;OUT["cpu_seconds"]=time.process_time()-CPU
    (ROOT/"results.json").write_text(json.dumps(OUT,indent=2)+"\n")


def begin(name):
    OUT["configurations_started"]+=1;OUT["actual_work"]["new_exact_arithmetic_configurations"]+=1
    OUT["active_configuration"]=name;save()


def finish(name,data):
    OUT["configurations"].append({"id":name,"status":"PASS",**data})
    OUT["configurations_completed"]+=1;OUT["active_configuration"]=None;save()


def add(name,coefficients):
    p={int(k):F(v) for k,v in coefficients.items() if v}
    assert p and all(v>0 for v in p.values()) and name not in COMP
    COMP[name]=p
    OUT["actual_work"]["Hessian_component_polynomials"]+=1
    OUT["actual_work"]["nonzero_component_coefficients"]+=len(p)


def value(p,t):return sum((v*t**k for k,v in p.items()),F(0))
def encoded(p):return {str(k):str(v) for k,v in sorted(p.items())}
def displayed(v):return format(float(v),".14e")


try:
    reg=json.loads((ROOT/"registration.json").read_text())
    for n,hsh in reg["frozen_local_sha256"].items():
        assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==hsh,n
    inp=json.loads((ROOT/"input_manifest.json").read_text())
    data=json.loads((ROOT/"upstream_constants.json").read_text())
    assert len(inp["configurations"])==8
    T,Tmax,rho,h,R,g0=map(F,[inp["T"],inp["Tmax_for_reused_r_bounds"],inp["rho"],inp["face_halfwidth"],inp["scalar_cutoff_R"],inp["candidate_g0"]])
    A=[F(v)/27 for v in data["pi3_natural_A_order0to4"]]
    E=[F(v)/27 for v in data["pi3_extended_A_order0to4"]]
    rjet=[F(v)/27 for v in data["pi3_extended_r_order0to2"]]
    m0,m1,m2=[F(data[n])/27 for n in ["M0_coefficient_pi3","M1_sum_coefficient_pi3","M2_sum_coefficient_pi3"]]
    cq,cf=F(22),F(13)
    OUT["runtime"]={"python":sys.version,"executable":sys.executable,"arithmetic":"fractions.Fraction, integers only"}
    OUT["registration_sha256"]=hashlib.sha256((ROOT/"registration.json").read_bytes()).hexdigest()

    begin("G001-uniform-constants")
    assert 5*h*h<=rho*rho and F(1,3)-rho>F(1,4)
    assert 5**5<64**2 and 48==3*16
    assert F(48,32)==F(3,2)
    assert 32**2>=6*13**2 and 16**2*2>22**2
    n4,n3=comb(32,4),comb(32,3)
    assert n4==35960 and n3==4960
    assert F(factorial(7),64*6*7)==F(15,8)
    assert F(4*factorial(4),3*8)==4
    DS=F(12,8)+F(16,4)
    DSL=F(12,8)+F(16,4)*F(9,4)
    DSL3=F(12,8)+F(16,4)*F(9,4)**3
    assert (DS,DSL,DSL3)==(F(11,2),F(21,2),F(753,16))
    finish("G001-uniform-constants",{
        "determinant_lower":32,"cofactor_upper":48,"inverse_entry_absolute_upper":"3/2",
        "four_tail_inverse_row_squares_upper":"54","two_completion_inverse_row_squares_upper":"27",
        "four_subset_count_upper":n4,"three_subset_count_upper":n3,
        "full_coercivity_lower":str(cf),"quotient_coercivity_lower":str(cq),
        "sum_D":str(DS),"sum_D_Lambda":str(DSL),"sum_D_Lambda_cubed":str(DSL3),
        "full_face_box_volume":str((2*h)**5),"positive_unit_face_intersection_volume":str(2*h**5),
        "geometry":"face x1=1, other5 coordinates centered at (1,1,1,1/2,0); full symmetric box lies in rho ball"})

    begin("G002-rank4")
    a1,a2=4/R,20/R**2
    b1,b2=1+4/R,2+8/R+12/R**2
    pref=F(n4,32)/(729*R**8)
    add("rank4_aggregate",{1:pref*27*(b2+b1*b1),2:pref*108*a1*b1,3:pref*18*(a2+3*a1*a1)})
    finish("G002-rank4",{"coefficients":encoded(COMP["rank4_aggregate"]),"at_fixed_T":str(value(COMP["rank4_aggregate"],T)),
        "dependency":"086 author theorem; independent acceptance status tracked separately", "actual_independent_four_enumeration":0})

    begin("G003-cutoff-removal")
    s1,s2=1+2/R,2+4/R+6/R**2
    H1=max(4/R,s1);H2=max(20/R**2,s2,H1*H1)
    lam,V=F(9,4),F(5)
    pcut=DS/(27*R**2)
    add("unexpanded_cutoff_removal",{1:pcut*lam**2*H2*E[0],2:pcut*V*lam*H1*E[1],3:pcut*V**2*E[2]/3})
    finish("G003-cutoff-removal",{"coefficients":encoded(COMP["unexpanded_cutoff_removal"]),"Lambda":str(lam),"amplitude_transport":str(V),
        "H1":str(H1),"H2":str(H2),"weighted_columns_merged":False,"amplitude_expanded_before_cutoff_removal":False})

    begin("G004-quadratic-and-t5")
    beta1,beta2=F(1),F(5,4)
    add("weighted_quadratic_remainder",{
        1:E[2]*m0+12*beta1*E[2]*m1+14*beta2*E[2]*m2,
        2:(12*E[3]*m1+28*beta1*E[3]*m2)/2,
        3:14*E[4]*m2/3})
    hs2=max(s2,s1*s1)
    add("t5_inplane_correction",{1:m0*rjet[0]*hs2,2:m0*2*rjet[1]*s1,3:m0*4*rjet[2]/3})
    finish("G004-quadratic-and-t5",{"quadratic_coefficients":encoded(COMP["weighted_quadratic_remainder"]),
        "t5_coefficients":encoded(COMP["t5_inplane_correction"]),"fixed_fiber_score_bounds":[str(beta1),str(beta2)],
        "C_norm_squared":28,"C_norm_rational_upper":6,"fixed_tangent_transport_norm_upper":1,
        "frozen_r_bounds_valid_up_to":str(Tmax),"old085_tables_recomputed":False})

    begin("G005-signed-and-normal")
    add("negative_cusp",{0:F(data["cusp_semiconcavity_coefficient_pi4"])/81})
    F2=6/(3*cq**5)
    add("signed_renormalized_F_curvature",{0:2*A[1]*28*F2})
    add("signed_upper_tail_curvature",{-1:2*A[1]*rho*48*DSL3/27})
    add("signed_coefficient_cross_including_linear",{1:4*A[2]*3*DSL/(2*27)})
    add("signed_coefficient_second_including_linear",{2:2*A[3]*DS/(4*27)})
    Q0,Q2=1/(9*cq**3),4/(9*cq**5)
    add("normal_layer_constant",{1:A[2]*Q0})
    add("normal_deficit_upper_tail",{-1:Q2*(A[2]*rho**2/2+2*A[1]*rho+A[0])})
    finish("G005-signed-and-normal",{"signed_F_Hessian_majorant":str(F2),"q1_peak_majorant":str(Q0),"q_Hessian_scale_majorant":str(Q2),
        "signed_layer_identity":"L_j=F_j+T Dm_j(1,0).z-Q_j(T,z)",
        "linear_term_paid":"L_j and DL_j direct finite-time bounds pay all coefficient cross/second terms; linear term has zero pure z Hessian only",
        "normal_layer_identity":"a(y)[Tq1(0)+integral_T^infty t^3(q_t(0)-q_t(z))dt]",
        "coefficients":{name:encoded(COMP[name]) for name in list(COMP)[4:]}})

    begin("G006-original-tail-and-value")
    tail7=F(15,8)/(27*cf**8)
    add("original_full6D_t_ge_T",{-7:tail7})
    Rvalue=F(1,640)
    Mvalue=F(n3)/(243*cf**3*Rvalue**6)
    Pvalue=F(15,8)/(27*cf**6)
    Gbar=Mvalue*Tmax+Pvalue/(5*Tmax**5)
    finish("G006-original-tail-and-value",{"original_tail_coefficient_Tminus7":str(tail7),
        "value_qualification_distance_lower":"1/4","value_large_scalar_threshold":str(Rvalue),
        "small_time_density_upper":str(Mvalue),"large_time_density_scale_upper":str(Pvalue),
        "true_G_uniform_upper_on_rho_ball":str(Gbar),"display_G_uniform_upper":displayed(Gbar),
        "actual_density_evaluations":0,"scope":"analytic positive-three-tail value bound plus original6D Fourier tail, no Hessian conclusion from value bound"})

    begin("G007-complete-ledger")
    total={}
    for p in COMP.values():
        for k,v in p.items():total[k]=total.get(k,F(0))+v
    assert set(total)=={-7,-1,0,1,2,3}
    Gamma=value(total,T);penalty=Gamma*5*h*h/2;required=g0+penalty
    rows=[]
    for name,p in COMP.items():
        v=value(p,T)
        rows.append({"component":name,"coefficients_in_T":encoded(p),"Gamma_at_T":str(v),"display_Gamma":displayed(v),
                     "fraction_of_total_exact":str(v/Gamma),"display_fraction_of_total":float(v/Gamma)})
    assert sum(F(r["Gamma_at_T"]) for r in rows)==Gamma
    assert penalty>Gbar
    lower_gamma_power8=F(8**8,7**7)*total[1]**7*total[-7]
    lower_penalty_power8=lower_gamma_power8*(5*h*h/2)**8
    assert lower_penalty_power8>Gbar**8
    # Decimal-power bracket for an analytic envelope lower bound, no T sampling.
    exponent=0
    while F(10)**(8*(exponent+1))<=lower_penalty_power8:exponent+=1
    assert F(10)**(8*exponent)<=lower_penalty_power8<F(10)**(8*(exponent+1))
    finish("G007-complete-ledger",{"components":rows,"total_coefficients_in_T":encoded(total),
        "Gamma_rational_upper":str(Gamma),"display_Gamma":displayed(Gamma),
        "box_penalty_Gamma_times_5h2_over2":str(penalty),"display_box_penalty":displayed(penalty),
        "required_every_vertex_lower_at_least":str(required),"candidate_g0":str(g0),
        "true_G_uniform_upper":str(Gbar),"fixed_ledger_vertex_test_impossible_even_with_exact_vertices":True,
        "all_time_envelope_lower_bound_eighth_power":str(lower_gamma_power8),
        "all_time_penalty_lower_bound_eighth_power":str(lower_penalty_power8),
        "all_time_penalty_lower_bound_decimal_order":exponent,
        "all_time_penalty_lower_bound_exceeds_true_G_upper":True,
        "all_time_scope":"Every 0<T<=1/64 with the same frozen085 remainder bounds and this rational majorant expression; not actual minimal curvature or model exclusion",
        "unique_stationary_point_proof":"Gamma''(T)=2C2+6C3*T+2Cminus1/T^3+56Cminus7/T^9>0; no candidate T evaluated",
        "full_G_analytic_terms_remaining_unpaid":[],"dependency_086_final_acceptance":"pending status update from parent"})

    begin("G008-same-box-and-S")
    a,b=map(F,[inp["comparison_positive_window"]["a"],inp["comparison_positive_window"]["b"]])
    radius_upper=F(9,4)*h
    score_radius=radius_upper/a
    assert score_radius<1 and F(9,4)**2>5
    Fbox_multiplier=1/(1-score_radius)
    positive_window_penalty_multiplier=F(129,256)/a**2*5*h*h*Fbox_multiplier
    log_penalty=5*h*h/(32*a*a)
    assert 0<log_penalty<1
    Bcritical=F(32,6)*g0
    assert Bcritical==F(4,48750)
    Icritical_power5=64**5*(g0/6)**6
    finish("G008-same-box-and-S",{
        "comparison_window":{"a":str(a),"b":str(b)},
        "positive_window_box_upper_multiplier_of_center_Fupper":str(Fbox_multiplier),
        "positive_window_vertex_penalty_multiplier_of_center_Fupper":str(positive_window_penalty_multiplier),
        "display_positive_window_vertex_penalty_multiplier":float(positive_window_penalty_multiplier),
        "positive_window_condition":"min32 F_vertex_lower >= g0 + multiplier*F_center_upper; all matching scalar endpoints UNPAID",
        "original_mixture_log_curvature_penalty":str(log_penalty),"display_original_mixture_log_curvature_penalty":float(log_penalty),
        "zero_gradient_required_factor_lower":str(1+log_penalty),"zero_gradient_sufficient_factor_upper":str(1/(1-log_penalty)),
        "zero_gradient_hypothetical_sufficient_center_lower":str(g0/(1-log_penalty)),
        "original_mixture_condition":"log F_lower - h*gradient_l1_upper/F_lower - log_penalty >= log g0; matching center value and5 face gradients UNPAID",
        "kernel_decision":"Do not release32 true-G vertices for this excluded ledger. A new matching positive-window center+gradient pilot is a conditional candidate only, requiring separate parent release and production qualification.",
        "S_definition":"64*(g_L/6)^(6/5)-I_U",
        "candidate_Bcritical":str(Bcritical),"candidate_Icritical_fifth_power":str(Icritical_power5),
        "global_g_L":"UNPAID","complete_I_U":"UNPAID","certified_S":"UNPAID","new_global_g":False,"new_I":False,"new_Cn":False})
    assert OUT["configurations_completed"]==8
    OUT["status"]="PASS_EIGHT_CONFIGURATIONS_COMPLETE_RATIONAL_GAMMA_LEDGER_ENVELOPE_NONCOMPETITIVE"
    OUT["complete_Gamma_formula_computed"]=True
    OUT["numeric_vertex_certificate"]=False
    OUT["Gamma_dependency_086_final_acceptance"]="pending parent update"
    OUT["actual_G_or_F_evaluated"]=False
    OUT["global_g_or_I_or_S_certificate"]=False
    save()
    print(json.dumps({"status":OUT["status"],"pid":os.getpid(),"wall_seconds":OUT["wall_seconds"],"cpu_seconds":OUT["cpu_seconds"],
                      "Gamma_display":displayed(Gamma),"penalty_display":displayed(penalty),"true_G_upper_display":displayed(Gbar),"actual_work":OUT["actual_work"]}))
except BaseException as exc:
    OUT["status"]="FAILED_FIRST_STOP_NO_AUTOMATIC_RETRY"
    OUT["failure"]={"type":type(exc).__name__,"message":str(exc),"traceback":traceback.format_exc()}
    save();print(json.dumps({"status":OUT["status"],"pid":os.getpid(),"failure":str(exc)}));raise
