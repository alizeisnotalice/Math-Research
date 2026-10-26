#!/usr/bin/env python3
"""One fixed relational audit of saved certificates, not an A093 replay."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys


def require(ok,message):
    if not ok:raise AssertionError(message)


def main():
    root=Path(__file__).resolve().parent
    reg=json.loads((root/'registration.json').read_text())
    require(len(reg['configurations'])==1,'One relational certificate configuration')
    author=json.loads(Path(reg['author_results_path']).read_text())
    prior=json.loads(Path(reg['upstream085_results_path']).read_text())
    extracted=json.loads(Path(reg['author_extraction_path']).read_text())
    ac={row['id']:row for row in author['configurations']}
    pc={row['id']:row for row in prior['configurations']}
    extraction_pairs=[
        ('pi3_natural_A_order0to4','Q003-amplitude-C4','operator_norm_majorants_order0to4'),
        ('pi3_extended_A_order0to4','Q003-amplitude-C4','extended_amplitude_operator_majorants_order0to4'),
        ('pi3_extended_r_order0to2','Q004-time-remainder-C2','extended_remainder_operator_majorants_order0to2')]
    for field,config,key in extraction_pairs:
        require(extracted[field]==pc[config][key],'Saved upstream table copied exactly: '+field)
    require(extracted['cusp_semiconcavity_coefficient_pi4']==pc['Q007-single-cusp-budget']['single_negative_cusp_semiconcavity_bound_coefficient_pi4'],
            'Saved cusp coefficient copied exactly')
    for field,key in [('M0_coefficient_pi3','pi3_M0_majorant'),
                      ('M1_sum_coefficient_pi3','pi3_M1_sum_majorant'),
                      ('M2_sum_coefficient_pi3','pi3_M2_sum_majorant')]:
        require(extracted[field]==pc['Q005-labelled-quotient-bases'][key],
                'Original28-variable moment coefficient copied exactly: '+field)
    ledger=ac['G007-complete-ledger']
    rows=ledger['components']
    required=reg['required_component_names']
    require(len(rows)==12 and {row['component'] for row in rows}==set(required),'Exact12 component names, no duplication')
    coeff={int(k):F(v) for k,v in ledger['total_coefficients_in_T'].items()}
    require(set(coeff)=={-7,-1,0,1,2,3} and all(v>0 for v in coeff.values()),'Six positive envelope coefficients')
    residual=coeff.copy()
    nonzero=0
    stored_gamma=F(ledger['Gamma_rational_upper'])
    stored_value_residual=stored_gamma
    for row in rows:
        part={int(k):F(v) for k,v in row['coefficients_in_T'].items()}
        require(all(v>0 for v in part.values()),'Positive component coefficients')
        nonzero+=len(part)
        for power,value in part.items():residual[power]-=value
        saved=F(row['Gamma_at_T'])
        stored_value_residual-=saved
        require(F(row['fraction_of_total_exact'])*stored_gamma==saved,'Saved component fraction matches value')
    require(nonzero==20 and all(v==0 for v in residual.values()),'Stored coefficient ledger closes with20 entries')
    require(stored_value_residual==0,'Stored evaluated component amounts sum to saved total; no T reevaluation')
    h=F(reg['face_halfwidth']); penalty_scale=5*h*h/2
    penalty=F(ledger['box_penalty_Gamma_times_5h2_over2'])
    upper=F(ac['G006-original-tail-and-value']['true_G_uniform_upper_on_rho_ball'])
    require(penalty==penalty_scale*stored_gamma and penalty>upper,'Stored penalty exceeds separate true-G upper certificate')
    require(upper==F(ledger['true_G_uniform_upper']),'Same-ball true-G upper amount')
    # Verify the supplied exact all-T witness algebraically. No optimizer,
    # alternate T, time grid, or decimal-order search is evaluated.
    witness=F(ledger['all_time_envelope_lower_bound_eighth_power'])
    require(witness*7**7==8**8*coeff[1]**7*coeff[-7],'Weighted AM-GM eighth-power witness')
    penalty_witness=F(ledger['all_time_penalty_lower_bound_eighth_power'])
    require(penalty_witness==penalty_scale**8*witness,'Penalty witness normalization')
    require(F(10)**280<=penalty_witness<F(10)**288 and penalty_witness>upper**8,'Saved uniform lower witness proves10^35 bracket and exclusion')
    require(ledger['all_time_penalty_lower_bound_decimal_order']==35,'Author claimed only the verified decimal bracket')
    compare=ac['G008-same-box-and-S']
    a=F(compare['comparison_window']['a'])
    box_multiplier=F(compare['positive_window_box_upper_multiplier_of_center_Fupper'])
    window_penalty=F(compare['positive_window_vertex_penalty_multiplier_of_center_Fupper'])
    q=F(compare['original_mixture_log_curvature_penalty'])
    require(box_multiplier*(1-F(9,4)*h/a)==1,'Stored exponential majorant uses same box/window')
    require(window_penalty/q==F(129,8)*box_multiplier,'Independent dimensionless relation between window and log penalties')
    gate=F(compare['candidate_Bcritical'])
    gate5=F(compare['candidate_Icritical_fifth_power'])
    require(gate**6==gate5,'Equivalent B and I conditional S thresholds')
    require(compare['global_g_L']==compare['complete_I_U']==compare['certified_S']=='UNPAID','No fabricated global endpoints or S')
    output={'schema':'B093-saved-relational-certificate-results-v1','status':'PASS',
            'new_configurations':1,'source_tables_and_cusp_matched':4,'moment_coefficients_matched':3,'component_names_verified':required,
            'component_count':12,'nonzero_component_coefficients':20,
            'coefficient_checksum_residual':{str(k):str(v) for k,v in sorted(residual.items())},
            'saved_value_checksum_residual':str(stored_value_residual),
            'fixed_penalty_exceeds_true_G_upper':True,'all_T_penalty_eighth_power_witness_verified':True,
            'all_T_penalty_lower_at_least':'10^35','dimensionless_window_comparison_identity_verified':True,
            'conditional_B_I_gate_identity_verified':True,
            'new_time_candidates':0,'old_author_configuration_replays':0,
            'component_formulas_reconstructed':0,'density_or_kernel_calls':0,
            'actual_curvature_lower_bound_proved':False,'model_excluded':False,
            'dependencies':['Python standard library Fraction','Read-only saved093 and085 certificates'],
            'python':sys.version}
    (root/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','new_configurations':1,'components':12,'new_time_candidates':0,'kernel_calls':0}))


if __name__=='__main__':main()
