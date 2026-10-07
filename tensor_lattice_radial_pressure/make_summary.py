"""Summarize saved immutable rounds; no estimator or parameter reselection."""
import hashlib
import json
import math
from pathlib import Path

root=Path(__file__).resolve().parent
tags=['v2_round1','v2_round2','v2_round3','v2_L65537_control1']
data=[json.loads((root/(tag+'.json')).read_text()) for tag in tags]
lines=['# Tensor lattice radial pressure: completed numerical receipt','',
       'All three main rounds and the single pre-registered L control completed. This is a floating hardband g1 relaxation of a source-centred radial energy; it does not certify original softFIRST geom, an endpoint estimate, a supplied complete lower-bound family or a dimension-uniform energy order.','',
       'Old round1/round2 are superseded because global x subtraction distorted exact cone-face r=R events. V2 uses local dx and resolves identifiable c=0 band threshold equalities by integer arithmetic. Numerical winner ties take smaller R only when computed floating scores are exactly equal; arbitrary mathematical/near score ties remain uncertified. Details and remaining limits are in README.md.','',
       '| Tag | n | L | c | Energy mean | Approx. marginal 95% MC interval | Max nested-grid difference |','|---|---:|---:|---:|---:|---|---:|']
receipts=[]
for d in data:
    for row in d['rows']:
        lo,hi=row['approx95_MC_t_interval']
        lines.append(f"| {d['round']} | {row['n']} | {row['L']} | {row['c']:g} | {row['energy_mean']:.8g} | [{lo:.7g}, {hi:.7g}] | {row['nested_grid_max_abs_difference']:.5g} |")
    for row in d['rows'][::4]:
        reps=row['replicate_data']; samples=d['protocol']['samples_per_side'];n=row['n'];L=row['L']
        boundary=sum(q[s]['boundary_sources'] for q in reps for s in ['screens_a','screens_b'])/(2*len(reps)*samples)
        total=lambda key:sum(q[s][key] for q in reps for s in ['screens_a','screens_b'])
        receipts.append({'tag':d['round'],'n':n,'L':L,'samples_per_side':samples,'nodes':d['protocol']['nodes'],'replicate_pairs':len(reps),'log10_total_atoms':n*math.log10(L),'boundary_exact_probability':-math.expm1(n*math.log1p(-2/L)),'boundary_observed_fraction':boundary,'symbolic_c0_ties_resolved':total('symbolic_c0_ties_resolved'),'floating_band_classifications_changed':total('c0_band_values_changed_by_exact_ties'),'unresolved_near_hardband':total('unresolved_near_hardband')})
lines+=['','| Tag | n | Boundary source probability | Observed fraction | c=0 floating classifications repaired | Unresolved near-hardband screens |','|---|---:|---:|---:|---:|---:|']
for r in receipts:
    lines.append(f"| {r['tag']} | {r['n']} | {r['boundary_exact_probability']:.7g} | {r['boundary_observed_fraction']:.7g} | {r['floating_band_classifications_changed']} | {r['unresolved_near_hardband']} |")
main=data[2];control=data[3]
lines+=['','The n=128 large-L control retains the same normalized thresholds and window while changing only the finite tensor lattice support size; it is a different complete source. Any energy change measures boundary/support sensitivity as well as MC/grid error. It does not establish an asymptotic law.','',
       'Each run has eight independent replicate pairs and A/B sides are independent within each pair. The control tag suffix 1 shares its seed namespace with v2_round1; these different-L/sample runs may be coupled and all runs are not mutually independent. The final v2_round3 has suffix 3 and different seeds from the control. Old/new repair comparisons intentionally share seeds.','',
       'The approximate t intervals have only eight independent replicate pairs and exclude quadrature and floating error. A zero estimate/interval is not proof of zero. Separate marginal bounded Hoeffding grid intervals, usually vacuous here, are saved in every result JSON. Coarse/fine differences are diagnostics for discontinuous gates, not certified continuum integration bounds. Every threshold is pre-set; all c values are retained. No favorable threshold is selected and no order in n is fitted.','',
       'Verification: 480 explicit all-atom low-dimensional cross-check cases, nine local source-face checks, and 600 independent interior second-arrival formula cases passed. The latter does not replace the full finite-lattice boundary computation. No positive-width L1 realization was executed.']
(root/'summary.md').write_text('\n'.join(lines)+'\n')
files=['pressure.py','pressure_v1_superseded.py','internal_formula_check.py','validation_v2.json','internal_formula_validation.json','control_preregistration.json','audit_clarifications.json','repair_receipt.json','README.md','summary.md']+[tag+'.json' for tag in tags]
manifest={'status':'COMPLETE_FLOAT_PRESSURE_THREE_MAIN_ROUNDS_PLUS_FIXED_L_CONTROL','source_relation':'Basic finite uniform tensor lattice only; not full A13/A14 or B60-B63 reproduction.','independence_scope':'A/B sides and replicate pairs are independent within each run; control1 and main round1 share seed namespace and may couple; final main round3 has different seeds from control.','numeric_winner_tie_scope':'Smaller R is chosen only for exactly equal computed floating scores; arbitrary mathematical and near score ties remain unverified.','main_dimensions':[8,32,128],'main_L':257,'control_L':65537,'superseded':['round1.json','round2.json'],'receipts':receipts,'sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files},'limits':['no original softFIRST geom certification','no endpoint or general energy theorem','no continuum quadrature bound','no interval-arithmetic floating certificate','no dimension-uniform energy order fit','no positive-width L1 transfer']}
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'status':manifest['status'],'receipts':receipts},indent=2))
