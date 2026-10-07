from pathlib import Path
import json,math,hashlib
import numpy as np
D=Path(__file__).resolve().parent
runs=[json.loads((D/f'round{j}.json').read_text()) for j in [1,2,3]]
err=json.loads((D/'root_error_brackets.json').read_text())['records'];cross=json.loads((D/'future_average_crosscheck.json').read_text())['records']
lines=['# Frozen tensor soft future-cap diagnostics','',
'Four selected receivers reuse the complete old tensor measure, old lambda and old A-stream seeds. They are selected from hard-band/nonconcentrated short-shell rows before looking at soft responses. No population frequency or dimensional order is inferred.', '',
'| n / saved pair | round | all L nodes | visible raw sigma | K/q at finite candidate | K/q scalar bounds* |',
'|---|---:|---:|---:|---:|---:|']
for run in runs:
    for row in run['rows']:
        k=row['future_average'];rt=row['finite_visible_FIRST'];pair=row['point_id'].split('_')[1]
        lines.append(f"| {row['n']} / {pair} | {run['round']} | {row['L_nodes']} | {rt['sigma_raw']:.9g} | {k['normalized_K_mid']:.9g} | [{k['normalized_K_lower']:.9g}, {k['normalized_K_upper']:.9g}] |")
lines+=['', '*The scalar quadrature remainder is analytic, but the floating function evaluation and guard are not outward certified. These are numerical bounds, not formal machine intervals.', '',
'All original coordinate arrival radii remain closed-cube cuts. Every jump has a saved left/right response trace. Targeted rounds retain prior cuts and add up to64 geometric L midpoints per point, besides the nested logL and sigma grids. The raw bisection root is only the last visible crossing of the finite physical-node maximum, not a true last-crossing boundary.','',
'Possible full last-crossing sigma corridors (round3):','',
'| n / pair | last clear finite-node violation lower | full continuous numerical cap start upper | finite K / theory epsilon_n | root-L hard jump |',
'|---|---:|---:|---:|---:|']
receipts=[]
for row in runs[-1]['rows']:
    pair=row['point_id'].split('_')[1];passed=next((q for q in row['future_cap_trials'] if q['status'].startswith('NUMERICALLY_ENCLOSED')),None)
    low=row['last_clear_finite_sigma_violation'];hi=passed['sigma_start'] if passed else None
    k=row['future_average'];jump=next(z for z in cross if z['point_id']==row['point_id'])
    lines.append(f"| {row['n']} / {pair} | {low:.9g} | {hi:.9g} | {k['ratio_K_mid_to_epsilon_n']:.6g} | {jump['root_L_jump_size']:.6g} |")
    receipts.append(dict(point_id=row['point_id'],n=row['n'],possible_full_last_crossing_sigma_lower=low,possible_full_last_crossing_sigma_upper=hi,scope='Numerical corridor only: lower from a finite node with response scalar-lower>q, upper from exhaustive continuous sigma/L cell numerical enclosure. Floating functions are not outward certified. Raw bisection root is display-only.',K_full_sigma_L_candidate_box='UNKNOWN; actual soft physical winner Ls and actual FIRST parameter are not certified. The candidate-pair K above epsilon does not delete actual original output.',future_average=k,root_error_brackets=[z for z in err if z['point_id']==row['point_id']],final_numeric_cap_trial=passed))
lines+=['', 'The final numerical corridors lie above1/n and S0=1/(16sqrt n). This supplies middle-softness candidates for the original shared-input necessary cap layer, without proving actual FIRST/GOOD/score/CP/GP/continuation qualification. The numerical fullcap upper is at a sigma greater than the finite visible root; it is not equality P=q.', '',
'K is the original weighted future average with alpha=ceil sqrt n. Positive Bernstein integration eliminates v quadrature: the field is a polynomial of degree n, and weights are alpha/(n+alpha) product_{j=1}^{alpha-1}(k+j)/(n+j). Independent Gauss integration at a polynomial-exact order differs by at most4.94e-15 on the final four candidate pairs.', '',
'At the four finite candidate pairs, scalar-lower K/q exceeds the prechosen theory threshold epsilon_n=log(n+2)^(-4). The n128 values also exceed1/16; n512 values do not. This is not evidence of a remaining low-future-average branch. K over the entire possible sigma/L winner box remains UNKNOWN: actual Ls and actual FIRST are not certified, so no actual original output is declared deleted.', '',
'Finite-grid future maxima alone only exclude starts when a node is above q; all nodes below q cannot establish continuous cap. The separate cell oracle uses a positive coordinate endpoint envelope in sigma and (u/l)^n P_{sigma,u} in L, retains every hard boundary, and subdivides unresolved sigma cells. Failure or budget remainder stays UNKNOWN. Completed numerical envelopes cover every sigma/L cell, but machine rounding prevents promoting them to continuous certificates.', '',
'No old hard energy simulation was rerun. All old source/round/NPZ files were read-only. Three soft rounds, amendments, source coordinates/seeds, complete response grids and jump traces are saved. There is no new weak-type fee or fitted exponent.']
(D/'summary.md').write_text('\n'.join(lines)+'\n')
(D/'completion_receipt.json').write_text(json.dumps(dict(status='THREE_ROUNDS_COMPLETE_SOFT_NECESSARY_CAP_NUMERICAL_DIAGNOSTICS_ONLY',elapsed_seconds=[z['elapsed_seconds'] for z in runs],final_receivers=receipts,validation=json.loads((D/'validation.json').read_text())['status'],endpoint_guard=json.loads((D/'cap_endpoint_validation.json').read_text())['status'],future_average_crosscheck_max_error=max(z['abs_difference'] for z in cross),qualification='No certified actual FIRST, no actual Ls selection, no outward floating intervals. Full candidate-box K remains unknown, so candidate-pair K does not certify original output deletion.',artifacts=[dict(name=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(D.glob('*')) if p.is_file() and p.name!='completion_receipt.json']),indent=2))
print('\n'.join(lines))
