from pathlib import Path
import json,hashlib,math
import numpy as np
D=Path(__file__).resolve().parent
rounds=[json.loads((D/f'round{r}.json').read_text()) for r in [1,2,3]]
lines=['# Implicit overlapping tensor mixture: three frozen numerical rounds','',
'B41 shared-resolution mechanism on the full finite measure. No parameter or lambda selection; no actualgeom or FIRST gates. Values below are independent A/B cross-product finite-node energies. `cert` means the screened, count-bound-certified residual subset; `upper` retains every unresolved residual query and floating boundary ambiguity. These are numerical screens, not interval certificates.','',
'| round | n | samples / side | nodes | Gamma1 cert | Theta hard cert | Gamma vstar cert | Gamma vstar raw cert | Gamma vstar upper |',
'|---|---:|---:|---:|---:|---:|---:|---:|---:|']
receipt=[]
for rd in rounds:
    for row in rd['rows']:
        e=np.array(row['energy_mean']);reps=row['replicates'];M=rd['samples'];K=rd['nodes']
        screens=[z for r in reps for k in ['screen_A','screen_B'] for z in r[k]]
        cert=sum(z['certified'] for z in screens)/(len(screens)*M)
        unk=1-cert;fixed=sum(z['cert_fixed'] for z in screens)/(len(screens)*M)
        gap=max(abs(np.array(r['energy'])-np.array(r['coarse_energy'])).max() for r in reps)
        lines.append(f"| {rd['round']} | {row['n']} | {M} | {K} | {e[0,0]:.6g} | {e[1,0]:.6g} | {e[3,0]:.6g} | {e[3,1]:.6g} | {e[3,2]:.6g} |")
        receipt.append(dict(round=rd['round'],n=row['n'],eta=row['eta'],samples=M,nodes=K,cert_fraction_all_queries=cert,unknown_fraction_all_queries=unk,fixed_eta0_cert_fraction_all_queries=fixed,
            max_nested_difference_all_quantities_and_axes=float(gap),near_band_count=sum(z['near_band'] for z in screens),near_delay_count=sum(z['near_delay'] for z in screens),near_score_tie_count=sum(z['near_score_tie'] for z in screens),energy_mean=row['energy_mean'],energy_se=(np.array(row['energy_sd'])/2).tolist(),approx95_t3_interval=row['approx95_t3_interval']))
lines+=['', 'The final n=16/32/64 values are three separate fixed-dimensional inputs; no growth law is fitted. All rounds use the same source formula and frozen log(lambda)=-n log64, but different seeds, sample sizes and nested time grids. The upper column omits the microbox gate: it is a full hard-relaxation upper envelope that retains unknown queries, and does not prove a high residual energy. The certified residual is the exact complement of the box-existence event Eexists, conditional on computed winners/counts; no finite greedy-cover complement is used.', '',
'Approximate t3 intervals and bounded Hoeffding intervals for each finite grid statistic are stored in `round*.json`; four replicate pairs give little precision, and the bounded intervals are wide. Nested coarse/fine differences are discontinuous-integrand resolution diagnostics, without a deterministic continuum quadrature error bound. Atomic-to-L1 transfer and arbitrary allowed-scale or actualgeom qualification remain open.','',
'All source/cone draws use the complete probability measure. Distinct atom counts are 256^n, with log10 counts approximately 38.53,77.06,154.13. Overlapping source sites carry summed component masses. No sparse finite subsample replaces the measure when computing maximal averages.']
(D/'summary.md').write_text('\n'.join(lines)+'\n')
(D/'completion_receipt.json').write_text(json.dumps(dict(status='THREE_ROUNDS_COMPLETE_FLOAT_NUMERICAL_PRESSURE_ONLY',round_elapsed_seconds=[r['elapsed_seconds'] for r in rounds],protocol_sha256=hashlib.sha256((D/'preregistration.json').read_bytes()).hexdigest(),uncertainty_scope='Mean averages four independent A/B cross-product pair finite-grid statistics. SE=replicate SD/2; approximate t3 CI and marginal bounded Hoeffding cover finite-node MC targets only, excluding floating error, unknown residual, quadrature, actualgeom and L1 transfer.',validation=json.loads((D/'validation.json').read_text())['status'],microbox_validation=json.loads((D/'microbox_validation.json').read_text())['status'],screens=receipt,artifacts=[dict(name=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(D.glob('*')) if p.is_file() and p.name!='completion_receipt.json']),indent=2))
print('\n'.join(lines))
