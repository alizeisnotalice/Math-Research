"""New independent controls; imports frozen engine without changing main rounds."""
from pathlib import Path
import importlib.util,json,hashlib,math,time
import numpy as np
D=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('probe',D/'probe.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
PLAN=[(128,32,513),(512,16,1025)]
PRE=D/'highdim_preregistration.json'
assert not PRE.exists(),'Control preregistration already exists; refuse rerun.'
engine_sha=hashlib.sha256((D/'probe.py').read_bytes()).hexdigest()
assert engine_sha==json.loads((D/'preregistration.json').read_text())['script_sha256']
PRE.write_text(json.dumps(dict(status='FROZEN_BEFORE_HIGH_DIMENSION_CONTROL',controls=[dict(n=n,M=M,K=K,pairs=2,delta_t=n*math.log(2)/(K-1),seeds=[[202610087000+n*100+rep*2,202610087001+n*100+rep*2] for rep in range(2)]) for n,M,K in PLAN],full_source_formula='Unchanged complete three-component mu with L64,q=[1,2,4],p=[.25,.5,.25],W=1',lambda_formula='loglambda=-nlog64',eta_formula='(49/65536)/sqrt(n)',vstar=p.VSTAR,engine_sha256=engine_sha,qualification='No FIRST/actualgeom. Independent additional control tags; no original round rerun; no optimization or exponent fit.',uncertainty='Two independent pair finite-node cross energies. Approx t1=12.706204736 intervals and marginal Hoeffding bounds are broad. Gamma1 short-shell width v1 and hardband discontinuities may alias time grids. Coarse-fine difference is diagnostic, not quadrature bound. At n512 diagnostic exp(logm) may underflow, while scoring and microbox ratios stay logarithmic.'),indent=2))
start=time.time();rows=[]
for n,M,K in PLAN:
    tag=f'highdim_n{n}';assert not (D/f'{tag}.json').exists()
    reps=[];energies=[]
    for rep in range(2):
        seed=202610087000+n*100+rep*2
        t,a,fa,sa,da=p.profile(n,64,M,K,seed)
        _,b,fb,sb,db=p.profile(n,64,M,K,seed+1)
        e=np.trapezoid(a*b,t,axis=-1);coarse=np.trapezoid((a*b)[:,:,::2],t[::2],axis=-1)
        np.savez_compressed(D/f'{tag}_pair{rep}.npz',t=t,profile_A=a,profile_B=b,fixed_eta0_A=fa,fixed_eta0_B=fb)
        reps.append(dict(rep=rep,draw_A=da,draw_B=db,energy=e.tolist(),coarse_energy=coarse.tolist(),screen_A=sa,screen_B=sb,fixed_eta0_diagnostic_energy=np.trapezoid(fa*fb,t,axis=-1).tolist()))
        energies.append(e)
        (D/f'{tag}_partial.json').write_text(json.dumps(reps,indent=2))
        print(json.dumps(dict(tag=tag,rep=rep,cert_energies=e[:,0].tolist(),fullhard_upper=e[:,2].tolist(),elapsed_seconds=time.time()-start)),flush=True)
    ee=np.array(energies);mean=ee.mean(axis=0);sd=ee.std(axis=0,ddof=1);se=sd/math.sqrt(2);T=n*math.log(2);ho=T*math.sqrt(math.log(40)/4)
    screens=[s for rr in reps for key in ['screen_A','screen_B'] for s in rr[key]]
    row=dict(tag=tag,n=n,L=64,q=p.Q.tolist(),p=p.P.tolist(),M=M,K=K,pairs=2,delta_t=T/(K-1),log_lambda=-n*math.log(64),eta=p.ETA0/math.sqrt(n),eta0=p.ETA0,vstar=p.VSTAR,
       distinct_atom_log10_count=n*math.log10(256),quantity_order=['Gamma_v1','Theta_hard','Theta_hard_v1','Gamma_vstar'],profile_axis=['screened_lower_cert_res','raw_cert_res','floating_upper_including_unknown'],energy_mean=mean.tolist(),energy_sd=sd.tolist(),energy_se=se.tolist(),approx95_t1_interval=np.stack([np.maximum(0,mean-12.706204736*se),mean+12.706204736*se],axis=-1).tolist(),hoeffding95_grid_interval=np.stack([np.maximum(0,mean-ho),np.minimum(T,mean+ho)],axis=-1).tolist(),
       cert_fraction_all_queries=sum(s['certified'] for s in screens)/(len(screens)*M),unknown_fraction_all_queries=sum(s['unknown'] for s in screens)/(len(screens)*M),fixed_eta0_cert_fraction_all_queries=sum(s['cert_fixed'] for s in screens)/(len(screens)*M),
       max_nested_difference_per_quantity_axis=np.max([abs(np.array(rr['energy'])-np.array(rr['coarse_energy'])) for rr in reps],axis=0).tolist(),near_band_count=sum(s['near_band'] for s in screens),near_score_tie_count=sum(s['near_score_tie'] for s in screens),mass_diagnostic_zero_sample_nodes=sum(s['mean_mass']==0 for s in screens),
       replicates=reps,elapsed_seconds=time.time()-start,uncertainty_scope='Finite-node statistic only. Two pairs; approximate t1 intervals do not guarantee coverage. Coarse/fine differences do not bound continuum quadrature or Gamma1 aliasing. No outward floating intervals or actualgeom gates.')
    (D/f'{tag}.json').write_text(json.dumps(row,indent=2));rows.append(row)
lines=['# Independent high-dimensional controls','', 'The original three rounds were not rerun. Source formula, all three overlapping tensor components, fixed lambda, eta_n and vstar are unchanged; controls use independent registered seeds. No dimension order is fitted.','', '| n | M / side | nodes | Δt | Γ1 screened cert | Θhard screened cert | Γvstar screened cert | Γvstar fullhard upper | cert fraction all queries |', '|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for row in rows:
    e=np.array(row['energy_mean']);lines.append(f"| {row['n']} | {row['M']} | {row['K']} | {row['delta_t']:.6g} | {e[0,0]:.6g} | {e[1,0]:.6g} | {e[3,0]:.6g} | {e[3,2]:.6g} | {row['cert_fraction_all_queries']:.6g} |")
lines+=['', 'Screened cert is a subset passing the sufficient microbox count bound for the exact nonconcentrated residual, conditional on computed winner/counts. Failed certificates remain unknown. Fullhard upper omits the microbox gate and is not cooperative residual energy.','',
'Each energy mean averages two independent A/B cross-product pair estimates. SE is the sample SD divided by sqrt2. Stored approximate t1 intervals and marginal bounded Hoeffding intervals address finite-node MC statistics only; two pairs give broad, fragile intervals. Γ1 and hardband discontinuities may alias the finite time grids, especially Δt=.34657 at n512. Stored coarse/fine differences diagnose sensitivity without giving a deterministic quadrature error bound. General floating comparisons are screened, without outward interval arithmetic.','',
'At n512, diagnostic exponentiation of winner log mass can underflow to zero; winner scores, lambda comparisons and microbox share computations remain in log space. That mean-mass diagnostic is not a numerical zero-mass claim.','',
'Full control JSON contains per-node screens, sources/seeds/hashes and all quantity/axis energies. Four NPZ files preserve full A/B profiles.']
(D/'highdim_summary.md').write_text('\n'.join(lines)+'\n')
# Independent readback of saved profiles.
maxerr=0.
for row in rows:
    for rr in row['replicates']:
        data=np.load(D/f"{row['tag']}_pair{rr['rep']}.npz")
        calc=np.trapezoid(data['profile_A']*data['profile_B'],data['t'],axis=-1)
        maxerr=max(maxerr,float(abs(calc-np.array(rr['energy'])).max()))
(D/'highdim_completion_receipt.json').write_text(json.dumps(dict(status='INDEPENDENT_HIGH_DIMENSION_CONTROLS_COMPLETE',elapsed_seconds=time.time()-start,NPZ_readback_max_energy_error=maxerr,preregistration_sha256=hashlib.sha256(PRE.read_bytes()).hexdigest(),artifacts=[dict(name=x.name,sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in sorted(D.glob('highdim*')) if x.is_file() and x.name!='highdim_completion_receipt.json']),indent=2))
print('\n'.join(lines),flush=True)
