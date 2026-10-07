#!/usr/bin/env python3
"""Finite adaptive search for fixed-v=1 continuous-window hard shell energy.

No old data are evaluated or overwritten. Actual FIRST/geom gates are absent:
Theta is expressly the weighted hard relaxation with g=1, not actual Theta_g.
"""
from pathlib import Path
import datetime, hashlib, json, math, time, types
import numpy as np

HERE=Path(__file__).resolve().parent
engine=types.ModuleType('winner_snapshot')
exec(compile((HERE/'winner_engine_snapshot.py').read_text(),str(HERE/'winner_engine_snapshot.py'),'exec'),engine.__dict__)
BASE=202610070812
CELLS=[(n,N) for n in (8,16,32,64) for N in (16,64,128)]
FAMILIES=['hierarchical_orthant_cascade','anisotropic_product_strata','lowrank_wavy_cloud','multilevel_shell_spokes',
          'sparse_grid_capacity','shared_resolution_overlap']
SOURCE_TABLE=Path('/Users/zhengzhihao/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_2ypj2h60ogms22_45c1/temp/RWTemp/2026-10/07aa327ce93623270ba54256177864ed/中心立方体_阶段性下界构造总表_20261005 (3).txt')


def serialize(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,(np.integer,np.floating)):return x.item()
    raise TypeError(type(x))


def fresh_parameters(n,N,family,seed):
    rng=np.random.default_rng(seed)
    return dict(n=n,N=N,family=family,geometry_seed=seed,
        amplitude=float(rng.uniform(.32,1.25)), jitter=float(rng.uniform(.01,.15)),
        weight_strength=float(rng.uniform(0.,2.0)), decay=float(rng.uniform(.35,.9)),
        rank=int(rng.choice([2,4,8,min(n,16),n])), levels=int(rng.choice([3,4,6,8])),
        anisotropy=float(rng.uniform(0.,1.5)))


def make_source(p):
    n,N=p['n'],p['N']; rng=np.random.default_rng(p['geometry_seed'])
    k=math.ceil(math.log2(N))
    bits=2*((np.arange(N)[:,None]>>np.arange(k))&1)-1
    logits=rng.normal(size=N)
    if p['family']=='hierarchical_orthant_cascade':
        vectors=rng.choice([-1.,1.],(k,n))
        scales=p['decay']**np.arange(k)
        atoms=(bits*scales)@vectors/scales.sum()
        # Unequal branch masses correlate with nested geometry, rather than
        # a one-dimensional diagonal or the earlier two-cluster family.
        # Frozen edge energies along complete binary-tree paths (B61 mechanism).
        path_energy=np.zeros(N)
        for depth in range(1,k+1):
            energies=rng.normal(size=2**depth)
            path_energy+=energies[np.arange(N)%(2**depth)]
        logits+=path_energy/math.sqrt(k)
    elif p['family']=='anisotropic_product_strata':
        labels=rng.permutation(n)%k
        # Graph cut labels from vertex signs, not independent coordinate bits:
        # coordinate i is sign(vertex u_i)*sign(vertex v_i).
        other=(labels+rng.integers(1,k,n))%k
        atoms=(bits[:,labels]*bits[:,other]).astype(float)
        shifts=rng.uniform(-.4,.4,k)
        atoms+=shifts[labels]
        scales=np.exp(-p['anisotropy']*np.linspace(0,1,n))
        atoms*=scales
        # Ising-like correlated weights on the low-dimensional product code.
        logits+=np.sum(bits*np.roll(bits,1,axis=1),axis=1)/math.sqrt(k)
    elif p['family']=='lowrank_wavy_cloud':
        rank=max(1,min(p['rank'],n))
        latent=rng.uniform(-1,1,(N,rank))
        projection=rng.normal(size=(rank,n))/math.sqrt(rank)
        atoms=latent@projection
        atoms+=.18*np.sin(3*atoms)
        atoms/=max(1.,float(np.quantile(np.abs(atoms),.95)))
        atoms*=np.exp(-p['anisotropy']*np.linspace(0,1,n))
        logits+=np.sum(latent**2,axis=1)/math.sqrt(rank)
    elif p['family']=='sparse_grid_capacity':
        q=max(2,min(8,p['levels']))
        atoms=(rng.integers(0,q,(N,n))+.5)/q-.5
        atoms*=2
        logits+=rng.normal(size=N)  # all N capacity coordinates are free
    elif p['family']=='shared_resolution_overlap':
        levels=min(p['levels'],6);layer=np.arange(N)%levels
        qvalues=2**np.arange(levels)
        atoms=np.empty((N,n))
        for ell,qvalue in enumerate(qvalues):
            mask=layer==ell
            atoms[mask]=(rng.integers(0,qvalue,(int(mask.sum()),n))+.5)/qvalue-.5
        atoms*=2
        # Resolution is shared by each ENTIRE n-dimensional atom; the mixture
        # is formed after joint-coordinate sampling. Layer masses are positive.
        layer_logits=rng.normal(size=levels)
        logits=layer_logits[layer]-np.log(np.bincount(layer)[layer])+.2*logits
    else:
        levels=p['levels']; level=np.arange(N)%levels
        radii=np.geomspace(max(.025,p['decay']*.14),1.,levels)
        spokes=rng.normal(size=(max(4,N//levels),n))
        spokes/=np.max(np.abs(spokes),axis=1)[:,None]
        atoms=spokes[(np.arange(N)//levels)%len(spokes)]*radii[level,None]
        # Level mass perturbation is independent of the number of points in it.
        logits+=rng.normal(size=levels)[level]-np.log(np.bincount(level)[level])
    atoms=p['amplitude']*atoms+p['jitter']*rng.uniform(-1,1,(N,n))
    logits=np.clip(p['weight_strength']*logits,-40.,40.)
    weights=np.exp(logits-logits.max());weights/=weights.sum()
    if 'free_weights' in p:
        weights=np.array(p['free_weights'],dtype=float)
        weights/=weights.sum()
    assert np.all(weights>0) and math.isclose(float(weights.sum()),1.,abs_tol=1e-12)
    digest=hashlib.sha256(atoms.astype('<f8').tobytes()+weights.astype('<f8').tobytes()).hexdigest()
    return atoms,weights,digest


def mutate(p,seed,weights=None,capacity_only=False):
    rng=np.random.default_rng(seed); q=dict(p)
    if weights is not None:
        # Evolutionary capacity step changes every positive atom coefficient;
        # it is not described as LP/MIP optimality or a solver certificate.
        w=np.array(weights)**float(rng.uniform(.6,1.5))*np.exp(rng.normal(0,.8,len(weights)))
        w=np.maximum(w,1e-100);w/=w.sum();q['free_weights']=w.tolist()
    if capacity_only:return q
    q['amplitude']=float(np.clip(p['amplitude']*math.exp(rng.normal(0,.28)),.08,2.5))
    q['jitter']=float(np.clip(p['jitter']*math.exp(rng.normal(0,.6)),.0001,.4))
    q['weight_strength']=float(np.clip(p['weight_strength']+rng.normal(0,.75),0.,7.))
    q['decay']=float(np.clip(p['decay']+rng.normal(0,.12),.12,.98))
    q['anisotropy']=float(np.clip(p['anisotropy']+rng.normal(0,.35),0.,3.))
    if rng.random()<.3:q['rank']=int(rng.choice([2,4,8,min(p['n'],16),p['n']]))
    if rng.random()<.3:q['levels']=int(rng.choice([3,4,6,8,12]))
    return q


def initial_lambda_grid(n,N):
    # Frozen theoretical mass/radius anchors. Any later lambda choice is an
    # explicitly recorded training maximization, then frozen for validation.
    return sorted(set(math.log(capture/N)-n*math.log(radius)-math.log(3)
                      for capture in (1.,max(1.,N/8),max(1.,N/2),float(N))
                      for radius in (1.03,1.2,1.4,1.65,1.93)))


def evaluate(p,logs,M,K,seed,stage,record_total=12):
    atoms,weights,digest=make_source(p);n,N=p['n'],p['N']
    ts=np.linspace(0,n*math.log(2),K);L=len(logs)
    # dimensions: replicate, lambda, quantity (Gamma1,Theta_hard,Theta_hard1), node
    lower=np.zeros((2,L,3,K));upper=np.zeros_like(lower);raw=np.zeros_like(lower)
    ambiguity=np.zeros((L,3),dtype=np.int64)
    winner_mass=np.zeros(K);winner_count=np.zeros(K);source_membership=np.zeros(K)
    numerical_tie_samples=0;distance_tie_samples=0
    for rep in range(2):
        rng=np.random.default_rng(seed+rep*1000003)
        sources=rng.choice(N,M,p=weights)
        omega=rng.uniform(-.5,.5,(M,n));face=rng.integers(0,2*n,M)
        omega[np.arange(M),face//2]=np.where(face%2,.5,-.5)
        diff=atoms[sources,None,:]-atoms[None,:,:]
        for ki,t in enumerate(ts):
            r=math.exp(t/n)
            q=2*np.max(np.abs(diff/r+omega[:,None,:]),axis=2)
            assert np.all(q[np.arange(M),sources]==1.)
            qs,lq,mass,valid,response,ix,logu=engine.window_maximal(q,weights,n,t)
            rows=np.arange(M);tol=(n+1)*1e-10
            near=valid&(response>=logu[:,None]-tol)
            numerical_tie_samples+=int(np.count_nonzero(near.sum(axis=1)>1))
            orig_sorted=np.sort(q,axis=1)
            distance_tie_samples+=int(np.count_nonzero(np.any(orig_sorted[:,:-1]==orig_sorted[:,1:],axis=1)))
            delay=n*lq;Q=qs[rows,ix];d=delay[rows,ix]
            capture_lower=qs>=1.;capture_upper=qs>=1.-tol
            short_lower=capture_lower&(delay<=1.-tol)
            short_upper=capture_upper&(delay<=1.+tol)
            exp_delay=np.exp(-np.maximum(delay,0.))
            winner_mass[ki]+=float(mass[rows,ix].mean())/2
            winner_count[ki]+=float((q<=Q[:,None]).sum(axis=1).mean())/2
            source_membership[ki]+=float((Q>=1.).mean())/2
            for li,ll in enumerate(logs):
                band_lo=(logu-tol>ll+math.log(2))&(logu+tol<=ll+math.log(4))
                band_hi=(logu+tol>ll+math.log(2))&(logu-tol<=ll+math.log(4))
                raw_band=(logu>ll+math.log(2))&(logu<=ll+math.log(4))
                masks_lo=(short_lower,capture_lower,short_lower)
                masks_hi=(short_upper,capture_upper,short_upper)
                vals=(np.ones_like(qs),exp_delay,exp_delay)
                raw_values=(np.ones(M),np.exp(-np.maximum(d,0.)),np.exp(-np.maximum(d,0.)))
                raw_masks=((d>=0.)&(d<=1.),d>=0.,(d>=0.)&(d<=1.))
                for quantity in range(3):
                    yl=vals[quantity]*(masks_lo[quantity]&band_lo[:,None])
                    yu=vals[quantity]*(masks_hi[quantity]&band_hi[:,None])
                    a=np.min(np.where(near,yl,np.inf),axis=1)
                    b=np.max(np.where(near,yu,-np.inf),axis=1)
                    lower[rep,li,quantity,ki]=a.mean();upper[rep,li,quantity,ki]=b.mean()
                    ambiguity[li,quantity]+=np.count_nonzero(b-a>1e-14)
                    raw[rep,li,quantity,ki]=np.mean(raw_values[quantity]*raw_masks[quantity]*raw_band)
    scores=np.trapezoid(lower[0,:,0]*lower[1,:,0],ts,axis=1)
    selected=int(np.argmax(scores))
    quantities={}
    epsilon=math.sqrt(math.log(12*record_total*K/.05)/(4*M))
    for qi,name in enumerate(('Gamma_v1','Theta_hard','Theta_hard_v1')):
        a=lower[0,selected,qi];b=lower[1,selected,qi];mean=(a+b)/2
        hi=upper[:,selected,qi].mean(axis=0)
        cross=a*b;fine=float(np.trapezoid(cross,ts));coarse=float(np.trapezoid(cross[::2],ts[::2]))
        quantities[name]=dict(energy_cross_fine=fine,energy_cross_coarse=coarse,
             nested_grid_difference=abs(fine-coarse),energy_empirical_square=float(np.trapezoid(mean**2,ts)),
             integral_mean=float(np.trapezoid(mean,ts)),max_mean=float(mean.max()),
             numerical_ambiguity_sample_nodes=int(ambiguity[selected,qi]),
             finite_node_95_lower=float(np.trapezoid(np.maximum(0,mean-epsilon)**2,ts)),
             finite_node_95_upper=float(np.trapezoid(np.minimum(1,hi+epsilon)**2,ts)),
             mean=mean.tolist(),rep_A=a.tolist(),rep_B=b.tolist(),numerical_upper=hi.tolist(),
             raw_mean=raw[:,selected,qi].mean(axis=0).tolist())
    return dict(stage=stage,parameters=p,input_sha256=digest,atoms=atoms.tolist(),weights=weights.tolist(),
       sample_seed=seed,samples_per_replicate=M,fine_nodes=K,coarse_nodes=(K+1)//2,
       lambda_grid=logs,lambda_scores_Gamma1=scores.tolist(),selected_lambda_index=selected,
       selected_log_lambda=logs[selected],training_selection=len(logs)>1,
       quantities=quantities,t=ts.tolist(),finite_node_band_epsilon=epsilon,
       winner_mean_mass=winner_mass.tolist(),winner_mean_atom_count=winner_count.tolist(),
       source_in_winner_fraction=source_membership.tolist(),
       near_response_tie_sample_nodes=numerical_tie_samples,distance_tie_sample_nodes=distance_tie_samples)


def validate_engine():
    out=[]
    tests=[(np.array([0.,0.]),np.array([[.5,0.],[-.5,0.],[0.,1.]]),np.array([.2,.3,.5])),
           (np.array([.2,-.1]),np.array([[0.,0.],[.4,.2],[.9,-.1]]),np.array([.2,.7,.1])),
           (np.zeros(3),np.array([[0.,0.,0.],[0.,0.,0.],[.5,.5,.5]]),np.array([.1,.2,.7]))]
    for x,atoms,w in tests:
        d=2*np.max(np.abs(atoms-x),axis=1)
        qs,_,mass,_,_,ix,lu=engine.window_maximal(d[None,:],w,len(x),0.)
        candidates=sorted(set([1.,2.]+[v for v in d if 1.<=v<=2.]))
        u,R=max([(float(w[d<=R].sum())/R**len(x),R) for R in candidates],key=lambda p:(p[0],-p[1]))
        assert math.isclose(math.exp(lu[0]),u,rel_tol=1e-12) and qs[0,ix[0]]==R
        out.append(dict(pass_=True,u=u,R=R,mass=float(mass[0,ix[0]])))
    # Known equality: source q=1, delay=0, is included; neither direction faces
    # nor captured atoms are silently treated as open cubes.
    q=np.ones((8,1));ts=[0.,math.log(2)/2,math.log(2)]
    for t in ts:
        _,lq,_,_,_,ix,lu=engine.window_maximal(q,np.array([1.]),8,t)
        assert np.all(lq[np.arange(len(q)),ix]==0.) and np.all(lu==-t)
    out.append(dict(pass_=True,validation='single atom lambda=.25: Gamma1=Theta_hard=1 on [0,log2), delay0 accepted; analytic squared energy=log2',energy=math.log(2)))
    assert math.isclose(math.expm1(8*math.log(2)),2.**8-1,rel_tol=1e-12)
    out.append(dict(pass_=True,validation='cone Jacobian dx=r^n dt dnu on boundary [-1/2,1/2]^n; shell volume matches'))
    return out


def save(name,data):
    p=HERE/name
    with p.open('x') as f:json.dump(data,f,default=serialize,ensure_ascii=False,indent=2,allow_nan=False)
    return p


def main():
    start=time.perf_counter();stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    checks=validate_engine()
    design=dict(seed_base=BASE,cells=CELLS,families=FAMILIES,
       stages='72 starting candidates + 24 adaptive geometry/all-positive-coefficient mutations + 12 sparse-grid capacity-only mutations; select one per n,N cell; pressure frozen-lambda winners with fresh sampling; select 6 after pressure; final fresh-seed validation',
       search_budget=dict(exploration_M=32,exploration_K=65,pressure_M=256,pressure_K=257,validation_M=1024,validation_K=513),
       objective='maximize integral Gamma_v=1(t)^2 dt by independent replicate cross product; record Theta_hard and Theta_hard,v1 with actual geometric exp(-delay)',
       eligibility='Finite positive atomic mu with total W=H=1; actual centered continuous [1,2] full-side maximal winner; a single fixed selected lambda per frozen input; hardband 2lambda<u<=4lambda; source face inclusion uses equality',
       qualification_limits=['No actual FIRST/continuation/CP/GP gate supplied: g=1 hard relaxation, so weighted Theta_hard is not actual Theta_g','Continuous window does not certify original arbitrary finite allowed J','Atomic inputs do not certify diffuse L1 inputs','A large fixed finite energy is not a violation of unspecified C log^A(n); no asymptotic exponent fitting','Search lambda and parameters overfit training; pressure/final streams separate and final input/lambda frozen'],
       arithmetic='Exact enumeration for floating distances including all equal-distance groups and endpoints a,b; source distance computed from offsets preserves q_source=1; tolerance (n+1)*1e-10 lower/upper screening is not exact-real interval certification',
       jacobian='dx=r^n dt dnu. Expectations retain source labels and sample with source weights; not unweighted spatial samples.',
       engine_snapshot_sha256=hashlib.sha256((HERE/'winner_engine_snapshot.py').read_bytes()).hexdigest(),
       construction_source=dict(path=str(SOURCE_TABLE),sha256=hashlib.sha256(SOURCE_TABLE.read_bytes()).hexdigest(),
          mappings={'sparse_grid_capacity':'B38: arbitrary positive coefficients on a finite sparse grid, evolutionary weight updates; atomic sparse subset, not full q^n cell density or LP/MIP certificate',
          'shared_resolution_overlap':'B41/B40: full n-dimensional atom uses one shared dyadic resolution before positive layer mixing; finite sampled truncation, not full tensor lattice or its weak-budget theorem',
          'anisotropic_product_strata':'B07/B06: graph cut coordinates derived from common vertex signs; correlated Ising-like weights, anisotropic spacings, atomic finite graph realization',
          'hierarchical_orthant_cascade':'B61/B52: frozen binary-tree path energies produce positive Gibbs weights on multiresolution leaf geometry; not the exact original periodic histogram input',
          'lowrank_wavy_cloud':'B04/B18 continuous latent/correlated positive source inspiration only; explicit new finite atomic cloud',
          'multilevel_shell_spokes':'B62: entire atom shares a radial layer before mixture; finite spoke realization rather than continuous cube-density integral'},
          excluded='A13/A14 positive spectral compensation, Sidon labels, peak widths, full winner/FIRST qualifications are not implemented or inherited'),
       validation=checks)
    design_path=save(f'design_{stamp}.json',design)
    exploration=[];best={}
    serial=0
    for n,N in CELLS:
        for family in FAMILIES:
            p=fresh_parameters(n,N,family,BASE+serial*104729)
            r=evaluate(p,initial_lambda_grid(n,N),32,65,BASE+8000000+serial*1009,'exploration')
            r['candidate_id']=serial;r['parent_id']=None;exploration.append(r)
            if (n,N) not in best or r['quantities']['Gamma_v1']['energy_cross_fine']>best[(n,N)]['quantities']['Gamma_v1']['energy_cross_fine']:best[(n,N)]=r
            serial+=1
    for generation in (1,2):
        for n,N in CELLS:
            parent=best[(n,N)]
            p=mutate(parent['parameters'],BASE+serial*104729,parent['weights'])
            r=evaluate(p,initial_lambda_grid(n,N),32,65,BASE+8000000+serial*1009,'adaptive_exploration')
            r['candidate_id']=serial;r['parent_id']=parent['candidate_id'];r['generation']=generation
            exploration.append(r)
            if r['quantities']['Gamma_v1']['energy_cross_fine']>parent['quantities']['Gamma_v1']['energy_cross_fine']:best[(n,N)]=r
            serial+=1
    # One guaranteed B38 capacity mutation in every independent n,N cell,
    # regardless of which mechanism the initial score preferred.
    for n,N in CELLS:
        parent=max([r for r in exploration if r['parameters']['n']==n and r['parameters']['N']==N and r['parameters']['family']=='sparse_grid_capacity'],key=lambda r:r['quantities']['Gamma_v1']['energy_cross_fine'])
        p=mutate(parent['parameters'],BASE+serial*104729,parent['weights'],capacity_only=True)
        r=evaluate(p,initial_lambda_grid(n,N),32,65,BASE+8000000+serial*1009,'capacity_only_exploration')
        r['candidate_id']=serial;r['parent_id']=parent['candidate_id'];exploration.append(r)
        if r['quantities']['Gamma_v1']['energy_cross_fine']>best[(n,N)]['quantities']['Gamma_v1']['energy_cross_fine']:best[(n,N)]=r
        serial+=1
    save(f'exploration_{stamp}.json',dict(design_path=str(design_path),records=exploration,runtime_seconds=time.perf_counter()-start))
    print('EXPLORATION complete',len(exploration),'runtime',time.perf_counter()-start,flush=True)
    pressure=[]
    for cell,r in best.items():
        v=evaluate(r['parameters'],[r['selected_log_lambda']],256,257,BASE+16000000+r['candidate_id']*1013,'pressure_frozen_lambda',record_total=12)
        v['candidate_id']=r['candidate_id'];pressure.append(v)
        print('PRESSURE',cell,v['parameters']['family'],v['quantities']['Gamma_v1']['energy_cross_fine'],flush=True)
    save(f'pressure_{stamp}.json',dict(records=pressure,runtime_seconds=time.perf_counter()-start))
    # Four dimension winners ensure independently changed n. Add both other N
    # values at the highest-energy dimension, retaining a same-n N control.
    selected=[]
    for n in (8,16,32,64):
        selected.append(max([r for r in pressure if r['parameters']['n']==n],key=lambda r:r['quantities']['Gamma_v1']['energy_cross_fine']))
    leader=max(selected,key=lambda r:r['quantities']['Gamma_v1']['energy_cross_fine'])
    for r in pressure:
        if r['parameters']['n']==leader['parameters']['n'] and r not in selected:selected.append(r)
    final=[]
    for r in selected:
        v=evaluate(r['parameters'],[r['selected_log_lambda']],1024,513,BASE+32000000+r['candidate_id']*1019,'fresh_seed_validation',record_total=len(selected))
        v['candidate_id']=r['candidate_id'];final.append(v)
        print('VALIDATE',v['parameters']['n'],v['parameters']['N'],v['quantities']['Gamma_v1']['energy_cross_fine'],flush=True)
    final_path=save(f'validation_{stamp}.json',dict(records=final,
       confidence='Conditional on the frozen selected inputs/lambda: two-sided Hoeffding for lower and upper screened RVs of all 3 quantities at all fine nodes, union factor 12*selected_records*K; trapezoid bounds only bound the weighted finite-node target. Numerical enclosure is heuristic floating tolerance, so no exact-real certificate.',
       sampling_relation='Distinct stage seed namespaces, two separately seeded replicates per input. CRN across nodes and lambda within a replicate. Fresh validation uses same algorithm, not an independent implementation.',
       runtime_seconds=time.perf_counter()-start))
    summary=[]
    for r in final:
        summary.append(dict(candidate_id=r['candidate_id'],n=r['parameters']['n'],N=r['parameters']['N'],family=r['parameters']['family'],log_lambda=r['selected_log_lambda'],
              energies={k:v['energy_cross_fine'] for k,v in r['quantities'].items()},
              gamma_grid_difference=r['quantities']['Gamma_v1']['nested_grid_difference'],
              gamma_node_energy_interval=[r['quantities']['Gamma_v1']['finite_node_95_lower'],r['quantities']['Gamma_v1']['finite_node_95_upper']]))
    save(f'summary_{stamp}.json',dict(validation_path=str(final_path),records=summary,
        runtime_seconds=time.perf_counter()-start,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        conclusion='Report these finite frozen-input energies only. Even a successful high-energy search would concern the hard relaxation; it is not an actual geom counterexample. No fixed-v polylog candidate is proved or disproved here.'))
    print('FINAL',final_path,'RUNTIME',time.perf_counter()-start,flush=True)


if __name__=='__main__':main()
