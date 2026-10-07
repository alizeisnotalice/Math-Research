#!/usr/bin/env python3
"""Real finite-J centered cubes, full positive sources, deterministic S pressure.
Source-cone conditional IS restores receiver Lebesgue measure. No old imports.
"""
from pathlib import Path
import hashlib,json,time,sys
import numpy as np

BASE=Path(__file__).resolve().parent
PREFIX='source_occupation_tail_geometry_20261007'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_input(n,layers,family):
    weights=np.arange(1,layers+1,dtype=float)
    weights/=weights.sum()
    if family=='concentric_box_mixture':
        widths=np.exp(np.linspace(-.5*np.log(n),.5*np.log(n),layers))
        weights=1/np.arange(1,layers+1,dtype=float);weights/=weights.sum()
        components=[dict(kind='box',width=float(d)) for d in widths]
    else:
        A=16*n
        components=[dict(kind='grid',k=l+1,A=A,q=2*A*(l+1)+1,
             phase=.25*(l%2) if family=='capacity_weighted_grid' else 0.,
             alternating_positive_weights=family=='capacity_weighted_grid') for l in range(layers)]
    data=dict(n=n,layers=layers,family=family,components=components,weights=weights.tolist(),W=1.)
    data['input_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()
    return data


def sample_source(data,count,rng):
    n=data['n'];labels=rng.choice(data['layers'],size=count,p=data['weights'])
    y=np.empty((count,n))
    for l,c in enumerate(data['components']):
        at=np.flatnonzero(labels==l)
        if not len(at):continue
        if c['kind']=='box':y[at]=(rng.random((len(at),n))-.5)*c['width']
        else:
            if c['alternating_positive_weights']:
                p=1+2*(np.arange(c['q'])%2);p=p/p.sum()
                ix=rng.choice(c['q'],size=(len(at),n),p=p)
            else:ix=rng.integers(c['q'],size=(len(at),n))
            y[at]=(ix-c['A']*c['k']+c['phase'])/c['k']
    return y,labels


def responses(data,x,L):
    n=data['n'];logu=np.full((len(x),len(L)),-np.inf)
    arithmetic_boundary_count=0
    xx=x[:,None,:];R=L[None,:,None]
    for weight,c in zip(data['weights'],data['components']):
        if c['kind']=='box':
            d=c['width']
            left=xx-R/2;right=xx+R/2
            overlap=np.maximum(0,np.minimum(right,d/2)-np.maximum(left,-d/2))
            # Exact geometric plateau, avoiding cancellation R/R at fully interior cubes.
            fully_inside=(left>=-d/2)&(right<=d/2)
            axis_density=np.where(fully_inside,1/d,overlap/(d*R))
            arithmetic_boundary_count+=int(np.count_nonzero((np.abs(left+d/2)<1e-11)|(np.abs(right-d/2)<1e-11)))
        else:
            low=(xx-R/2)*c['k']-c['phase'];high=(xx+R/2)*c['k']-c['phase']
            arithmetic_boundary_count+=int(np.count_nonzero((np.abs(low-np.rint(low))<1e-11)|(np.abs(high-np.rint(high))<1e-11)))
            ilo=np.maximum(np.ceil(low)+c['A']*c['k'],0)
            ihi=np.minimum(np.floor(high)+c['A']*c['k'],c['q']-1)
            count=np.maximum(ihi-ilo+1,0)
            if c['alternating_positive_weights']:
                odd=np.where(count>0,np.floor((ihi+1)/2)-np.floor(ilo/2),0)
                mass=(count+2*odd)/(2*c['q']-1)
            else:mass=count/c['q']
            axis_density=mass/R
        with np.errstate(divide='ignore'):
            component=np.log(axis_density).sum(axis=-1)+np.log(weight)
        logu=np.logaddexp(logu,component)
    return logu,arithmetic_boundary_count


def receiver_samples(y,total,n,rng):
    C=1+n*np.log(2)
    core=rng.random(total)<1/C
    omega=rng.random((total,n))-.5
    radii=np.exp(rng.random(total)*np.log(2))
    face=rng.integers(n,size=total);sign=2*rng.integers(2,size=total)-1
    omega[np.arange(total),face]=sign*.5
    offsets=omega*radii[:,None]
    # Core is genuinely uniform cube1, its actual radius can be <1.
    offsets[core]=rng.random((int(core.sum()),n))-.5
    radial=2*np.max(np.abs(offsets),axis=-1)
    log_sup=np.where(core,0.,-n*np.log(radii))
    return y+offsets,radial,log_sup,core


def scalar_validation():
    import itertools
    data=dict(n=2,layers=1,family='validation',components=[dict(kind='grid',k=1,A=1,q=3,phase=0.,alternating_positive_weights=True)],weights=[1.])
    atoms=np.array(list(itertools.product([-1,0,1],repeat=2)))
    one=np.array([1,3,1])/5
    weights=np.array([one[i]*one[j] for i,j in itertools.product(range(3),repeat=2)])
    checks=0
    for x in [np.array([0.,0.]),np.array([.5,.5]),np.array([1.,0.]),np.array([.13,-.47])]:
        L=np.array([1.,1.5,2.])
        logu,_=responses(data,x[None,:],L)
        for j,R in enumerate(L):
            # Closed boundaries, same distances and full cumulative mass.
            distances=2*np.max(np.abs(atoms-x),axis=1)
            direct=weights[distances<=R].sum()/R**2
            assert np.isclose(np.exp(logu[0,j]),direct,atol=1e-15,rtol=2e-14)
            checks+=1
    box=dict(n=2,layers=1,family='validation',components=[dict(kind='box',width=2.)],weights=[1.])
    z,_=responses(box,np.array([[0.,0.]]),np.array([1.,1.5,2.]))
    assert np.all(z[0]==z[0,0]) and np.argmax(z[0])==0
    return dict(scalar_closed_atom_response_checks=checks,exact_box_plateau_tie_check=True)


def evaluate_case(plan,family,confidence_records):
    n=plan['n'];T=plan['inner_receivers_per_replicate'];H=plan['outer_sources']
    C=1+n*np.log(2);rng=np.random.default_rng(plan['seed']+{'capacity_weighted_grid':1000,'overlap_resolution_mixture':2000,'concentric_box_mixture':3000}[family])
    data=make_input(n,plan['layers'],family);y,labels=sample_source(data,H,rng)
    L=np.exp(np.linspace(0,np.log(2),plan['J']))
    alpha=np.array([.1,.5,2.]);logbase=-n*np.log(2) if family=='concentric_box_mixture' else -n*np.log(32*n)
    logtau=logbase+np.log(alpha)
    values=np.empty((2,H,3));lower=np.empty_like(values);upper=np.empty_like(values)
    all_logM=[];all_winner=[];all_events=[]
    threshold_near=np.zeros(3,dtype=int);response_tie_near=0;boundary_count=0;construction_boundary=0
    for rep in range(2):
        yy=np.repeat(y,T,axis=0)
        x,radial,logsup,core=receiver_samples(yy,len(yy),n,rng)
        scores=np.empty((len(x),3));lo=np.empty_like(scores);hi=np.empty_like(scores)
        logM=np.empty(len(x));winner=np.empty(len(x),dtype=int);event=np.empty((len(x),3),dtype=bool)
        for begin in range(0,len(x),512):
            end=min(len(x),begin+512)
            logu,bd=responses(data,x[begin:end],L);boundary_count+=bd
            wi=np.argmax(logu,axis=1);m=logu[np.arange(len(wi)),wi]
            logM[begin:end]=m;winner[begin:end]=wi
            near=logu>=m[:,None]-1e-11
            response_tie_near+=int(np.count_nonzero(near.sum(axis=1)>1))
            rr=radial[begin:end];ls=logsup[begin:end]
            construction_boundary+=int(np.count_nonzero(np.abs(rr[:,None]-L[None,:])<=1e-12))
            capture=rr[:,None]<=L[None,:]
            kernel_IS=C*np.exp(-n*np.log(L)[None,:]-ls[:,None])*capture
            assert np.max(kernel_IS)<=C*(1+1e-12)
            fixed=kernel_IS[np.arange(len(wi)),wi]
            low=np.min(np.where(near,kernel_IS,np.inf),axis=1)
            high=np.max(np.where(near,kernel_IS,-np.inf),axis=1)
            for a,lt in enumerate(logtau):
                diff=m-lt;ev=diff>0;amb=np.abs(diff)<=1e-11
                threshold_near[a]+=int(amb.sum());event[begin:end,a]=ev
                with np.errstate(over='ignore'):
                    g=np.exp(-np.maximum(diff,0))*ev
                scores[begin:end,a]=g*fixed
                lo[begin:end,a]=np.where(amb,0,g*low)
                hi[begin:end,a]=np.where(amb,high,g*high)
        values[rep]=scores.reshape(H,T,3).mean(axis=1)
        lower[rep]=lo.reshape(H,T,3).mean(axis=1);upper[rep]=hi.reshape(H,T,3).mean(axis=1)
        all_logM.append(logM);all_winner.append(winner);all_events.append(event)
    # Global finite-record concentration (arithmetic screen remains heuristic).
    source_band_count=confidence_records*max(256,H)*2
    e_inner=C*np.sqrt(np.log(2*source_band_count/.02)/(2*T))
    e_outer=C*np.sqrt(np.log(2*confidence_records*8/.03)/(2*H))
    e_square=C*C*np.sqrt(np.log(2*confidence_records*8/.03)/(2*H))
    records=[]
    for a,av in enumerate(alpha):
        sh=values[:,:,a].mean(axis=0)
        sl=np.maximum(0,lower[:,:,a].mean(axis=0)-e_inner)
        su=np.minimum(C,upper[:,:,a].mean(axis=0)+e_inner)
        I1=float(sh.mean());I2=float(np.mean(values[0,:,a]*values[1,:,a]))
        I1lo=max(0,float(sl.mean())-e_outer);I1hi=min(C,float(su.mean())+e_outer)
        I2lo=max(0,I2-e_square);I2hi=min(C*C,I2+e_square)
        tails=[]
        for multiplier in [.5,1.,2.]:
            K=multiplier*np.sqrt(n)
            plugin=float(np.mean(np.maximum(sh-K,0)))
            tail_lo=max(0,float(np.mean(np.maximum(sl-K,0)))-e_outer)
            tail_hi=min(C,float(np.mean(np.maximum(su-K,0)))+e_outer)
            tails.append(dict(K=K,plugin_biased_integral=plugin,plugin_ratio=plugin/I1 if I1 else None,
                finite_sample_bound=[tail_lo,tail_hi],ratio_bound=[tail_lo/I1hi if I1hi else 0,tail_hi/I1lo if I1lo else None]))
        records.append(dict(alpha=float(av),log_tau=float(logtau[a]),I1=I1,I2_cross_replicate=I2,
            I2_over_I1=I2/I1 if I1 else None,I2_over_I1_sqrtn=I2/I1/np.sqrt(n) if I1 else None,
            I1_finite_sample_bound=[I1lo,I1hi],I2_finite_sample_bound=[I2lo,I2hi],
            max_observed_S_estimate=float(sh.max()),nonzero_receiver_events=int(sum(z[:,a].sum() for z in all_events)),
            near_threshold_events=int(threshold_near[a]),tails=tails))
    profile=BASE/(PREFIX+'_'+plan['name']+'_'+family+'_profiles.npz')
    assert not profile.exists()
    np.savez_compressed(profile,source_y=y,source_component=labels,L_nodes=L,
        S_replicates=values,S_arithmetic_lower=lower,S_arithmetic_upper=upper,
        receiver_logM=np.stack(all_logM),receiver_winner=np.stack(all_winner),receiver_event=np.stack(all_events))
    return dict(plan=plan,input=data,records=records,profile_path=str(profile),profile_sha256=sha(profile),
        envelope_C=C,inner_radius=e_inner,outer_radius=e_outer,square_radius=e_square,
        near_max_response_ties=response_tie_near,near_atom_or_box_boundaries=boundary_count,
        near_constructed_source_capture_boundaries=construction_boundary,
        arithmetic_screen='near-response/threshold/boundary checks are heuristic, not interval arithmetic',
        confidence_scope='Ideal exact finite-source model; floating-point evaluations not certified; global union .02 inner+.03 outer')


def main():
    regpath=BASE/(PREFIX+'_registration.json');reg=json.loads(regpath.read_text())
    resultpath=BASE/(PREFIX+'_results.json');assert not resultpath.exists()
    start=time.perf_counter();checks=scalar_validation();cases=[]
    for plan in reg['rounds']:
        for family in reg['families']:
            case=evaluate_case(plan,family,33);cases.append(case)
            receipt=BASE/(PREFIX+'_progress.json')
            receipt.write_text(json.dumps(dict(status='running',completed_cases=len(cases),last_plan=plan['name'],last_family=family,profiles=[dict(path=c['profile_path'],sha256=c['profile_sha256']) for c in cases]),indent=2)+'\n')
            print(json.dumps(dict(case=plan['name'],family=family,nonempty_records=sum(r['nonzero_receiver_events']>0 for r in case['records']),ratios=[r['I2_over_I1'] for r in case['records']])),flush=True)
    for plan in reg['same_n_dependency_checks']:
        case=evaluate_case(plan,reg['same_n_checks_family'],33);cases.append(case)
        print(json.dumps(dict(case=plan['name'],family=reg['same_n_checks_family'],ratios=[r['I2_over_I1'] for r in case['records']])),flush=True)
    payload=dict(status='complete',registration_sha256=sha(regpath),script_sha256=sha(Path(__file__)),
        cases=cases,validation=checks,elapsed_seconds=time.perf_counter()-start,numpy_version=np.__version__,python_version=sys.version,
        no_asymptotic_claim=True,no_actual_geom_certificate=True,
        square='unbiased independent conditional replicate product before arithmetic error',
        confidence='33 frozen records; conservative global Hoeffding; bounds often vacuous; not arithmetic certification')
    resultpath.write_text(json.dumps(payload,indent=2)+'\n')
    (BASE/(PREFIX+'_progress.json')).write_text(json.dumps(dict(status='complete',result_path=str(resultpath),result_sha256=sha(resultpath),elapsed_seconds=payload['elapsed_seconds']),indent=2)+'\n')
    print(json.dumps(dict(status='complete',records=sum(len(c['records']) for c in cases),elapsed_seconds=payload['elapsed_seconds'])),flush=True)


if __name__=='__main__':main()
