"""Finite tensor lattice: exact event candidates, floating arithmetic, fixed MC protocol."""
import argparse
import itertools
import json
import math
from pathlib import Path
import time
import numpy as np

OUT = Path(__file__).resolve().parent
CS = [-0.5, 0.0, 0.5, 1.0]


def tensor_maximal(x, z, L):
    """Arrays (...,n); returns log(u L^n), smallest maximizing R and screens.

    Closed cubes with side R in [1,2]. Each coordinate uses z_i-2,...,z_i+2.
    x is restricted to z+[-1,1]^n, so these cover every relevant integer.
    Equal arrival values are merged before evaluating the product count.
    """
    n = x.shape[-1]
    k = z[..., :, None] + np.arange(-2, 3)
    d = np.where((k >= 0) & (k < L), 2 * np.abs(x[..., :, None] - k), np.inf)
    d.sort(axis=-1)
    c1 = np.sum(d <= 1., axis=-1)
    c2 = np.sum(d <= 2., axis=-1)
    logbase = np.log(np.maximum(c1, 1)).sum(axis=-1)
    zeros = np.sum(c1 == 0, axis=-1)
    previous = np.arange(5)
    inc = np.zeros(5)
    inc[1:] = np.log(previous[1:] + 1) - np.log(previous[1:])
    pending = (d > 1.) & (d <= 2.)
    eventr = np.where(pending, d, np.inf).reshape(x.shape[:-1] + (5*n,))
    loginc = np.where(pending, inc, 0.).reshape(eventr.shape)
    first = (pending & (previous == 0)).reshape(eventr.shape)
    order = np.argsort(eventr, axis=-1, kind='stable')
    eventr = np.take_along_axis(eventr, order, axis=-1)
    loginc = np.take_along_axis(loginc, order, axis=-1)
    first = np.take_along_axis(first, order, axis=-1)
    remain = zeros[..., None] - np.cumsum(first, axis=-1)
    logprod = logbase[..., None] + np.cumsum(loginc, axis=-1)
    last_tie = np.ones(eventr.shape, dtype=bool)
    last_tie[..., :-1] = eventr[..., :-1] != eventr[..., 1:]
    with np.errstate(invalid='ignore'):
        val = logprod - n*np.log(eventr)
    val = np.where((remain == 0) & last_tie & np.isfinite(eventr), val, -np.inf)
    val1 = np.where(zeros == 0, logbase, -np.inf)
    val2 = np.where(np.all(c2 > 0, axis=-1), np.log(np.maximum(c2, 1)).sum(axis=-1)-n*math.log(2), -np.inf)
    vals = np.concatenate((val1[..., None], val, val2[..., None]), axis=-1)
    radii = np.concatenate((np.ones(val1.shape+(1,)), eventr, np.full(val1.shape+(1,), 2.)), axis=-1)
    arg = np.argmax(vals, axis=-1)
    best = np.take_along_axis(vals, arg[..., None], axis=-1)[..., 0]
    winner = np.take_along_axis(radii, arg[..., None], axis=-1)[..., 0]
    # Detection only, not interval certification or an outward correction.
    finite = np.isfinite(d)
    screen = {'arrivals_near_window_endpoint': int(np.sum(finite & ((np.abs(d-1)<1e-12)|(np.abs(d-2)<1e-12))))}
    return best, winner, screen


def explicit_maximal(x, L):
    atoms = np.array(list(itertools.product(range(L), repeat=len(x))), dtype=float)
    radii = 2*np.max(np.abs(atoms-x), axis=-1)
    candidates = sorted(set([1., 2.] + radii[(radii > 1)&(radii <= 2)].tolist()))
    vals = [math.log(int(np.sum(radii <= R)))-len(x)*math.log(R) if np.any(radii<=R) else -math.inf for R in candidates]
    arg = int(np.argmax(vals))
    return vals[arg], candidates[arg]


def validate():
    rng = np.random.default_rng(837103)
    records=[]
    for n in [1,2,3,4]:
        for L in [2,3,4]:
            for _ in range(30):
                z=rng.integers(L,size=n)
                x=z+rng.uniform(-1,1,size=n)
                v,R,_=tensor_maximal(x,z,L)
                e,er=explicit_maximal(x,L)
                assert abs(v-e)<2e-12 and abs(R-er)<2e-12, (x,v,e,R,er)
                records.append((n,L))
    special=[]
    for n in [1,2,3,4]:
        for L in [2,3,4]:
            for zp in [0,L-1]:
                for offset in [0.,.5,1.,-.5,-1.]:
                    z=np.full(n,zp);x=z+offset
                    v,R,s= tensor_maximal(x,z,L)
                    e,er=explicit_maximal(x,L)
                    assert abs(v-e)<2e-12 and abs(R-er)<2e-12, (x,v,e,R,er)
                    special.append({'n':n,'L':L,'z':zp,'offset':offset,'log_u_normalized':float(v),'R':float(R),'screen':s})
    result={'status':'PASS_FLOAT_EVENT_VS_EXPLICIT','random_cases':len(records),'ties_boundary_cases':len(special),'special_cases':special,'tolerance':2e-12,'claim':'Implementation cross-check, not interval proof.'}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2))
    return result


def draw(n,L,samples,seed):
    rng=np.random.default_rng(seed)
    z=rng.integers(0,L,size=(samples,n))
    omega=rng.uniform(-.5,.5,size=(samples,n))
    face=rng.integers(0,n,size=samples)
    sign=rng.choice([-1.,1.],size=samples)
    omega[np.arange(samples),face]=.5*sign
    return z,omega


def profiles(n,L,samples,nodes,seed):
    z,omega=draw(n,L,samples,seed)
    t=np.linspace(0,n*math.log(2),nodes)
    r=np.exp(t/n)
    all_y=np.empty((len(CS),samples,nodes))
    screens={'arrivals_near_window_endpoint':0,'winner_near_r':0,'hardband_near_threshold':0,'boundary_sources':int(np.sum(np.any((z==0)|(z==L-1),axis=-1)))}
    for start in range(0,samples,16):
        zz=z[start:start+16,None,:]
        xx=zz+r[None,:,None]*omega[start:start+16,None,:]
        values,winner,scr=tensor_maximal(xx,zz,L)
        screens['arrivals_near_window_endpoint']+=scr['arrivals_near_window_endpoint']
        screens['winner_near_r']+=int(np.sum(np.abs(winner-r)<1e-12))
        kernel=np.exp(n*(np.log(r)[None,:]-np.log(winner)))
        kernel=np.where(r[None,:]<=winner,kernel,0.)
        for j,c in enumerate(CS):
            shifted=values-c*math.sqrt(n)
            band=(shifted>math.log(2)) & (shifted<=math.log(4))
            screens['hardband_near_threshold']+=int(np.sum((np.abs(shifted-math.log(2))<1e-10)|(np.abs(shifted-math.log(4))<1e-10)))
            all_y[j,start:start+16]=band*kernel
    return t,all_y,screens


def run_round(tag,samples,nodes,replicates=8):
    path=OUT/f'{tag}.json'
    assert not path.exists(), 'Round outputs are immutable; select a new tag.'
    rows=[]; tic=time.time()
    for n in [8,32,128]:
        energy=np.empty((replicates,len(CS)))
        coarse=np.empty_like(energy)
        simple=np.empty_like(energy)
        summaries=[]
        for rep in range(replicates):
            # Independent A/B samples; within-side common paths across t retain correlations.
            seed=817209+int(tag[-1])*1000000+n*1000+rep*10
            t,ya,sa=profiles(n,257,samples,nodes,seed)
            _,yb,sb=profiles(n,257,samples,nodes,seed+1)
            pa=ya.mean(axis=1);pb=yb.mean(axis=1)
            energy[rep]=np.trapezoid(pa*pb,t,axis=-1)
            coarse[rep]=np.trapezoid((pa*pb)[:,::2],t[::2],axis=-1)
            # Biased same-sample square is reported solely as a diagnostic.
            simple[rep]=np.trapezoid(pa*pa,t,axis=-1)
            summaries.append({'replicate':rep,'seeds':[seed,seed+1],'energy':energy[rep].tolist(),'coarse_energy':coarse[rep].tolist(),'same_sample_square_diagnostic':simple[rep].tolist(),'screens_a':sa,'screens_b':sb})
            np.savez_compressed(OUT/f'{tag}_n{n}_rep{rep}.npz',t=t,phi_a=pa,phi_b=pb)
        for j,c in enumerate(CS):
            mean=float(energy[:,j].mean());sd=float(energy[:,j].std(ddof=1));se=sd/math.sqrt(replicates)
            # t(7)=2.3646 for the fixed eight replicate protocol, approximate MC CI only.
            ci=[max(0.,mean-2.364624251*se),mean+2.364624251*se]
            # Truly bounded concentration for grid statistic Z in [0,T], generally vacuous.
            T=n*math.log(2);hoe=T*math.sqrt(math.log(40)/(2*replicates))
            rows.append({'n':n,'L':257,'log10_total_atoms':n*math.log10(257),'c':c,'log_lambda':-n*math.log(257)+c*math.sqrt(n),'energy_mean':mean,'energy_replicate_sd':sd,'approx95_MC_t_interval':ci,'hoeffding95_grid_interval':[max(0.,mean-hoe),min(T,mean+hoe)],'nested_grid_mean_abs_difference':float(np.mean(np.abs(energy[:,j]-coarse[:,j]))),'nested_grid_max_abs_difference':float(np.max(np.abs(energy[:,j]-coarse[:,j]))),'same_sample_square_diagnostic_mean':float(simple[:,j].mean()),'replicate_data':summaries})
        print(json.dumps({'round':tag,'n':n,'samples_per_side':samples,'nodes':nodes,'means':energy.mean(axis=0).tolist(),'elapsed_seconds':time.time()-tic}),flush=True)
        # Intermediate snapshots recover completed dimensions if a later one fails.
        (OUT/f'{tag}_partial.json').write_text(json.dumps(rows,indent=2))
    result={'round':tag,'protocol':{'samples_per_side':samples,'nodes':nodes,'coarse_nodes':(nodes+1)//2,'replicates':replicates,'n':[8,32,128],'L':257,'c':CS,'window_side_length':[1,2],'winner_tie':'smallest R, exact floating equal scores','cone':'uniform 2n faces, uniform free coordinates','jacobian':'dx = r^n dt d(normalized cone); no 2^-n','integrand':'Y=1_{2lambda<u<=4lambda}(r/R)^n 1_{r<=R}; energy=int(EY)^2 dt','status':'Hardband g1 relaxation only; no original softFIRST geom qualification.'},'elapsed_seconds':time.time()-tic,'rows':rows}
    path.write_text(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--round',type=int);p.add_argument('--validate',action='store_true');a=p.parse_args()
    if a.validate: print(json.dumps({'validation':validate()['status']}))
    if a.round:
        protocol={1:(128,129),2:(512,257),3:(1024,513)}
        run_round('round'+str(a.round),*protocol[a.round])
