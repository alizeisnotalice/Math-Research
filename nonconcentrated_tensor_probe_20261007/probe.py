#!/usr/bin/env python3
"""B41-inspired finite overlapping tensor mixture, floating pressure only."""
from pathlib import Path
import argparse, hashlib, itertools, json, math, time
import numpy as np
OUT=Path(__file__).resolve().parent
Q=np.array([1,2,4]); P=np.array([.25,.5,.25]); ETA0=49/65536
DIMS=[16,32,64]; REPS=4; VSTAR=math.log(65536/49)


def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2))


def draws(n,L,M,seed):
    rng=np.random.default_rng(seed)
    # ONE shared resolution for the whole n-dimensional source atom.
    labels=rng.choice(3,M,p=P)
    Z=np.empty((M,n),dtype=np.int64)
    for j,q in enumerate(Q):
        rows=np.where(labels==j)[0]
        Z[rows]=rng.integers(0,L*q,size=(len(rows),n))*(4//q)
    omega=rng.uniform(-1,1,(M,n));face=rng.integers(0,2*n,M)
    omega[np.arange(M),face//2]=np.where(face%2,1.,-1.)
    return Z,omega,labels


def candidates(dx,Z,L):
    """All closed-cube coordinate arrivals, R in [1,2], local coordinates.

    Mixture component log masses are updated separately; all identical floating
    arrivals across coordinates/components are completed before scoring.
    Mathematical arbitrary score ties remain numerically screened.
    """
    M,n=dx.shape
    event=[];incs=[];first=[];which=[];logbase=[];zero=[];c2all=[]
    for j,q in enumerate(Q):
        step=4//q; residue=Z%step; base=Z//step
        offsets=np.arange(-2*q,2*q+2)
        indices=base[:,:,None]+offsets
        local=(step*offsets[None,None,:]-residue[:,:,None])/4
        dist=np.where((indices>=0)&(indices<L*q),2*np.abs(dx[:,:,None]-local),np.inf)
        dist.sort(axis=-1)
        c1=np.sum(dist<=1,axis=-1);c2=np.sum(dist<=2,axis=-1)
        logbase.append(np.log(np.maximum(c1,1)).sum(axis=1));zero.append((c1==0).sum(axis=1));c2all.append(c2)
        k=np.arange(len(offsets));inc=np.zeros(len(offsets));inc[1:]=np.log(k[1:]+1)-np.log(k[1:])
        pending=(dist>1)&(dist<=2)
        event.append(np.where(pending,dist,np.inf).reshape(M,-1))
        incs.append(np.where(pending,inc,0).reshape(M,-1))
        first.append((pending&(k==0)).reshape(M,-1))
        which.append(np.full(n*len(offsets),j))
    er=np.concatenate(event,axis=1);li=np.concatenate(incs,axis=1);fi=np.concatenate(first,axis=1);wh=np.concatenate(which)
    order=np.argsort(er,axis=1,kind='stable')
    er=np.take_along_axis(er,order,axis=1);li=np.take_along_axis(li,order,axis=1);fi=np.take_along_axis(fi,order,axis=1);wh=wh[order]
    logs=[];alive=[]
    for j in range(3):
        zz=zero[j][:,None]-np.cumsum(fi&(wh==j),axis=1)
        prod=logbase[j][:,None]+np.cumsum(np.where(wh==j,li,0),axis=1)
        # c2 log product independently recomputed; guards accumulated drift.
        l2=np.log(np.maximum(c2all[j],1)).sum(axis=1)
        lp=np.concatenate((logbase[j][:,None],prod,l2[:,None]),axis=1)
        ok=np.concatenate(((zero[j]==0)[:,None],zz==0,np.all(c2all[j]>0,axis=1)[:,None]),axis=1)
        logs.append(np.where(ok,math.log(P[j])-n*math.log(L*Q[j])+lp,-np.inf));alive.append(ok)
    lm=np.logaddexp.reduce(np.stack(logs),axis=0)
    radii=np.concatenate((np.ones((M,1)),er,np.full((M,1),2.)),axis=1)
    last=np.ones(er.shape,dtype=bool);last[:,:-1]=er[:,:-1]!=er[:,1:]
    valid=np.concatenate((np.ones((M,1),dtype=bool),last&np.isfinite(er),np.ones((M,1),dtype=bool)),axis=1)
    scores=np.where(valid,lm-n*np.log(radii),-np.inf)
    ix=np.argmax(scores,axis=1);rows=np.arange(M)
    # In production h=2/n is smaller than every lattice spacing. For general
    # validation use global safe floor(h*q)+1 per coordinate and zero if empty.
    h=2/n
    caplogs=np.stack([np.where(alive[j],math.log(P[j])+n*(math.log(math.floor(h*Q[j])+1)-math.log(L*Q[j])),-np.inf) for j in range(3)])
    lc=np.logaddexp.reduce(caplogs,axis=0)
    with np.errstate(invalid="ignore"):
        logratio=lc-lm
    return scores,radii,lm,logratio,ix,valid


def explicit(n,L,x):
    # Independent construction: enumerate complete grids, merge identical atoms.
    masses={}
    for q,p in zip(Q,P):
        for k in itertools.product(range(L*q),repeat=n):
            z=tuple(v*(4//q) for v in k)
            masses[z]=masses.get(z,0.)+float(p)/(L*q)**n
    atoms=np.array(list(masses))/4;w=np.array(list(masses.values()))
    d=2*np.max(np.abs(atoms-x),axis=1)
    rr=sorted(set([1.,2.]+d[(d>1)&(d<=2)].tolist()))
    mm=np.array([w[d<=r].sum() for r in rr]);vals=np.where(mm>0,np.log(np.maximum(mm,1e-300))-n*np.log(rr),-np.inf)
    ix=int(np.argmax(vals))
    return vals[ix],rr[ix],mm[ix],float(w.sum()),len(masses)


def validate():
    rng=np.random.default_rng(202610076100);maxerr=0.;cases=0;winnererr=0.
    for n in [1,2,3,4]:
        L=2
        for k in range(12):
            Z,omega,labels=draws(n,L,1,202610070000+n*100+k)
            r=float(rng.uniform(1,2));dx=r*omega/2
            sc,rr,lm,lr,ix,v=candidates(dx,Z,L)
            ev,er,em,total,count=explicit(n,L,Z[0]/4+dx[0])
            err=abs(sc[0,ix[0]]-ev);maxerr=max(maxerr,err);winnererr=max(winnererr,abs(rr[0,ix[0]]-er))
            assert err<5e-12 and abs(rr[0,ix[0]]-er)<5e-12
            assert abs(total-1)<5e-12;cases+=1
    # Formula mixture priors; empirical sampling jointly checks resolution and
    # all-coordinate divisibility, with binomial deviation recorded, not asserted exact.
    Z,om,label=draws(3,2,100000,202610076111)
    freq=np.bincount(label,minlength=3)/len(label)
    for j,q in enumerate(Q):assert np.all(Z[label==j]%(4//q)==0)
    # Deterministic overlap mass from all component labels at the origin.
    origin_mass=sum(float(p)/(2*q)**3 for q,p in zip(Q,P))
    save('validation.json',dict(status='PASS_FLOAT_SMALL_EXPLICIT',cases=cases,max_log_score_error=maxerr,max_winner_error=winnererr,
        exact_priors=P.tolist(),empirical_priors=freq.tolist(),empirical_prior_errors=(freq-P).tolist(),sample_count=100000,
        whole_atom_shared_resolution=True,origin_merged_mass=origin_mass,scope='Independent enumerated full source/winner cross-check; floating arithmetic, not interval proof.'))
    print('VALIDATION PASS',maxerr,flush=True)


def profile(n,L,M,K,seed):
    Z,omega,labels=draws(n,L,M,seed)
    t=np.linspace(0,n*math.log(2),K);out=np.zeros((4,3,K));fixed=np.zeros((2,K));stats=[]
    ll=-n*math.log(L);eta=ETA0/math.sqrt(n)
    for ki,tt in enumerate(t):
        r=math.exp(tt/n);dx=r*omega/2
        sums=np.zeros((4,3));fixedsum=np.zeros(2);st=dict(certified=0,unknown=0,cert_fixed=0,winner_capture=0,near_score_tie=0,near_band=0,near_delay=0,min_log_ratio=math.inf,max_log_ratio=-math.inf,mean_mass=0.,mean_R=0.,winner_equals_r=0)
        for start in range(0,M,16):
            sc,rr,lm,lr,ix,valid=candidates(dx[start:start+16],Z[start:start+16],L);rows=np.arange(len(ix))
            logu=sc[rows,ix];R=rr[rows,ix];delay=n*(np.log(rr)-math.log(r));tol=(n+1)*1e-10
            near=valid&(sc>=logu[:,None]-tol)
            # Count-based microbox bound proves branch implication analytically;
            # tol merely screens floating log comparisons, not an interval certificate.
            cert=lr<math.log(eta)-tol
            band=(logu>ll+math.log(2))&(logu<=ll+math.log(4))
            blo=(logu>ll+math.log(2)+tol)&(logu<=ll+math.log(4)-tol)
            bhi=(logu>ll+math.log(2)-tol)&(logu<=ll+math.log(4)+tol)
            cap=rr>=r;short=cap&(delay<=1)
            caplo=rr>r;caplo|=rr==r
            shortlo=caplo&(delay<=1-tol)
            caphi=rr>=r*(1-tol/n);shorthi=caphi&(delay<=1+tol)
            value=np.exp(-np.maximum(delay,0.))
            masks=[short,cap,short,cap&(delay<=VSTAR)];mlo=[shortlo,caplo,shortlo,caplo&(delay<=VSTAR-tol)];mhi=[shorthi,caphi,shorthi,caphi&(delay<=VSTAR+tol)]
            w=[np.ones_like(value),value,value,np.ones_like(value)]
            rawcert=cert[rows,ix]
            for qi in range(4):
                raw=w[qi][rows,ix]*masks[qi][rows,ix]*band*rawcert
                lo=np.min(np.where(near,w[qi]*mlo[qi]*cert*blo[:,None],np.inf),axis=1)
                # Unknown residual is retained in this upper envelope.
                hi=np.max(np.where(near,w[qi]*mhi[qi]*bhi[:,None],-np.inf),axis=1)
                sums[qi]+=np.array([lo.sum(),raw.sum(),hi.sum()])
            fcert=lr[rows,ix]<math.log(ETA0)-tol
            fixedsum+=np.array([np.sum(short[rows,ix]*band*fcert),np.sum(value[rows,ix]*cap[rows,ix]*band*fcert)])
            st['certified']+=int(rawcert.sum());st['unknown']+=int((~rawcert).sum());st['cert_fixed']+=int(fcert.sum())
            st['winner_capture']+=int((R>=r).sum());st['near_score_tie']+=int((near.sum(axis=1)>1).sum())
            st['near_band']+=int(((abs(logu-ll-math.log(2))<tol)|(abs(logu-ll-math.log(4))<tol)).sum())
            d=n*(np.log(R)-math.log(r));st['near_delay']+=int(((abs(d)<tol)|(abs(d-1)<tol)).sum())
            st['winner_equals_r']+=int((R==r).sum());st['min_log_ratio']=min(st['min_log_ratio'],float(lr[rows,ix].min()));st['max_log_ratio']=max(st['max_log_ratio'],float(lr[rows,ix].max()))
            st['mean_mass']+=float(np.exp(lm[rows,ix]).sum())/M;st['mean_R']+=float(R.sum())/M
        out[:,:,ki]=sums/M;fixed[:,ki]=fixedsum/M;stats.append(st)
    return t,out,fixed,stats,dict(seed=seed,component_label_counts=np.bincount(labels,minlength=3).tolist(),source_sha256=hashlib.sha256(Z.astype('<i8').tobytes()+omega.astype('<f8').tobytes()).hexdigest(),boundary_sources=int(np.any((Z==0)|(Z>=4*L-4),axis=1).sum()))


def run_round(j):
    M,K={1:(32,65),2:(64,129),3:(128,257)}[j];path=OUT/f'round{j}.json'
    assert not path.exists(),'Frozen output refuses overwrite'
    results=[];tic=time.time()
    for n in DIMS:
        records=[];energies=[]
        for rep in range(REPS):
            seed=202610076000+j*100000+n*100+rep*2
            t,a,fa,sa,da=profile(n,64,M,K,seed)
            _,b,fb,sb,db=profile(n,64,M,K,seed+1)
            cross=np.trapezoid(a*b,t,axis=-1);coarse=np.trapezoid((a*b)[:,:,::2],t[::2],axis=-1)
            rec=dict(rep=rep,draw_A=da,draw_B=db,energy=cross.tolist(),coarse_energy=coarse.tolist(),screen_A=sa,screen_B=sb,
                fixed_eta0_diagnostic_energy=np.trapezoid(fa*fb,t,axis=-1).tolist())
            records.append(rec);energies.append(cross)
            np.savez_compressed(OUT/f'round{j}_n{n}_pair{rep}.npz',t=t,profile_A=a,profile_B=b,fixed_eta0_A=fa,fixed_eta0_B=fb)
            print(json.dumps(dict(round=j,n=n,pair=rep,energy_cert=cross[:,0].tolist(),energy_upper_unknown=cross[:,2].tolist(),seconds=round(time.time()-tic,2))),flush=True)
        ee=np.array(energies);T=n*math.log(2);ho=T*math.sqrt(math.log(40)/(2*REPS));mean=ee.mean(axis=0);sd=ee.std(axis=0,ddof=1)
        results.append(dict(n=n,L=64,q=Q.tolist(),p=P.tolist(),log_lambda=-n*math.log(64),eta=ETA0/math.sqrt(n),eta0=ETA0,
            labeled_atom_log10_counts=[n*math.log10(64*q) for q in Q],distinct_atom_log10_count=n*math.log10(256),
            quantity_order=['Gamma_v1','Theta_hard','Theta_hard_v1','Gamma_vstar'],vstar=VSTAR,profile_axis=['screened_lower_cert_res','raw_cert_res','floating_upper_including_unknown'],
            energy_mean=mean.tolist(),energy_sd=sd.tolist(),approx95_t3_interval=np.stack([np.maximum(0,mean-3.182446305*sd/2),mean+3.182446305*sd/2],axis=-1).tolist(),
            hoeffding95_grid_intervals=np.stack([np.maximum(0,mean-ho),np.minimum(T,mean+ho)],axis=-1).tolist(),replicates=records))
        save(f'round{j}_partial.json',results)
    save(f'round{j}.json',dict(round=j,samples=M,nodes=K,coarse_nodes=(K+1)//2,pairs=REPS,elapsed_seconds=time.time()-tic,rows=results))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');p.add_argument('--round',type=int);args=p.parse_args()
    if args.validate:validate()
    if args.round:run_round(args.round)
