#!/usr/bin/env python3
"""Two fixed positive inner rules; one outer top2 midpoint partition."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from itertools import product,combinations,permutations
from collections import defaultdict
import hashlib,importlib.util,json,sys,time
import numpy as np
from scipy.special import gamma
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
HELPER=HERE.parent/'upper_order002'/'top2_certificate.py'
OLD_RESULT=HERE.parent/'upper_order002'/'result.json'
spec=importlib.util.spec_from_file_location('locked_moment_helper',HELPER)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def main():
    started=time.monotonic();ctx.prec=256;pid=os.getpid()
    hashes={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'moment_source_sha256':hashlib.sha256(HELPER.read_bytes()).hexdigest(),'tail_certificate_sha256':hashlib.sha256(OLD_RESULT.read_bytes()).hexdigest()}
    state={'status':'STARTED','pid':pid,**hashes,'inner_rows':0,'orbit_index':0,'next_batch':0,'partial_baseline':0.,'partial_values':[[0.,0.],[0.,0.]],'completed_orbits':[],
      'min_schur':1.,'max_Kopt_over_Kfixed':0.,'min_aeff_over_a':1e99,'max_aeff_over_m':0.,'started_unix':time.time()}
    grid=np.full((256,65),np.nan);cdf=np.full((256,65),np.nan)
    if '--resume' in sys.argv:
        saved=json.loads((HERE/'checkpoint.json').read_text());assert all(saved[k]==v for k,v in hashes.items())
        for k in state:
            if k not in ('pid','status','started_unix',*hashes.keys()):state[k]=saved[k]
        cache=np.load(HERE/'inner_checkpoint.npz');grid=cache['grid'];cdf=cache['cdf']
    def checkpoint(status,save_inner=False):
        state.update(status=status,seconds_this_checkpoint=time.monotonic()-started)
        if save_inner:np.savez_compressed(HERE/'inner_checkpoint.npz',grid=grid,cdf=cdf)
        (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
        if time.monotonic()-started>120:
            state['status']='INCOMPLETE_SAVED_RESUMABLE'
            np.savez_compressed(HERE/'inner_checkpoint.npz',grid=grid,cdf=cdf)
            (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
            raise SystemExit(75)
    checkpoint('STARTED')
    print(json.dumps({'pid':pid,'status':'started',**hashes}),flush=True)
    directions=[(1,)+t for t in product((-1,1),repeat=4)];lookup={v:i for i,v in enumerate(directions)}
    V=np.array(directions,dtype=float).T;outer=np.einsum('ai,bi->iab',V,V)
    P=np.eye(16)-V.T@V/16;off=P*P;np.fill_diagonal(off,0);p2=(11/16)**2
    actions=[]
    for p in permutations(range(4)):
        for sign in (-1,1):
            action=[]
            for v in directions:
                out=tuple(v[j] for j in p)+(sign*v[4],)
                if out[0]<0:out=tuple(-a for a in out)
                action.append(lookup[out])
            actions.append(action)
    unseen=set(combinations(range(16),2));orbits=[]
    while unseen:
        pair=min(unseen);orbit={tuple(sorted((a[pair[0]],a[pair[1]]))) for a in actions};unseen-=orbit;orbits.append((pair,len(orbit)))
    assert len(orbits)==11 and sum(v for p,v in orbits)==120
    endpoints=[arb(1)/4]
    for _ in range(256):endpoints.append(endpoints[-1]*17/16)
    estat=[old.stats(t) for t in endpoints];mids_arb=[(l*r).sqrt() for l,r in zip(endpoints,endpoints[1:])]
    mstat=[old.stats(t) for t in mids_arb];mids=np.array([float(t) for t in mids_arb]);ms=np.array([float(v[1]) for v in mstat]);ns=np.array([float(1/v[2]) for v in mstat])
    F=[v[0] for v in estat];delta=[r-l for l,r in zip(F,F[1:])];mw=[(r**15-l**15)/15 for l,r in zip(F,F[1:])]
    weights=np.empty((256,256))
    for i in range(256):
        l,r=F[i:i+2];weights[i,i]=float(2*(r*(r**15-l**15)/15-(r**16-l**16)/16))
        for j in range(i):weights[i,j]=weights[j,i]=float(delta[i]*mw[j])
    for row in range(state['inner_rows'],256):
        ratio=mids_arb[row]*4
        for j in range(65):
            s=arb(1)/4 if j==0 else (mids_arb[row] if j==64 else ratio**(arb(j)/64)/4)
            grid[row,j]=float(s);cdf[row,j]=float(old.stats(s)[0])
        state['inner_rows']=row+1
        if row%16==0:checkpoint('BUILDING_INNER_TABLE',True)
    checkpoint('INNER_TABLE_COMPLETE',True)
    low_probability=float(estat[0][0]);low_mean=float(estat[0][1]);rules={}
    for level in (32,64):
        stride=64//level;nodes=grid[:,::stride];Fs=cdf[:,::stride]
        probs=np.column_stack((np.full(256,low_probability),np.diff(Fs,axis=1)))
        reps=np.column_stack((np.full(256,low_mean),np.sqrt(nodes[:,:-1]*nodes[:,1:])))
        assert np.all(probs>0) and np.max(np.abs(probs.sum(axis=1)-np.array([float(v[0]) for v in mstat])))<1e-12
        rules[level]=(reps,probs)
    first,second=np.indices((256,256));ii=first.ravel();jj=second.ravel();ww=weights.ravel();I=np.eye(5);x=np.array([1.,1.,1.,1.,0.]);batchsize=2048
    kappa=3*gamma(1.75)/(2**3.5*np.pi**2.5)
    for orbit_index in range(state['orbit_index'],len(orbits)):
        pair,count=orbits[orbit_index];r1,r2=pair;groups=defaultdict(list);ret_inner=int(V[:,r1]@V[:,r2])
        for l in range(16):
            if l in pair:continue
            b=int(V[:,r1]@V[:,l]);d=int(V[:,r2]@V[:,l]);groups[(b*b,d*d,ret_inner*b*d)].append(l)
        group_indices=list(groups.values());representatives=[g[0] for g in group_indices];group_outer=np.array([outer[g].sum(axis=0) for g in group_indices]);Vg=V[:,representatives];ng=len(representatives)
        for beg in range(state['next_batch'],len(ii),batchsize):
            sl=slice(beg,min(beg+batchsize,len(ii)));imin=np.minimum(ii[sl],jj[sl]);z1=mids[ii[sl]];z2=mids[jj[sl]];m=ms[imin];n=ns[imin];a=1/n
            Bm=16*m[:,None,None]*I+(z1-m)[:,None,None]*outer[r1]+(z2-m)[:,None,None]*outer[r2]
            Ba=16*a[:,None,None]*I+(z1-a)[:,None,None]*outer[r1]+(z2-a)[:,None,None]*outer[r2]
            invBa=np.linalg.inv(Ba);h=np.einsum('ag,nab,bg->ng',Vg,invBa,Vg,optimize=True);schur=1-a[:,None]*h
            assert np.min(schur)>11/16-1e-8
            nvec=np.repeat(n[:,None],16,axis=1);nvec[:,r1]=1/z1;nvec[:,r2]=1/z2
            Kfixed=(nvec@off)[:,representatives]/p2;Kopt=h/schur
            assert np.all(Kopt>0) and np.all(Kopt<=Kfixed*(1+1e-8))
            state['min_schur']=min(state['min_schur'],float(schur.min()));state['max_Kopt_over_Kfixed']=max(state['max_Kopt_over_Kfixed'],float(np.max(Kopt/Kfixed)))
            K=np.concatenate((Kfixed,Kopt),axis=1)
            sign,logdet=np.linalg.slogdet(Bm);assert np.all(sign==1)
            qbase=np.einsum('a,nab,b->n',x,invBa,x);base=np.exp(-.5*logdet-1.75*np.log(qbase))
            state['partial_baseline']+=float(ww[sl]@base)
            for lev_index,level in enumerate((32,64)):
                reps,probs=rules[level];nodes=reps[imin];ps=probs[imin];A0=np.zeros_like(K);A1=np.zeros_like(K)
                for node in range(level+1):
                    s=nodes[:,node,None];weighted=ps[:,node,None]/(1+K*s);A0+=weighted;A1+=s*weighted
                eff=A1/A0
                assert np.all(eff>0)
                state['min_aeff_over_a']=min(state['min_aeff_over_a'],float(np.min(eff/a[:,None])))
                state['max_aeff_over_m']=max(state['max_aeff_over_m'],float(np.max(eff/m[:,None])))
                for variant in range(2):
                    ae=eff[:,variant*ng:(variant+1)*ng]
                    Beff=np.einsum('ng,gab->nab',ae,group_outer,optimize=True)+z1[:,None,None]*outer[r1]+z2[:,None,None]*outer[r2]
                    rhs=np.broadcast_to(x,(len(m),5))[:,:,None];sol=np.linalg.solve(Beff,rhs)[:,:,0];qeff=sol@x
                    assert np.all(qeff>0)
                    kernel=np.exp(-.5*logdet-1.75*np.log(qeff))
                    state['partial_values'][lev_index][variant]+=float(ww[sl]@kernel)
            state['next_batch']=sl.stop;checkpoint('ACCUMULATING_ORBIT')
        state['completed_orbits'].append({'representative_indices':list(pair),'multiplicity':count,'random_column_groups':ng,'baseline_midpoint':kappa*count*state['partial_baseline'],
          'Green_values_by_inner_level_and_column':[[kappa*count*v for v in row] for row in state['partial_values']]})
        state.update(orbit_index=orbit_index+1,next_batch=0,partial_baseline=0.,partial_values=[[0.,0.],[0.,0.]])
        checkpoint('ORBIT_COMPLETED')
    baseline=sum(v['baseline_midpoint'] for v in state['completed_orbits']);values=np.array([v['Green_values_by_inner_level_and_column'] for v in state['completed_orbits']]).sum(axis=0)
    prior=json.loads(OLD_RESULT.read_text())
    result={'status':'NUMERICAL_RESOLVENT_TWO_INNER_LEVEL_DIAGNOSTIC_ONLY','point':[1,1,1,1,0],'pid':pid,'baseline_Green_midpoint':baseline,
      'inner_levels':[32,64],'column_variants':['normalized_fixed_projection','optimal_null_columns'],'Green_values':values.tolist(),'level64_gain_over_baseline':(values[1]/baseline-1).tolist(),'level64_optimal_over_fixed':float(values[1,1]/values[1,0]-1),'level64_optimal_ratio_to1e_minus4':float(values[1,1]/1e-4),
      'relative_level_changes':(values[1]/values[0]-1).tolist(),'guard_diagnostics':{k:state[k] for k in ('min_schur','max_Kopt_over_Kfixed','min_aeff_over_a','max_aeff_over_m')},
      'low_cell':{'delta':'1/4','probability':low_probability,'conditional_mean_atom':low_mean,'scope':'Approximate positive mean atom; analytic convexity error described inproof.md, not propagated to an integral enclosure.'},
      'inherited_Green_low_omission_upper':prior['Green_low_omission_upper'],'inherited_Green_high_omission_upper':prior['Green_high_omission_upper'],
      'orbit_results':state['completed_orbits'],'scope':'Outer midpoint plus positive inner quadrature, not a strict full-integral certificate. Level differences are not error bounds or confidence intervals.',
      'seconds':time.monotonic()-started,'python':sys.executable,**hashes}
    (HERE/'pilot_result.json').write_text(json.dumps(result,indent=2)+'\n')
    state.update(status='COMPLETE_EXIT0',exit_code=0,seconds_this_checkpoint=time.monotonic()-started)
    (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='orbit_results'},indent=2),flush=True)

if __name__=='__main__':main()
