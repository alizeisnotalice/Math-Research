#!/usr/bin/env python3
"""One registered midpoint matrix pilot; no strict interior quadrature claim."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
from itertools import product,combinations,permutations
import hashlib,importlib.util,json,sys,time
import numpy as np
import sympy as sp
from scipy.special import gamma
from flint import ctx,arb

HERE=Path(__file__).resolve().parent
MOMENT_SOURCE=HERE.parent/'upper_order002'/'top2_certificate.py'
OLD_RESULT=HERE.parent/'upper_order002'/'result.json'
spec=importlib.util.spec_from_file_location('locked_top2_moments',MOMENT_SOURCE)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def main():
    started=time.monotonic();ctx.prec=256;pid=os.getpid()
    hashes={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'moment_source_sha256':hashlib.sha256(MOMENT_SOURCE.read_bytes()).hexdigest(),
            'tail_certificate_sha256':hashlib.sha256(OLD_RESULT.read_bytes()).hexdigest()}
    state={'status':'RUNNING_NUMERICAL_MATRIX_PILOT','pid':pid,**hashes,'orbit_index':0,'next_batch':0,'partial_baseline':0.0,'partial_adaptive':0.0,'max_kernel_ratio':1.0,'completed_orbits':[],'started_unix':time.time()}
    if '--resume' in sys.argv:
        saved=json.loads((HERE/'checkpoint.json').read_text())
        assert all(saved[k]==v for k,v in hashes.items())
        state.update({k:saved[k] for k in ('orbit_index','next_batch','partial_baseline','partial_adaptive','max_kernel_ratio','completed_orbits')})
    def checkpoint(status):
        state.update(status=status,seconds_this_checkpoint=time.monotonic()-started)
        (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
        if time.monotonic()-started>120:
            state['status']='INCOMPLETE_SAVED_RESUMABLE'
            (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
            raise SystemExit(75)
    checkpoint('STARTED')
    print(json.dumps({'pid':pid,'status':'started',**hashes}),flush=True)
    directions=[(1,)+t for t in product((-1,1),repeat=4)];lookup={v:i for i,v in enumerate(directions)}
    Vsp=sp.Matrix(5,16,lambda i,j:directions[j][i]);Psp=sp.eye(16)-Vsp.T*Vsp/16
    assert Vsp*Psp==sp.zeros(5,16) and Psp*Psp==Psp and all(Psp[i,i]==sp.Rational(11,16) for i in range(16))
    actions=[]
    for p in permutations(range(4)):
        for sign in (-1,1):
            action=[]
            for v in directions:
                out=tuple(v[j] for j in p)+(sign*v[4],)
                if out[0]<0:out=tuple(-a for a in out)
                action.append(lookup[out])
            actions.append(action)
    pairs=set(combinations(range(16),2));unseen=pairs.copy();orbits=[]
    while unseen:
        pair=min(unseen);orbit={tuple(sorted((a[pair[0]],a[pair[1]]))) for a in actions}
        assert orbit<=pairs
        unseen-=orbit;orbits.append((pair,len(orbit)))
    assert len(orbits)==11 and sum(n for pair,n in orbits)==120
    V=np.array(directions,dtype=float).T;outer=np.einsum('ai,bi->iab',V,V)
    P=np.eye(16)-V.T@V/16;off=P*P;np.fill_diagonal(off,0)
    diagonal2=(11/16)**2
    endpoints=[arb(1)/4]
    for j in range(256):endpoints.append(endpoints[-1]*17/16)
    endpoint_stats=[old.stats(z) for z in endpoints]
    midpoints=[(l*r).sqrt() for l,r in zip(endpoints,endpoints[1:])]
    midpoint_stats=[old.stats(z) for z in midpoints]
    F=[v[0] for v in endpoint_stats]
    delta=[r-l for l,r in zip(F,F[1:])]
    minweight=[(r**15-l**15)/15 for l,r in zip(F,F[1:])]
    diag=[2*(r*(r**15-l**15)/15-(r**16-l**16)/16) for l,r in zip(F,F[1:])]
    weights=np.empty((256,256))
    for i in range(256):
        weights[i,i]=float(diag[i])
        for j in range(i):weights[i,j]=weights[j,i]=float(delta[i]*minweight[j])
    exact_mass=F[-1]**16-16*F[-1]*F[0]**15+15*F[0]**16
    assert abs(weights.sum()*120-float(exact_mass))<1e-12
    first,second=np.indices((256,256));ii=first.ravel();jj=second.ravel();ww=weights.ravel()
    mids=np.array([float(z) for z in midpoints]);ms=np.array([float(v[1]) for v in midpoint_stats]);ns=np.array([float(1/v[2]) for v in midpoint_stats])
    x=np.array([1.,1.,1.,1.,0.]);I=np.eye(5);batchsize=2048
    kappa=3*gamma(1.75)/(2**3.5*np.pi**2.5)
    for orbit_index in range(state['orbit_index'],len(orbits)):
        pair,count=orbits[orbit_index];ret1,ret2=pair
        random=np.array([j for j in range(16) if j not in pair])
        start_batch=state['next_batch'] if orbit_index==state['orbit_index'] else 0
        for beg in range(start_batch,len(ii),batchsize):
            sl=slice(beg,min(beg+batchsize,len(ii)));imin=np.minimum(ii[sl],jj[sl])
            z1=mids[ii[sl]];z2=mids[jj[sl]];t=mids[imin];m=ms[imin];n=ns[imin];a=1/n;d=m*n-1
            variance_upper=m*(t-m)
            assert np.all(d>0) and np.all(variance_upper>0)
            nvec=np.repeat(n[:,None],16,axis=1);nvec[:,ret1]=1/z1;nvec[:,ret2]=1/z2
            csum=nvec@off
            Hbar=diagonal2*(m*d)[:,None]+variance_upper[:,None]*csum
            correction=diagonal2*d[:,None]**2/Hbar
            correction[:,ret1]=0;correction[:,ret2]=0
            neff=nvec-correction;aeff=1/neff
            assert np.all(neff>0)
            assert np.all(aeff[:,random]>=a[:,None]-1e-12) and np.all(aeff[:,random]<=m[:,None]+1e-12)
            assert np.all(neff[:,ret1]==1/z1) and np.all(neff[:,ret2]==1/z2)
            Bm=16*m[:,None,None]*I+(z1-m)[:,None,None]*outer[ret1]+(z2-m)[:,None,None]*outer[ret2]
            Ba=16*a[:,None,None]*I+(z1-a)[:,None,None]*outer[ret1]+(z2-a)[:,None,None]*outer[ret2]
            Beff=np.einsum('nj,jab->nab',aeff,outer,optimize=True)
            sign,logdet=np.linalg.slogdet(Bm);assert np.all(sign==1)
            rhs=np.broadcast_to(x,(len(m),5))[:,:,None]
            qbase=np.linalg.solve(Ba,rhs)[:,:,0]@x
            qadaptive=np.linalg.solve(Beff,rhs)[:,:,0]@x
            assert np.all(qadaptive>0) and np.all(qadaptive<=qbase*(1+1e-8))
            Kbase=np.exp(-.5*logdet-1.75*np.log(qbase))
            Kadaptive=np.exp(-.5*logdet-1.75*np.log(qadaptive))
            state['partial_baseline']+=float(ww[sl]@Kbase)
            state['partial_adaptive']+=float(ww[sl]@Kadaptive)
            state['max_kernel_ratio']=max(state['max_kernel_ratio'],float(np.max(Kadaptive/Kbase)))
            state['next_batch']=sl.stop
            checkpoint('ACCUMULATING_ONE_ORBIT')
        state['completed_orbits'].append({'representative_indices':list(pair),'multiplicity':count,
          'baseline_Green_midpoint':kappa*count*state['partial_baseline'],'adaptive_Green_midpoint':kappa*count*state['partial_adaptive']})
        state.update(orbit_index=orbit_index+1,next_batch=0,partial_baseline=0.0,partial_adaptive=0.0)
        checkpoint('ORBIT_COMPLETED')
    baseline=sum(r['baseline_Green_midpoint'] for r in state['completed_orbits'])
    adaptive=sum(r['adaptive_Green_midpoint'] for r in state['completed_orbits'])
    prior=json.loads(OLD_RESULT.read_text())
    result={'status':'NUMERICAL_FIXED_TOP2_BOX_MATRIX_PILOT_NOT_STRICT_INTERIOR_CERTIFICATE','point':[1,1,1,1,0],'pid':pid,'partition':{'left':'1/4','ratio':'17/16','bins':256,'full_square_midpoints':65536,'actual_stabilizer_orbits':11,'unordered_pairs':120},
      'baseline_Green_midpoint':baseline,'adaptive_Green_midpoint':adaptive,'relative_gain':adaptive/baseline-1,'ratio_to1e_minus4':adaptive/1e-4,'ratio_to128e_minus6':adaptive/128e-6,'largest_pointwise_kernel_ratio':state['max_kernel_ratio'],
      'orbit_contributions':state['completed_orbits'],'inherited_strict_Green_low_omission_upper':prior['Green_low_omission_upper'],'inherited_strict_Green_high_omission_upper':prior['Green_high_omission_upper'],
      'checks':{'VP_zero_P_projection_exact':True,'actual48_element_stabilizer':True,'probability_mass_closure_numeric':True,'retained_inverse_coefficients_exactly1_over_z':True,'sampled_effective_weights_between_a_and_m':True,'sampled_Qadaptive_at_most_Qbase':True},
      'scope':'Analytic adaptive LS improvement and full-support tail control are proved; interior midpoint sum has no certified quadrature error and is not a rigorous bound on the full formula oractualG.',
      'seconds':time.monotonic()-started,'python':sys.executable,**hashes}
    (HERE/'pilot_result.json').write_text(json.dumps(result,indent=2)+'\n')
    state.update(status='COMPLETE_EXIT0',exit_code=0,seconds_this_checkpoint=time.monotonic()-started)
    (HERE/'checkpoint.json').write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='orbit_contributions'},indent=2),flush=True)

if __name__=='__main__':main()
