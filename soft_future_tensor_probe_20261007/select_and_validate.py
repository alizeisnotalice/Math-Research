from pathlib import Path
import importlib.util,json,math,hashlib
import numpy as np
D=Path(__file__).resolve().parent;OLD=D.parent/'nonconcentrated_tensor_probe_20261007'
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
p=module('hard_engine',OLD/'probe.py');soft=module('soft_engine',D/'model.py')
points=[]
for n,M,K in [(128,32,513),(512,16,1025)]:
 for rep in [0,1]:
  original=json.loads((OLD/f'highdim_n{n}.json').read_text())
  seed=original['replicates'][rep]['draw_A']['seed'];Z,om,label=p.draws(n,64,M,seed)
  data=np.load(OLD/f'highdim_n{n}_pair{rep}.npz');times=data['t'];found=None
  for ki in np.where(data['profile_A'][3,0]>0)[0]:
   tt=float(times[ki]);r=math.exp(tt/n);dx=r*om/2
   sc,rr,lm,lr,ix,valid=p.candidates(dx,Z,64);rows=np.arange(M)
   u=sc[rows,ix]+n*math.log(64);R=rr[rows,ix];delay=n*(np.log(R)-math.log(r));tol=(n+1)*1e-10
   good=(u>math.log(2)+tol)&(u<=math.log(4)-tol)&(lr[rows,ix]<math.log(p.ETA0/math.sqrt(n))-tol)&(R>=r)&(delay<=p.VSTAR-tol)&(n*np.log(R)>=p.VSTAR)
   if good.any():
    si=int(np.where(good)[0][0]);zz=Z[si];dd=dx[si];arrivals=[]
    for q in p.Q:
     dist=2*np.abs(dd[:,None]+(zz[:,None]-(4//q)*np.arange(64*q)[None,:])/4)
     arrivals.extend(dist[(dist>1)&(dist<=2)].tolist())
    breaks=np.unique([1.,2.]+arrivals)
    found=dict(id=f'n{n}_pair{rep}_A_firsteligible',n=n,seed=seed,sample_count=M,old_nodes=K,old_pair=rep,sample_index=si,node_index=int(ki),t=tt,r=r,Z=zz.tolist(),omega=om[si].tolist(),dx=dd.tolist(),source_component_label=int(label[si]),hard_R=float(R[si]),hard_beta=float(n*math.log(R[si])),hard_delay=float(delay[si]),normalized_log_u=float(u[si]),microbox_log_share_upper=float(lr[si,ix[si]]),eta=p.ETA0/math.sqrt(n),L_breakpoints=breaks.tolist())
    break
  assert found is not None,(n,rep)
  points.append(found)
  print('SELECTED',found['id'],found['normalized_log_u'],found['hard_beta'],len(found['L_breakpoints']),flush=True)
(D/'frozen_points.json').write_text(json.dumps(dict(rule='For each n128/n512 pair0/pair1 saved side-A stream, first saved time with lower-cert Gamma_vstar>0 yielding first source row with hardband, microbox-cert, d<=vstar, beta>=vstar. No soft response inspected before selection.',points=points),indent=2))
# Independent scalar positive Gauss-Legendre integration, explicit coordinate grid.
x,w=np.polynomial.legendre.leggauss(512);v=(x+1)/2;weight=w/2;maxdiff=0.;records=[]
for n in [2,3]:
 Z=np.arange(n)*17;dx=np.linspace(-.7,.7,n);Ls=[1.,1.13,1.57,2.]
 h,b,e,c=soft.coordinates(Z,dx,Ls,128)
 for j,q in enumerate(p.Q):
  for li,L in enumerate(Ls):
   for i in range(n):
    pos=np.abs(dx[i]+(Z[i]-(4//q)*np.arange(64*q))/4)/L
    a=np.abs(pos-.5);bb=pos+.5
    inside=pos<=.5
    ff=np.where(inside[:,None],v[None,:]*(2-np.exp(-a[:,None]/v)-np.exp(-bb[:,None]/v)),v[None,:]*(np.exp(-a[:,None]/v)-np.exp(-bb[:,None]/v)))
    ref=float((ff@weight).sum()/(q*L));diff=abs(ref-b[j,li,i]);maxdiff=max(maxdiff,diff)
    assert diff<=e[j,li,i]+1e-11,(j,li,i,diff,e[j,li,i])
    records.append(dict(n=n,q=int(q),L=L,coord=i,geometric_trapezoid_value=float(b[j,li,i]),independent_direct_GL512=ref,difference=diff,analytic_trapezoid_error_plus_unverified_floating_guard=float(e[j,li,i])))
# Hard faces: forced tagged R=r retained for every source component to which z belongs.
faces=[]
for n in [128,512]:
 for r in [1.,math.sqrt(2),2.]:
  Z=np.full(n,252);dx=np.zeros(n);dx[0]=r/2
  h,b,e,c=soft.coordinates(Z,dx,[r],128)
  faces.append(dict(n=n,r=r,tagged_source_arrival=2*abs(dx[0]),hard_coordinate_count=[float(h[j,0,0]*p.Q[j]*r) for j in range(3)]))
(D/'validation.json').write_text(json.dumps(dict(status='PASS_FLOAT_EXACT_GRID_AGGREGATION_VS_DIRECT_GAUSS',cases=len(records),max_abs_difference=maxdiff,records=records,local_source_face_checks=faces,scope='Independent floating check. Analytic trapezoid error is rigorous in exact arithmetic; NumPy exp/log/sum rounding not outward certified.'),indent=2))
(D/'preregistration.json').write_text(json.dumps(dict(status='FROZEN_BEFORE_SOFT_PRESSURE',point_file_sha256=hashlib.sha256((D/'frozen_points.json').read_bytes()).hexdigest(),frozen_source='Original full three tensor-component mu, W1, qres=[1,2,4],p=[1/4,1/2,1/4],support[0,64)^n',lambda_q='Original loglambda=-nlog64, common q=lambda fixed, no optimization',rounds=[dict(round=r,spectral_trapezoid_N=N,uniform_L_nodes=LN,sigma_k_nodes=SN,broad_sigma_nodes=BN) for r,N,LN,SN,BN in [(1,128,17,65,17),(2,256,33,129,33),(3,512,65,257,65)]],L_policy='EVERY coordinate arrival of all full components plus1/2 endpoints and nested logL grid; all are closed cube evaluations. Left/right traces of hard jumps saved separately.',sigma_policy='k=n sigma uniform[0,32] plus uniform[0,1] sigma and dyadic k powers beyond32; common full future grid, not monotone assumption. Candidate roots are finite-L last-visible root brackets, not actual continuous FIRST.',future_average='At candidate sigma and its finite-L winner compute exact positive Bernstein integration for integer alpha=ceil(sqrt n): K/q=alpha integral v^(alpha-1) F_{sigma+(1-sigma)v,L}/q dv. Polynomial identity, floating arithmetic plus propagated scalar integral error; no new theory fee claimed.',certification='Analytic scalar trapezoid error plus documented unverified floating guard. Finite node violations with margin above analytic error are robust numerical exclusions, not formal outward-rounded certificates. Finite-grid success cannot prove full future cap. L sandwich oracle may numerically enclose cells; unresolved cells remain unknown. No FIRST/GOOD/CPGP/fullactualgeom assertion.',model_sha256=hashlib.sha256((D/'model.py').read_bytes()).hexdigest()),indent=2))
