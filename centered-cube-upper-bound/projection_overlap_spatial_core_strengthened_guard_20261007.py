"""Original-symbol periodic spatial diagnostics. No continuum/interval certification."""
from pathlib import Path
from fractions import Fraction
import numpy as np
import json, hashlib, time, math

HERE=Path(__file__).resolve().parent
PREFIX='projection_overlap_spatial_core_strengthened'
DATE='20261007'
REG=HERE/f'{PREFIX}_registration_{DATE}.json'
OUT=HERE/f'{PREFIX}_results_{DATE}.json'
assert REG.exists() and not OUT.exists(), 'Never overwrite registered data'
registration=json.loads(REG.read_text())
started=time.time()
records=[]; sources=[]; files=[]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def frac(v):return float(Fraction(str(v)))
def conv(v, symbol):return np.fft.ifft2(np.fft.fft2(v,axes=(-2,-1))*symbol,axes=(-2,-1)).real
def save(name, **arrays):
    p=HERE/f'{PREFIX}_{name}_{DATE}.npz'
    assert not p.exists()
    np.savez_compressed(p,**arrays)
    files.append({'file':p.name,'sha256':sha(p)})
    return p.name
def bump(t,width):return np.maximum(1-(2*t/width)**2,0)**3
def wrapped(v):return (v+np.pi)%(2*np.pi)-np.pi
def source(n,N,spec):
    phase=2*np.pi*np.arange(N)/N
    U=np.zeros((N,N));Deff=0.
    for center,cw in zip(spec['centers_pi_units'],spec['center_weights']):
        cw=frac(cw)
        dx=wrapped(phase-np.pi*center[0]);dy=wrapped(phase-np.pi*center[1])
        for width,lw in zip(spec['widths'],spec['layer_weights']):
            weight=cw*frac(lw);Deff+=weight*width
            U+=weight*bump(dx,width)[:,None]*bump(dy,width)[None,:]
    freq=np.fft.fftfreq(N,d=1/N)
    v=freq**2;g=np.ones(N);positive=v>0;g[positive]=np.log1p(v[positive])/v[positive]
    S=n/2*((1-g)[:,None]+(1-g)[None,:])
    C0=101/100
    amp=35/(8*C0*n*Deff);u=amp*U;Su=conv(u,S)
    omega=U>0
    f=np.where(omega,1+Su,0)
    mu=np.where(omega,1.,-Su)
    filename=save(f'source_n{n}_N{N}_{spec["name"]}',U=U,u=u,f=f,mu=mu,Omega=omega,
                  phase=phase,g1_coordinate_symbol=g,S_symbol=S)
    screens={'min_nu':float(f.min()),'min_mu':float(mu.min()),'max_mu':float(mu.max()),
             'source_mu_mass_error':float(abs(f.mean()-mu.mean())),
             'obstacle_identity_error':float(np.max(np.abs(f-mu-Su))),
             'Gi_coordinate_kernel_min':float(np.fft.ifft(g).real.min())}
    source_record={'n':n,'N':N,'name':spec['name'],'A':amp,'D_eff':Deff,
                   'W_phase_average':float(f.mean()),'Omega_volume_fraction':float(omega.mean()),
                   'nu_peak':float(f.max()),'file':filename,'screens':screens,
                   'formula_input_frozen':True}
    sources.append(source_record)
    return {'f':f,'mu':mu,'u':u,'Omega':omega,'g':g,'S':S,'spec':spec,'file':filename,
            'source_record':source_record,'FFT_f':np.fft.fft2(f)}

def symbols(n,g,J):
    r=(np.arange(J+1)/J)**2
    logs=np.log1p(r[:,None]*(g[None,:]-1))
    P1=np.exp((n/4)*logs)
    V1=np.exp((n/4)*(logs[1:]-logs[:-1]))
    kernelP=np.fft.ifft(P1[1:],axis=1).real
    kernelV=np.fft.ifft(V1,axis=1).real
    screen={'P_coordinate_kernel_min':float(kernelP.min()),
            'V_coordinate_kernel_min':float(kernelV.min()),
            'P_coordinate_mass_error':float(np.max(np.abs(kernelP.sum(axis=1)-1))),
            'V_coordinate_mass_error':float(np.max(np.abs(kernelV.sum(axis=1)-1)))}
    return r,P1,V1,screen

def evaluate(n,N,inp,J,alpha,extra=False):
    f,mu,u=inp['f'],inp['mu'],inp['u'];omega=inp['Omega'];W=float(f.mean())
    r,P1,V1,kernel_screen=symbols(n,inp['g'],J)
    maximum=f.copy();second=np.full_like(f,-np.inf);labels=np.zeros((N,N),np.uint8)
    for j in range(1,J+1):
        p=P1[j];k=np.outer(p*p,p*p)
        response=np.fft.ifft2(inp['FFT_f']*k).real
        new=response>maximum
        second=np.where(new,maximum,np.maximum(second,response))
        maximum=np.where(new,response,maximum);labels=np.where(new,j,labels).astype(np.uint8)
    region=(~omega)&(maximum>alpha)
    theta=np.zeros_like(f);theta[region]=alpha/maximum[region]
    weak_value=float(alpha*region.mean())
    near_threshold=(~omega)&(np.abs(maximum-alpha)<=1e-10*(1+np.abs(maximum)))
    near_tie=region&((maximum-second)<=1e-10*(1+np.abs(maximum)))
    weights=np.stack([theta,region.astype(float)])
    # All recursions below use two selector versions but the same original nodes.
    c=np.zeros((2,N,N));d=np.zeros_like(c);h=np.zeros_like(c);pair=np.zeros_like(c)
    direct_def=np.zeros(2);direct_pair=np.zeros(2);min_pi=math.inf;max_pi=-math.inf
    for j in range(J,0,-1):
        pmat=np.outer(P1[j],P1[j])
        q=conv(weights*(labels==j)[None,:,:],pmat)
        min_pi=min(min_pi,float(q.min()));max_pi=max(max_pi,float(q.max()))
        if j==J:
            c=q.copy();h=q.copy()
        else:
            vmat=np.outer(V1[j],V1[j]) # Increment from j to j+1.
            stack=conv(np.concatenate([c,d,h,pair],axis=0),vmat)
            vc,vd,vh,ve=stack[:2],stack[2:4],stack[4:6],stack[6:8]
            Pf=np.fft.ifft2(inp['FFT_f']*pmat).real
            direct_def+=np.mean(Pf[None,:,:]*q*vc,axis=(1,2))
            direct_pair+=np.mean(Pf[None,:,:]*q*vh,axis=(1,2))
            c=q+(1-q)*vc
            d=vd+q*vc
            h=q+vh
            pair=ve+q*vh
    pmat=np.outer(P1[1],P1[1])
    projected=conv(np.concatenate([c,d,h,pair],axis=0),pmat)
    Cstop,Cdef,Ctot,Cpair=projected[:2],projected[2:4],projected[4:6],projected[6:8]
    survival=1-Cstop
    signed_SCdef=conv(Cdef,inp['S'])
    metrics=[]
    for v,kind in enumerate(['weak_normalized','uncapped_comparison']):
        T=float(np.mean(f*Ctot[v]));stop=float(np.mean(f*Cstop[v]));D=float(np.mean(f*Cdef[v]));P=float(np.mean(f*Cpair[v]))
        target=weak_value if v==0 else float(np.mean(region*maximum))
        cap_mass=float(np.mean(mu*Cdef[v]));potential=float(np.mean(u*signed_SCdef[v]))
        metrics.append({'kind':kind,'selected_traffic':target,'source_projected_traffic':T,
                        'T_over_W':target/W,'Tstop_over_W':stop/W,'D_over_W':D/W,'pair_over_W':P/W,
                        'mu_Cdef_over_W':cap_mass/W,'signed_u_S_Cdef_over_W':potential/W,
                        'mu_Cdef_cap_area_bound_over_W':float(region.mean())/W,
                        'mu_Cdef_cap_exact_column_bound_over_W':float(Cdef[v].mean())/W,
                        'source_projection_residual':abs(T-target),'stop_plus_D_residual':abs(stop+D-target),
                        'direct_future_union_D_residual':abs(D-direct_def[v]),
                        'direct_future_pair_residual':abs(P-direct_pair[v]),
                        'obstacle_D_identity_residual':abs(D-cap_mass-potential),
                        'min_D_column':float(Cdef[v].min()),'min_pair_minus_D_column':float((Cpair[v]-Cdef[v]).min()),
                        'min_stop_column':float(Cstop[v].min()),'max_stop_column':float(Cstop[v].max()),
                        'min_signed_S_Cdef':float(signed_SCdef[v].min()),
                        'Ctot_stop_D_max_residual':float(np.max(np.abs(Ctot[v]-Cstop[v]-Cdef[v])))})
    fn=save(f'n{n}_N{N}_{inp["spec"]["name"]}_J{J}_alpha{alpha:g}'+('_spatial' if extra else ''),
            r=r,P_coordinate_symbols=P1,V_coordinate_symbols=V1,M_grid=maximum,
            second_response=second,winner_labels=labels,region=region,theta_weight=theta,
            Cstop=Cstop,Cdef=Cdef,Ctotal=Ctot,Cpair=Cpair,survival=survival,signed_S_Cdef=signed_SCdef)
    record={'n':n,'N':N,'source':inp['spec']['name'],'J':J,'alpha':alpha,
            'source_file':inp['file'],'profiles_file':fn,'extra_spatial_stability':extra,
            'W_phase_average':W,'receiver_region_fraction':float(region.mean()),
            'receiver_weak_value':weak_value,'near_threshold_receivers':int(near_threshold.sum()),
            'near_top_two_receivers_on_region':int(near_tie.sum()),
            'minimum_pi':min_pi,'maximum_pi':max_pi,'kernel_screens':kernel_screen,'metrics':metrics,
            'r_grid_is_continuous_winner':False,'numeric_zero_is_analytic_zero':False}
    records.append(record)
    return record,maximum,labels,region

stability=[]
for round_spec in registration['rounds']:
    n,N=round_spec['n'],round_spec['phase_grid']
    for spec in registration['sources']:
        inp=source(n,N,spec)
        previous={}
        for J in registration['parameter_panels']:
            for alpha in registration['threshold_alpha']:
                record,mx,labels,region=evaluate(n,N,inp,J,alpha)
                if alpha in previous:
                    old_record,old_mx,old_labels,old_region=previous[alpha]
                    old_r=(old_labels/old_record['J'])**2;new_r=(labels/J)**2
                    union=old_region|region
                    record['nested_r_diagnostics']={
                        'coarse_J':old_record['J'],'max_response_increase':float(np.max(mx-old_mx)),
                        'max_response_decrease_roundoff':float(np.min(mx-old_mx)),
                        'receiver_region_disagreement_fraction':float(np.mean(old_region!=region)),
                        'winner_r_disagreement_on_union_fraction':float(np.mean((old_r!=new_r)[union])) if union.any() else 0.,
                        'weak_D_over_W_difference':record['metrics'][0]['D_over_W']-old_record['metrics'][0]['D_over_W'],
                        'weak_pair_over_W_difference':record['metrics'][0]['pair_over_W']-old_record['metrics'][0]['pair_over_W']}
                previous[alpha]=(record,mx,labels,region)
        if n==32 and spec['name']=='nested_shared_resolution':
            stability.append(previous[1])
        print(json.dumps({'n':n,'N':N,'source':spec['name'],'elapsed_seconds':time.time()-started}),flush=True)

for N in [128,512]:
    spec=registration['sources'][1]
    inp=source(32,N,spec)
    stability.append(evaluate(32,N,inp,32,1,extra=True))
stability.sort(key=lambda t:t[0]['N'])
spatial=[]
for old,new in zip(stability,stability[1:]):
    orow,om,ol,oa=old;nrow,nm,nl,na=new
    step=nrow['N']//orow['N'];sample_m=nm[::step,::step];sample_l=nl[::step,::step];sample_a=na[::step,::step]
    union=oa|sample_a
    spatial.append({'coarse_N':orow['N'],'fine_N':nrow['N'],
                    'max_response_difference_on_common_nodes':float(np.max(np.abs(sample_m-om))),
                    'receiver_region_disagreement_common_nodes':float(np.mean(oa!=sample_a)),
                    'winner_label_disagreement_common_union':float(np.mean((ol!=sample_l)[union])) if union.any() else 0.,
                    'W_difference':nrow['W_phase_average']-orow['W_phase_average'],
                    'weak_D_over_W_difference':nrow['metrics'][0]['D_over_W']-orow['metrics'][0]['D_over_W'],
                    'weak_pair_over_W_difference':nrow['metrics'][0]['pair_over_W']-orow['metrics'][0]['pair_over_W']})

payload={'status':'COMPLETED_NUMERICAL_DIAGNOSTIC_NOT_INTERVAL_CERTIFIED',
         'runtime_seconds':time.time()-started,'registration_sha256':sha(REG),'script_sha256':sha(Path(__file__)),
         'sources':sources,'records':records,'spatial_stability':spatial,'saved_arrays':files,
         'scope':registration['limitations'],'confidence':'No MC confidence claim; grids/FFT float errors are not certified.'}
OUT.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps({'status':payload['status'],'records':len(records),'runtime_seconds':payload['runtime_seconds'],
                  'max_weak_D_over_W':max(r['metrics'][0]['D_over_W'] for r in records),
                  'max_weak_pair_over_W':max(r['metrics'][0]['pair_over_W'] for r in records),
                  'max_obstacle_identity_residual':max(m['obstacle_D_identity_residual'] for r in records for m in r['metrics'])}))
