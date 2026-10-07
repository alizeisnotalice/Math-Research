"""Original finite-Fourier correlated sources; nonlinear integrals are numerical screens."""
from pathlib import Path
from fractions import Fraction as Q
from math import comb, pi
from hashlib import sha256
import json,time
import numpy as np

H=Path(__file__).resolve().parent
REG=H/'ordered_source_weighted_dual_registration_20261007.json'
OUT=H/'ordered_source_weighted_dual_results_20261007.json'
assert REG.exists() and not OUT.exists(), 'Preserve every registered result'
started=time.time()
lams=np.asarray([0.125,0.5,1.,2.,8.])
grids=(64,128,256); tol=1e-7
checks=0; frac_checks=0; records=[]; saved=[]
max_grid_diff=0.; max_tel=0.; min_integrated=0.; min_pointwise=0.

def g(c,z):
    z=np.asarray(z,dtype=float)
    b=np.zeros_like(z); good=z>0
    b[good]=z[good]/np.log1p(z[good])-1
    return 1/(1+c*b)
def phi(u,lam):
    u=np.maximum(u,0)
    x=u/lam
    ans=u-lam*np.log1p(x)
    small=x<1e-4
    z=x[small]
    ans[small]=lam*(z*z/2-z**3/3+z**4/4-z**5/5+z**6/6)
    return ans
def values(coeff,modes,N):
    spec=np.zeros((N,N),dtype=np.complex128)
    spec[modes[:,0]%N,modes[:,1]%N]=coeff*N*N
    return np.fft.ifft2(spec).real

for n,d in ((8,3),(32,5),(128,7)):
    vw=[(1,0),(0,1)]
    cyc=[(1,1),(1,-1),(2,1),(-1,2),(0,1),(1,0)]
    vw+= [cyc[(j-2)%len(cyc)] for j in range(2,n)]
    vw=np.asarray(vw)
    modes=np.asarray([(k,l) for k in range(-d,d+1) for l in range(-d,d+1)])
    k,l=modes[:,0],modes[:,1]
    base=(1-np.abs(k)/(d+1))*(1-np.abs(l)/(d+1))
    one=(1-np.abs(k)/(d+1))*(l==0)
    source=[base.astype(complex),
            0.4*base+0.4*base*np.exp(1j*(k*pi/3+l*pi/2))+0.2*one]
    paths=[np.arange(n),np.arange(n)[::-1],
           np.random.default_rng(202610071200+n).permutation(n),
           np.random.default_rng(202610071201+n).permutation(n)]
    # Exact random-permutation edge weights, not an exponential enumeration.
    for j in range(n):
        beta=Q(1,comb(n,j)*(n-j))
        assert beta==Q(1,n*comb(n-1,j))
        frac_checks+=1
    for c in (1.,0.75,0.5):
        mult=np.stack([g(c,(k*a+l*b)**2) for a,b in vw])
        full_multiplier=mult.prod(axis=0)
        for si,coef0 in enumerate(source):
            all_entropies=[]
            for pathid,path in enumerate(paths):
                coeff=coef0.copy()
                path_entropy=np.empty((n+1,len(grids),len(lams)))
                local_point=[]
                for step in range(n+1):
                    previous=coeff.copy()
                    if step: coeff=coeff*mult[path[step-1]]
                    for gi,N in enumerate(grids):
                        u=values(coeff,modes,N)
                        assert u.min()>-tol
                        assert abs(u.mean()-1)<tol
                        checks+=2
                        for li,lam in enumerate(lams):
                            path_entropy[step,gi,li]=phi(u,lam).mean()
                        if step and step in (1,n//2,n) and N in (128,256):
                            old=values(previous,modes,N)
                            ix=np.fft.fftfreq(N,d=1/N)
                            a,b=vw[path[step-1]]
                            action=g(c,(ix[:,None]*a+ix[None,:]*b)**2)
                            for li,lam in enumerate(lams):
                                oldphi=phi(old,lam)
                                gp=np.fft.ifft2(np.fft.fft2(oldphi)*action).real
                                defect=gp-phi(u,lam)
                                minimum=float(defect.min())
                                min_pointwise=min(min_pointwise,minimum)
                                assert minimum>=-tol
                                checks+=1
                                local_point.append({'step':step,'grid':N,'lambda':float(lam),
                                                    'min_pointwise_defect':minimum,
                                                    'defect_mean':float(defect.mean())})
                assert np.max(np.abs(coeff-coef0*full_multiplier))<1e-12
                checks+=1
                edges=path_entropy[:-1]-path_entropy[1:]
                minedge=float(edges.min())
                min_integrated=min(min_integrated,minedge)
                assert minedge>=-tol
                assert path_entropy.min()>=-tol and path_entropy.max()<=1+tol
                checks+=2
                discrepancy=edges.sum(axis=0)-(path_entropy[0]-path_entropy[-1])
                tel=float(np.abs(discrepancy).max())
                max_tel=max(max_tel,tel)
                assert tel<1e-12
                checks+=1
                dif=float(np.max(np.abs(path_entropy[:,1:,:]-path_entropy[:,:-1,:])))
                max_grid_diff=max(max_grid_diff,dif)
                npz=H/f'ordered_source_weighted_dual_n{n}_c{c:g}_s{si}_p{pathid}_20261007.npz'
                assert not npz.exists()
                np.savez_compressed(npz,entropy=path_entropy,edge_defect=edges,
                                    path=path,phase_vectors=vw,lambda_values=lams,
                                    grids=np.asarray(grids),source_coeff=coef0,modes=modes)
                saved.append({'file':npz.name,'sha256':sha256(npz.read_bytes()).hexdigest()})
                records.append({'n':n,'degree':d,'c':c,'source':si,'path':pathid,
                                'min_integrated_edge_defect':minedge,
                                'max_nested_grid_difference':dif,
                                'max_telescoping_discrepancy':tel,
                                'finest_initial_entropy':path_entropy[0,-1,:].tolist(),
                                'finest_final_entropy':path_entropy[-1,-1,:].tolist(),
                                'pointwise_selected_edges':local_point})
                all_entropies.append(path_entropy)
            avg=np.mean(all_entropies,axis=0)
            avgfee=np.mean([v[0]-v[-1] for v in all_entropies],axis=0)
            assert np.max(np.abs(avgfee-(avg[0]-avg[-1])))<1e-12
            assert avgfee.min()>=-tol and avgfee.max()<=1+tol
            checks+=2
    print(json.dumps({'n':n,'completed':True,'seconds':time.time()-started}),flush=True)
res={'status':'PASS_EXACT_AND_NUMERICAL_SCREEN','fraction_checks':frac_checks,
     'numerical_checks':checks,'elapsed_seconds':time.time()-started,
     'max_nested_grid_difference':max_grid_diff,'max_telescoping_discrepancy':max_tel,
     'min_integrated_edge_defect':min_integrated,'min_pointwise_defect':min_pointwise,
     'registration_sha256':sha256(REG.read_bytes()).hexdigest(),
     'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
     'records':records,'saved_data':saved,
     'scope':'Original G_c correlated two-phase torus sources; no interval or graph-sampling/weak endpoint certification.'}
OUT.write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({key:res[key] for key in ('status','fraction_checks','numerical_checks',
             'elapsed_seconds','max_nested_grid_difference','max_telescoping_discrepancy',
             'min_integrated_edge_defect','min_pointwise_defect')}))
