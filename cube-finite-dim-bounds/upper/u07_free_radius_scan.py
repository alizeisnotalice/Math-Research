"""Free-R H46 discovery; complete algebraic boundary-face candidate search.
Floating polynomial roots and sphere integrals are not certified bounds.
"""
import sys
sys.dont_write_bytecode=True
import hashlib,json,math
from pathlib import Path
from functools import lru_cache
import numpy as np
from numpy.polynomial import Polynomial as P
from numpy.polynomial.legendre import leggauss
from scipy.special import roots_jacobi
from scipy.optimize import minimize,minimize_scalar
from u06_h46_scan import sphere_extrema,hs,mass as allroot_mass

class Found(Exception):
    def __init__(self,row):self.row=row

def stationary_roots(poly,n,k,m):
    S=P([0.,1.]);p=(n-2)/2
    equation=poly**2*S**3-p*p if n==3 else poly*S**2-p
    roots=[]
    for z in equation.roots():
        if abs(z.imag)<1e-7 and k<z.real<k+m:
            s=float(z.real)
            if poly(s)>0 and abs(poly(s)-p*s**(-n/2))<1e-6:roots.append(s)
    return roots

@lru_cache(maxsize=4096)
def cube_candidates(n,B4,B6):
    S=P([0.,1.]);a=3/(n+2);b=15/(n+8);g=30/((n+8)*(n+4));p=(n-2)/2;candidates=[]
    def add(k,free,kind):
        y=np.r_[np.ones(k),np.asarray(free),np.zeros(n-k-len(free))]
        if np.any(y<0) or np.any(y>1):return
        h4,h6=hs(y,n);value=y.sum()**((2-n)/2)-B4*h4-B6*h6
        candidates.append((float(value),y.tolist(),kind))
    def grad(u,p2):return 2*a*B4*S+b*B6*p2-3*g*B6*S*S+(-2*B4+2*b*B6*S)*u-3*B6*u*u
    for k in range(1,n+1):
        add(k,[],'vertex')
        for m in range(1,n-k+1):
            u=(S-k)/m;p2=k+m*u*u
            for ss in stationary_roots(grad(u,p2),n,k,m):
                uu=float(u(ss))
                if 0<uu<1:add(k,[uu]*m,'one_free_value')
            if not B6 or m<2:continue
            T=10*S/(n+8)-2*B4/(3*B6)
            for pp in range(1,m):
                qq=m-pp
                if pp!=qq:
                    u=(S-k-qq*T)/(pp-qq);v=T-u;p2=k+pp*u*u+qq*v*v
                    for ss in stationary_roots(grad(u,p2),n,k,m):
                        uu,vv=float(u(ss)),float(v(ss))
                        if 0<uu<1 and 0<vv<1:add(k,[uu]*pp+[vv]*qq,'two_free_values')
                else:
                    ell=pp;delta=1-10*ell/(n+8);assert abs(delta)>1e-12
                    ss=(k-2*ell*B4/(3*B6))/delta
                    if not k<ss<k+m:continue
                    tt=10*ss/(n+8)-2*B4/(3*B6)
                    uv=(p*ss**(-n/2)-2*a*B4*ss-b*B6*(k+ell*tt*tt)+3*g*B6*ss*ss)/(B6*(3-2*b*ell))
                    disc=tt*tt-4*uv
                    if disc>=0:
                        uu=(tt+math.sqrt(disc))/2;vv=(tt-math.sqrt(disc))/2
                        if 0<uu<1 and 0<vv<1:add(k,[uu]*ell+[vv]*ell,'equal_multiplicity')
    value,y,kind=min(candidates,key=lambda z:z[0])
    return value,{'squared_point':y,'point':np.sqrt(y).tolist(),'kind':kind,'candidate_count':len(candidates)}

@lru_cache(maxsize=4096)
def radius_oracle(n,B4,B6):
    lo=math.sqrt(n)+.01;hi=8.
    def env(r):return r**(2-n)-sphere_extrema(n,B4,B6,r)[0][0]
    opt=minimize_scalar(env,bounds=(lo,hi),method='bounded',options={'maxiter':70,'xatol':1e-8})
    level,r=min([(float(opt.fun),float(opt.x)),(env(lo),lo),(env(hi),hi)])
    return r,level,{'radius_nfev':int(opt.nfev)+2,'radius_success':bool(opt.success),'radius_interval':[lo,hi]}

@lru_cache(maxsize=4)
def directions(n,fine=False):
    if n==3:
        nu,nf=(48,96) if fine else (24,48)
        z,w=leggauss(nu);u=(z+1)/2;wu=w/2
        z,w=leggauss(nf);f=(z+1)*np.pi/4;wf=w*np.pi/4
        U,F=np.meshgrid(u,f,indexing='ij')
        y=np.array([(1-U*U)*np.cos(F)**2,(1-U*U)*np.sin(F)**2,U*U]).reshape(n,-1).T
        weights=(8*wu[:,None]*wf[None,:]).ravel()
    else:
        nu,nv,nf=(16,16,32) if fine else (12,12,24)
        z,w=roots_jacobi(nu,.5,0);u=(z+1)/2;wu=w/(2*np.sqrt(2))*np.sqrt(1+u)
        z,w=leggauss(nv);v=(z+1)/2;wv=w/2
        z,w=leggauss(nf);f=(z+1)*np.pi/4;wf=w*np.pi/4
        U,V,F=np.meshgrid(u,v,f,indexing='ij')
        y=np.array([(1-U*U)*(1-V*V)*np.cos(F)**2,(1-U*U)*(1-V*V)*np.sin(F)**2,(1-U*U)*V*V,U*U]).reshape(n,-1).T
        weights=(16*wu[:,None,None]*wv[None,:,None]*wf[None,None,:]).ravel()
    h4,h6=hs(y,n)
    return h4,h6,weights

def integral(n,A,B4,B6,R,fine=False):
    h4,h6,w=directions(n,fine);lo=np.zeros_like(h4);hi=np.full_like(h4,R)
    for _ in range(60):
        mid=(lo+hi)/2;positive=mid**(2-n)-A-B4*h4*mid**4-B6*h6*mid**6>0
        lo=np.where(positive,mid,lo);hi=np.where(positive,hi,mid)
    r=(lo+hi)/2
    return float(w@(r*r/2-A*r**n/n-B4*h4*r**(n+4)/(n+4)-B6*h6*r**(n+6)/(n+6)))

def evaluate(n,v,fine=False):
    A,B4,B6=map(float,v)
    if not(.02<=A<=1.25 and -.05<=B4<=.05 and -.02<=B6<=.02):return None
    R,outer,oracle=radius_oracle(n,B4,B6);m,loc=cube_candidates(n,B4,B6)
    if A<=outer+1e-6 or A>=m-1e-6:return None
    M=integral(n,A,B4,B6,R,fine)
    row={'n':n,'A':A,'B4':B4,'B6':B6,'R':R,'outer_negative_margin':A-outer,
         'cube_min_before_A_numeric':m,'cube_positive_margin_numeric':m-A,'cube_min_candidate':loc,
         'mass_dimensionless':M,'cost':M/(2**n*(m-A)),'sphere_nodes':len(directions(n,fine)[0]),**oracle}
    if row['cost']<(8 if n==3 else 16):raise Found(row)
    return row

def run_start(n,seed,budget):
    v=np.array(seed,dtype=float);R,outer,_=radius_oracle(n,v[1],v[2]);m,_=cube_candidates(n,v[1],v[2])
    if m-outer<=4e-6:raise ValueError('registered seed has no feasible radius')
    if v[0]<=outer+1e-6 or v[0]>=m-1e-6:v[0]=outer+.2*(m-outer)
    projected=v.tolist();best=None;valid=0
    def objective(x):
        nonlocal best,valid
        row=evaluate(n,x)
        if row:
            valid+=1
            if best is None or row['cost']<best['cost']:best=row
            return row['cost']
        return 1e4+float(np.sum(np.square(x)))
    simplex=np.vstack([v,v+[.005,0,0],v+[0,.0005,0],v+[0,0,.0001]])
    opt=minimize(objective,v,method='Nelder-Mead',options={'maxfev':budget,'initial_simplex':simplex,'xatol':2e-8,'fatol':1e-8})
    return {'n':n,'seed':list(seed),'projected_seed':projected,'nfev':int(opt.nfev),'success':bool(opt.success),'valid_evaluations':valid,'best':best}

def main():
    starts={3:[(.4089479272201547,.005004477192161995,.0009386738283819751),(.4012454957860509,.004,.0005)],4:[(.1289749174365545,.001,0.),(.1289749174365545,.001,.0001)]}
    runs=[];refined=[];checks=[];early=None
    try:
        for n in (3,4):
            for seed in starts[n]:runs.append(run_start(n,seed,360 if n==3 else 300))
            best=min((z['best'] for z in runs if z['n']==n and z['best']),key=lambda z:z['cost'])
            refined.append(evaluate(n,[best['A'],best['B4'],best['B6']],True))
            M,info=allroot_mass(n,best['A'],best['B4'],best['B6'],best['R'],True)
            checks.append({'n':n,'allroot_mass_same_coarse_directions':M,'bisection_mass':best['mass_dimensionless'],**info})
    except Found as e:early=e.row
    src=Path(__file__)
    data={'batch':'U07','status':'floating free-radius discovery, no continuous upper certificate','runs':runs,'refinements':refined,'allroot_crosschecks':checks,'threshold_candidate':early,'source_sha256':{src.name:hashlib.sha256(src.read_bytes()).hexdigest(),'u06_h46_scan.py':hashlib.sha256(src.with_name('u06_h46_scan.py').read_bytes()).hexdigest()}}
    src.with_name('u07_results.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2))

if __name__=='__main__':main()
