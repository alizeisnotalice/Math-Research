"""Finite discovery only: true jump-rate killed stable generator, not spectral power.

Generator normalization: jump density |z|^(-1-alpha), global c_alpha omitted.
Cell source/density approximation; no claim of a continuous Green certificate.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh, toeplitz
from scipy.integrate import quad_vec

ALPHAS=(.5,.8,1.,1.2,1.5,1.8)
GRIDS=(65,129)

def generator(alpha,ncells):
    h=2/ncells
    x=(np.arange(ncells)-(ncells-1)/2)*h
    k=np.arange(1,ncells,dtype=float)
    rates=h**(-alpha)/alpha*((k-.5)**(-alpha)-(k+.5)**(-alpha))
    first=np.r_[-2**(alpha+1)*h**(-alpha)/alpha,rates]
    q=toeplitz(first)
    killing=((1-x)**(-alpha)+(1+x)**(-alpha))/alpha
    defect=float(np.max(np.abs(q.sum(axis=1)+killing)))
    assert np.all(q[np.triu_indices(ncells,1)]>0)
    assert defect<1e-8
    return q,x,h,killing,defect

def scan(alpha,ncells):
    q,x,h,killing,defect=generator(alpha,ncells)
    lam,v=eigh(-q)
    assert lam[0]>0
    src=ncells//2
    targets=np.flatnonzero((x>=.25)&(x<=.75))
    coeff=(v[src,:][None,:]*v[targets,:])
    scoeff=v[src,:]*v.sum(axis=0)
    scaled=lam/lam[0]
    def integrand(u):
        decay=np.exp(-scaled*u)
        p=(coeff@decay)/h
        s=float(scoeff@decay)
        return np.r_[s**3,s**4,p**3,p**4]
    values,error=quad_vec(integrand,0,np.inf,epsabs=2e-10,epsrel=2e-9,limit=300)
    nt=len(targets)
    rows=[]
    for n,offset in ((3,2),(4,2+nt)):
        mass=values[n-3]/lam[0]
        green=values[offset:offset+nt]/lam[0]
        costs=mass/(2**n*x[targets]**n*green)
        best=int(np.argmin(costs))
        rows.append({'n':n,'alpha':alpha,'N':ncells,'h':h,
                     'a':float(x[targets[best]]),'R':float(1/x[targets[best]]),
                     'cost':float(costs[best]),'unit_mass':float(mass),
                     'unit_corner_green':float(green[best]),
                     'a_count':nt,'a_first':float(x[targets[0]]),'a_last':float(x[targets[-1]]),
                     'quad_vec_error_estimate_scaled':float(error),
                     'target_improved_discretely':bool(costs[best]<(8 if n==3 else 16))})
    return rows,{'alpha':alpha,'N':ncells,'lambda_min':float(lam[0]),
                 'row_killing_defect':defect,'source_reconstruction_error':float(np.max(np.abs(v[src]@v.T-np.eye(ncells)[src]))),
                 'killing_min':float(killing.min()),'killing_max':float(killing.max())}

def main():
    rows=[];checks=[]
    for alpha in ALPHAS:
        for ncells in GRIDS:
            r,c=scan(alpha,ncells);rows.extend(r);checks.append(c)
    result={'batch':'U03','status':'finite discrete discovery, not continuous bounds',
            'generator':'cell-integrated full-space |z|^(-1-alpha) jumps with outside killing',
            'density_normalization':'p=P/h; dt scale lambda_min cancels in ratio',
            'rows':rows,'matrix_checks':checks,
            'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
