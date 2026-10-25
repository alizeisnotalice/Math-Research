#!/usr/bin/env python3
"""Finite diagnostics for the r=2b trace tilt and exact positive series."""
import json
import math
import os
import time
import numpy as np
from scipy.special import gamma, gammaln
from scipy.stats import qmc
from rank2_angle import ROOT, TRIPLES, WEIGHTS, geometry, digest, dump

START = time.monotonic()
BETA = .75
B = .15
RHO = 247/256

def classes(triple):
    exact, gram = geometry(triple)
    grouped = {}
    import sympy as sp
    for i in range(16):
        x, y = exact.row(i)
        norm = sp.simplify(x*x + y*y)
        if norm == 0:
            continue
        z = sp.simplify((x+sp.I*y)**2/norm)
        coefficient = (float(norm)/16)**BETA
        grouped[z] = grouped.get(z, 0) + coefficient
    z = np.array([complex(a) for a in grouped])
    c = np.array(list(grouped.values()))
    return z, c, gram

def sample(z, c):
    out = []
    for seed in (1025, 2025):
        u = qmc.Sobol(2*len(z), scramble=True, seed=seed).random_base2(16)
        angle = math.pi*u[:, :len(z)]
        expvar = -np.log(u[:, len(z):])
        log_s = np.log(np.sin(BETA*angle)) - np.log(np.sin(angle))/BETA
        log_s += (1-BETA)/BETA*(np.log(np.sin((1-BETA)*angle))-np.log(expvar))
        x = np.exp(log_s)*c[None,:]**(1/BETA)
        det = np.zeros(len(u))
        for i in range(len(z)):
            for j in range(i):
                det += (1-(z[i]*z[j].conjugate()).real)/2*x[:,i]*x[:,j]
        out.append({'seed': seed, 'sample_count':len(u),
                    'E_D_minus_b':float(np.mean(det**(-B))),
                    'E_D_minus_b_rho':float(np.mean(det**(-B*RHO))),
                    'E_D_minus_half':float(np.mean(det**(-.5)))})
    return out

def series(z,c,N,s):
    # r=2s is the exact trace tilt, not beta.
    r=2*s
    p=c/sum(c)
    phi=2*np.pi*np.arange(2*N+1)/(2*N+1)
    direction=(z[:,None]*np.exp(-1j*phi)[None,:]).real.astype(np.longdouble)
    h=np.zeros((2*N+1,2*N+1),dtype=np.longdouble)
    h[0]=1
    powa=np.ones_like(direction)
    binom_signed=np.longdouble(1)
    for n in range(1,2*N+1):
        powa*=direction
        binom_signed*=np.longdouble(n-1-BETA)/n
        h[n]=binom_signed*np.sum(p[:,None]*powa,axis=0)
    g=np.zeros_like(h)
    g[0]=1
    q=r/BETA
    for n in range(1,2*N+1):
        j=np.arange(1,n+1,dtype=np.longdouble)
        g[n]=-np.sum((n-j+q*j)[:,None]*h[1:n+1]*g[n-1::-1],axis=0)/n
    mu=np.ones(N+1)
    for k in range(1,N+1):
        logfactor=k*math.log(4)+2*gammaln(k+1)-gammaln(r+2*k)+gammaln(r)
        mu[k]=float(np.mean(g[2*k]))*math.exp(logfactor)
    coeff=np.ones(N+1)
    for k in range(1,N+1):
        coeff[k]=coeff[k-1]*(s+k-1)/k
    factor=2**r*gamma(1+r/BETA)/gamma(1+r)*sum(c)**(-r/BETA)
    return {'N':N,'s':s,'trace_tilt':r,'partial_sum':float(factor*np.dot(coeff,mu)),
            'prefactor':factor,'mu_first_6':mu[:6].tolist(),'mu_last':float(mu[-1]),
            'mu_monotone_diagnostic':bool(np.all(np.diff(mu)<=1e-14)),
            'mu_nonnegative_diagnostic':bool(np.min(mu)>=0),
            'last_increment':float(factor*coeff[-1]*mu[-1])}

def main():
    reg=ROOT/'joint_registration.json'
    assert json.loads(reg.read_text())['status']=='registered_before_computation'
    rows=[]
    script=digest(__file__)
    checkpoint=ROOT/'joint_checkpoint.json'
    if checkpoint.exists():
        prior=json.loads(checkpoint.read_text())
        if prior.get('script_sha256')==script['sha256']:
            rows=prior['representatives']
    for triple in TRIPLES[len(rows):]:
        if time.monotonic()-START>120:
            return
        z,c,gram=classes(triple)
        item={'triple':triple,'gram_determinant':gram,'line_classes':len(z),
              'spectral_weights':(c/sum(c)).tolist(),'stable_total_scale_power':float(sum(c)),
              'QMC_nonrigorous':sample(z,c),'series':[]}
        for s in (B,B*RHO):
            for N in (64,128):
                item['series'].append(series(z,c,N,s))
        rows.append(item)
        dump('joint_checkpoint.json',{'script_sha256':script['sha256'],'pid':os.getpid(),
                                    'representatives':rows,'next_representative':len(rows)})
        print(json.dumps(item),flush=True)
    m=lambda s: gamma(1+s/BETA)/gamma(1+s)
    qmc_means=[np.mean([x['E_D_minus_b_rho'] for x in r['QMC_nonrigorous']]) for r in rows]
    j=2*gamma(1.75)*m(B*RHO)**(3/RHO)*math.exp(sum(t/RHO*math.log(x) for t,x in zip(WEIGHTS,qmc_means)))
    c0=3/(2**7.5*math.pi**2.5)
    result={'status':'nonrigorous_joint_diagnostic_complete','script':script,'registration':digest(reg),
            'angle_script':digest(ROOT/'rank2_angle.py'),'representatives':rows,
            'joint_J_fractional_diagnostic':j,'joint_I_fractional_diagnostic':c0**(10/7)*math.pi**2/60*j**(10/7),
            'seconds':time.monotonic()-START,'pid':os.getpid(),'new_C5_upper_bound':False,
            'warning':'QMC estimates are not bounds. Finite positive series needs a certified singular tail for an upper bound.'}
    dump('joint_diagnostic.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='representatives'},indent=2),flush=True)

if __name__=='__main__':
    main()
