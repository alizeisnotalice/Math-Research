#!/usr/bin/env python3
"""One registered MC probe and an independent locked-moment Gauss upper screen."""
import hashlib
import itertools
import json
import math
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import arb,ctx

P=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_upper016')
START=time.monotonic()

def digest(path):
    p=Path(path)
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}

def dump(name,data):
    tmp=P/(name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    tmp.replace(P/name)

def ar(q):
    q=F(q)
    return arb(q.numerator)/q.denominator

def gauss():
    manifest=json.loads((SRC/'manifest.json').read_text())
    assert manifest['status']=='locked_two_necessary_screens_and_selected_m3_cubic_minorant_ceiling'
    for item in manifest['white_list']:
        assert digest(SRC/item['path'])['sha256']==item['sha256'],item['path']
    screen=json.loads((SRC/'corner_screen.json').read_text())
    herm=json.loads((SRC/'hermite_moments.json').read_text())
    assert screen['beta']=='1/2' and screen['completed_W_power']==9
    assert herm['completed_W_power']==6
    vals=[screen['moments']['one']]+[herm['moments'][str(k)] for k in (1,2,3)]
    mu=[sum((F(math.comb(6,n)*2**(6-n))*F(v[n]) for n in range(7)),F(0))/12**6 for v in vals]
    rows=json.loads((SRC/'screen_verification.json').read_text())['rows']
    row=next(v for v in rows if v['m']==3)
    assert mu[0]==F(row['E_u_2m']) and mu[1]==F(row['E_d_u_2m'])
    v=F(row['E_B_u_3m'])
    assert v==F(965344619,31225862828160)
    variance=mu[0]*mu[2]-mu[1]**2
    assert variance>0
    A=(mu[0]*mu[3]-mu[1]*mu[2])/variance
    B=(mu[1]*mu[3]-mu[2]**2)/variance
    delta=A*A-4*B
    assert delta>0
    ctx.prec=256
    nodes=[(ar(A)-ar(delta).sqrt())/2,(ar(A)+ar(delta).sqrt())/2]
    z1,z2=nodes
    weights=[(ar(mu[0])*z2-ar(mu[1]))/(z2-z1),(ar(mu[1])-ar(mu[0])*z1)/(z2-z1)]
    assert 0<z1<z2<1 and all(w>0 for w in weights)
    for k in range(4):
        assert sum((w*z**k for w,z in zip(weights,nodes)),arb(0)).contains(ar(mu[k]))
    Ucap=sum((w*z.sqrt() for w,z in zip(weights,nodes)),arb(0))
    c1=1/(32*arb.pi()**3)
    required=ar(F(15198363,100000000000))
    Ucrit=(required*ar(v)**2/c1)**(arb(1)/3)
    gcap=c1*Ucap**3/ar(v)**2
    return {'status':'locked_moment_Gauss_upper_screen','source_lock':digest(SRC/'manifest.json'),
            'all_11_lock_hashes_verified':True,'mu_0_to_3':[str(x) for x in mu],
            'V':str(v),'gauss_polynomial':['1',str(-A),str(B)],
            'nodes':[str(x) for x in nodes],'weights':[str(x) for x in weights],
            'U_upper':str(Ucap),'Ucrit':str(Ucrit),'all_a_g_upper':str(gcap),
            'strictly_excludes_required_g':bool(gcap<required),'U_upper_float':float(Ucap),
            'Ucrit_float':float(Ucrit),'g_upper_float':float(gcap),
            'precision_bits':ctx.prec,'seconds':time.monotonic()-START}

def monte_carlo():
    vectors=np.array([(1,)+v for v in itertools.product((-1,1),repeat=4)],dtype=float)
    coefficient=((vectors@np.ones(5))**2-1)/24
    assert np.count_nonzero(coefficient==1)==1 and np.count_nonzero(coefficient==1/3)==5
    assert np.count_nonzero(coefficient==0)==10
    rng=np.random.default_rng(26026)
    n=32768
    values=np.empty(n);weights=np.empty(n);conds=np.empty(n);mineigs=np.empty(n);maxeigs=np.empty(n)
    logdets=np.empty(n);us=np.empty(n)
    ztilt=1/(16*math.gamma(1.5))
    fallback=0;nonpositive_float=0;nonfinite=0
    done=0
    if (P/'checkpoint_arrays.npz').exists():
        cp=json.loads((P/'checkpoint.json').read_text())
        assert cp['script_sha256']==digest(__file__)['sha256']
        saved=np.load(P/'checkpoint_arrays.npz')
        done=cp['completed_samples'];rng.bit_generator.state=cp['rng_state']
        for name,dest in [('values',values),('weights',weights),('conds',conds),('mineigs',mineigs),('maxeigs',maxeigs),('logdets',logdets),('us',us)]:
            dest[:done]=saved[name]
        fallback=cp['fallback'];nonpositive_float=cp['nonpositive_float'];nonfinite=cp['nonfinite']
    for start in range(done,n,4096):
        if time.monotonic()-START>120:
            return None
        stop=min(start+4096,n)
        y=rng.gamma(shape=.5,scale=1,size=(stop-start,16))
        assert np.all(y>0) and np.all(np.isfinite(y))
        S=1/(4*y)
        T=np.sum(S,axis=1)
        M=np.einsum('ni,ia,ib->nab',S/T[:,None],vectors,vectors,optimize=True)
        eig=np.linalg.eigvalsh(M)
        signs,logdet=np.linalg.slogdet(M)
        condition=eig[:,-1]/eig[:,0]
        bad=(signs<=0)|(eig[:,0]<=0)|(~np.isfinite(logdet))|(condition>1e12)
        nonpositive_float+=int(np.count_nonzero((signs<=0)|(eig[:,0]<=0)))
        nonfinite+=int(np.count_nonzero(~np.isfinite(logdet)))
        if np.any(bad):
            import mpmath as mp
            mp.mp.dps=60
            for j in np.flatnonzero(bad):
                ss=[mp.mpf(float(v)) for v in S[j]];total=sum(ss)
                mat=mp.matrix([[sum(ss[k]*int(vectors[k,a])*int(vectors[k,b]) for k in range(16))/total for b in range(5)] for a in range(5)])
                dd=mp.det(mat)
                assert dd>0
                ee=mp.eigsy(mat,eigvals_only=True)
                assert ee[0]>0
                logdet[j]=float(mp.log(dd));eig[j,0]=float(ee[0]);eig[j,-1]=float(ee[-1])
                condition[j]=float(ee[-1]/ee[0]);fallback+=1
        u=(S@coefficient)/T
        weight=T**(-.5)/ztilt
        value=weight*np.exp(logdet/2)*u**6
        assert np.all(np.isfinite(value)) and np.all(value>=0)
        values[start:stop]=value;weights[start:stop]=weight;conds[start:stop]=condition
        mineigs[start:stop]=eig[:,0];maxeigs[start:stop]=eig[:,-1];logdets[start:stop]=logdet;us[start:stop]=u
        np.savez(P/'checkpoint_arrays.npz',values=values[:stop],weights=weights[:stop],conds=conds[:stop],mineigs=mineigs[:stop],maxeigs=maxeigs[:stop],logdets=logdets[:stop],us=us[:stop])
        dump('checkpoint.json',{'completed_samples':stop,'pid':os.getpid(),'script_sha256':digest(__file__)['sha256'],
                               'rng_state':rng.bit_generator.state,'fallback':fallback,'nonpositive_float':nonpositive_float,'nonfinite':nonfinite})
    v=float(F(965344619,31225862828160));c1=1/(32*math.pi**3)
    mean=float(np.mean(values));se=float(np.std(values,ddof=1)/math.sqrt(n))
    ucrit=(float(F(15198363,100000000000))*v*v/c1)**(1/3)
    return {'status':'one_seed_IID_MC_diagnostic_not_certificate','seed':26026,'sample_count':n,
            'Ztilt':ztilt,'U':mean,'sample_standard_error_not_a_rigorous_interval':se,'Ucrit':ucrit,'U_over_Ucrit':mean/ucrit,
            'all_a_g_from_sample_U_not_bound':c1*mean**3/v**2,'tilt_weight_mean':float(np.mean(weights)),
            'tilt_weight_SE':float(np.std(weights,ddof=1)/math.sqrt(n)),
            'tilt_effective_sample_size':float(np.sum(weights)**2/np.dot(weights,weights)),
            'eigenvalue_min':float(np.min(mineigs)),'eigenvalue_max':float(np.max(maxeigs)),
            'condition_quantiles':{str(q):float(np.quantile(conds,q)) for q in (0,.5,.9,.99,.999,1)},
            'determinant_range':list(map(float,np.exp([np.min(logdets),np.max(logdets)]))),
            'u_range':[float(np.min(us)),float(np.max(us))],
            'high_precision_fallback_count':fallback,'nonpositive_float_eigen_or_det_count':nonpositive_float,
            'nonfinite_float_det_count':nonfinite,'clipping':False,'self_normalized':False,
            'seconds':time.monotonic()-START}

def main():
    assert json.loads((P/'registration.json').read_text())['status']=='registered_before_computation'
    mc=monte_carlo()
    if mc is None:
        print('120 second checkpoint; rerun to resume',flush=True)
        return
    dump('mc_result.json',mc)
    print(json.dumps(mc,indent=2),flush=True)
    exact=gauss()
    dump('gauss_upper.json',exact)
    print(json.dumps(exact,indent=2),flush=True)
    dump('run_manifest.json',{'script':digest(__file__),'registration':digest(P/'registration.json'),
             'source_files':[digest(SRC/x) for x in ('manifest.json','corner_screen.json','hermite_moments.json','screen_verification.json')],
             'python':sys.executable,'pid':os.getpid(),'seconds':time.monotonic()-START,'new_global_g':False,'new_C5_upper_bound':False})

if __name__=='__main__':main()
