#!/usr/bin/env python3
"""Two predeclared rho values: exact V and one shared finite MC sample for U."""
import hashlib,itertools,json,math,os,sys,time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import sympy as sp
from flint import arb,ctx

P=Path(__file__).resolve().parent
START=time.monotonic()
RHOS=(F(7,8),F(31,32))

def digest(path):
    p=Path(path)
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}

def dump(name,obj):
    tmp=P/(name+'.tmp');tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n');tmp.replace(P/name)

def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator

def exact():
    ctx.prec=256
    vec=np.array([(1,)+v for v in itertools.product((-1,1),repeat=4)],dtype=np.int64)
    dot=np.abs(vec@np.ones(5,dtype=np.int64))
    I=np.eye(5,dtype=np.int64);J=np.ones((5,5),dtype=np.int64)
    assert np.array_equal(vec[dot==5].T@vec[dot==5],J)
    assert np.array_equal(vec[dot==3].T@vec[dot==3],4*I+J)
    assert np.array_equal(vec[dot==1].T@vec[dot==1],12*I-2*J)
    # Exact Mellin/Frullani normalization checks, before floating diagnostics.
    zinv=16*sp.gamma(sp.Rational(3,2))
    vpref=sp.simplify(zinv*sp.Rational(1,16)*2*sp.gamma(5)/sp.gamma(sp.Rational(9,2)))
    upref=sp.simplify(zinv*16**sp.Rational(5,2)/(sp.gamma(3)*2*sp.sqrt(sp.pi)*sp.pi**2))
    assert vpref==sp.Rational(128,35)
    assert upref==2048/sp.pi**2
    assert sp.Rational(128,35)*5*16**4/16**5==sp.Rational(8,7)
    rows=[]
    c1=1/(32*arb.pi()**3)
    target=ar(F(15198363,100000000000))
    for rho in RHOS:
        small=1-rho;middle=1-rho/3
        C=ar(small).sqrt()+5*ar(middle).sqrt()+10
        Dperp=12+4/ar(middle).sqrt()
        adjx=5*Dperp**4
        V=ar(F(128,35))*adjx/C**5
        Ucrit=(target*V**2/c1)**(arb(1)/3)
        rows.append({'rho':str(rho),'lambda_special':str(small),'lambda_middle':str(middle),'lambda_other':'1',
                     'C':str(C),'D_perp':str(Dperp),'x_adj_D_x':str(adjx),'V':str(V),'V_float':float(V),
                     'Ucrit':str(Ucrit),'Ucrit_float':float(Ucrit)})
    return {'status':'exact_normalizations_and_Arb_V','rho_zero_V_exact':'8/7',
            'V_prefactor_exact':str(vpref),'U_CB_prefactor_exact':str(upref),
            'category_Gram_sums_exactly_verified':True,'rows':rows,'Arb_bits':ctx.prec}

def sample(exact_result):
    vec=np.array([(1,)+v for v in itertools.product((-1,1),repeat=4)],dtype=float)
    rr=((vec@np.ones(5))**2-1)/24
    rho=np.array([float(x) for x in RHOS])
    n=32768
    values=np.empty((n,2));weights=np.empty(n);conds=np.empty(n);mins=np.empty(n);maxs=np.empty(n);dets=np.empty(n);us=np.empty(n)
    rng=np.random.default_rng(30030);done=0;fallback=0;float_nonpositive=0
    if (P/'checkpoint_arrays.npz').exists():
        cp=json.loads((P/'checkpoint.json').read_text());assert cp['script_sha256']==digest(__file__)['sha256']
        done=cp['completed'];rng.bit_generator.state=cp['rng_state'];fallback=cp['fallback'];float_nonpositive=cp['float_nonpositive']
        a=np.load(P/'checkpoint_arrays.npz')
        for key,dest in [('values',values),('weights',weights),('conds',conds),('mins',mins),('maxs',maxs),('dets',dets),('us',us)]:dest[:done]=a[key]
    ztilt=1/(16*math.gamma(1.5))
    for lo in range(done,n,4096):
        if time.monotonic()-START>120:return None
        hi=min(lo+4096,n)
        Y=rng.gamma(.5,1,size=(hi-lo,16));assert np.all(Y>0)
        S=1/(4*Y);T=np.sum(S,axis=1)
        M=np.einsum('ni,ia,ib->nab',S/T[:,None],vec,vec,optimize=True)
        ee=np.linalg.eigvalsh(M);sg,ld=np.linalg.slogdet(M);cond=ee[:,-1]/ee[:,0]
        bad=(ee[:,0]<=0)|(sg<=0)|(~np.isfinite(ld))|(cond>1e12)
        float_nonpositive+=int(np.count_nonzero((ee[:,0]<=0)|(sg<=0)))
        if np.any(bad):
            import mpmath as mp
            mp.mp.dps=60
            for j in np.flatnonzero(bad):
                ss=[mp.mpf(float(q)) for q in S[j]];total=sum(ss)
                mat=mp.matrix([[sum(ss[k]*int(vec[k,a])*int(vec[k,b]) for k in range(16))/total for b in range(5)] for a in range(5)])
                eig=mp.eigsy(mat,eigvals_only=True);det=mp.det(mat);assert det>0 and eig[0]>0
                ld[j]=float(mp.log(det));ee[j,0]=float(eig[0]);ee[j,-1]=float(eig[4]);cond[j]=float(eig[4]/eig[0]);fallback+=1
        u=(S@rr)/T;L=1-u[:,None]*rho[None,:];assert np.all(L>0)
        weight=T**(-.5)/ztilt
        value=weight[:,None]*np.exp(ld/2)[:,None]/L**3
        assert np.all(np.isfinite(value)) and np.all(value>=0)
        values[lo:hi]=value;weights[lo:hi]=weight;conds[lo:hi]=cond;mins[lo:hi]=ee[:,0];maxs[lo:hi]=ee[:,-1];dets[lo:hi]=np.exp(ld);us[lo:hi]=u
        np.savez(P/'checkpoint_arrays.npz',values=values[:hi],weights=weights[:hi],conds=conds[:hi],mins=mins[:hi],maxs=maxs[:hi],dets=dets[:hi],us=us[:hi])
        dump('checkpoint.json',{'completed':hi,'pid':os.getpid(),'script_sha256':digest(__file__)['sha256'],
                              'rng_state':rng.bit_generator.state,'fallback':fallback,'float_nonpositive':float_nonpositive})
    rows=[]
    c1=1/(32*math.pi**3)
    for j,ex in enumerate(exact_result['rows']):
        U=float(np.mean(values[:,j]));se=float(np.std(values[:,j],ddof=1)/math.sqrt(n));V=ex['V_float']
        rows.append({'rho':str(RHOS[j]),'U_diagnostic':U,'sample_SE_not_CI':se,'V_certified_expression_float':V,
                     'Ucrit':ex['Ucrit_float'],'U_over_Ucrit':U/ex['Ucrit_float'],
                     'g_opt_diagnostic_not_bound':c1*U**3/V**2,'sample_max_U_integrand':float(np.max(values[:,j]))})
    return {'status':'one_seed_two_rho_MC_diagnostic_not_certificate','seed':30030,'samples':n,'rows':rows,
            'tilt_mean':float(np.mean(weights)),'tilt_SE':float(np.std(weights,ddof=1)/math.sqrt(n)),
            'tilt_ESS':float(np.sum(weights)**2/np.dot(weights,weights)),
            'condition_quantiles':{str(q):float(np.quantile(conds,q)) for q in (0,.5,.9,.99,.999,1)},
            'eigenvalue_min':float(np.min(mins)),'eigenvalue_max':float(np.max(maxs)),
            'det_range':[float(np.min(dets)),float(np.max(dets))],'u_range':[float(np.min(us)),float(np.max(us))],
            'fallback_count':fallback,'nonpositive_float_count':float_nonpositive,'clipping':False,'self_normalized':False}

def main():
    assert json.loads((P/'registration.json').read_text())['status']=='registered_before_computation'
    ex=exact();dump('exact_V.json',ex);print(json.dumps(ex,indent=2),flush=True)
    mc=sample(ex)
    if mc is None:
        print('120 second checkpoint; rerun to resume',flush=True);return
    dump('mc_result.json',mc);print(json.dumps(mc,indent=2),flush=True)
    dump('run_manifest.json',{'script':digest(__file__),'registration':digest(P/'registration.json'),
                             'pid':os.getpid(),'python':sys.executable,'seconds':time.monotonic()-START,
                             'new_global_g':False,'new_C5_upper_bound':False})

if __name__=='__main__':main()
