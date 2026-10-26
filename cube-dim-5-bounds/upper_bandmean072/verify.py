#!/usr/bin/env python3
"""Readback with independent box moments and algebraic R, no new H calls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import hashlib,json,os,signal,sys,time
from flint import arb,ctx,fmpq

HERE=Path(__file__).resolve().parent
if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ar(x):
    x=F(x)
    return arb(fmpq(x.numerator,x.denominator))


def ep(x,upper=False):
    return F(str((x.upper() if upper else x.lower()).fmpq()))


def main():
    start=time.monotonic()
    def stop(signum,frame):
        raise SystemExit(75)
    signal.signal(signal.SIGALRM,stop)
    signal.setitimer(signal.ITIMER_REAL,175)
    ctx.prec=384
    reg=json.loads((HERE/'registration.json').read_text())
    result=json.loads((HERE/'results.json').read_text())
    receipt=json.loads((HERE/'run_receipt.json').read_text())
    cp=json.loads((HERE/'checkpoint.json').read_text())
    assert sha(HERE/'mechanism.md')==reg['mechanism_sha256']
    for src in reg['sources'].values():
        assert sha(src['path'])==src['sha256']
    assert sha(HERE/'certify.py')==receipt['source_sha256']
    assert sha(HERE/'registration.json')==receipt['registration_sha256']
    assert sha(HERE/'results.json')==receipt['results_sha256']
    assert sha(HERE/'checkpoint.json')==receipt['checkpoint_sha256']
    assert result['fingerprint']==cp['fingerprint']==receipt['fingerprint']
    olddir=Path(reg['sources']['manifest.json']['path']).parent
    for row in json.loads((olddir/'manifest.json').read_text())['files']:
        assert sha(olddir/row['name'])==row['sha256']
    old=json.loads(Path(reg['sources']['results.json']['path']).read_text())
    oldcp=json.loads(Path(reg['sources']['checkpoint.json']['path']).read_text())
    Hlo=F(oldcp['rows'][0]['H']['lower'])
    center=list(map(F,reg['center']))
    Q=[[F(i==j,4)-F(1,36) for j in range(4)] for i in range(4)]
    Mtx=[[1,1,1,1],[-1,-1,-1,2],[-1,-1,2,-1],[-1,2,-1,-1],[2,-1,-1,-1]]
    caps=[('6731619/100000000','178899/100000000000000'),
          ('2576953/10000000000','6133/50000000000000000'),
          ('248023/2000000000','305013/100000000000000000000')]
    outputs=[]
    for ci,(registered,r) in enumerate(zip(reg['cells'],result['cells'])):
        h=F(r['h'])
        assert r['h']==registered['h']
        lower=[max(F(0),c-h)-c for c in center[1:]]
        upper=[min(F(1),c+h)-c for c in center[1:]]
        assert lower==list(map(F,r['u_lower'])) and upper==list(map(F,r['u_upper']))
        mu=[(l+u)/2 for l,u in zip(lower,upper)]
        variance=[(u-l)**2/12 for l,u in zip(lower,upper)]
        eta=sum(Q[i][i]*variance[i] for i in range(4))
        eta+=sum(Q[i][j]*mu[i]*mu[j] for i in range(4) for j in range(4))
        lam=sum(abs(sum(F(z,6)*v for z,v in zip(row,mu))) for row in Mtx)
        assert mu==list(map(F,r['mu'])) and eta==F(r['eta']) and lam==F(r['lambda'])
        assert all(variance[i]+mu[i]*mu[i]==F(r['second_moments'][i]) for i in range(4))
        volume=prod(u-l for l,u in zip(lower,upper))
        assert volume==F(r['volume'])
        cornermax=max(sum(Q[i][j]*u[i]*u[j] for i in range(4) for j in range(4))
                      for u in product(*zip(lower,upper)))
        assert cornermax==h*h==F(r['maximum_metric_corner'])
        if ci==0:
            assert mu==[-F(26137,26810624),0,0,0] and lam==-mu[0]
        else:
            assert all(v==0 for v in mu) and eta==F(8,27)*h*h and lam==0
        mean=mean2=L=L2=U=U2=F(0)
        octave_mean=[F(0) for _ in range(12)]
        first_taylor=None
        for bi,(b,v) in enumerate(zip(old['bands'],r['bands'])):
            assert v==cp['cells'][ci]['bands'][bi]
            t=F(b['start'])
            R=F(v['R_upper'])
            coef=min(R,1+lam/t+F(21,8)*R*eta/(t*t))
            assert coef==F(v['chosen_mean_coefficient'])
            assert F(v['Taylor_mean_coefficient'])==1+lam/t+F(21,8)*R*eta/(t*t)
            assert v['chosen']==('Taylor_mean' if coef<R else 'range')
            if v['chosen']=='Taylor_mean' and first_taylor is None:
                first_taylor=b['start']
            ml,mu_mass=F(b['mass_lower']),F(b['mass_upper'])
            assert 0<ml<=mu_mass and 1<=coef<=R
            mean+=mu_mass*coef
            octave_mean[bi//8]+=mu_mass*coef
            if registered['reuse_old_range']:
                assert R==F(b['R_upper']) and v['R_source']=='reused067'
                R2=R
            else:
                x=1+h*h/(10*t*t)
                P=16*x**5-20*x**3+5*x
                R2=ep(ar(P)+ar(P*P-1).sqrt(),True)
                L+=ml/R
                U+=mu_mass*R
                L2+=ml/R2
                U2+=mu_mass*R2
            coef2=min(R2,1+lam/t+F(21,8)*R2*eta/(t*t))
            mean2+=mu_mass*coef2
        if registered['reuse_old_range']:
            L=L2=F(old['whole_box_H_lower'])
            U=U2=F(old['whole_box_H_upper'])
        assert L==F(r['H_range_lower']) and U==F(r['H_range_upper'])
        assert mean==F(r['mean_H_upper']) and L<=mean<=U and L2<=mean2<=U2
        q=arb(5)/4
        chord2=ar(volume)*(ar(L2)**q+(ar(U2)**q-ar(L2)**q)*ar(mean2-L2)/ar(U2-L2))
        cu=ep(chord2,True)
        mean_cap,power_cap=map(F,caps[ci])
        assert max(mean,mean2)<=mean_cap
        assert max(cu,F(r['chord_integral_upper']))<=power_cap
        expected=ci==2
        assert r['meaningful_gate_passed']==expected and r['ideal_gate_passed']==expected
        # Both failed cases lie strictly above the gate even allowing outward radius.
        ratio=chord2/(ar(volume)*ar(Hlo)**q)
        assert (ep(ratio,True)<F(5,4)) if expected else (ep(ratio)>2)
        outputs.append({'h':str(h),'mean_upper_cap':str(mean_cap),'power_integral_upper_cap':str(power_cap),
                        'independent_mean_upper':str(mean2),'independent_power_upper':str(cu),
                        'first_Taylor_band_start':first_taylor,
                        'octave_mean_contribution_upper':list(map(str,octave_mean)),
                        'meaningful_and_ideal_gate_passed':expected})
    signal.setitimer(signal.ITIMER_REAL,0)
    out={'status':'passed','pid':os.getpid(),'source_sha256':sha(__file__),
         'elapsed_seconds':time.monotonic()-start,'new_H_stencil_jet_cones':0,
         'old067_manifest_intact':True,'cells':outputs}
    assert out['elapsed_seconds']<180
    (HERE/'readback.json').write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({'status':'passed','seconds':out['elapsed_seconds'],
                     'caps':caps,'old067_intact':True},indent=2))


if __name__=='__main__':
    main()
