#!/usr/bin/env python3
"""Independent retained-time, box moment, min-ratio, mean and chord readback."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import hashlib,json,os,signal,sys,time
from flint import arb,ctx,fmpq

HERE=Path(__file__).resolve().parent
if hasattr(sys,'set_int_max_str_digits'): sys.set_int_max_str_digits(0)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ar(x):
    x=F(x)
    return arb(fmpq(x.numerator,x.denominator))


def ep(x,upper=False): return F(str((x.upper() if upper else x.lower()).fmpq()))


def main():
    began=time.monotonic()
    def stop(signum,frame): raise SystemExit(75)
    signal.signal(signal.SIGALRM,stop)
    signal.setitimer(signal.ITIMER_REAL,175)
    ctx.prec=384
    reg=json.loads((HERE/'registration.json').read_text())
    r=json.loads((HERE/'results.json').read_text())
    cp=json.loads((HERE/'checkpoint.json').read_text())
    receipt=json.loads((HERE/'run_receipt.json').read_text())
    assert sha(HERE/'mechanism.md')==reg['mechanism_sha256']
    for src in reg['sources'].values(): assert sha(src['path'])==src['sha256']
    for key in ['067_manifest','072_manifest','075_manifest']:
        p=Path(reg['sources'][key]['path'])
        for row in json.loads(p.read_text())['files']:
            assert sha(p.parent/row.get('name',row.get('path')))==row['sha256']
    assert sha(HERE/'certify.py')==receipt['source_sha256']
    assert sha(HERE/'registration.json')==receipt['registration_sha256']
    assert sha(HERE/'results.json')==receipt['results_sha256']
    assert sha(HERE/'checkpoint.json')==receipt['checkpoint_sha256']
    assert r['fingerprint']==cp['fingerprint']==receipt['fingerprint']
    assert receipt['new_H_stencil_jet_cone_evaluations']==0 and not receipt['076_consumed']
    old=json.loads(Path(reg['sources']['067_results']['path']).read_text())
    oldpoints=json.loads(Path(reg['sources']['067_endpoints']['path']).read_text())['rows']
    bands=old['bands'][16:]
    assert len(bands)==73 and r['positive_bands']==73
    assert bands[0]['start']=='1/512' and bands[-1]['stop']=='infinity'
    assert all(bands[i]['stop']==bands[i+1]['start'] for i in range(72))
    assert r['center_H_eps']==oldpoints[16]['H'] and oldpoints[16]['t']=='1/512'
    Hlo=F(r['center_H_eps']['lower'])
    Hhi=F(r['center_H_eps']['upper'])
    assert 0<sum(F(b['mass_lower']) for b in bands)<=Hlo<=Hhi<=sum(F(b['mass_upper']) for b in bands)
    c=list(map(F,reg['center']))
    Q=[[F(i==j,4)-F(1,36) for j in range(4)] for i in range(4)]
    ntail=[[1,1,1,1],[-1,-1,-1,2],[-1,-1,2,-1],[-1,2,-1,-1],[2,-1,-1,-1]]
    caps=[('1966477/10000000000','141449/10000000000000000','1821499/1000000'),
          ('1308479/10000000000','731983/1000000000000000000','1090249/1000000'),
          ('1237487/10000000000','486579/10000000000000000000','1014877/1000000')]
    outputs=[]
    for ci,(cell,spec) in enumerate(zip(r['cells'],reg['cells'])):
        h=F(cell['h'])
        assert cell['h']==spec['h'] and cell['eps']=='1/512'
        ul=[max(F(0),x-h)-x for x in c[1:]]
        uu=[min(F(1),x+h)-x for x in c[1:]]
        assert ul==list(map(F,cell['u_lower'])) and uu==list(map(F,cell['u_upper']))
        mu=[(l+u)/2 for l,u in zip(ul,uu)]
        var=[(u-l)**2/12 for l,u in zip(ul,uu)]
        eta=sum(Q[i][i]*var[i] for i in range(4))+sum(Q[i][j]*mu[i]*mu[j] for i in range(4) for j in range(4))
        oldlambda=sum(abs(sum(F(z,6)*m for z,m in zip(row,mu))) for row in ntail)
        liplambda=max(map(abs,mu))
        lam=min(oldlambda,liplambda)
        heff=max(abs(x) for pair in zip(ul,uu) for x in pair)
        assert eta==F(cell['eta']) and mu==list(map(F,cell['mu']))
        assert oldlambda==F(cell['lambda_old']) and liplambda==F(cell['lambda_lip'])
        assert lam==F(cell['lambda_effective']) and heff==h==F(cell['h_effective'])
        assert all(var[i]+mu[i]**2==F(cell['second_moments'][i]) for i in range(4))
        vmax=max(sum(Q[i][j]*u[i]*u[j] for i in range(4) for j in range(4)) for u in product(*zip(ul,uu)))
        assert vmax==h*h
        volume=prod(u-l for l,u in zip(ul,uu))
        assert volume==F(cell['volume'])
        L=U=M=L2=U2=M2=F(0)
        octave=[F(0) for _ in range(10)]
        for bi,(v,b) in enumerate(zip(cell['bands'],bands)):
            assert v==cp['cells'][ci]['bands'][bi] and v['old_band_index']==16+bi
            assert v['start']==b['start'] and v['stop']==b['stop']
            t=F(v['start'])
            R=min(F(v['R_old_upper']),F(v['R_lip_upper']))
            assert R==F(v['R_upper'])
            assert min(F(v['R_old_lower']),F(v['R_lip_lower']))==F(v['R_lower'])
            # This finite batch selects lip strictly; no all-parameter assertion.
            assert F(v['R_lip_upper'])<F(v['R_old_lower']) and v['ratio_choice']=='lip'
            coeff=min(R,1+lam/t+F(21,8)*R*eta/(t*t))
            assert coeff==F(v['chosen_mean_coefficient'])
            assert 1+lam/t+F(21,8)*R*eta/(t*t)==F(v['Taylor_mean_coefficient'])
            assert v['mean_choice']==('Taylor' if coeff<R else 'range')
            ml,mu_mass=F(b['mass_lower']),F(b['mass_upper'])
            assert 0<ml<=mu_mass
            L+=ml/R;U+=mu_mass*R;M+=mu_mass*coeff
            # Independent old expression, and split-exponential for new expression.
            x=1+h*h/(10*t*t)
            P=16*x**5-20*x**3+5*x
            Rold2=ep(ar(P)+ar(P*P-1).sqrt(),True)
            Rlip2=ep(ar(heff/(2*t)).exp()**2,True)
            R2=min(Rold2,Rlip2)
            C2=min(R2,1+lam/t+F(21,8)*R2*eta/(t*t))
            L2+=ml/R2;U2+=mu_mass*R2;M2+=mu_mass*C2
            octave[bi//8]+=mu_mass*C2
        assert L==F(cell['H_range_lower']) and U==F(cell['H_range_upper']) and M==F(cell['mean_H_upper'])
        assert 0<L<=M<=U and 0<L2<=M2<=U2
        q=arb(5)/4
        chord=ar(volume)*(ar(L2)**q+(ar(U2)**q-ar(L2)**q)*ar(M2-L2)/ar(U2-L2))
        ratio=chord/(ar(volume)*ar(Hlo)**q)
        mean_cap,power_cap,ratio_cap=map(F,caps[ci])
        assert max(M,M2)<=mean_cap
        assert max(F(cell['chord_integral_upper']),ep(chord,True))<=power_cap
        assert max(F(cell['normalized_power_upper']),ep(ratio,True))<=ratio_cap
        pass11=ci>0;pass105=ci==2
        assert cell['summary']['passes_1_1']==pass11 and cell['summary']['passes_1_05']==pass105
        assert (ep(ratio,True)<F(11,10)) if pass11 else (ep(ratio)>F(11,10))
        assert (ep(ratio,True)<F(21,20)) if pass105 else (ep(ratio)>F(21,20))
        outputs.append({'h':str(h),'mean_cap':str(mean_cap),'power_cap':str(power_cap),
                        'normalized_cap':str(ratio_cap),'independent_mean_upper':str(M2),
                        'independent_power_upper':str(ep(chord,True)),
                        'independent_normalized_upper':str(ep(ratio,True)),
                        'octave_mean_contribution_upper':list(map(str,octave)),
                        'pass_1_1':pass11,'pass_1_05':pass105})
    signal.setitimer(signal.ITIMER_REAL,0)
    out={'status':'passed','pid':os.getpid(),'source_sha256':sha(__file__),
         'elapsed_seconds':time.monotonic()-began,'precision_bits':384,
         'old_067_072_075_manifests_intact':True,'new_H_stencil_jet_cones':0,
         '076_consumed':False,'cells':outputs}
    assert out['elapsed_seconds']<180
    (HERE/'readback.json').write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({'status':'passed','seconds':out['elapsed_seconds'],'caps':caps},indent=2))


if __name__=='__main__': main()
