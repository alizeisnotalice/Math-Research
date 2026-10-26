#!/usr/bin/env python3
"""Three registered band-mean/chord objects, zero H/stencil/jet/cone calls."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import argparse, hashlib, json, os, resource, signal, sys, time
import flint
from flint import arb, ctx, fmpq

HERE=Path(__file__).resolve().parent
if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atomic(p,value):
    t=p.with_suffix('.tmp')
    t.write_text(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n')
    os.replace(t,p)


def ar(x):
    x=F(x)
    return arb(fmpq(x.numerator,x.denominator))


def ep(x, upper=False):
    return str((x.upper() if upper else x.lower()).fmpq())


def main():
    started=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    reg=json.loads((HERE/'registration.json').read_text())
    assert reg['max_program_wall_seconds']==180
    receipt={'pid':os.getpid(),'python':sys.executable,'flint':flint.__version__,
             'source_sha256':sha(__file__),'registration_sha256':sha(HERE/'registration.json'),
             'new_H_stencil_jet_cone_evaluations':0,'old_ranges_recomputed':0}
    def stop(signum=None,frame=None):
        receipt.update({'status':'checkpointed_exit75','elapsed_seconds':time.monotonic()-started})
        atomic(HERE/'run_receipt.json',receipt)
        raise SystemExit(75)
    signal.signal(signal.SIGALRM,stop)
    signal.setitimer(signal.ITIMER_REAL,reg['hard_alarm_seconds'])
    assert sha(HERE/'mechanism.md')==reg['mechanism_sha256']
    for src in reg['sources'].values():
        assert sha(src['path'])==src['sha256']
    olddir=Path(reg['sources']['manifest.json']['path']).parent
    oldmanifest=json.loads((olddir/'manifest.json').read_text())
    for row in oldmanifest['files']:
        assert sha(olddir/row['name'])==row['sha256']
    old=json.loads(Path(reg['sources']['results.json']['path']).read_text())
    oldcp=json.loads(Path(reg['sources']['checkpoint.json']['path']).read_text())
    center=list(map(F,reg['center']))
    assert center==list(map(F,old['center']))
    ctx.prec=reg['precision_bits']
    assert ctx.prec==384
    fingerprint=hashlib.sha256((sha(__file__)+sha(HERE/'registration.json')).encode()).hexdigest()
    receipt['fingerprint']=fingerprint
    state={'fingerprint':fingerprint,'cells':[],'segments':[]}
    cppath=HERE/'checkpoint.json'
    if cppath.exists():
        assert args.resume,'Existing checkpoint requires --resume.'
        state=json.loads(cppath.read_text())
        assert state['fingerprint']==fingerprint
    Ninv=[[F(z,6) for z in row] for row in
          [[2,1,1,1,1],[1,-1,-1,-1,2],[1,-1,-1,2,-1],
           [1,-1,2,-1,-1],[1,2,-1,-1,-1]]]
    Q=[[F(i==j,4)-F(1,36) for j in range(4)] for i in range(4)]
    results=[]
    summaries=[]
    for ci,cell in enumerate(reg['cells']):
        h=F(cell['h'])
        lows=list(map(F,cell['box_lower']))
        highs=list(map(F,cell['box_upper']))
        assert lows==[max(F(0),x-h) for x in center[1:]]
        assert highs==[min(F(1),x+h) for x in center[1:]]
        ul=[l-c for l,c in zip(lows,center[1:])]
        uu=[u-c for u,c in zip(highs,center[1:])]
        mu=[(l+u)/2 for l,u in zip(ul,uu)]
        second=[(l*l+l*u+u*u)/3 for l,u in zip(ul,uu)]
        eta=sum(Q[i][j]*(second[i] if i==j else mu[i]*mu[j])
                for i in range(4) for j in range(4))
        lam=sum(abs(sum(Ninv[i][j+1]*mu[j] for j in range(4))) for i in range(5))
        corners=[sum(Q[i][j]*u[i]*u[j] for i in range(4) for j in range(4))
                 for u in product(*zip(ul,uu))]
        assert all(l<=0<=u for l,u in zip(ul,uu)) and max(corners)==h*h
        if not any(mu):
            assert lam==0 and eta==F(8,27)*h*h
        volume=prod(u-l for l,u in zip(lows,highs))
        if ci==len(state['cells']):
            state['cells'].append({'h':str(h),'bands':[]})
        data=state['cells'][ci]
        assert data['h']==str(h)
        for bi,b in enumerate(old['bands'][len(data['bands']):],len(data['bands'])):
            t=F(b['start'])
            if cell['reuse_old_range']:
                rl,ru=F(b['R_lower']),F(b['R_upper'])
                Rsource='reused067'
            else:
                R=( (1+ar(h*h/(20*t*t))).sqrt()+ar(h)/(2*arb(5).sqrt()*ar(t)) )**10
                rl,ru=F(ep(R)),F(ep(R,True))
                Rsource='new_two_widths_only'
            taylor=1+lam/t+F(21,8)*ru*eta/(t*t)
            coeff=min(ru,taylor)
            assert 1<=coeff<=ru
            data['bands'].append({'index':bi,'start':str(t),'R_lower':str(rl),'R_upper':str(ru),
                                  'R_source':Rsource,'Taylor_mean_coefficient':str(taylor),
                                  'chosen_mean_coefficient':str(coeff),
                                  'chosen':'Taylor_mean' if taylor<ru else 'range'})
            state['segments']=[{'pid':os.getpid(),'completed_cell':ci,
                               'completed_band_count':len(data['bands']),
                               'elapsed_seconds':time.monotonic()-started}]
            atomic(cppath,state)
            if time.monotonic()-started>=reg['soft_checkpoint_seconds']:
                stop()
        assert len(data['bands'])==89
        M=sum(F(b['mass_upper'])*F(r['chosen_mean_coefficient'])
              for b,r in zip(old['bands'],data['bands']))
        if cell['reuse_old_range']:
            L,U=F(old['whole_box_H_lower']),F(old['whole_box_H_upper'])
        else:
            L=sum(F(b['mass_lower'])/F(r['R_upper']) for b,r in zip(old['bands'],data['bands']))
            U=sum(F(b['mass_upper'])*F(r['R_upper']) for b,r in zip(old['bands'],data['bands']))
        assert 0<L<=M<=U
        q=arb(5)/4
        chord=ar(volume)*(ar(L)**q+(ar(U)**q-ar(L)**q)*ar(M-L)/ar(U-L))
        constant=ar(volume)*ar(U)**q
        Hlo=F(oldcp['rows'][0]['H']['lower'])
        baseline_lower=ar(volume)*ar(Hlo)**q
        ratio=chord/baseline_lower
        summary={'h':str(h),'clipped':any(uu[i]!=h or ul[i]!=-h for i in range(4)),
                 'volume':str(volume),'lambda':str(lam),'eta':str(eta),
                 'H_lower':float(L),'H_upper':float(U),'mean_upper':float(M),
                 'mean_over_Hcenter':float(M/Hlo),'chord_upper':float(F(ep(chord,True))),
                 'constant_upper':float(F(ep(constant,True))),
                 'chord_over_center_baseline_upper':float(F(ep(ratio,True))),
                 'constant_over_chord':float(F(ep(constant))/F(ep(chord,True))),
                 'Taylor_selected_bands':sum(r['chosen']=='Taylor_mean' for r in data['bands'])}
        result={'h':str(h),'box_lower':list(map(str,lows)),'box_upper':list(map(str,highs)),
                'u_lower':list(map(str,ul)),'u_upper':list(map(str,uu)),
                'mu':list(map(str,mu)),'second_moments':list(map(str,second)),
                'eta':str(eta),'lambda':str(lam),'volume':str(volume),
                'maximum_metric_corner':str(max(corners)),
                'H_range_lower':str(L),'H_range_upper':str(U),
                'mean_H_upper':str(M),'chord_integral_upper':ep(chord,True),
                'chord_integral_expression':str(chord),
                'constant_integral_lower':ep(constant),'constant_integral_upper':ep(constant,True),
                'chord_over_center_baseline_upper':ep(ratio,True),
                'meaningful_gate_passed':F(ep(ratio,True))<=F(reg['meaningful_power_ratio_gate']),
                'ideal_gate_passed':F(ep(ratio,True))<=F(reg['ideal_power_ratio_gate']),
                'old_range_reused':cell['reuse_old_range'],'bands':data['bands'],
                'summary':summary}
        results.append(result)
        summaries.append(summary)
        atomic(HERE/'results.json',{'fingerprint':fingerprint,'cells':results,'new_global_bound':False})
    receipt.update({'status':'complete','elapsed_seconds':time.monotonic()-started,
                    'new_range_widths':2,'new_mean_chord_cells':3,
                    'results_sha256':sha(HERE/'results.json'),'checkpoint_sha256':sha(cppath),
                    'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    assert receipt['elapsed_seconds']<reg['max_program_wall_seconds']
    atomic(HERE/'run_receipt.json',receipt)
    atomic(HERE/'summary.json',summaries)
    signal.setitimer(signal.ITIMER_REAL,0)
    print(json.dumps({'runtime':receipt,'cells':summaries},indent=2))


if __name__=='__main__':
    main()
