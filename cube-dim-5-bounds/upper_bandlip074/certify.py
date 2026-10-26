#!/usr/bin/env python3
"""Three new eps/width pairs, audited body-Lipschitz and positive band means."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import argparse,hashlib,json,os,signal,sys,time,resource
import flint
from flint import arb,ctx,fmpq

HERE=Path(__file__).resolve().parent
if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atomic(p,o):
    tmp=p.with_suffix('.tmp')
    tmp.write_text(json.dumps(o,sort_keys=True,separators=(',',':'))+'\n')
    os.replace(tmp,p)


def ar(x):
    x=F(x)
    return arb(fmpq(x.numerator,x.denominator))


def ep(x,upper=False):
    return str((x.upper() if upper else x.lower()).fmpq())


def main():
    began=time.monotonic()
    parser=argparse.ArgumentParser()
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    reg=json.loads((HERE/'registration.json').read_text())
    receipt={'pid':os.getpid(),'source_sha256':sha(__file__),
             'registration_sha256':sha(HERE/'registration.json'),'python':sys.executable,
             'flint':flint.__version__,'new_H_stencil_jet_cone_evaluations':0,
             '076_consumed':False}
    def stop(signum=None,frame=None):
        receipt.update({'status':'checkpointed_exit75','elapsed_seconds':time.monotonic()-began})
        atomic(HERE/'run_receipt.json',receipt)
        raise SystemExit(75)
    signal.signal(signal.SIGALRM,stop)
    signal.setitimer(signal.ITIMER_REAL,reg['hard_alarm_seconds'])
    assert sha(HERE/'mechanism.md')==reg['mechanism_sha256']
    for src in reg['sources'].values(): assert sha(src['path'])==src['sha256']
    for name in ['067_manifest','072_manifest','075_manifest']:
        p=Path(reg['sources'][name]['path'])
        for row in json.loads(p.read_text())['files']:
            assert sha(p.parent/row.get('name',row.get('path')))==row['sha256']
    old=json.loads(Path(reg['sources']['067_results']['path']).read_text())
    endpoints=json.loads(Path(reg['sources']['067_endpoints']['path']).read_text())['rows']
    first=reg['first_old_band_index']
    bands=old['bands'][first:]
    assert first==16 and len(bands)==73 and bands[0]['start']==reg['eps']=='1/512'
    assert bands[-1]['stop']=='infinity'
    assert all(b['stop']==bands[i+1]['start'] for i,b in enumerate(bands[:-1]))
    center_H=endpoints[first]['H']
    assert endpoints[first]['t']==reg['eps']
    assert sum(F(b['mass_lower']) for b in bands)<=F(center_H['lower'])
    assert sum(F(b['mass_upper']) for b in bands)>=F(center_H['upper'])
    center=list(map(F,reg['center']))
    assert center==list(map(F,old['center']))
    ctx.prec=reg['precision_bits']
    assert ctx.prec==384 and reg['max_program_wall_seconds']==180
    fingerprint=hashlib.sha256((sha(__file__)+sha(HERE/'registration.json')).encode()).hexdigest()
    receipt['fingerprint']=fingerprint
    state={'fingerprint':fingerprint,'cells':[]}
    cp=HERE/'checkpoint.json'
    if cp.exists():
        assert args.resume,'Existing checkpoint requires --resume.'
        state=json.loads(cp.read_text())
        assert state['fingerprint']==fingerprint
    Q=[[F(i==j,4)-F(1,36) for j in range(4)] for i in range(4)]
    Ntail=[[1,1,1,1],[-1,-1,-1,2],[-1,-1,2,-1],[-1,2,-1,-1],[2,-1,-1,-1]]
    cells=[]
    summaries=[]
    for ci,cell in enumerate(reg['cells']):
        h=F(cell['h'])
        bl=list(map(F,cell['box_lower']))
        bu=list(map(F,cell['box_upper']))
        assert bl==[max(F(0),c-h) for c in center[1:]]
        assert bu==[min(F(1),c+h) for c in center[1:]]
        ul=[x-c for x,c in zip(bl,center[1:])]
        uu=[x-c for x,c in zip(bu,center[1:])]
        mu=[(l+u)/2 for l,u in zip(ul,uu)]
        second=[(l*l+l*u+u*u)/3 for l,u in zip(ul,uu)]
        eta=sum(Q[i][j]*(second[i] if i==j else mu[i]*mu[j]) for i in range(4) for j in range(4))
        lambda_old=sum(abs(sum(F(a,6)*m for a,m in zip(row,mu))) for row in Ntail)
        lambda_lip=max(map(abs,mu))
        lam=min(lambda_old,lambda_lip)
        heff=max(abs(x) for pair in zip(ul,uu) for x in pair)
        corners=[sum(Q[i][j]*u[i]*u[j] for i in range(4) for j in range(4)) for u in product(*zip(ul,uu))]
        assert max(corners)==h*h and heff==h and all(l<=0<=u for l,u in zip(ul,uu))
        volume=prod(u-l for l,u in zip(bl,bu))
        if ci==len(state['cells']): state['cells'].append({'h':str(h),'bands':[]})
        done=state['cells'][ci]
        assert done['h']==str(h)
        for bi,b in enumerate(bands[len(done['bands']):],len(done['bands'])):
            t=F(b['start'])
            oldR=((1+ar(h*h/(20*t*t))).sqrt()+ar(h)/(2*arb(5).sqrt()*ar(t)))**10
            lipR=ar(heff/t).exp()
            oldlo,oldhi=F(ep(oldR)),F(ep(oldR,True))
            liplo,liphi=F(ep(lipR)),F(ep(lipR,True))
            rl,ru=min(oldlo,liplo),min(oldhi,liphi)
            taylor=1+lam/t+F(21,8)*ru*eta/(t*t)
            coef=min(ru,taylor)
            assert 1<=rl<=ru and 1<=coef<=ru
            done['bands'].append({'old_band_index':first+bi,'start':str(t),'stop':b['stop'],
                                  'R_old_lower':str(oldlo),'R_old_upper':str(oldhi),
                                  'R_lip_lower':str(liplo),'R_lip_upper':str(liphi),
                                  'R_lower':str(rl),'R_upper':str(ru),
                                  'ratio_choice':'lip' if liphi<oldhi else 'old',
                                  'Taylor_mean_coefficient':str(taylor),
                                  'chosen_mean_coefficient':str(coef),
                                  'mean_choice':'Taylor' if taylor<ru else 'range'})
            state['progress']={'pid':os.getpid(),'cell':ci,'bands_complete':len(done['bands']),
                               'elapsed_seconds':time.monotonic()-began}
            atomic(cp,state)
            if time.monotonic()-began>=reg['soft_checkpoint_seconds']: stop()
        assert len(done['bands'])==73
        L=sum(F(b['mass_lower'])/F(r['R_upper']) for b,r in zip(bands,done['bands']))
        U=sum(F(b['mass_upper'])*F(r['R_upper']) for b,r in zip(bands,done['bands']))
        M=sum(F(b['mass_upper'])*F(r['chosen_mean_coefficient']) for b,r in zip(bands,done['bands']))
        assert 0<L<=M<=U
        q=arb(5)/4
        chord=ar(volume)*(ar(L)**q+(ar(U)**q-ar(L)**q)*ar(M-L)/ar(U-L))
        constant=ar(volume)*ar(U)**q
        baseline=ar(volume)*ar(F(center_H['lower']))**q
        ratio=chord/baseline
        ratiohi=F(ep(ratio,True))
        first_taylor=next((r['start'] for r in done['bands'] if r['mean_choice']=='Taylor'),None)
        summary={'h':str(h),'eps':reg['eps'],'clipped_coordinates':sum(u!=h or l!=-h for l,u in zip(ul,uu)),
                 'volume':str(volume),'lambda':str(lam),'eta':str(eta),
                 'H_center':float(F(center_H['lower'])),'L':float(L),'U':float(U),'M':float(M),
                 'mean_over_center':float(M/F(center_H['lower'])),
                 'chord_upper':float(F(ep(chord,True))),'constant_upper':float(F(ep(constant,True))),
                 'normalized_power_upper':float(ratiohi),
                 'passes_1_1':ratiohi<=F(11,10),'passes_1_05':ratiohi<=F(21,20),
                 'lip_selected_bands':sum(r['ratio_choice']=='lip' for r in done['bands']),
                 'Taylor_selected_bands':sum(r['mean_choice']=='Taylor' for r in done['bands']),
                 'first_Taylor_band':first_taylor}
        cells.append({'h':str(h),'eps':reg['eps'],'box_lower':list(map(str,bl)),'box_upper':list(map(str,bu)),
                      'u_lower':list(map(str,ul)),'u_upper':list(map(str,uu)),
                      'mu':list(map(str,mu)),'second_moments':list(map(str,second)),
                      'eta':str(eta),'lambda_old':str(lambda_old),'lambda_lip':str(lambda_lip),
                      'lambda_effective':str(lam),'h_effective':str(heff),'volume':str(volume),
                      'H_range_lower':str(L),'H_range_upper':str(U),'mean_H_upper':str(M),
                      'chord_integral_upper':ep(chord,True),'constant_integral_upper':ep(constant,True),
                      'normalized_power_upper':ep(ratio,True),'gates':summary['passes_1_1'],
                      'bands':done['bands'],'summary':summary})
        summaries.append(summary)
        atomic(HERE/'results.json',{'fingerprint':fingerprint,'center_H_eps':center_H,
                                   'time_domain':'[1/512,infinity)','center':reg['center'],
                                   'first_old_band':16,'positive_bands':73,'cells':cells,
                                   '076_consumed':False,'new_global_I_g_C5':False})
    receipt.update({'status':'complete','elapsed_seconds':time.monotonic()-began,
                    'results_sha256':sha(HERE/'results.json'),'checkpoint_sha256':sha(cp),
                    'new_eps_width_pairs':3,'positive_bands_per_pair':73,
                    'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    assert receipt['elapsed_seconds']<180
    atomic(HERE/'run_receipt.json',receipt)
    atomic(HERE/'summary.json',summaries)
    signal.setitimer(signal.ITIMER_REAL,0)
    print(json.dumps({'runtime':receipt,'cells':summaries},indent=2))


if __name__=='__main__': main()
