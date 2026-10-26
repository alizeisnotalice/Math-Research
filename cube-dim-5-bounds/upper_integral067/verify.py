#!/usr/bin/env python3
"""Independent readback: no H, log, DD, jet, cone, or time quadrature rerun."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib, json, os, sys, time
from flint import arb, ctx, fmpq

HERE=Path(__file__).resolve().parent
if hasattr(sys,'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ar(x):
    x=F(x)
    return arb(fmpq(x.numerator,x.denominator))


def endpoint(x, upper=False):
    return F(str((x.upper() if upper else x.lower()).fmpq()))


def main():
    start=time.monotonic()
    ctx.prec=384
    reg=json.loads((HERE/'registration.json').read_text())
    result=json.loads((HERE/'results.json').read_text())
    cp=json.loads((HERE/'checkpoint.json').read_text())
    receipt=json.loads((HERE/'run_receipt.json').read_text())
    assert sha(HERE/'mechanism.md')==reg['mechanism_sha256']
    for s in reg['sources'].values():
        assert sha(s['path'])==s['sha256']
    assert sha(HERE/'certify.py')==receipt['source_sha256']
    assert sha(HERE/'registration.json')==receipt['registration_sha256']
    assert sha(HERE/'results.json')==receipt['results_sha256']
    assert sha(HERE/'checkpoint.json')==receipt['endpoint_checkpoint_sha256']
    assert sha(HERE/'zero_stencil.json')==receipt['zero_stencil_sha256']
    assert result['fingerprint']==cp['fingerprint']==receipt['fingerprint']
    low=json.loads(Path(reg['sources']['low063']['path']).read_text())
    assert cp['rows'][0]['H']==low['center_H']
    h=F(reg['halfwidth'])
    c=list(map(F,result['center']))
    assert c==list(map(F,low['selected_center']['point']))
    assert list(map(F,result['box_lower']))==[x-h for x in c[1:]]
    assert list(map(F,result['box_upper']))==[x+h for x in c[1:]]
    assert all(0<x-h<x+h<1 for x in c[1:])
    assert F(result['volume'])==(2*h)**4==F(1,4294967296)
    corners=[sum(F(x*x,4) for x in s)-F(sum(s)**2,36)
             for s in product([-1,1],repeat=4)]
    assert corners==list(map(F,result['metric_corner_values'])) and max(corners)==1
    rows=cp['rows']
    times=list(map(F,reg['time_endpoints']))
    assert [F(r['t']) for r in rows]==times
    assert times==[F(2**i,2048)*(1+F(j,8)) for i in range(11) for j in range(8)]+[F(1)]
    weights=json.loads((HERE/'zero_stencil.json').read_text())['weights']
    assert len(weights)==140 and sum(map(F,weights.values()))==0
    L=U=L2=U2=F(0)
    ml_sum=mu_sum=F(0)
    groups=[]
    for i,(row,b) in enumerate(zip(rows,result['bands'])):
        lo,hi=F(row['H']['lower']),F(row['H']['upper'])
        assert 0<lo<=hi
        if i+1<len(rows):
            ml=lo-F(rows[i+1]['H']['upper'])
            mu=hi-F(rows[i+1]['H']['lower'])
            assert b['stop']==rows[i+1]['t']
        else:
            ml,mu=lo,hi
            assert b['stop']=='infinity'
        assert 0<ml<=mu
        assert ml==F(b['mass_lower']) and mu==F(b['mass_upper'])
        assert b['start']==row['t']
        R=F(b['R_upper'])
        assert 1<=F(b['R_lower'])<=R
        assert ml/R==F(b['contribution_lower']) and mu*R==F(b['contribution_upper'])
        L+=ml/R
        U+=mu*R
        ml_sum+=ml
        mu_sum+=mu
        # Independent algebraic formula:
        # cosh(10 asinh(s)) = T5(1+2s^2), R = P + sqrt(P^2-1).
        x=1+h*h/(10*times[i]*times[i])
        P=16*x**5-20*x**3+5*x
        assert P>1
        alternate=ar(P)+ar(P*P-1).sqrt()
        R2=endpoint(alternate,True)
        if i==0:
            assert 2*P*3125==3125**2+1
            R2=F(3125)
        L2+=ml/R2
        U2+=mu*R2
        if i%8==0:
            groups.append({'start':row['t'],'mass_upper':F(0),'contribution_upper':F(0)})
        groups[-1]['stop']=b['stop']
        groups[-1]['mass_upper']+=mu
        groups[-1]['contribution_upper']+=mu*R2
    assert L==F(result['whole_box_H_lower']) and U==F(result['whole_box_H_upper'])
    assert ml_sum<=F(rows[0]['H']['lower'])<=F(rows[0]['H']['upper'])<=mu_sum
    q=arb(5)/4
    vol=F(result['volume'])
    powlo=endpoint(ar(vol)*ar(L2)**q)
    powhi=endpoint(ar(vol)*ar(U2)**q,True)
    # Short rational certificate, verified independently by both evaluations.
    caps={'H_lower':'1136663/10000000000','H_upper':'167327/1250000000',
          'power_integral_lower':'136631/50000000000000000000',
          'power_integral_upper':'83811/25000000000000000000'}
    assert F(caps['H_lower'])<=min(L,L2)
    assert max(U,U2)<=F(caps['H_upper'])
    assert F(caps['power_integral_lower'])<=min(powlo,F(result['power_integral_lower']))
    assert max(powhi,F(result['power_integral_upper']))<=F(caps['power_integral_upper'])
    assert max(U,U2)<F(5,4)*F(rows[0]['H']['lower'])
    assert result['ideal_gate_passed'] and result['meaningful_gate_passed']
    assert not result['new_global_I_or_g_or_C5']
    out={'status':'passed','pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,
         'source_sha256':sha(__file__),'precision_bits':384,
         'independent_R_formula':'P+sqrt(P^2-1), P=T5(1+h^2/(10t^2))',
         'new_H_or_logs_or_jet_or_cones':0,'simple_rational_certificate':caps,
         'independent_H_lower':str(L2),'independent_H_upper':str(U2),
         'independent_power_lower':str(powlo),'independent_power_upper':str(powhi),
         'groups':[{k:str(v) for k,v in g.items()} for g in groups]}
    (HERE/'readback.json').write_text(json.dumps(out,sort_keys=True,separators=(',',':'))+'\n')
    print(json.dumps({'status':out['status'],'seconds':out['elapsed_seconds'],
                     'caps':caps,'time_bands':len(rows),'new_H':0},indent=2))


if __name__=='__main__':
    main()
