#!/usr/bin/env python3
"""Standard-library exact replay of explicit rational posterior-space gains.
Uses same-directory replay_hull.py and r732/replay_exact.py; no optimization.
"""
import hashlib
import importlib.util
import json
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
helper=HERE/'replay_hull.py'
spec=importlib.util.spec_from_file_location('hull_replay',helper)
hull=importlib.util.module_from_spec(spec)
spec.loader.exec_module(hull)
exact=hull.exact


def check(base,r,scale):
    prior=hull.check(base,r)
    mass=list(map(Q,base['atom_masses']))
    vol=list(map(Q,base['row_volumes']))
    masks=base['atom_masks']
    fullpi=list(map(Q,base['pi']))
    active=base['active_rows']
    pi=[fullpi[i] for i in active]
    q=[m/sum(mass) for m in mass]
    ker=[[m/vol[i] if mask&(1<<i) else Q(0)
          for m,mask in zip(mass,masks)] for i in active]
    pm=[sum(pi[i]*ker[i][j] for i in range(len(pi))) for j in range(len(mass))]
    alpha=[[pi[i]*ker[i][j]/pm[j] for i in range(len(pi))] for j in range(len(mass))]
    rounded=list(map(Q,r['rounded_B']))
    phi=list(map(Q,r['phi']))
    step=Q(r['step_t'])
    z=list(map(Q,r['space_weights']))
    nu=list(map(Q,r['induced_nu']))
    assert len(rounded)==len(phi)==len(pm)
    assert all((b*scale).denominator==1 for b in rounded)
    center=sum(p*b for p,b in zip(pm,rounded))
    assert phi==[b-center for b in rounded]
    assert sum(p*x for p,x in zip(pm,phi))==0
    maximum=max(map(abs,phi))
    assert step==(min(Q(1,2),1/(2*maximum)) if maximum else Q(1,2))
    assert 0<step<=Q(1,2) and all(abs(step*x)<=Q(1,2) for x in phi)
    assert z==[p*(1+step*x) for p,x in zip(pm,phi)]
    assert sum(z)==1 and min(z)>0
    assert nu==[sum(z[j]*alpha[j][i] for j in range(len(pm))) for i in range(len(pi))]
    w=[x/y for x,y in zip(pm,q)]
    grad0=[exact.lin([(ker[i][j],w[j]) for j in range(len(pm))]) for i in range(len(pi))]
    B=[(sum(alpha[j][i]*grad0[i][0] for i in range(len(pi))),
        sum(alpha[j][i]*grad0[i][1] for i in range(len(pi)))) for j in range(len(pm))]
    half=Q(1,2*scale)
    assert all(b-half<=interval[0]<=interval[1]<=b+half for b,interval in zip(rounded,B))
    coeff=[sum((nu[i]-pi[i])*ker[i][j] for i in range(len(pi))) for j in range(len(pm))]
    linear=exact.lin(list(zip(coeff,w)))
    penalty=step*step*sum(p*x*x for p,x in zip(pm,phi))
    analytic=linear[0]-penalty
    assert Q(r['analytic_direction_gain_lower'])<=analytic
    pn=[sum(nu[i]*ker[i][j] for i in range(len(pi))) for j in range(len(pm))]
    fl,fh=exact.lin([(p,p/qq) for p,qq in zip(pn,q)]+[(-v,v/p) for v,p in zip(nu,pi)])
    kl=exact.lin([(p,p/qq) for p,qq in zip(pm,q)])
    direct=fl-Q(base['KL_pi_upper'])
    assert Q(r['direct_gain_lower'])<=direct
    if r['certified_strict_gain']:
        assert Q(r['direct_gain_lower'])>0 and direct>0
    return dict(case=r['case'],dimension=base['dimension'],target=r['target'],
                prior_hull_replay_pass=True,exact_centering=True,
                exact_step_and_positivity=True,exact_space_tilt=True,
                exact_induced_nu=True,rounded_B_error_window_verified=True,
                direct_gain_lower_verified=True,analytic_gain_lower_verified=True,
                direct_gain_lower_64=hull.decimal(direct),
                direct_gain_upper_64=hull.decimal(fh-kl[0],True),
                analytic_gain_lower_64=hull.decimal(analytic),
                analytic_penalty_exact=str(penalty),
                strict_direct_gain=direct>0,strict_analytic_gain=analytic>0,
                stored_strict_gain_verified=bool(r['certified_strict_gain']))


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled; run without -O/PYTHONOPTIMIZE.')
    start=time.monotonic()
    raw=(HERE.parent/'r732'/'capacity.json').read_bytes()
    data=json.loads((HERE/'tilt_gain.json').read_text())
    digest=hashlib.sha256(raw).hexdigest()
    assert digest==data['source_sha256']
    base=json.loads(raw)['results']
    byid={r['case']:r for r in base}
    rows=data['results']
    assert len(byid)==len(base)==len(rows) and {r['case'] for r in rows}==set(byid)
    result=[check(byid[r['case']],r,data['rounding_scale']) for r in rows]
    out=dict(status='complete',all_pass=True,cases=len(result),
             source_sha256=digest,log_series_terms=64,
             strict_direct_gains=sum(r['strict_direct_gain'] for r in result),
             strict_analytic_gains=sum(r['strict_analytic_gain'] for r in result),
             maximum_direct_gain_lower=max((Q(r['direct_gain_lower_64']) for r in result)),
             python=sys.version.split()[0],elapsed_seconds=time.monotonic()-start,
             scope='Exact rational rounded-direction gains at fixed LF pi; posterior hull certificate replay included; no optimizer, no LF-face optimization, no central asymptotic claim.',
             results=result)
    out['maximum_direct_gain_lower']=hull.decimal(out['maximum_direct_gain_lower'])
    (HERE/'replay_gain.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','all_pass','cases','strict_direct_gains','strict_analytic_gains','maximum_direct_gain_lower','elapsed_seconds']}))


if __name__=='__main__':
    main()
