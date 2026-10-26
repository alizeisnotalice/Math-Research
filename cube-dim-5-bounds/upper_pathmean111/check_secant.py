#!/usr/bin/env python3
"""Six preregistered same-catalog cases, with independent linear-weight check."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,time
import pathmean as P
import secant as S
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pack(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return list(map(pack,x))
    return x
start=time.monotonic();print('secant111_pid',os.getpid(),flush=True)
registration=json.loads((HERE/'secant_registration.json').read_text())
for p,h in registration['sources'].items():assert sha(HERE/p)==h
reg=json.loads((HERE/'registration.json').read_text());cat=json.loads((HERE/'ratio_catalog.json').read_text());rows=[]
for case in reg['finite_cases']:
    plan=P.darboux_plan(case['offsets'],case['a'],case['panels_per_piece']);got=S.secant_envelope(plan,cat)
    assert got['J_upper']<=got['old_J_upper']
    independent=F(0)
    for cell,weight in zip(plan['cells'],got['positive_endpoint_weights']):
        l,r=cell['left'],cell['right'];W=next(p['W'] for p in plan['profile']['pieces'] if p['lo']<=l<r<=p['hi'])
        # Multiply W by the two barycentric linear polynomials, rather than
        # using the first-moment identity in the implementation.
        alpha=P.integral(P.mul(W,[r/(r-l),-1/(r-l)]),l,r)
        beta=P.integral(P.mul(W,[-l/(r-l),1/(r-l)]),l,r)
        assert alpha==weight['alpha'] and beta==weight['beta'] and alpha>=0 and beta>=0
        Rr=F(cat[str(r/plan['a'])]['upper']);Rl=F(cat[str(l/plan['a'])]['upper'])
        independent+=alpha*min(Rl,Rr)+beta*Rr
    if plan['profile']['h']:assert independent==got['secant_J_upper']
    rows.append({'name':case['name'],'result':got,
                 'strict_mean_factor_improvement':got['factor_upper']<got['old085_factor_upper']})
out={'status':'passed','pid':os.getpid(),'wall_seconds':time.monotonic()-start,
     'registration_sha256':sha(HERE/'secant_registration.json'),'secant_sha256':sha(HERE/'secant.py'),
     'code_sha256':sha(__file__),'cases':rows,'new_ratio_H_jet_log_calls':0,
     'strict_mean_factor_improvements':sum(row['strict_mean_factor_improvement'] for row in rows)}
(HERE/'secant_results.json').write_text(json.dumps(pack(out),sort_keys=True,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ('status','pid','wall_seconds','strict_mean_factor_improvements','new_ratio_H_jet_log_calls')}))
for row in rows:
    g=row['result'];print(row['name'],'factor',float(g['factor_upper']),'old',float(g['old085_factor_upper']),'curvature_fraction',float(g['curvature_ratio_upper']))
