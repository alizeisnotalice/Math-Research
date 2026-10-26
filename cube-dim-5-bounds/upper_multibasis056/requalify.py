#!/usr/bin/env python3
"""One frozen 162-box Fraction-only requalification; no kernel imports."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import os
import resource
import signal
import sys
import time
import sympy as sp

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--output-dir',type=Path,default=HERE)
parser.add_argument('--resume',action='store_true')
args=parser.parse_args();OUT=args.output_dir.resolve();OUT.mkdir(parents=True,exist_ok=True)
start=time.monotonic();rows=[];segments=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def digest(x):return hashlib.sha256(canon(x)).hexdigest()
reg=json.loads((HERE/'registration.json').read_text())
code_hash=sha(Path(__file__))
assert code_hash==reg['audit_code_sha256']
fp=digest({'code':code_hash,'registration':sha(HERE/'registration.json'),'state':reg['input_state_sha256']})
def save(stage,include_rows=True):
    x={'pid':os.getpid(),'stage':stage,'elapsed_seconds':time.monotonic()-start,'deadline_seconds':120,
       'fingerprint':fp,'completed_boxes':len(rows),'segments':segments,
       'recoverable_command':f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {__file__} --output-dir {OUT} --resume'}
    if include_rows:x['rows']=rows
    (OUT/'checkpoint.json').write_text(json.dumps(x,indent=2)+'\n')
def timeout(_sig,_frame):
    save('timeout_recoverable')
    raise SystemExit(75)
signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
if args.resume:
    old=json.loads((OUT/'checkpoint.json').read_text());assert old['fingerprint']==fp
    assert old['stage']!='complete'
    rows=old['rows'];segments=old.get('segments',[])
save('hashes_and_input')
statepath=Path(reg['input_state']);regpath=Path(reg['input_registration'])
assert sha(statepath)==reg['input_state_sha256']
assert sha(regpath)==reg['input_registration_sha256']
assert sha(Path(reg['source046']))==reg['source046_sha256']
state=json.loads(statepath.read_text());oldreg=json.loads(regpath.read_text())
assert state['fingerprint']==reg['input_fingerprint']
assert state['source052_hashes']==oldreg['sources']==reg['source052_hashes']
assert state['source_hashes']['arb_cover.py']==reg['source046_sha256']
assert state['epsilon']==oldreg['epsilon']==reg['model']['epsilon']=='1/64'
assert state['target']==oldreg['target']==reg['model']['target']=='1/10000'
derived=digest({'registration_sha256':reg['input_registration_sha256'],'sources':oldreg['sources'],
    'parent':oldreg['parent_state_sha256'],'selection':oldreg['selected_sha256'],'locks':oldreg['locks']})
assert derived==state['fingerprint']
assert digest(state['roots'])==state['initial_root_sha256']
assert len(state['splits'])==162 and oldreg['parent_splits']==0
assert state['batch_number']==2 and not state['inflight'] and state['qualification']['passed']
cases=[]
for split in state['splits']:
    node=split['node'];assert node[1]==''
    cert=split['failed_certificate'];assert cert['passed'] is False
    root=state['roots'][node[0]];assert len(root)==4 and all(0<=x<32 for x in root)
    record=state['cache'][cert['anchor']]
    cases.append({'node':node,'root':root,'old_certificate':cert,
                  'point':record['point'],'values':record['values'],'bits':record['bits']})
assert len({tuple(c['node']) for c in cases})==162
del state

save('signed_basis_algebra')
N=sp.Matrix.hstack(*(sp.Matrix(v) for v in [
    [1,1,1,1,1],[1,-1,-1,-1,1],[1,-1,-1,1,-1],[1,-1,1,-1,-1],[1,1,-1,-1,-1]]))
assert N.det()==-48
for sign in product((-1,1),repeat=4):
    Ds=sp.diag(1,*sign);Ns=Ds*N;Ni=Ns.inv();ss=sp.Matrix(sign)
    Qs=(Ni.T*Ni)[1:5,1:5]
    assert Qs==sp.eye(4)/4-ss*ss.T/36
    assert Qs.eigenvals()=={sp.Rational(1,4):3,sp.Rational(5,36):1}
    assert all(Ns[0,j]==1 and all(abs(Ns[i,j])==1 for i in range(5)) for j in range(5))
def baseq(v):return sum(x*x for x in v)/F(4)-sum(v)**2/F(36)
def psi(v):return sum(x*x for x in v)/F(4)-sum(abs(x) for x in v)**2/F(36)
global_concavity_gap=-psi((F(0),F(5),F(0),F(0)))-(-psi((F(1),F(5),F(0),F(0))))
assert global_concavity_gap==F(-1,18)
algebra={'signed_bases_checked':16,'Q_s':'I4/4-ss^T/36','eigenvalues':{'1/4':3,'5/36':1},
         'global_concavity_counterexample_midpoint_gap':str(global_concavity_gap),
         'global_curvature_bounds':'(5/36)||v||2^2 <= psi(v) <= (2/9)||v||2^2'}

eps=F(1,64);target=F(1,10000)
for index in range(len(rows),len(cases)):
    case=cases[index];cert=case['old_certificate']
    point=list(map(F,case['point']));assert point[0]==1
    assert ','.join(map(str,point))==cert['anchor']
    lo=[F(x,32) for x in case['root']];hi=[x+F(1,32) for x in lo]
    assert all(a<=p<=b for a,p,b in zip(lo,point[1:],hi))
    values=[tuple(map(F,pair)) for pair in case['values']]
    assert len(values)==5 and all(a<=b for a,b in values)
    Hlo,Hhi=values[0];assert 0<target<Hlo<=Hhi
    gradient=[tuple(map(F,pair)) for pair in cert['loggradient_directed']]
    assert len(gradient)==4
    for (a,b),(gl,gu) in zip(values[1:],gradient):
        ratios=[a/Hlo,a/Hhi,b/Hlo,b/Hhi]
        assert gl<=min(ratios)<=max(ratios)<=gu
    Delta=F(cert['Delta_lower']);assert Delta>0
    lower=[a-p for a,p in zip(lo,point[1:])]
    upper=[b-p for b,p in zip(hi,point[1:])]
    def linear(v):return sum(min(gl*x,gu*x) for x,(gl,gu) in zip(v,gradient))
    old=[]
    for bits in product((0,1),repeat=4):
        v=tuple((upper[j] if bits[j] else lower[j]) for j in range(4))
        old.append((Delta+linear(v)-baseq(v)/(eps*eps),bits))
    oldmin,oldbits=min(old)
    assert oldmin==F(cert['minimum_log_margin'])<=0
    assert list(oldbits)==cert['worst_corner_bits']
    sets=[sorted(set((a,F(0),b))) for a,b in zip(lower,upper)]
    new=[]
    for v in product(*sets):
        curvature=psi(v);norm2=sum(x*x for x in v)
        assert F(5,36)*norm2<=curvature<=F(2,9)*norm2
        assert curvature<=baseq(v)
        new.append((Delta+linear(v)-curvature/(eps*eps),v))
    newmin,worst=min(new);assert newmin>=oldmin
    rows.append({'node':case['node'],'root':case['root'],'anchor':case['point'],
      'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),
      'anchor_inside':True,'anchor_strictly_inside':all(a<p<b for a,p,b in zip(lo,point[1:],hi)),
      'cached_bits':case['bits'],'H_lower':str(Hlo),'Delta_lower':str(Delta),
      'loggradient_directed':cert['loggradient_directed'],
      'old_16_minimum':str(oldmin),'old_worst_corner_bits':list(oldbits),
      'new_point_count':len(new),'new_minimum':str(newmin),'new_worst_offset':list(map(str,worst)),
      'new_worst_on_zero_coordinate':any(x==0 for x in worst),'new_certificate_passed':newmin>0,
      'exact_improvement':str(newmin-oldmin)})
    if len(rows)%16==0:save('requalification_prefix')

signal.alarm(0)
assert sha(statepath)==reg['input_state_sha256'] and sha(regpath)==reg['input_registration_sha256']
passed=[r for r in rows if r['new_certificate_passed']]
failed=[r for r in rows if not r['new_certificate_passed']]
minnew=min(F(r['new_minimum']) for r in rows)
def floor_decimal(x,den=1000000):return F((x.numerator*den)//x.denominator,den)
summary={'boxes':len(rows),'original_failed_roots':len(rows),'newly_certified_parents':len(passed),
    'still_uncertified':len(failed),'all_anchors_inside':all(r['anchor_inside'] for r in rows),
    'strictly_inside_anchors':sum(r['anchor_strictly_inside'] for r in rows),
    'old_points_checked':16*len(rows),'new_points_checked':sum(r['new_point_count'] for r in rows),
    'new_worst_points_with_zero_coordinate':sum(r['new_worst_on_zero_coordinate'] for r in rows),
    'minimum_new_log_margin':str(minnew),'minimum_new_log_margin_rational_floor':str(floor_decimal(minnew)),
    'minimum_exact_improvement':str(min(F(r['exact_improvement']) for r in rows)),
    'passed_node_ids':[r['node'][0] for r in passed],'uncertified_node_ids':[r['node'][0] for r in failed],
    'H_gradient_Green_kernel_calls':0,'log_exp_Gamma_Arb_calls':0,
    'source_state_modified':False,'live_queue_or_controller_touched':False,'global_g_or_Cn_claim':False}
result={'status':'complete_frozen_162_root_multibasis_requalification','input_state_sha256':reg['input_state_sha256'],
        'input_fingerprint':reg['input_fingerprint'],'algebra':algebra,'summary':summary,'certificates':rows}
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
segments.append({'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'completed_boxes':len(rows)})
save('complete',include_rows=False)
receipt={'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'exit_code':0,'python':sys.version.split()[0],
         'sympy':sp.__version__,'threads':1,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
         'source_state_sha256_after':sha(statepath)}
(OUT/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'summary':{k:v for k,v in summary.items() if not k.endswith('_node_ids')},'receipt':receipt},indent=2))
