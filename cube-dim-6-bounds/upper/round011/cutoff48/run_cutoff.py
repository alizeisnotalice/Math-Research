#!/usr/bin/env python3
"""Exactly one new cutoff; reuse two saved exact Green values and stop as registered."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from time import perf_counter
import hashlib,importlib.util,json,os,signal
import flint
from flint import arb,fmpq,fmpz,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
G2=PROJECT/'upper/round010/short_minimum_diagnostic/batch2_results.json'
POINT2=PROJECT/'upper/round010/short_minimum_diagnostic/point2.json'
G0=PROJECT/'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json'
SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py'
FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
T=arb(6)/125;target64=arb(3)/256000;target625=arb(3)/250000

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(v):return {'ball':v.str(100,more=True),'lower':v.lower().str(100,more=True),'upper':v.upper().str(100,more=True)}
def exact_g(n,d):return arb(fmpq(fmpz(n),fmpz(d)))/arb.pi()**5
start=perf_counter();b2=json.loads(G2.read_text());p2=json.loads(POINT2.read_text());g0=json.loads(G0.read_text())
assert b2['point_numerators']==p2['candidate_numerators'] and b2['denominator']==p2['denominator']
assert b2['G']['complete_group_terms']==967680 and b2['exit_code']==0
v2=exact_g(b2['G']['G_times_pi5_numerator'],b2['G']['G_times_pi5_denominator'])
v0=exact_g(g0['Q_over_pi5_numerator'],g0['Q_over_pi5_denominator'])
assert (v2-arb(b2['G']['G_arb384']['ball'])).contains(0) and (v0-arb(g0['G_arb_384'])).contains(0)
points=[{'label':'B010_point2','numerators':b2['point_numerators'],'denominator':b2['denominator'],'G':v2,'G_source':str(G2)},
 {'label':'A007_reference','numerators':g0['point_numerators'],'denominator':g0['denominator'],'G':v0,'G_source':str(G0)}]
# Time-tag scan is deliberately broader than a point-only scan; no prior T=.048 record is accepted silently.
roots=[PROJECT/'upper',Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper'),Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper')]
scanned=0;matches=[]
for root in roots:
 for p in root.rglob('*.json'):
  if HERE in p.parents:continue
  scanned+=1;txt=p.read_text()
  if '"6/125"' in txt:matches.append({'path':str(p),'sha256':sha(p)})
assert not matches,'Prior exact T tag found; inspect before any kernel'
assert sha(FAN)=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
manifest={'status':'FROZEN_BEFORE_KERNEL','T':'6/125','point_order':[{'label':p['label'],'numerators':p['numerators'],'denominator':p['denominator'],'G_reconstructed':ball(p['G'])}for p in points],
 'duplicate_scan':{'roots':list(map(str,roots)),'json_files_scanned':scanned,'exact_T_tag_matches':matches},
 'max_H_passes':2,'max_seconds':1200,'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in [G2,POINT2,G0,SOURCE,FAN,HERE/'registration.md',Path(__file__)]]}
(HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
spec=importlib.util.spec_from_file_location('b011_existing_H',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.FAN=FAN
out={'status':'RUNNING_ONE_CUTOFF','pid':os.getpid(),'T':'6/125','new_H_passes':0,'new_G_passes':0,'new_gradients':0,'new_spatial_points':0,'rows':[],'source_sha256':sha(__file__),'manifest_sha256':sha(HERE/'input_manifest.json'),
 'python':sys.version,'python_executable':sys.executable,'python_flint':flint.__version__,'arb_bits':384,'global_m_T_certified':False,'new_C6_bound':False}
def write():
 out['elapsed_seconds']=perf_counter()-start;(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
def timeout(a,b):raise TimeoutError('Registered 20 minute cap')
signal.signal(signal.SIGALRM,timeout);signal.alarm(1200)
try:
 cones=m.fan();write()
 for p in points:
  hr=m.evaluate(p['numerators'],p['denominator'],T,cones,start,gradient=False)
  H=arb(hr['H_arb384']);K=p['G']-H;assert K>0 and hr['complete_group_terms']==967680
  row={'label':p['label'],'point_numerators':p['numerators'],'denominator':p['denominator'],'G_reconstructed_no_kernel':ball(p['G']),'H_arb384':hr['H_arb384'],'K_T':ball(K),'formal_point_cost':ball(T/(64*K)),
   'threshold64':'3/256000','target625':'3/250000','strictly_refutes64':bool(K<=target64),'necessary64_condition_passes':bool(K>target64),'strictly_refutes625':bool(K<target625),'necessary625_condition_passes':bool(K>=target625),
   'complete_group_terms':hr['complete_group_terms'],'computed_sign_paired_terms':hr['computed_sign_paired_terms'],'elapsed_seconds':hr['elapsed_seconds']}
  out['rows'].append(row);out['new_H_passes']+=1;write();print(json.dumps(row),flush=True)
  if row['strictly_refutes64']:
   out['status']='STOP_FIXED_T_REFUTED_FOR_STRICT64';break
  if not row['necessary64_condition_passes']:
   out['status']='STOP_PRECISION_UNRESOLVED';break
 else:out['status']='COMPLETE_TWO_NECESSARY_CONDITIONS_ONLY_STOPPED'
 out['exit_code']=0;write();print(json.dumps({k:v for k,v in out.items()if k!='rows'},indent=2),flush=True)
except BaseException as e:out.update({'status':'ERROR_STOPPED','error':repr(e),'exit_code':1});write();raise
finally:signal.alarm(0)
