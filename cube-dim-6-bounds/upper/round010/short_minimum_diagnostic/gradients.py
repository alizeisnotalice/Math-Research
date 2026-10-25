#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter,time
import hashlib,importlib.util,json,os,signal
import flint
from flint import arb,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
OLD=PROJECT/'upper/round007/g_time_moments/results.json'
SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py'
FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
start=perf_counter();epoch=time();old=json.loads(OLD.read_text());Ts=[F(4,125),F(8,125)]
assert all(F(r['tau'])not in Ts for r in old['gradients'])
# Scan only complete gradient records with the same exact center; old null-valued H rows are not gradients.
scan=[];roots=[PROJECT/'upper',Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper'),Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')]
for root in roots:
 for p in root.rglob('*.json'):
  if HERE in p.parents:continue
  txt=p.read_text()
  if '998836847' not in txt or 'gradient_x2_to_x6' not in txt:continue
  z=json.loads(txt)
  if z.get('center_numerators')!=old['center_numerators']:continue
  grads=z.get('gradients',[])
  hits=[r for r in grads if F(r.get('tau','0'))in Ts and r.get('gradient_x2_to_x6')is not None]
  scan.append({'path':str(p),'sha256':sha(p),'new_target_gradient_matches':len(hits)})
  assert not hits,'Existing target gradient '+str(p)
spec=importlib.util.spec_from_file_location('audited_b010',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.FAN=FAN
inputs=[OLD,SOURCE,FAN,HERE/'registration.md',Path(__file__)]
out={'status':'RUNNING_TWO_OLD_CENTER_GRADIENTS','pid':os.getpid(),'batch_start_epoch':epoch,'max_total_seconds':1200,'center_numerators':old['center_numerators'],'denominator':old['denominator'],'gradient_passes':[],'new_full_passes':0,'duplicate_scan':scan,'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in inputs],'python':sys.version,'python_flint':flint.__version__,'arb_bits':ctx.prec}
def write():
 out['elapsed_seconds']=perf_counter()-start
 (HERE/'gradient_results.json').write_text(json.dumps(out,indent=2)+'\n')
def timeout(a,b):raise TimeoutError('B010 total 20 minute cap')
signal.signal(signal.SIGALRM,timeout);signal.alarm(1200)
try:
 cones=m.fan();write()
 for T in Ts:
  r=m.evaluate(old['center_numerators'],old['denominator'],arb(T.numerator)/T.denominator,cones,start,gradient=True)
  prior=next(v for v in old['time_nodes']if F(v['tau'])==T)
  assert (arb(r['H_arb384'])-arb(prior['H_arb384'])).contains(0)
  r['tau']=str(T);r['old_H_value_crosscheck']=True
  out['gradient_passes'].append(r);out['new_full_passes']+=1;write()
  print(json.dumps({'T':str(T),'H':r['H_arb384'],'gradient':r['gradient_x2_to_x6'],'seconds':r['elapsed_seconds']}),flush=True)
 out.update({'status':'COMPLETE_TWO_NEW_GRADIENT_PASSES','exit_code':0});write()
 print(json.dumps({'status':out['status'],'pid':out['pid'],'elapsed_seconds':out['elapsed_seconds']}),flush=True)
except BaseException as e:out.update({'status':'ERROR_STOPPED','error':repr(e),'exit_code':1});write();raise
finally:signal.alarm(0)
