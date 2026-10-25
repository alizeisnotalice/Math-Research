#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from time import perf_counter,time
import hashlib,importlib.util,json,os,signal
import flint
from flint import arb,fmpz,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py';FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(v):return {'ball':v.str(100,more=True),'lower':v.lower().str(100,more=True),'upper':v.upper().str(100,more=True)}
r=json.loads((HERE/'point1_results.json').read_text());gr=json.loads((HERE/'gradient_results.json').read_text())
assert r['status']=='STOP_FIRST_TRUE_TARGET_COUNTEREXAMPLE' and r['H_rows'][0]['T']=='8/125' and r['H_rows'][0]['refutes_strict_64_target']
remaining=1200-(time()-gr['batch_start_epoch']);assert remaining>0
G=arb(fmpq(fmpz(r['G']['G_times_pi5_numerator']),fmpz(r['G']['G_times_pi5_denominator'])))/arb.pi()**5
assert (G-arb(r['G']['G_arb384']['ball'])).contains(0)
spec=importlib.util.spec_from_file_location('b010_extra_H',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.FAN=FAN
start=perf_counter()
out={'status':'RUNNING_ONE_ADDITIONAL_H_SAME_POINT','pid':os.getpid(),'T':'4/125','point_numerators':r['point_numerators'],'denominator':r['denominator'],'G_reconstructed_no_replay':ball(G),'new_H_passes':0,'new_G_passes':0,'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in [HERE/'point1_results.json',HERE/'continuation_registration.md',SOURCE,FAN,Path(__file__)]],'python':sys.version,'python_flint':flint.__version__,'arb_bits':384}
def write():
 out['elapsed_seconds']=perf_counter()-start;out['batch_elapsed_seconds']=time()-gr['batch_start_epoch'];(HERE/'other_time_results.json').write_text(json.dumps(out,indent=2)+'\n')
def timeout(a,b):raise TimeoutError('Original total batch deadline')
signal.signal(signal.SIGALRM,timeout);signal.alarm(int(remaining))
try:
 write();h=m.evaluate(r['point_numerators'],r['denominator'],arb(4)/125,m.fan(),start,gradient=False)
 H=arb(h['H_arb384']);K=G-H;assert K>0
 out.update({'status':'COMPLETE_OTHER_TIME_STOPPED','H_arb384':h['H_arb384'],'K_T':ball(K),'formal_point_cost':ball((arb(4)/125)/(64*K)),
 'refutes_m_T_ge_T_over_4000':bool(K<arb(1)/125000),'refutes_strict_64_target':bool(K<=arb(1)/128000),
 'complete_group_terms':h['complete_group_terms'],'computed_sign_paired_terms':h['computed_sign_paired_terms'],'new_H_passes':1,'cumulative_new_full_passes':5,'global_minimum_claimed':False,'new_C6_bound':False,'exit_code':0})
 write();print(json.dumps({k:v for k,v in out.items()if k!='inputs'},indent=2),flush=True)
except BaseException as e:out.update({'status':'ERROR_STOPPED','error':repr(e),'exit_code':1});write();raise
finally:signal.alarm(0)
