#!/usr/bin/env python3
"""Separate authorized batch: one exact rational G and one H_.032, then stop."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter,time
import ast,hashlib,importlib.util,itertools,json,math,os,signal
import flint
from flint import arb,fmpq,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
SOURCE=PROJECT/'upper/round007/g_time_moments/run_moments.py';GSOURCE=HERE/'evaluate_point.py'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ball(v):return {'ball':v.str(100,more=True),'lower':v.lower().str(100,more=True),'upper':v.upper().str(100,more=True)}
p=json.loads((HERE/'point2.json').read_text());assert p['max_new_passes']==2 and p['T']=='4/125'
# Load only the exact arithmetic function definitions; never execute the old script driver.
names={'generic_dd','laurent_dd','anchors','parse_fan','exact_green'}
tree=ast.parse(GSOURCE.read_text());selected=[n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name in names]
assert {n.name for n in selected}==names
exec(compile(ast.Module(body=selected,type_ignores=[]),str(GSOURCE),'exec'),globals())
spec=importlib.util.spec_from_file_location('b010_batch2_H',SOURCE);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.FAN=FAN
start=perf_counter();remaining=1200-(time()-p['batch2_start_epoch']);assert remaining>0
out={'status':'RUNNING_SEPARATE_BATCH2','pid':os.getpid(),'point_numerators':p['candidate_numerators'],'denominator':p['denominator'],'T':'4/125','new_full_passes':0,'original_batch_new_passes':5,'global_minimum_claimed':False,'new_C6_bound':False,'python':sys.version,'python_flint':flint.__version__,'arb_bits':384,'inputs':[{'path':str(f.resolve()),'sha256':sha(f)}for f in [HERE/'point2.json',HERE/'batch1_receipt.json',HERE/'batch2_registration.md',SOURCE,GSOURCE,FAN,Path(__file__)]]}
def write():
 out['elapsed_seconds']=perf_counter()-start;out['batch2_wall_seconds']=time()-p['batch2_start_epoch'];(HERE/'batch2_results.json').write_text(json.dumps(out,indent=2)+'\n')
def timeout(a,b):raise TimeoutError('Separate batch 2 deadline')
signal.signal(signal.SIGALRM,timeout);signal.alarm(int(remaining))
try:
 write();out['G']=exact_green(p['candidate_numerators'],p['denominator'],parse_fan());out['new_full_passes']+=1;write()
 print(json.dumps({'G':out['G']['G_arb384']['ball'],'seconds':out['G']['elapsed_seconds'],'zero_terms':out['G']['zero_node_terms'],'repeat_terms':out['G']['repeated_node_terms']}),flush=True)
 G=arb(out['G']['G_arb384']['ball']);h=m.evaluate(p['candidate_numerators'],p['denominator'],arb(4)/125,m.fan(),start,gradient=False)
 H=arb(h['H_arb384']);K=G-H;assert K>0
 out.update({'status':'COMPLETE_BATCH2_ALL_SEARCH_STOPPED','H_arb384':h['H_arb384'],'K_T':ball(K),'formal_point_cost':ball((arb(4)/125)/(64*K)),
 'refutes_strict_64_target':bool(K<=arb(1)/128000),'refutes_m_T_ge_T_over_4000':bool(K<arb(1)/125000),
 'H_complete_group_terms':h['complete_group_terms'],'H_computed_sign_paired_terms':h['computed_sign_paired_terms'],'H_elapsed_seconds':h['elapsed_seconds'],'new_full_passes':2,'cumulative_new_full_passes':7,'exit_code':0})
 write();print(json.dumps({k:v for k,v in out.items()if k not in ['G','inputs']},indent=2),flush=True)
except BaseException as e:out.update({'status':'ERROR_STOPPED','error':repr(e),'exit_code':1});write();raise
finally:signal.alarm(0)
