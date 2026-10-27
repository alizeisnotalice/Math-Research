import json,hashlib,time,resource,subprocess,os,sys,itertools,math
from pathlib import Path
from fractions import Fraction as Q
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rd=lambda p:json.loads(p.read_text())
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def pins():
 for x in rd(P/'source_pins.json')['files']:
  q=Path(x['path']);assert sha(q)==x['sha256'] and q.stat().st_size==x['bytes']
r=rd(P/'registration.json');assert sha(P/'check_saved_768.py')==r['script_sha256'] and sha(P/'source_pins.json')==r['source_pins_sha256']
with(P/'started_review.json').open('x')as f:json.dump({'observer_pid':os.getpid(),'unix_time':time.time(),'registration_sha256':sha(P/'registration.json')},f)
pins();a=Path(r['author_directory']);rc=rd(a/'run_receipt.json');obs=rd(a/'observed_tool_receipt.json');rel=rd(a/'parent_release.json')
assert rc['status']=='COMPLETE_FIXED_FINITE_BOX_PILOT' and rc['failure']is None and obs['exit_code']==0 and obs['child_reaped'] and obs['release_sha256']==sha(a/'parent_release.json') and rel['script_sha256']==sha(a/'run_prepared.py')
spec=rd(P/'inputs.json');old=rd(Path(spec['old_caps_path']));up=rd(Path(spec['round128_input_path']));table=rd(P/'cone_table.json')
assert list(map(Q,spec['old_beta2_interval']))==[Q(old['caps']['b_lower'][2]),Q(old['caps']['beta_upper'][2])]
assert Q(old['h'])==Q(spec['h'])==Q(up['h_values'][5])==Q(1,128) and Q(spec['old_beta2_interval'][1])==Q(up['beta_by_h'][5][1]) and spec['p_interval']==up['p_intervals'][5]
V=[(1,)+x for x in itertools.product((-1,1),repeat=5)];tokens=list(map(int,Path(spec['fan_path']).read_text().split()));cursor=1;assert tokens[0]==len(table['cones'])==21
for cid,cone in enumerate(table['cones']):
 det=tokens[cursor];cursor+=1;rays=[];rates=[]
 for _ in range(6):rates.append(tokens[cursor]);rays.append(tokens[cursor+1:cursor+7]);cursor+=7
 assert cone['id']==cid and cone['rays']==rays and cone['rates']==rates and cone['determinant']==det
 assert det>0 and min(rates)>0 and Q(cone['alpha'])==Q(det,math.prod(rates))
 nodes=[]
 for v in V:
  ds=[sum(a*b for a,b in zip(v,ray))for ray in rays];assert min(ds)>=0 or max(ds)<=0
  nodes.append([str(Q(abs(d),rate))for d,rate in zip(ds,rates)])
 assert cone['nodes_by_body']==nodes
assert cursor==len(tokens) and 46080*sum(Q(c['alpha'])for c in table['cones'])==Q(spec['R0'])
save('input_chain_review.json',{'status':'PASS_EXACT_SOURCE_TABLE_AND_OLD_CAP_CHAIN','cones':21,'body_labels':32,'pairs':496,'source_fan_producer_replayed':False,'old_baseline_recomputed':False,'producer_terminal_exit':obs['exit_code'],'source_files_hashed':len(rd(P/'source_pins.json')['files'])})
t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with(P/'checker_stdout.txt').open('x')as out,(P/'checker_stderr.txt').open('x')as err:
 proc=subprocess.Popen(['/usr/bin/nice','-n','10','/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(P/'check_saved_768.py')],cwd=P,stdout=out,stderr=err)
 save('checker_process.json',{'pid':proc.pid,'observer_pid':os.getpid(),'unix_time':time.time(),'nice_minimum':10,'time_limit':None,'script_sha256':sha(P/'check_saved_768.py')})
 code=proc.wait()
after=resource.getrusage(resource.RUSAGE_CHILDREN);unchanged=True;failure=None
try:pins()
except BaseException as e:unchanged=False;failure=repr(e)
receipt={'exit_code':code,'complete_wall_seconds':time.monotonic()-t,'complete_checker_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'checker_pid':proc.pid,'observer_pid':os.getpid(),'child_reaped':True,'source_pins_unchanged_after':unchanged,'source_pin_failure':failure,'scope':'sole checker spawn through exit/reap and final sourcepin check; observer receipt write excluded','children':1,'time_limit':None,'stdout_sha256':sha(P/'checker_stdout.txt'),'stderr_sha256':sha(P/'checker_stderr.txt')}
save('actual_tool_receipt.json',receipt)
raise SystemExit(code if unchanged else 1)
