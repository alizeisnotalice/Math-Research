import os,sys,json,time,resource,hashlib,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def pins():
 for row in rd(P/'source_pins.json')['files']:
  p=Path(row['path']);assert sha(p)==row['sha256']and p.stat().st_size==row['bytes']
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));reg=rd(P/'registration.json');A=Path(reg['author_directory']);assert sha(P/'check_saved_768.py')==reg['checker_sha256']and sha(P/'source_pins.json')==reg['source_pins_sha256']
with(P/'started_review.json').open('x')as f:json.dump({'observer_pid':os.getpid(),'registration_sha256':sha(P/'registration.json'),'unix_time':time.time()},f)
pins();terminal=rd(A/'exit_receipt.json');release=rd(A/'parent_release.json');started=rd(A/'started.json');primary=rd('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round134/primary_gate/terminal_receipt.json')
assert terminal['status']=='FINITE_FIXED_BATCH_COMPLETE'and terminal['native_calls']==0 and terminal['result_sha256']==sha(A/'result.json')and terminal['point_cache_sha256']==sha(A/'point_cache.json')
assert primary['exit_code']==0 and primary['child_reaped']and release['main_approved']is True and release['independent_review_approved']is True and started['release_sha256']==sha(A/'parent_release.json')
before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.monotonic()
with(P/'checker_stdout.txt').open('x')as out,(P/'checker_stderr.txt').open('x')as err:
 child=subprocess.Popen(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(P/'check_saved_768.py')],cwd=P,stdout=out,stderr=err)
 save('checker_process.json',{'pid':child.pid,'observer_pid':os.getpid(),'unix_time':time.time(),'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'script_sha256':sha(P/'check_saved_768.py')});code=child.wait()
after=resource.getrusage(resource.RUSAGE_CHILDREN);unchanged=True;failure=None
try:pins()
except BaseException as e:unchanged=False;failure=repr(e)
save('actual_tool_receipt.json',{'exit_code':code,'complete_wall_seconds':time.monotonic()-t,'complete_child_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'child_reaped':True,'checker_pid':child.pid,'observer_pid':os.getpid(),'source_pins_unchanged_after':unchanged,'source_pin_failure':failure,'scope':'sole saved checker spawn through exit/reap and final source hashes; observer receipt write excluded','children':1,'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'stdout_sha256':sha(P/'checker_stdout.txt'),'stderr_sha256':sha(P/'checker_stderr.txt')})
raise SystemExit(code if unchanged else 1)
