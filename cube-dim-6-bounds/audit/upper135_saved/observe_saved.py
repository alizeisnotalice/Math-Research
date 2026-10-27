import sys,os,json,hashlib,time,resource,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def pins():
 for row in rd(P/'source_pins.json')['files']:
  q=Path(row['path']);assert sha(q)==row['sha256'] and q.stat().st_size==row['bytes']
 for row in rd(A/'source_manifest.json')['sources']:
  q=Path(row['path']);assert sha(q)==row['sha256'] and q.stat().st_size==row['bytes']
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
r=rd(P/'registration.json');A=Path(r['author_directory']);assert sha(P/'check_new_point.py')==r['checker_sha256'] and sha(P/'source_pins.json')==r['source_pins_sha256']
with(P/'started_review.json').open('x')as f:json.dump({'observer_pid':os.getpid(),'unix_time':time.time(),'registration_sha256':sha(P/'registration.json')},f)
pins();obs=rd(A/'observed_tool_receipt.json');run=rd(A/'run_receipt.json');release=rd(A/'parent_release.json');start=rd(A/'started.json')
assert obs['exit_code']==0 and obs['child_reaped'] and run['status']=='COMPLETE_ONE_FULL_GREEN_POINT'and run['failure']is None
assert obs['release_sha256']==start['release_sha256']==sha(A/'parent_release.json') and release['script_sha256']==sha(A/'run_point.py')and release['registration_sha256']==sha(A/'registration.json')
t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with(P/'checker_stdout.txt').open('x')as out,(P/'checker_stderr.txt').open('x')as err:
 child=subprocess.Popen(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(P/'check_new_point.py')],cwd=P,stdout=out,stderr=err)
 save('checker_process.json',{'pid':child.pid,'observer_pid':os.getpid(),'unix_time':time.time(),'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'script_sha256':sha(P/'check_new_point.py')});code=child.wait()
after=resource.getrusage(resource.RUSAGE_CHILDREN);unchanged=True;failure=None
try:pins()
except BaseException as e:unchanged=False;failure=repr(e)
save('actual_tool_receipt.json',{'exit_code':code,'complete_wall_seconds':time.monotonic()-t,'complete_child_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'child_reaped':True,'checker_pid':child.pid,'observer_pid':os.getpid(),'all_new_and_archived_source_pins_unchanged_after':unchanged,'source_pin_failure':failure,'scope':'sole independent new-point checker spawn through exit/reap and final source hash check; observer receipt writing excluded','children':1,'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'stdout_sha256':sha(P/'checker_stdout.txt'),'stderr_sha256':sha(P/'checker_stderr.txt')})
raise SystemExit(code if unchanged else 1)
