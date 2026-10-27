import sys,os,json,time,resource,hashlib,subprocess
from pathlib import Path
P=Path(__file__).resolve().parent
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
def pins():
 for row in rd(P/'source_pins.json')['files']:
  q=Path(row['path']);assert sha(q)==row['sha256'] and q.stat().st_size==row['bytes']
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
r=rd(P/'registration.json');assert sha(P/'check_saved_768.py')==r['checker_sha256'] and sha(P/'source_pins.json')==r['source_pins_sha256']
with(P/'started_review.json').open('x')as f:json.dump({'observer_pid':os.getpid(),'unix_time':time.time(),'registration_sha256':sha(P/'registration.json')},f)
pins();A=Path(r['author_directory']);release=rd(A/'parent_release.json');started=rd(A/'started.json');terminal=rd(A/'exit_receipt.json')
assert release['main_approved']is True and release['independent_review_approved']is True and sha(A/'parent_release.json')==started['release_sha256']
assert release['independent_manifest_sha256']==sha(Path(r['preflight_review_path']).parent/'manifest.json')
assert terminal['status']=='LOCAL_BATCH_COMPLETE' and terminal['corner_exp_used']==95232 and started['pid']==30702 and terminal['native_calls']==0
before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.monotonic()
with(P/'checker_stdout.txt').open('x')as out,(P/'checker_stderr.txt').open('x')as err:
 child=subprocess.Popen(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-B',str(P/'check_saved_768.py')],cwd=P,stdout=out,stderr=err)
 save('checker_process.json',{'pid':child.pid,'observer_pid':os.getpid(),'unix_time':time.time(),'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'script_sha256':sha(P/'check_saved_768.py')});code=child.wait()
after=resource.getrusage(resource.RUSAGE_CHILDREN);unchanged=True;failure=None
try:pins()
except BaseException as e:unchanged=False;failure=repr(e)
receipt={'exit_code':code,'complete_wall_seconds':time.monotonic()-t,'complete_child_CPU_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,'child_reaped':True,'checker_pid':child.pid,'observer_pid':os.getpid(),'source_pins_unchanged_after':unchanged,'source_pin_failure':failure,'scope':'sole checker spawn through exit/reap plus final source hash check; observer receipt write excluded','children':1,'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None,'stdout_sha256':sha(P/'checker_stdout.txt'),'stderr_sha256':sha(P/'checker_stderr.txt')};save('actual_tool_receipt.json',receipt)
raise SystemExit(code if unchanged else 1)
