from pathlib import Path
import os,json,subprocess,time,resource
P=Path(__file__).resolve().parent
with(P/'observer_started.json').open('x')as f:json.dump({'pid':os.getpid()},f)
t=time.monotonic();r=resource.getrusage(resource.RUSAGE_CHILDREN)
with(P/'stdout.txt').open('x')as out,(P/'stderr.txt').open('x')as err:
 child=subprocess.Popen(['/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python','-E','-B','/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round132/cross_center_band_envelopes/run_batch.py'],stdin=subprocess.DEVNULL,stdout=out,stderr=err)
 with(P/'child_started.json').open('x')as f:json.dump({'pid':child.pid},f)
 code=child.wait()
a=resource.getrusage(resource.RUSAGE_CHILDREN)
d={'exit_code':code,'child_reaped':True,'producer_pid':child.pid,'observer_pid':os.getpid(),'complete_child_wall_seconds':time.monotonic()-t,'complete_child_CPU_seconds':a.ru_utime+a.ru_stime-r.ru_utime-r.ru_stime,'scope':'child spawn through exit/reap and stdio close; observer own receipt excluded','time_limit':None}
with(P/'terminal_receipt.json').open('x')as f:json.dump(d,f,indent=2)
raise SystemExit(code)
