import subprocess,time,json,os,signal,sys
from pathlib import Path
base=Path(__file__).resolve().parent
out=Path(sys.argv[2])
out.mkdir(exist_ok=True)
cmd=[sys.executable,str(base/'variance72.py'),sys.argv[1],str(out)]
start=time.monotonic();peak=0.;rss=0
with (out/'native.log').open('w') as log:
 p=subprocess.Popen(cmd,stdout=log,stderr=log,start_new_session=True)
 while p.poll() is None:
  rows=[]
  for line in subprocess.check_output(['ps','-axo','pid=,ppid=,pcpu=,rss=,comm='],text=True).splitlines():
   a=line.split(None,4)
   if len(a)==5: rows.append((int(a[0]),int(a[1]),float(a[2]),int(a[3]),a[4]))
  owned={p.pid}
  for _ in range(6): owned.update(x[0] for x in rows if x[1] in owned)
  active=[x for x in rows if x[0] in owned and 'python' in x[4].lower()]
  peak=max(peak,sum(x[2] for x in active));rss=max(rss,sum(x[3] for x in active))
  if time.monotonic()-start>300:
   os.killpg(p.pid,signal.SIGKILL);p.wait();break
  time.sleep(.2)
 code=p.wait()
receipt={'exit_code':code,'wall_seconds':time.monotonic()-start,'peak_CPU_percent':peak,'peak_RSS_KiB':rss,'sample_interval_seconds':.2}
(out/'runtime72.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt));print((out/'native.log').read_text()[-3000:]);sys.exit(code)
