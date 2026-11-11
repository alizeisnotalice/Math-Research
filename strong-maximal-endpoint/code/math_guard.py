"""Run a local math batch while preserving available memory; only kill its own tree."""
import argparse,subprocess,time,os,signal,re,json
from pathlib import Path

def available_bytes(total):
    try:
        s=subprocess.check_output(['memory_pressure'],timeout=4,text=True)
        z=re.search(r'System-wide memory free percentage:\s*(\d+)%',s)
        return total*int(z[1])/100 if z else None
    except (subprocess.SubprocessError,OSError):return None

def process_tree(root):
    rows=[]
    for line in subprocess.check_output(['ps','-Ao','pid=,ppid=,rss=']).decode().splitlines():
        p=line.split()
        if len(p)==3:rows.append(tuple(map(int,p)))
    own={root}
    for _ in range(8):own|={p for p,pp,r in rows if pp in own}
    return [(p,r*1024) for p,pp,r in rows if p in own]
def stop_tree(root):
    rows=process_tree(root)
    for p,r in reversed(rows):
        try:os.kill(p,signal.SIGTERM)
        except ProcessLookupError:pass

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reserve-gib',type=float,default=1);ap.add_argument('--max-rss-gib',type=float,default=6.5);ap.add_argument('--log',required=True);ap.add_argument('--report',required=True);ap.add_argument('command',nargs=argparse.REMAINDER);a=ap.parse_args();cmd=a.command[1:] if a.command[:1]==['--'] else a.command;assert cmd
    total=int(subprocess.check_output(['sysctl','-n','hw.memsize']));reserve=a.reserve_gib*2**30;av=available_bytes(total)
    if av is None or av<=reserve:raise SystemExit('Insufficient measurable available memory; batch not started.')
    budget=min(a.max_rss_gib*2**30,av-reserve);Path(a.log).parent.mkdir(parents=True,exist_ok=True);peak=0;low=av;reason=None;bad=0;count=0;start=time.monotonic()
    with open(a.log,'w') as log:
        p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
        def interrupted(sig,frame):
            stop_tree(p.pid);raise KeyboardInterrupt
        signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
        try:
            while p.poll() is None:
                rss=sum(r for pid,r in process_tree(p.pid));peak=max(peak,rss);now=available_bytes(total);count+=1
                if now is None:bad+=1
                else:bad=0;low=min(low,now)
                if rss>budget:reason='owned_process_RSS_budget';stop_tree(p.pid);break
                if now is not None and now<reserve:reason='system_available_memory_reserve';stop_tree(p.pid);break
                if bad>=3:reason='available_memory_probe_failed';stop_tree(p.pid);break
                time.sleep(1)
            try:code=p.wait(timeout=10)
            except subprocess.TimeoutExpired:stop_tree(p.pid);code=p.wait(timeout=10)
        finally:
            if p.poll() is None:stop_tree(p.pid)
    z=dict(exit_code=code,stop_reason=reason,reserve_GiB=a.reserve_gib,owned_RSS_budget_GiB=budget/2**30,peak_owned_RSS_MiB=peak/2**20,min_system_available_MiB=low/2**20,samples=count,seconds=time.monotonic()-start,
      scope='macOS memory_pressure available-percentage estimate; monitor plus per-kernel MemoryConstrained, not a hard OS reservation.')
    Path(a.report).write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(z));raise SystemExit(0 if code==0 and reason is None else 3)
if __name__=='__main__':main()
