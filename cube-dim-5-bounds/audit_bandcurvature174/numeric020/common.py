"""Exact metadata IO; no production or H evaluator imports."""
from pathlib import Path
import hashlib,json,gzip,os,time
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def read(p):
    p=Path(p)
    with gzip.open(p,'rt') if p.suffix=='.gz' else p.open() as f:return json.load(f)
def lines(p):
    p=Path(p)
    with gzip.open(p,'rt') if p.suffix=='.gz' else p.open() as f:
        for raw in f:
            assert raw.endswith('\n'),'truncated prefix requires manual review'
            yield json.loads(raw)
def atomic(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    with tmp.open('w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def chain(p):
    previous='0'*64;out=[]
    for r in lines(p):
        body={k:v for k,v in r.items() if k!='hash'}
        assert r['seq']==len(out)+1 and r['previous']==previous and digest(body)==r['hash']
        previous=r['hash'];out.append(r)
    return out
def append(f,body):
    h=digest(body);f.write((json.dumps({**body,'hash':h},sort_keys=True,separators=(',',':'))+'\n').encode());f.flush();os.fsync(f.fileno());return h
def locks(d):
    for p,h in d.items():assert sha(p)==h,p
def fault_cleanup(ps,stop,start,release,path):
    stop.set();start.set();release.set();began=time.monotonic();actions=[]
    atomic(path,{'status':'fault_observed_before_bounded_cleanup','worker_pids':[p.pid for p in ps],'normal_wall_cap':None})
    def wait(seconds):
        end=time.monotonic()+seconds
        while any(p.is_alive() for p in ps) and time.monotonic()<end:
            for p in ps:p.join(timeout=0)
            time.sleep(min(.02,max(0,end-time.monotonic())))
    wait(5)
    for p in ps:
        if p.is_alive():actions.append({'action':'terminate','pid':p.pid});p.terminate()
    wait(2)
    for p in ps:
        if p.is_alive():actions.append({'action':'kill','pid':p.pid});p.kill()
    wait(2)
    for p in ps:p.join(timeout=0)
    r={'status':'fault_only_cleanup','actions':actions,'workers':[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive(),'join_called':True} for p in ps],'seconds_after_fault':time.monotonic()-began,'normal_wall_cap':None};atomic(path,r);return r
