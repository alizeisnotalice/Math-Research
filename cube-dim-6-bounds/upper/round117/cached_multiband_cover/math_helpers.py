"""C117 math/helpers extracted verbatim from released115; no launcher."""
import hashlib,importlib.util,json,os,subprocess
from pathlib import Path
from fractions import Fraction as Q

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def read(p):return json.loads(Path(p).read_text())

def canon(x):return(json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()

def digest(x):return hashlib.sha256(canon(x)).hexdigest()

def require(x,s):
    if not x:raise ValueError(s)

def new(p,x):
    with Path(p).open('xb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())

def atomic(p,x):
    p=Path(p);t=p.with_name(p.name+'.pending')
    with t.open('wb')as f:f.write(canon(x));f.flush();os.fsync(f.fileno())
    os.replace(t,p)

def load(name,p):
    s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def stop_owned(proc):
    if proc is not None and proc.poll()is None:
        proc.terminate()
        try:proc.wait(timeout=5)
        except subprocess.TimeoutExpired:proc.kill();proc.wait()

def add(i,j):return(i[0]+j[0],i[1]+j[1])

def scale(i,c):return tuple(sorted(c*x for x in i))

def euler(F,g,pa,pt,expected):
    a=Q(expected['window']['a']);T=Q(expected['window']['T'])
    terms=[scale(F,Q(5)),scale(pt,T),scale(pa,-a)]+[scale(v,Q(x))for v,x in zip(g,expected['center'])]
    out=tuple(sum(t[k]for t in terms)for k in[0,1]);require(out[0]<=0<=out[1],'Merged exact Euler');return out

def jet_from(expected,F,g,source,provenance):
    return{**{k:expected[k]for k in['model','model_definition','window','center']},'F_interval':list(map(str,F)),
        'face_gradient_intervals':[list(map(str,v))for v in g[1:]],'source_proof':source,'provenance':provenance}

def box_certificate(expected,jet):
    # Explicit general-window schema; sum jet is not forged native metadata.
    data={**{k:expected[k]for k in['model','model_definition','window','center','face_box','target','quadratic']},'jet':jet}
    with checker_lock:
        c=checker.certify(data,precision_bits=384);base.check_saved_box(jet,expected,c)
    return c
