"""Audit-owned utilities. Only the accepted081 JSON streaming reader is shared."""
import sys
sys.dont_write_bytecode = True
import hashlib, importlib.util, json, os
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict

HERE = Path(__file__).resolve().parent
DESK = Path('/Users/zhengzhihao/Desktop/T')
SCRATCH = DESK/'cube-finite-dim-scratch/audit_g123_sol'
PROD = DESK/'cube-finite-dim-scratch/upper_continue123'
AUTHOR = DESK/'Math-Research-cube5-upper/cube-dim-5-bounds/upper_continue123'
AUDITS = DESK/'Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds'
FAN = DESK/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz'
FAN_SHA = '8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
KERNEL_SHA = '128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
MODEL = 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
TARGET = F(11,100000)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8<<20), b''): h.update(block)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()).hexdigest()

def write(path,value):
    path=Path(path); tmp=Path(str(path)+'.tmp')
    with tmp.open('w') as f:
        json.dump(value,f,sort_keys=True,indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
    tmp.replace(path)

def nodup(pairs):
    out={}
    for key,value in pairs:
        assert key not in out, 'duplicate JSON key'; out[key]=value
    return out

def read(path):
    return json.loads(Path(path).read_text(),object_pairs_hook=nodup)

def lines(path):
    with Path(path).open() as f:
        for row in f: yield json.loads(row,object_pairs_hook=nodup)

def check_locks(locks):
    for path,expected in locks.items(): assert sha(path)==expected, ('source changed',path)

_reader_path=AUDITS/'audit_kernel081/common.py'
assert sha(_reader_path)=='acc65ff6b32761cc242099a5f77f501d8541ecfc22a765c2d89ae3f9459ec028'
_spec=importlib.util.spec_from_file_location('readonly081_jsonreader',_reader_path)
_reader=importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_reader)
Stream=_reader.Stream

def node(raw):
    assert type(raw) is list and len(raw)==2
    r,p=raw; assert type(r) is int and r>=0 and type(p) is str and set(p)<=set('01')
    return r,p

def box(roots,raw):
    r,path=node(list(raw)); assert r<len(roots)
    lo=[];hi=[]
    for j,c in enumerate(roots[r]):
        bits=path[j::4];den=32*(1<<len(bits));i=(c<<len(bits))+(int(bits,2) if bits else 0)
        lo.append(F(i,den));hi.append(F(i+1,den))
    return lo,hi

def tree(roots,leaves,splits):
    leaves=list(map(lambda n:node(list(n)),leaves)); splits=set(map(lambda n:node(list(n)),splits))
    assert len(set(leaves))==len(leaves)
    by=defaultdict(list); prefixes=set()
    for r,p in leaves:
        assert r<len(roots);by[r].append(p);prefixes.update((r,p[:j]) for j in range(len(p)))
    assert set(by)==set(range(len(roots))) and prefixes==splits
    depth=0
    for r,paths in by.items():
        d=max(map(len,paths));depth=max(depth,d);cursor=0
        for a,b in sorted(((int(p,2) if p else 0)<<(d-len(p)),((int(p,2) if p else 0)+1)<<(d-len(p))) for p in paths):
            assert a==cursor and b>a, ('prefix overlap/Kraft gap',r);cursor=b
        assert cursor==1<<d
    assert len(leaves)==len(roots)+len(splits)
    return dict(closed_Kraft_roots=len(roots),leaves=len(leaves),splits=len(splits),maximum_depth=depth,prefix_free=True,all_root_Kraft_exactly_one=True)

def geometry(point,fan):
    from math import lcm
    p=list(map(F,point));den=lcm(*(x.denominator for x in p));ints=[int(x*den) for x in p]
    k=[sum(a*b for a,b in zip(v,ints)) for v in fan['vertices']]
    return dict(common_denominator=str(den),distinct_count=len(set(k)),zero_count=k.count(0),eligible=len(set(k))==280 and 0 not in k)

def snapshot(path,needed=()):
    """No H/log/gradient evaluation: stream old objects and compare exact SHA."""
    needed=set(needed);out=dict(cache={},records={},passed=[],splits=[],meta={});r=Stream(path)
    while True:
        key=r.key()
        if key is None:break
        if key=='cache_by_epsilon':
            r.char('{')
            for i in range(2):
                if i:r.char(',')
                eps=r.value();r.char(':');assert eps in ('1/128','1/64') and eps not in out['cache']
                hashes={};kept={}
                for p,rec in r.items(True):
                    assert ','.join(map(str,map(F,rec['point'])))==p and len(rec['values'])==5
                    hashes[p]=digest(rec)
                    if p in needed:kept[p]=rec
                out['cache'][eps]=hashes;out['records'][eps]=kept
            r.char('}')
        elif key in ('passed','new_passed'):
            for row in r.items():out['passed'].append((node(row['node']),digest(row)))
        elif key in ('splits','new_splits'):
            for row in r.items():
                n=node(row['node']);assert row['axis']==len(n[1])%4 and row['children']==[[n[0],n[1]+'0'],[n[0],n[1]+'1']]
                out['splits'].append((n,digest(row)))
        else:out['meta'][key]=r.value()
    r.close();return out
