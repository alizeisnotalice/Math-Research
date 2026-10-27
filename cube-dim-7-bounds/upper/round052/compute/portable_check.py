#!/usr/bin/env python3
"""R052 retained certificates; stdlib, relative paths, no full-grid rerun."""
import json,hashlib,collections
from pathlib import Path
from fractions import Fraction as Q
D=Path(__file__).resolve().parent;sha=lambda b:hashlib.sha256(b).hexdigest()
rraw=(D/'registration.json').read_bytes();reg=json.loads(rraw);data=json.loads((D/'results.json').read_bytes());assert data['registration_sha256']==sha(rraw) and reg['code_sha256']==sha((D/'r052_integer.py').read_bytes()) and reg['equivalence_sha256']==sha((D/'equivalence.json').read_bytes())
assert len(data['top64'])==64 and len(data['workers'])==15 and data['chunks_completed']==1800 and data['totals']['tuples']==7372800
valid=8*sum(min(4096,(2**20-1)//(4*A)) for A in range(32,257));assert data['totals']['valid']==valid==3975632 and data['totals']['d_guard']==3397168 and data['totals']['sigma_guard']==data['totals']['positive_deficit_guard']==0
assert sum(w['chunks'] for w in data['workers'].values())==1800 and sum(w['tuples'] for w in data['workers'].values())==7372800

def verify(v):
 a=Q(v['a']);rho=Q(v['rho']);s=Q(v['s']);d=Q(v['d']);P=Q(3579,4096);e=Q(1,1024)
 A=8*a;k=(2**20)*d;t=16*rho;assert A.denominator==k.denominator==t.denominator==1 and 32<=A<=256 and 1<=k<=4096 and int(t) in (1,2,4,8,16,32,64,128)
 assert a>3 and s==2*a+rho and 0<d<=Q(1,16) and d<Q(1,32)/a
 ulo=P/(2*(1+e));uhi=1/(P*(1-d/2)*(1-e));mass=(1-Q(1,64))*d/(2*(1+3*d/4));M=max(Q(4),2*uhi/(1-3*d/4));br=(1-3/a)*ulo-2*e;sig=s*d*d;assert br>0 and sig<Q(1,2)
 lim=mass*br/M;delta=(1-e)*lim;c=delta-sig/(1-sig);calc={'u_lo':ulo,'u_hi':uhi,'mass':mass,'M':M,'bracket':br,'delta_lim':lim,'delta':delta,'sigma':sig,'penalty':sig/(1-sig),'c_raw':c}
 assert all(Q(v[k])==x for k,x in calc.items());return (-c,a,rho,d)
keys=[verify(v) for v in data['top64']];assert keys==sorted(keys) and len(set(keys))==64
for w in data['workers'].values():
 assert len(w['top'])==32
 for v in w['top']:verify(v)
occupancy=collections.Counter((int(8*Q(v['a'])),(int(2**20*Q(v['d']))-1)//512) for v in data['top64']);assert max(occupancy.values())<32
print(json.dumps({'status':'PASS','retained_global':64,'workers':15,'worker_retained_each':32,'max_top64_chunk_occupancy':max(occupancy.values()),'independent_list_count':valid,'scope':'retained arithmetic only; top64 completeness conditional on frozen chunk execution/retention; no full-grid rescan or independent operator proof'},indent=2))
