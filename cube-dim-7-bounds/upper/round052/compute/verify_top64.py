#!/usr/bin/env python3
"""Independent Fraction recomputation of retained candidates, not a grid rescan."""
import json,hashlib,os,time,resource,collections
from pathlib import Path
from fractions import Fraction as Q
D=Path(__file__).resolve().parent;start=time.monotonic();cpu=time.process_time();data=json.loads((D/'results.json').read_bytes());reg=json.loads((D/'registration.json').read_bytes());assert data['totals']['tuples']==7372800==reg['tuple_count'] and data['chunks_completed']==1800

def check(v):
 a=Q(v['a']);rho=Q(v['rho']);s=Q(v['s']);d=Q(v['d']);P=Q(3579,4096);e=Q(1,1024);assert a>3 and s==2*a+rho and rho>0 and 0<d<=Q(1,16) and d<Q(1,32)/a
 U=P/(2*(1+e));H=1/(P*(1-d/2)*(1-e));mass=(1-Q(1,64))*d/(2*(1+3*d/4));M=max(Q(4),2*H/(1-3*d/4));br=(1-3/a)*U-2*e;assert br>0
 sig=s*d*d;assert sig<Q(1,2);dl=mass*br/M;delta=(1-e)*dl;c=delta-sig/(1-sig)
 vals={'u_lo':U,'u_hi':H,'mass':mass,'M':M,'bracket':br,'delta_lim':dl,'delta':delta,'sigma':sig,'penalty':sig/(1-sig),'c_raw':c}
 assert all(Q(v[k])==value for k,value in vals.items());return c
T=data['top64'];assert len(T)==64
for z in T:check(z)
keys=[(-Q(v['c_raw']),Q(v['a']),Q(v['rho']),Q(v['d'])) for v in T];assert keys==sorted(keys) and len(set(keys))==64
for w in data['workers'].values():
 assert len(w['top'])==32
 for z in w['top']:check(z)
chunks=collections.Counter()
for v in T:
 A=Q(v['a'])*8;k=Q(v['d'])*2**20;assert A.denominator==k.denominator==1;chunks[(int(A),(int(k)-1)//512)]+=1
mx=max(chunks.values());assert mx<32
# Run's frozen source keeps top 32 per chunk before global merge. If a discarded candidate
# belonged to top 64, all 32 higher-ranked retained candidates of that same chunk would
# also belong to top 64. None of the published chunks contributes 32, contradiction.
out={'status':'PASS','method':'Fraction recomputation of 64 global + 32 per-worker candidates; retained chunk occupancy <32 supplies conditional completeness of top64 from frozen per-chunk top32 algorithm','top64_count':len(T),'worker_count':len(data['workers']),'max_global_top64_per_chunk':mx,'chunk_counts':{str(k):v for k,v in chunks.items()},'raw_sha256':hashlib.sha256((D/'results.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'checks retained candidates and data-dependent top64 completeness assuming recorded frozen exhaustive integer enumeration; not an independent rescan or operator theorem'};(D/'verify_top64.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
