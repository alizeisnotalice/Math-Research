#!/usr/bin/env python3
"""Replay recorded nonedge pairs, not the original all-pair search; optional terminal scan."""
import sys,os,json,gzip,hashlib,time,resource,concurrent.futures
from fractions import Fraction as Q
from pathlib import Path
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent
sha=lambda b:hashlib.sha256(b).hexdigest()
canonical=lambda obj:json.dumps(obj,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()
def bits(C):
 while C:
  bit=C&-C;C-=bit;yield bit.bit_length()-1

def step(q,h,C,i,j):
 pay=h[i]*q[i]+h[j]*q[j];assert pay>0
 weights=[sum((q[k] for k in bits(c)),Q()) for c in C]
 Ai=max((v-q[i] for c,v in zip(C,weights) if c>>i&1),default=Q());Aj=max((v-q[j] for c,v in zip(C,weights) if c>>j&1),default=Q());rest=max((v for c,v in zip(C,weights) if not c&((1<<i)|(1<<j))),default=Q())
 cross=max(Q(),min(pay/h[i],(pay+h[j]*(Aj-Ai))/(h[i]+h[j])));m=max(rest,Ai+cross,Aj+(pay-h[i]*cross)/h[j]);lo=max(Q(),(pay-h[j]*(m-Aj))/h[i]);hi=min(pay/h[i],m-Ai);assert lo<=hi
 x=max(lo,min(hi,pay/(h[i]+h[j])));y=(pay-h[i]*x)/h[j]
 delta=h[i]*(x*x-q[i]*q[i])+h[j]*(y*y-q[j]*q[j]);return x,y,m,delta

def verify(c,terminal):
 t=time.monotonic();cpu=time.process_time();h=list(map(Q,c['witness_budgets']));q=list(map(Q,c['q_before']));C=[sum((1<<i for i in ids),0) for ids in c['maximal_cliques']];n=len(q);U=Q(c['U']);adj=[0]*n
 for clique in C:
  for i in bits(clique):adj[i]|=clique&~(1<<i)
 L=max(sum((q[i] for i in bits(clique)),Q()) for clique in C);assert L==Q(c['L_before']) and sum((h[i]*q[i] for i in range(n)),Q())==U
 for i,j in c['accepted_pairs']:
  assert i!=j and not (adj[i]>>j)&1
  x,y,Ln,de=step(q,h,C,i,j);assert (Ln,de)<(L,Q());before=h[i]*q[i]+h[j]*q[j];q[i]=x;q[j]=y;assert before==h[i]*x+h[j]*y;L=Ln
 assert L==Q(c['L_after']) and sum((h[i]*q[i] for i in range(n)),Q())==U and sha(canonical(list(map(str,q))))==c['q_after_canonical_sha256']
 if c['stop']=='iteration_cap':assert len(c['accepted_pairs'])==128
 elif c['stop']=='stop_pair_coordinate' and terminal:
  for i in range(n):
   for j in range(i+1,n):
    if (adj[i]>>j)&1 or not (q[i] or q[j]):continue
    _,_,Ln,de=step(q,h,C,i,j);assert not (Ln,de)<(L,Q())
 else:assert c['stop'] in ('stop_pair_coordinate','iteration_cap')
 return {'case':c['raw_filename'],'steps':len(c['accepted_pairs']),'stop':c['stop'],'terminal_scan_done':bool(terminal and c['stop']=='stop_pair_coordinate'),'pid':os.getpid(),'wall_seconds':time.monotonic()-t,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
if __name__=='__main__':
 start=time.monotonic();compressed=(D/'phasea_compact.json.gz').read_bytes();man=json.loads((D/'compact_manifest.json').read_bytes());assert sha(compressed)==man['compressed_sha256'];body=gzip.decompress(compressed);assert sha(body)==man['canonical_content_sha256'];d=json.loads(body);assert len(d['cases'])==79
 case_map={c['raw_filename']:c for c in d['cases']}
 for a in man['three_lp_certificates']:
  raw=(D/a['filename']).read_bytes();assert sha(raw)==a['sha256'];cert=json.loads(raw);c=case_map[a['source_case']];assert cert['source']['sha256']==c['raw_sha256']
  h=list(map(Q,c['witness_budgets']));C=c['maximal_cliques'];z=list(map(Q,cert['primal_m']));y=list(map(Q,cert['dual_y']));U=Q(c['U'])
  assert all(x>=0 for x in z+y) and len(z)==len(C) and len(y)==len(h)
  assert all(sum((z[j] for j,V in enumerate(C) if i in V),Q())>=h[i] for i in range(len(h)))
  assert all(sum((y[i] for i in V),Q())<=1 for V in C)
  ml=sum((h[i]*y[i] for i in range(len(h))),Q());mu=sum(z,Q());assert ml==Q(cert['M_lower']) and mu==Q(cert['M_upper']) and 0<ml<=mu
  assert Q(cert['Lopt_lower'])==U/mu and Q(cert['Lopt_upper'])==U/ml and Q(cert['Lfinal'])==Q(c['L_after'])
  if cert['verdict']=='H_fail':assert Q(cert['L0'])-U/ml>=Q(cert['L0'])/10000 and 2*Q(cert['Lfinal'])>Q(cert['L0'])+U/ml
 terminal='--terminal' in sys.argv
 with concurrent.futures.ProcessPoolExecutor(max_workers=8) as ex:rows=list(ex.map(verify,d['cases'],[terminal]*len(d['cases'])))
 out={'status':'PASS','scope':'recorded pair steps only, not a replay of pair ranking or production search; terminal full scan only when requested','terminal_requested':terminal,'case_count':len(rows),'steps':sum(x['steps'] for x in rows),'lp_exact_certificates_checked':len(man['three_lp_certificates']),'bundle_sha256':man['compressed_sha256'],'script_sha256':sha(Path(__file__).read_bytes()),'wall_seconds':time.monotonic()-start,'results':rows}
 (D/'replay_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',len(rows),'cases',out['steps'],'steps; terminal scan:',terminal,'wall',out['wall_seconds'])
