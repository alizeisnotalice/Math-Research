#!/usr/bin/env python3
"""Portable exact R044/R038 bundle verifier, stdlib only; no optimization."""
import os,sys,json,hashlib,gzip,time,resource,collections
from fractions import Fraction as Q
from pathlib import Path
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
D=Path(__file__).resolve().parent;OLD=D.parent.parent/'round038/compute';sha=lambda b:hashlib.sha256(b).hexdigest();canon=lambda x:json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()
start=time.monotonic();cpu=time.process_time();man=json.loads((D/'archive_manifest.json').read_bytes());packed=(D/'r044_compact.json.gz').read_bytes();assert sha(packed)==man['bundle_sha256'];raw=gzip.decompress(packed);assert sha(raw)==man['content_sha256'];data=json.loads(raw);oldraw=(OLD/'phasea_compact.json.gz').read_bytes();assert sha(oldraw)==data['R038_compact_sha256']==man['source_R038_compact_sha256'];old=json.loads(gzip.decompress(oldraw));lookup={x['raw_filename']:x for x in old['cases']};assert len(data['cases'])==76
counts=collections.Counter()
for r in data['cases']:
 c=lookup[r['source_raw']];assert c['raw_sha256']==r['source_raw_sha256'] and c['q_after_canonical_sha256']==r['source_final_q_sha256']==sha(canon(r['q_initial']))
 h=list(map(Q,c['witness_budgets']));C=c['maximal_cliques'];q=list(map(Q,r['q_initial']));g=list(map(Q,r['global_proposed_q']));z=list(map(Q,r['mass_primal']));y=list(map(Q,r['mass_dual']));U=Q(r['U']);n=len(h)
 assert U==Q(c['U']) and r['n']==c['n'] and r['a']==c['a'] and len(q)==len(g)==len(y)==n and len(z)==len(C) and all(v>0 for v in h) and all(v>=0 for v in q+g+z+y)
 assert all(sum((z[j] for j,cl in enumerate(C) if i in cl),Q())>=h[i] for i in range(n));assert all(sum((y[i] for i in cl),Q())<=1 for cl in C)
 lo=sum((h[i]*y[i] for i in range(n)),Q());hi=sum(z,Q());assert lo==Q(r['M_lower']) and hi==Q(r['M_upper']) and 0<lo<=hi
 load=lambda v:max(sum((v[i] for i in cl),Q()) for cl in C);pay=lambda v:sum((h[i]*v[i] for i in range(n)),Q());assert pay(q)==U==pay(g) and sha(canon(list(map(str,q))))==r['source_final_q_sha256'];L=load(q);G=load(g);assert L==Q(r['L_initial'])==Q(c['L_after']) and G==Q(r['global_proposed_load']) and G<=U/lo and (G<L)==r['global_improves'];counts['global_improved']+=G<L
 rr=r['triple_retry'];v=list(map(Q,rr['best_q']));assert len(v)==n and all(x>=0 for x in v) and pay(v)==U;V=load(v);assert V==Q(rr['L_best']) and (V<L)==rr['strict_improvement'];tri=rr['best_triple'];assert (v==q if tri is None else {i for i in range(n) if v[i]!=q[i]}.issubset(set(tri)));counts['first32_improved']+=V<L
 b=r['B'];bc=b['candidate']
 if bc:
  w=list(map(Q,bc['q_new']));assert len(w)==n and all(x>=0 for x in w) and pay(w)==U;W=load(w);assert W==Q(bc['L_new'])<L and {i for i in range(n) if w[i]!=q[i]}.issubset(set(bc['triple']));counts['B_improved']+=1
 c0=r['C'];cand=c0['exact_candidate']
 if cand:
  w=list(map(Q,cand['q_feasible']));assert len(w)==n and all(x>=0 for x in w) and pay(w)==U;W=load(w);T=(L+G)/2;assert G<L and W==Q(cand['exact_load'])<=T==Q(cand['exact_target']);changed=[i for i in range(n) if w[i]!=q[i]];assert changed==cand['exact_changed_indices'] and len(changed)==cand['exact_support_count'];counts['C_feasible']+=1;counts['C_support3']+=len(changed)==3
 else:counts['C_no_exact_vector']+=1
 counts['cases']+=1
out={'status':'PASS','counts':dict(counts),'bundle_sha256':man['bundle_sha256'],'R038_compact_sha256':man['source_R038_compact_sha256'],'script_sha256':sha(Path(__file__).read_bytes()),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'verifies exact saved vectors and LP feasibility against audited cached cliques; producer-data replay, not independent geometry reconstruction or integer optimality proof'};(D/'archive_verify_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
