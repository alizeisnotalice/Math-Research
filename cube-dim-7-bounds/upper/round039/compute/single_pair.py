#!/usr/bin/env python3
"""R039: exact one-step scan with 4 clique branches and lex(L,S,E)."""
import os,sys,time,json,hashlib,resource,concurrent.futures
from pathlib import Path
from fractions import Fraction as Q
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
P=Path(__file__).resolve();D=P.parent;OLD=D.parent.parent/'round038/compute';PROOF=Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round039/joint_cover/PROOF.md')
FILES=('phasea_C_2_05_full_disjoint.json','phasea_C_1_06_full_disjoint.json','phasea_A_3_15_partial.json')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def bits(C):
 while C:
  b=C&-C;C-=b;yield b.bit_length()-1

def optimize(q,h,cliques,weights,i,j):
 hi=h[i];hj=h[j];pay=hi*q[i]+hj*q[j]
 if not pay:return None
 groups={}
 for C,w in zip(cliques,weights):
  bi=(C>>i)&1;bj=(C>>j)&1;tag=str(bi)+str(bj);outside=w-bi*q[i]-bj*q[j]
  if tag not in groups or outside>groups[tag]:groups[tag]=outside
 ratio=hi/hj;branches=[]
 for tag,A in sorted(groups.items()):
  bi=int(tag[0]);bj=int(tag[1]);branches.append((tag,Q(bi)-bj*ratio,A+bj*pay/hj))
 bound=pay/hi;xs={Q(),bound}
 for k,(_,ak,bk) in enumerate(branches):
  for _,al,bl in branches[k+1:]:
   if ak!=al:
    x=(bl-bk)/(ak-al)
    if 0<=x<=bound:xs.add(x)
 def phi(x):return max((a*x+b for _,a,b in branches))
 m=min(phi(x) for x in xs);lo=Q();up=bound
 for _,a,b in branches:
  if a>0:up=min(up,(m-b)/a)
  elif a<0:lo=max(lo,(m-b)/a)
  else:assert b<=m
 assert lo<=up
 x=up if hi>hj else lo if hi<hj else max(lo,min(up,pay/(hi+hj)))
 y=(pay-hi*x)/hj;assert x>=0 and y>=0 and hi*x+hj*y==pay
 ds=x+y-q[i]-q[j];de=hi*(x*x-q[i]*q[i])+hj*(y*y-q[j]*q[j]);assert phi(x)==m
 return {'i':i,'j':j,'P':pay,'groups':groups,'branches':branches,'candidates':sorted(xs),'min_load':m,'interval':[lo,up],'new_i':x,'new_j':y,'deltaS':ds,'deltaE':de}

def run(case,domain):
 start=time.monotonic();cpu=time.process_time();data=json.loads((OLD/case).read_text());h=list(map(Q,data['witness_budgets']));q=list(map(Q,data['q_after']));n=len(q);U=Q(data['U']);cliques=[sum((1<<i for i in c),0) for c in data['maximal_cliques']];weights=[sum((q[i] for i in bits(C)),Q()) for C in cliques];L=max(weights);S=sum(q,Q());E=sum((h[i]*q[i]*q[i] for i in range(n)),Q());assert sum((h[i]*q[i] for i in range(n)),Q())==U and L==Q(data['L_after'])
 adj=[0]*n
 for C in cliques:
  for i in bits(C):adj[i]|=C & ~(1<<i)
 best=None;scanned=0;zeropay=0
 for i in range(n):
  for j in range(i+1,n):
   if domain=='nonedge' and (adj[i]>>j)&1:continue
   scanned+=1;v=optimize(q,h,cliques,weights,i,j)
   if v is None:zeropay+=1;continue
   key=(v['min_load'],v['deltaS'],v['deltaE'],i,j)
   if best is None or key<best[0]:best=(key,v)
 improved=bool(best and best[0][:3]<(L,Q(),Q()))
 selected=best[1] if best else None;qp=q.copy()
 if improved:qp[selected['i']]=selected['new_i'];qp[selected['j']]=selected['new_j']
 Lnew=max(sum((qp[i] for i in bits(C)),Q()) for C in cliques);Snew=sum(qp,Q());Enew=sum((h[i]*qp[i]*qp[i] for i in range(n)),Q());assert sum((h[i]*qp[i] for i in range(n)),Q())==U and (Lnew,Snew,Enew)<=(L,S,E)
 if improved:assert (Lnew,Snew,Enew)==(selected['min_load'],S+selected['deltaS'],E+selected['deltaE'])
 detail=None
 if selected:
  detail={'i':selected['i'],'j':selected['j'],'P':str(selected['P']),'class_max_outside':{k:str(v) for k,v in selected['groups'].items()},'branches':[{'class':tag,'slope':str(a),'intercept':str(b)} for tag,a,b in selected['branches']],'candidates':[str(x) for x in selected['candidates']],'min_load':str(selected['min_load']),'optimal_interval':[str(x) for x in selected['interval']],'new_i':str(selected['new_i']),'new_j':str(selected['new_j']),'deltaS':str(selected['deltaS']),'deltaE':str(selected['deltaE'])}
 out={'case':case,'domain':domain,'input_sha256':sha(OLD/case),'scanned_pairs':scanned,'zero_payment_pairs':zeropay,'best_pair':detail,'strict_merit_decrease':improved,'change_class':'load' if Lnew<L else 'sum_q' if Snew<S else 'energy' if Enew<E else 'no_change','L_before':str(L),'L_after':str(Lnew),'S_before':str(S),'S_after':str(Snew),'E_before':str(E),'E_after':str(Enew),'q_after':[str(x) for x in qp],'payment_conserved':True,'no_iteration':True,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 name='scan_'+case.removesuffix('.json')+'_'+domain+'.json';(D/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'scan','case':case,'domain':domain,'pairs':scanned,'strict':improved,'class':out['change_class'],'pid':os.getpid(),'wall':out['wall_seconds']}),flush=True);return {'file':name,'sha256':sha(D/name),'case':case,'domain':domain,'change_class':out['change_class'],'strict':improved,'pid':os.getpid(),'wall_seconds':out['wall_seconds'],'cpu_seconds':out['cpu_seconds'],'peak_rss_bytes':out['peak_rss_bytes']}

def register():
 out={'batch':'R039','code_sha256':sha(P),'proof_sha256':sha(PROOF),'input_hashes':{f:sha(OLD/f) for f in FILES},'tasks':[{'case':f,'domain':d} for f in FILES for d in ('nonedge','all')],'max_scans':6,'max_steps_per_scan':1,'workers':6,'registered_unix_time':time.time(),'method':'four exact clique-class affine lines; load plateau then S then E; no new geometry/atoms/LP or iteration'};(D/'registration.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'registration_sha256':sha(D/'registration.json'),'code_sha256':out['code_sha256'],'inputs':out['input_hashes'],'registered_unix_time':out['registered_unix_time']},indent=2))
if __name__=='__main__':
 if sys.argv[1:]==['--register']:register();sys.exit()
 assert sys.argv[1:]==['--run'];reg=json.loads((D/'registration.json').read_text());assert reg['code_sha256']==sha(P) and reg['proof_sha256']==sha(PROOF) and all(sha(OLD/f)==x for f,x in reg['input_hashes'].items());start=time.monotonic();print(json.dumps({'event':'start','parent_pid':os.getpid(),'registration_sha256':sha(D/'registration.json'),'workers':6}),flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as ex:rows=list(ex.map(run,[x['case'] for x in reg['tasks']],[x['domain'] for x in reg['tasks']]))
 out={'batch':'R039','registration_sha256':sha(D/'registration.json'),'wall_seconds':time.monotonic()-start,'results':rows};(D/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'event':'complete','wall_seconds':out['wall_seconds'],'results':rows}),flush=True)
