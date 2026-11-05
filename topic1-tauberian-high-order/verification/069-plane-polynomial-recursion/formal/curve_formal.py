"""Exact Fraction bivariate (t,delta) polynomial curves modulo delta^(R+1)."""
from fractions import Fraction as F
from math import comb
from concurrent.futures import ProcessPoolExecutor
import argparse,json,os,sys,time,resource
sys.set_int_max_str_digits(0)
def add(*ps):
 r={}
 for p in ps:
  for k,v in p.items():r[k]=r.get(k,F(0))+v
 return {k:v for k,v in r.items() if v}
def sc(p,c):return {k:v*c for k,v in p.items() if v*c}
def mul(p,q,R):
 out={}
 for (i,j),x in p.items():
  for (k,l),y in q.items():
   if j+l<=R:out[(i+k,j+l)]=out.get((i+k,j+l),F(0))+x*y
 return {k:v for k,v in out.items() if v}
def diff(p):return {(i-1,j):i*x for (i,j),x in p.items() if i}
def integ(p):return {(i+1,j):x/F(i+1) for (i,j),x in p.items()}
def eval1(p):
 out={}
 for (i,j),v in p.items():out[(0,j)]=out.get((0,j),F(0))+v
 return {k:v for k,v in out.items() if v}
def rev(p):
 out={}
 for (i,j),x in p.items():
  for k in range(i+1):out[(k,j)]=out.get((k,j),F(0))+x*comb(i,k)*(-1)**k
 return {k:v for k,v in out.items() if v}
def const(x):return {(0,0):F(x)} if x else {}
def shift_delta(p,R):return {(i,j+1):x for (i,j),x in p.items() if j<R}
def summary(p,R):return [str(p.get((0,j),0)) for j in range(R+1)]
def moments(piece,R):
 r,h=piece['r'],piece['h'];z=add(mul(add(const(1),r),h,R),sc(mul(h,h,R),-1))
 return sc(eval1(integ(mul(z,diff(r),R))),2)
def area(piece,R):return sc(eval1(integ(mul(piece['h'],diff(piece['r']),R))),2)
def child(job):
 idx,p,total,R=job;start=time.monotonic();cpu=time.process_time()
 r,h,P=(p[x] for x in ('r','h','P'));q=h
 M=add(sc(mul(q,q,R),F(1,2)),mul(q,add(const(1),sc(r,-1)),R),P)
 low={'r':q,'h':shift_delta(add(q,sc(M,-1)),R)}
 q=add(const(1),sc(r,-1),h);M=add(sc(total,2),sc(P,-1),sc(mul(h,h,R),F(1,2)))
 high={'r':rev(q),'h':rev(shift_delta(add(q,sc(M,-1)),R))}
 info={'index':idx,'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'maxrss_native':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 return idx,low,high,info
def next_layer(pieces,R,pool):
 total=eval1(pieces[-1]['P']);jobs=[(i,p,total,R) for i,p in enumerate(pieces)]
 data=list(pool.map(child,jobs)) if pool else [child(j) for j in jobs]
 data.sort();new=[x[1] for x in data]+[x[2] for x in reversed(data)]
 prefix={};joins=[]
 for i,p in enumerate(new):
  if i:joins.append(eval1(new[i-1]['r'])=={k:v for k,v in p['r'].items() if k[0]==0} and eval1(new[i-1]['h'])=={k:v for k,v in p['h'].items() if k[0]==0})
  p['P']=add(integ(mul(p['h'],diff(p['r']),R)),prefix);prefix=eval1(p['P'])
 assert all(joins)
 return new,[x[3] for x in data]
def main():
 p=argparse.ArgumentParser();p.add_argument('--order',type=int,required=True);p.add_argument('--depth',type=int,required=True);p.add_argument('--workers',type=int,default=1);args=p.parse_args();R=args.order
 start=time.monotonic();mastercpu=time.process_time();pieces=[{'r':{(1,0):F(1)},'h':{(1,1):F(1)},'P':{(2,1):F(1,2)}}];rows=[];workerdata=[]
 pool=ProcessPoolExecutor(max_workers=args.workers) if args.workers>1 else None
 try:
  for j in range(args.depth+1):
   t=time.monotonic();m=add(*(moments(p,R) for p in pieces));ar=add(*(area(p,R) for p in pieces));
   rows.append({'layer':j,'pieces':len(pieces),'moment_I_coefficients':summary(m,R),'area_coefficients':summary(ar,R),'max_t_degree':max(max((i for i,k in p[x]),default=0) for p in pieces for x in ('r','h','P')),'terms':sum(len(p[x]) for p in pieces for x in ('r','h','P')),'seconds':time.monotonic()-t})
   print('layer',j,'pieces',len(pieces),'terms',rows[-1]['terms'],file=sys.stderr,flush=True)
   if j<args.depth:
    pieces,info=next_layer(pieces,R,pool);workerdata.extend(info)
 finally:
  if pool:pool.shutdown()
 result={'order':R,'depth':args.depth,'requested_workers':args.workers,'actual_worker_pids':sorted(set(x['pid'] for x in workerdata)),'elapsed_seconds':time.monotonic()-start,'master_cpu_seconds':time.process_time()-mastercpu,'worker_cpu_seconds_sum':sum(x['cpu_seconds'] for x in workerdata),'maxrss_native_master':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'maxrss_native_worker_max':max((x['maxrss_native'] for x in workerdata),default=0),'rows':rows,'worker_data':workerdata}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
