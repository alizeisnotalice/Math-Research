#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from flint import arb,acb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,hashlib,math,itertools,time,datetime,resource,sys,flint,subprocess,threading
D=Path(__file__).resolve().parent;ctx.prec=256;T=json.loads((D/'testcases.json').read_text());rays=json.loads((D/'r009-rays.json').read_text());inp=json.loads((D/'r009-inputs.json').read_text());cells=[json.loads(s) for s in (D/'cells.jsonl').read_text().splitlines()];old=json.loads((D/'r010-manifest.json').read_text())
for f in ('cells.jsonl','r009-rays.json','r009-inputs.json'):assert hashlib.sha256((D/f).read_bytes()).hexdigest()==old['files'][f]
W=[r['w_primitive'] for r in rays];XI=[r['xi_primitive'] for r in rays];BODY=inp['prefix_normals'];Q=[sum(abs(sum(a*b for a,b in zip(n,w))) for n in BODY) for w in W];L=math.lcm(*Q);ZETA=[[L//Q[i]*z for z in XI[i]] for i in range(46)];S=[];seen=set()
for c in cells:
 for s in c['simplices']:
  key=tuple(s['vertices']);assert key not in seen;seen.add(key);assert s['q']==[Q[i] for i in key];a=Fraction(s['alpha']);assert a==Fraction(abs(s['detW']),math.prod(s['q']));S.append((key,a))
assert len(cells)==135 and len(S)==212 and len(Q)==46
ASUM=sum((a for _,a in S),Fraction());assert 0<5040*ASUM<=Fraction(1,512)
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def A(v):
 if isinstance(v,Fraction):return arb(fmpq(v.numerator,v.denominator))
 return arb(v)
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(35)}
def width(x):return (x.upper()-x.lower()).upper()
def dd_real(nodes,beta,cache,counts):
 key=tuple(nodes)
 if key in cache:counts['dd_cache_hits']+=1;return cache[key]
 counts['unique_subdd']+=1;k=len(key)-1
 if k==0:
  counts['log_evaluations']+=1;assert key[0]!=0;out=A(abs(key[0])).log()
 elif key[0]==key[-1]:
  counts['hermite_evaluations']+=1;assert key[0]!=0;out=A(Fraction((-1)**(k-1),k*key[0]**k))
 else:out=(dd_real(key[1:],beta,cache,counts)-dd_real(key[:-1],beta,cache,counts))/(key[-1]-key[0])
 cache[key]=out;return out
def dd_complex(nodes,beta):
 cache={}
 def go(key):
  if key in cache:return cache[key]
  k=len(key)-1
  if k==0:out=acb(A(beta),-key[0]).log()
  elif key[0]==key[-1]:out=-(acb(0,1)**k)/(k*acb(A(beta),-key[0])**k)
  else:out=(go(key[1:])-go(key[:-1]))/(key[-1]-key[0])
  cache[key]=out;return out
 return go(tuple(nodes))
def counter():return dict(top_requests=0,unique_top=0,unique_subdd=0,dd_cache_hits=0,log_evaluations=0,hermite_evaluations=0,zero_ray_knots=0,repeated_top_nodes=0)
def phase(n,dx,eps,phase_id):
 assert eps==0;B=dx*L;beta=0;ar=[sum(n[j]*ZETA[i][j] for j in range(7)) for i in range(46)];cnt=counter();cnt['zero_ray_knots']=sum(a==0 for a in ar);assert cnt['zero_ray_knots']==0, f'zero exact knot in phase {phase_id}';top={}
 for ix,weight in S:
  key=tuple(sorted(ar[i] for i in ix));top[key]=top.get(key,Fraction())+weight;cnt['top_requests']+=1;cnt['repeated_top_nodes']+=7-len(set(key))
 cnt['unique_top']=len(top);cache={};val=arb(0)
 for key,weight in top.items():val+=A(weight*(B**6))*dd_real(key,beta,cache,cnt)
 return val,dict(phase_id=phase_id,counts=cnt,raw_phase=ep(val))
def phys(nums):
 d=math.lcm(*(x.denominator for x in nums));return [int(x*d) for x in nums],d
def symmetric(case,eps,precision):
 ctx.prec=precision;t=time.perf_counter();cpu=time.process_time();st=46080 if case=='e1' else 5040
 phases=[tuple(int(i==k) for i in range(7)) for k in range(7)] if case=='e1' else [(1,*x) for x in itertools.product((-1,1),repeat=6)]
 v=arb(0);records=[]
 for i,signs in enumerate(phases):
  z,r=phase(signs,1,eps,i);v+=z;records.append(r)
 norm=2*st/(2*arb.pi())**7;out=v*norm
 return dict(case=case,epsilon=str(eps),precision=precision,phase_count=len(phases),stabilizer=st,interval=ep(out),width=ep(width(out)),positive=bool(out>0),width_pass=bool(width(out)<=A(Fraction(1,10**12))),phase_records=records,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
def generic_chunk(chunk,precision):
 ctx.prec=precision;t=time.perf_counter();cpu=time.process_time();hexes=T['generic_case']['x_floathex'];vals=[Fraction.from_float(float.fromhex(v)) for v in hexes];assert all(vals) and len(set(abs(v) for v in vals))==7
 nums,dx=phys(vals);perm=next(itertools.islice(itertools.permutations(range(7)),chunk,chunk+1));total=arb(0);records=[]
 for signid,tail in enumerate(itertools.product((-1,1),repeat=6)):
  signs=(1,*tail);n=[signs[j]*nums[perm[j]] for j in range(7)];z,r=phase(n,dx,Fraction(0),64*chunk+signid);total+=z;records.append(r)
 out=2*total/(2*arb.pi())**7
 return dict(chunk=chunk,start_phase=64*chunk,end_phase=64*(chunk+1),precision=precision,interval=ep(out),phase_records=records,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
def ball(d):
 def dy(x):return arb(int(x[0]))*arb(2)**x[1]
 return dy(d['lower']).union(dy(d['upper']))
def scalar_cases():
 out=[];eps=Fraction(1,64)
 for c in T['scalar_DD_cases']:
  knots=[Fraction(x) for x in c['knots']];B=math.lcm(*(x.denominator for x in knots));a=tuple(sorted(int(v*B) for v in knots));beta=eps*B;cnt=counter();real=dd_real(a,beta,{},cnt)*B**6;complexval=dd_complex(a,beta)*B**6;assert real.overlaps(complexval.real)
  reference='real-complex overlap'
  if c['name']=='all_zero':ref=acb(A(Fraction(1,6)/eps**6));assert complexval.overlaps(ref);reference='exact 1/(6 eps^6)'
  if c['name']=='all_one':ref=1/(6*acb(A(eps),-1)**6);assert complexval.overlaps(ref);reference='complex 1/(6(eps-i)^6)'
  if c['name']=='near_distinct':
   M=32;delta=max(abs(v) for v in knots);theta=delta/eps;assert theta<1 and 1-theta*Fraction(M+7,M+2)>0;h=[Fraction(0)]*(M+1);h[0]=1
   for x in knots:
    for k in range(1,M+1):h[k]+=x*h[k-1]
   val=acb(0)
   for j in range(M+1):
    k=6+j;val+=-(acb(0,1)**k)*A(h[j]/(k*eps**k))
   tail=Fraction(1,6)/eps**6*math.comb(M+6,5)*theta**(M+1)/(1-theta*Fraction(M+7,M+2));err=A(tail);ref=acb(val.real+arb(0,err),val.imag+arb(0,err));assert complexval.overlaps(ref);reference='32 extra terms plus exact rational tail'
  out.append(dict(name=c['name'],knots=c['knots'],common_denominator=B,exact_distinct=len(set(a)),real=ep(real),complex_real=ep(complexval.real),complex_imag=ep(complexval.imag),reference=reference,counts=cnt))
 return out
def origin_cases():
 out=[]
 for text in T['origin_epsilons']:
  eps=Fraction(text);cnt=counter();B=L;val=dd_real((0,)*7,eps*B,{},cnt)*B**6;assert val.contains(A(Fraction(1,6)/eps**6));direct=840*A(ASUM)/(arb.pi()**7*A(eps)**6);bydd=645120*A(ASUM)*val/(2*arb.pi())**7;assert direct.overlaps(bydd) and direct>0
  out.append(dict(epsilon=text,interval=ep(direct),dd_identity_pass=True,opposite_factor=1,counts=cnt))
 return out
def monitor(stop,pid):
 while not stop.is_set():
  try:
   rows=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout.splitlines();mine=[]
   for line in rows:
    z=line.split(None,4)
    if len(z)==5 and (z[0]==str(pid) or z[1]==str(pid)):mine.append(line)
   sw=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),parent_pid=pid,processes=mine,swap=sw))+'\n')
  except Exception as e:
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(3)
def main():
 start=now();t=time.perf_counter();cpu=time.process_time();stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);thread.start();terms=[];extra=[];status='failed';error=None
 try:
  scalar=scalar_cases();origin=origin_cases();spec=[('e1',Fraction(1,64)),('e1',Fraction(1,256)),('ones',Fraction(1,64)),('ones',Fraction(1,256))]
  with ProcessPoolExecutor(max_workers=15) as ex:
   jobs={ex.submit(symmetric,c,e,256):('sym',c,str(e)) for c,e in spec};jobs.update({ex.submit(generic_chunk,k,256):('generic',k) for k in range(4)})
   for f in as_completed(jobs):
    try:r=f.result();terms.append(r);label=jobs[f]
    except BaseException as e:raise RuntimeError(f'task {jobs[f]} failed: {e!r}')
    with (D/'terminal.jsonl').open('a') as o:o.write(json.dumps(dict(label=label,status='completed',result=r))+'\n');o.flush();os.fsync(o.fileno())
  sym=[x for x in terms if 'case' in x];gens=sorted([x for x in terms if 'chunk' in x],key=lambda x:x['chunk']);assert len(sym)==4 and len(gens)==4
  for x in list(sym):
   if not x['width_pass']:
    z=symmetric(x['case'],Fraction(x['epsilon']),512);extra.append(dict(reason='width_only',initial=x,promoted=z));assert z['width_pass'];sym[sym.index(x)]=z
  for case in ('e1','ones'):
   x=next(z for z in sym if z['case']==case and z['epsilon']=='1/64');y=next(z for z in sym if z['case']==case and z['epsilon']=='1/256');assert x['positive'] and y['positive'] and ball(y['interval'])>ball(x['interval']);assert ball(origin[0]['interval'])>ball(x['interval']) and ball(origin[1]['interval'])>ball(y['interval'])
  def generic_sum(v):return sum((ball(x['interval']) for x in v),arb(0))
  gs=generic_sum(gens);target=Fraction(1,10**12)*Fraction(256,322560)
  if not width(gs)<=A(target):
   new=[generic_chunk(k,512) for k in range(4)];extra.append(dict(reason='generic_width_only',initial=gens,promoted=new));gens=new;gs=generic_sum(gens);assert width(gs)<=A(target)
  assert [r['phase_id'] for x in gens for r in x['phase_records']]==list(range(256));status='completed'
 except BaseException as e:error=repr(e)
 finally:stop.set();thread.join()
 out=dict(status=status,error=error,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,parent_cpu_seconds=time.process_time()-cpu,worker_cpu_seconds=sum(x.get('cpu_seconds',0) for x in terms),promotion_cpu_seconds=sum(y.get('cpu_seconds',0) for z in extra for y in ([z['promoted']] if isinstance(z['promoted'],dict) else z['promoted'])),peak_task_rss_bytes=max((x.get('peak_rss_bytes',0) for x in terms),default=0),python=sys.version,flint=flint.__version__,A_sum=str(ASUM),q_lcm=L,simplex_count=len(S),rays=len(Q),scalar_cases=scalar if 'scalar' in locals() else [],origin=origin if 'origin' in locals() else [],symmetry=sym if 'sym' in locals() else [],generic_chunks=gens if 'gens' in locals() else [],generic_partial=ep(gs) if 'gs' in locals() else None,generic_planning_width_pass=bool(width(gs)<=A(target)) if 'gs' in locals() else None,generic_partial_is_not_H=True,promotions=extra,worker_tasks=8,MC_samples=0)
 (D/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ('scalar_cases','origin','symmetry','generic_chunks','promotions')},indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
