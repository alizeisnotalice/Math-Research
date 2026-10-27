#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from fractions import Fraction
from pathlib import Path
from flint import arb,fmpq,ctx
from concurrent.futures import ProcessPoolExecutor,as_completed
import math,json,hashlib,time,datetime,resource,sys,flint,subprocess,threading
ctx.prec=512;D=Path(__file__).resolve().parent;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();m=json.loads((D/'r010-manifest.json').read_text())
for f in ('cells.jsonl','r009-rays.json','r009-inputs.json'):assert sha(D/f)==m['files'][f]
assert sha(D/'r016-result.json')==json.loads((D/'r016-manifest.json').read_text())['files']['result.json']
R=json.loads((D/'r009-rays.json').read_text());cells=[json.loads(x) for x in (D/'cells.jsonl').read_text().splitlines()];simp=[s for c in cells for s in c['simplices']];assert len(simp)==212 and len({tuple(s['vertices']) for s in simp})==212
Q={i:q for s in simp for i,q in zip(s['vertices'],s['q'])};assert len(Q)==46;L=math.lcm(*Q.values());XI=[r['xi_primitive'] for r in R];ZETA=[[L//Q[i]*v for v in XI[i]] for i in range(46)];FACS=[math.factorial(i) for i in range(7)];now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def ep(x):
 l=x.lower().man_exp();u=x.upper().man_exp();return {'lower':[str(l[0]),int(l[1])],'upper':[str(u[0]),int(u[1])],'display':x.str(45)}
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Fraction) else arb(x)
def dy(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def task(i):
 t=time.perf_counter();cpu=time.process_time();pairs=0;maxsupport=0;maxbits=0
 try:
  s=simp[i];idx=s['vertices'];zero=(0,)*7
  def mul(a,b):
   nonlocal pairs,maxsupport,maxbits
   pairs+=len(a)*len(b);c={}
   for x,v in a.items():
    for y,w in b.items():
     z=tuple(x[k]+y[k] for k in range(7));c[z]=c.get(z,0)+v*w
   maxsupport=max(maxsupport,len(c));maxbits=max(maxbits,max((abs(x).bit_length() for x in c.values()),default=0));return c
  def add(dest,src,scale=1):
   for k,v in src.items():dest[k]=dest.get(k,0)+scale*v
  def contract(a):return sum(v*math.prod(FACS[j] for j in ex) for ex,v in a.items())
  P2={};P4={};P6={}
  for coord in range(7):
   lin={tuple(int(k==j) for k in range(7)):ZETA[idx[j]][coord] for j in range(7) if ZETA[idx[j]][coord]}
   sq=mul(lin,lin);four=mul(sq,sq);six=mul(four,sq);add(P2,sq);add(P4,four);add(P6,six)
  P22=mul(P2,P2);P23=mul(P22,P2);P24=mul(P2,P4)
  C2=contract(P2);C4=contract(P4);C22=contract(P22);C6=contract(P6);C24=contract(P24);C23=contract(P23)
  v=[C2,C4,C22-C4,C6,C24-C6,C23-3*C24+2*C6];assert all(x>0 for x in v) and pairs<=59535 and maxsupport<=924
  return dict(id=i,status='completed',vertices=idx,alpha=s['alpha'],contractions=[str(x) for x in v],multiplication_pairs=pairs,max_support=maxsupport,max_coefficient_bits=maxbits,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except BaseException as e:return dict(id=i,status='failed',error=repr(e),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
def monitor(stop,pid):
 while not stop.is_set():
  try:
   rows=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout.splitlines();mine=[]
   for line in rows:
    b=line.split(None,4)
    if len(b)==5 and (b[0]==str(pid) or b[1]==str(pid)):mine.append(line)
   sw=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),parent_pid=pid,processes=mine,swap=sw))+'\n')
  except Exception as e:
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(5)
def main():
 start=now();t=time.perf_counter();cpu=time.process_time();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();rows=[];err=None;status='failed'
 try:
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={ex.submit(task,i):i for i in range(212)}
   for f in as_completed(fs):
    try:r=f.result()
    except BaseException as e:r=dict(id=fs[f],status='failed',error=repr(e))
    rows.append(r)
    with (D/'terminal.jsonl').open('a') as o:o.write(json.dumps(r)+'\n');o.flush();os.fsync(o.fileno())
  assert len(rows)==212 and {x['id'] for x in rows}==set(range(212)) and all(x['status']=='completed' for x in rows);rows.sort(key=lambda x:x['id'])
  S=sum((Fraction(x['alpha']) for x in rows),Fraction());R0=5040*S;assert S==Fraction(907157967345860903220972676779187,28127466839029475035590994122592223232000000000)
  factors=[720,720,120,720,120,24];powers=[2,4,4,6,6,6];names=['R1','R_part2','R11','R3','R21','R111'];Rs={'R0':R0}
  for k in range(6):Rs[names[k]]=Fraction(factors[k],L**powers[k])*sum((Fraction(x['alpha'])*int(x['contractions'][k]) for x in rows),Fraction())
  assert len(Rs)==7 and all(v>0 for v in Rs.values());rho={k:v/R0 for k,v in Rs.items() if k!='R0'};p0=A(R0)/arb.pi()**7;b=Fraction(65*67*69,8);a=Fraction(71,2)
  def val(x):
   s=[v*v for v in x];p1=sum(s);p2=sum(v*v for v in s);p3=sum(v**3 for v in s);e2=(p1*p1-p2)/2;e3=(p1**3-3*p1*p2+2*p3)/6
   J2=rho['R1']*p1;J4=rho['R_part2']*p2+6*rho['R11']*e2;J6=rho['R3']*p3+15*rho['R21']*(p1*p2-p3)+90*rho['R111']*e3
   m=J2/(2*a);mu2=J4/(12*a*(a+1));mu3=J6/(120*a*(a+1)*(a+2));V=mu2-m*m;W=mu3-2*m*mu2+m**3;H=m*mu3-mu2**2;assert m>0 and V>0 and W>0 and H>0
   g=(W**3+2*W*V**3+m*V**4)/H**3;return p0*A(g/b),dict(m=str(m),V=str(V),W=str(W),H=str(H),gauss_inverse3=str(g))
  toy_m=Fraction(3,2);toy2=Fraction(5,2);toy3=Fraction(9,2);tv=toy2-toy_m**2;tw=toy3-2*toy_m*toy2+toy_m**3;thh=toy_m*toy3-toy2**2;assert (tw**3+2*tw*tv**3+toy_m*tv**4)/thh**3==Fraction(9,16)
  TC=json.loads((D/'testcases.json').read_text());xx=[Fraction.from_float(float.fromhex(x)) for x in TC['generic_case']['x_floathex']];prior=json.loads((D/'r016-result.json').read_text());target=dy(prior['enclosures']['g_required_B_RB']['lower']);points=[('corner',[Fraction(1)]*7),('R015_index2710',xx),('e1',[Fraction(1)]+[Fraction(0)]*6)];gates=[]
  for name,x in points:
   y,meta=val(x);e=ep(y);upper=dy(e['upper']);gates.append(dict(name=name,point=[str(v) for v in x],interval=e,moments=meta,upper_below_target=upper<target))
  Rjson={k:str(v) for k,v in Rs.items()};out=dict(status='completed',R=Rjson,rho={k:str(v) for k,v in rho.items()},p0=ep(p0),L_common=L,gates=gates,corner_gate_failed=gates[0]['upper_below_target'],target_lower=str(target),source_target='R016 g_required_B_RB lower endpoint',global_cover_started=False,scalar_outputs=10,exact_constant_outputs=7,fixed_point_outputs=3,toy_law_pass=True,multiplication_pairs=sum(x['multiplication_pairs'] for x in rows),max_support=max(x['max_support'] for x in rows),max_coefficient_bits=max(x['max_coefficient_bits'] for x in rows))
  (D/'result.json').write_text(json.dumps(out,indent=2)+'\n');status='completed'
 except BaseException as e:err=repr(e)
 finally:stop.set();th.join()
 receipt=dict(status=status,error=err,start_utc=start,finish_utc=now(),wall_seconds=time.perf_counter()-t,parent_cpu_seconds=time.process_time()-cpu,sum_task_cpu_seconds=sum(x.get('cpu_seconds',0) for x in rows),peak_worker_rss_bytes=max((x.get('peak_rss_bytes',0) for x in rows),default=0),terminal_tasks=len(rows),processes=15,python=sys.version,flint=flint.__version__);(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
