#!/usr/bin/env python3
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from flint import fmpz_mat
from fractions import Fraction
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import itertools,hashlib,json,math,time,datetime,resource,subprocess,threading,sys,flint
D=Path(__file__).resolve().parent;old=json.loads((D/'r009-manifest.json').read_text())
for f,k in [('r009-rays.json','rays.json'),('r009-inputs.json','inputs.json')]:assert hashlib.sha256((D/f).read_bytes()).hexdigest()==old['files'][k]
R=json.loads((D/'r009-rays.json').read_text());inp=json.loads((D/'r009-inputs.json').read_text());W=[r['w_primitive'] for r in R];N=inp['arrangement_normals'];B=inp['prefix_normals'];CUT=inp['cut_body_indices'];assert len(R)==46 and len(N)==36 and len(CUT)==29 and all(r['id']==i for i,r in enumerate(R))
ALL=(1<<46)-1; dots=[[sum(n[k]*w[k] for k in range(7)) for w in W] for n in N];ZERO=[sum(1<<i for i,v in enumerate(row) if v==0) for row in dots];PLUS=[sum(1<<i for i,v in enumerate(row) if v>=0) for row in dots];MINUS=[sum(1<<i for i,v in enumerate(row) if v<=0) for row in dots];Ni=[sum((i+1)*v for i,v in enumerate(w)) for w in W];qi=[sum(abs(sum(a*b for a,b in zip(n,w))) for n in B) for w in W];assert min(qi)>0
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def ids(mask):return [i for i in range(46) if mask>>i&1]
def rank(mask):return fmpz_mat([W[i] for i in ids(mask)]).rank() if mask else 0
def monitor(stop,pid):
 while not stop.is_set():
  try:
   rows=subprocess.run(['ps','-axo','pid,ppid,pcpu,rss,command'],capture_output=True,text=True).stdout.splitlines();mine=[]
   for line in rows:
    bits=line.split(None,4)
    if len(bits)==5 and (bits[0]==str(pid) or bits[1]==str(pid)):mine.append(line)
   sw=subprocess.run(['sysctl','vm.swapusage'],capture_output=True,text=True).stdout.strip()
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),parent_pid=pid,processes=mine,swap=sw))+'\n')
  except Exception as e:
   with (D/'resources.jsonl').open('a') as f:f.write(json.dumps(dict(utc=now(),error=repr(e)))+'\n')
  stop.wait(5)
def cell_task(cell):
 t=time.perf_counter();c=time.process_time();key=cell['id'];mask=int(cell['mask']);sg=cell['signs'];rankcache={};nodes={};pullcache={}
 try:
  def rk(m):
   if m not in rankcache:rankcache[m]=rank(m)
   return rankcache[m]
  def pull(m):
   if m in pullcache:return pullcache[m]
   V=ids(m);r=rk(m);assert r>=1
   if r==1:
    assert len(V)==1;nodes[m]=dict(rank=1,anchor=V[0],facets=[]);pullcache[m]=[(V[0],)];return pullcache[m]
   facets=set()
   for zero in ZERO:
    f=m&zero
    if f and f!=m and rk(f)==r-1:facets.add(f)
   facets=sorted(facets);anchor=V[0];assert facets and any(not (f>>anchor&1) for f in facets)
   tris=[]
   for f in facets:
    if not (f>>anchor&1):
     for tri in pull(f):tris.append(tuple(sorted((anchor,)+tri)))
   assert len(tris)==len(set(tris)) and tris;tris=sorted(tris);nodes[m]=dict(rank=r,anchor=anchor,facets=[str(f) for f in facets]);pullcache[m]=tris;return tris
  tr=pull(mask);assert rk(mask)==7
  witness=[sum(W[i][j] for i in ids(mask)) for j in range(7)];assert min(witness)>0 and all(sg[j]*sum(N[j+7][k]*witness[k] for k in range(7))>0 for j in range(29))
  noncut=[B[i] for i in range(64) if i not in CUT];ell=[sum(n[k] for n in noncut)+sum(sg[j]*N[j+7][k] for j in range(29)) for k in range(7)]
  assert all(sum(ell[k]*W[i][k] for k in range(7))==qi[i] for i in ids(mask))
  simplices=[];vol=Fraction()
  for tt in tr:
   assert len(tt)==7 and len(set(tt))==7 and all(mask>>i&1 for i in tt)
   mat=fmpz_mat([[W[j][i] for j in tt] for i in range(7)]);det=int(mat.det());assert det
   nn=[Ni[i] for i in tt];qq=[qi[i] for i in tt];v=Fraction(math.factorial(7)*abs(det),math.prod(nn));alpha=Fraction(abs(det),math.prod(qq));vol+=v
   assert all(sg[j]*dots[j+7][i]>=0 for j in range(29) for i in tt)
   simplices.append(dict(vertices=list(tt),detW=det,N=nn,q=qq,alpha=str(alpha),volume=str(v)))
  return dict(id=key,status='completed',mask=str(mask),signs=sg,vertices=ids(mask),witness=witness,ell=ell,facets=nodes[mask]['facets'],simplices=simplices,volume=str(vol),face_dag={str(m):data for m,data in nodes.items()},wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 except BaseException as e:return dict(id=key,status='failed',mask=str(mask),error=repr(e),wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
def main():
 assert not (D/'levels.jsonl').exists();started=now();t=time.perf_counter();basecpu=time.process_time();stop=threading.Event();th=threading.Thread(target=monitor,args=(stop,os.getpid()),daemon=True);th.start();results=[];stagecpu=0.;faces={};seen=set();totalvol=Fraction();ntri=0;cells=[];status='failed'
 try:
  assert rank(ALL)==7;current={ALL:[]}
  for level in range(29):
   nxt={};parents=[]
   for mask,sg in sorted(current.items()):
    rec=dict(mask=str(mask),signs=sg,children=[])
    for sign,F in [(1,PLUS[level+7]),(-1,MINUS[level+7])]:
     child=mask&F;rr=rank(child);r=dict(sign=sign,mask=str(child),rank=rr,kept=(rr==7));rec['children'].append(r)
     if rr==7:
      path=sg+[sign]
      if child in nxt and nxt[child]!=path:raise ArithmeticError('conflicting signprefix same full-rank bitset')
      nxt[child]=path
    assert any(x['kept'] for x in rec['children']);parents.append(rec)
   assert nxt;row=dict(level=level+1,pilot=level<4,parent_count=len(current),child_count=len(nxt),parents=parents)
   with (D/'levels.jsonl').open('a') as f:f.write(json.dumps(row)+'\n');f.flush();os.fsync(f.fileno())
   current=nxt
  stagecpu=time.process_time()-basecpu;cells=[dict(id=i,mask=str(m),signs=sg) for i,(m,sg) in enumerate(sorted(current.items(),key=lambda x:tuple(x[1])))]
  def record(r):
   nonlocal totalvol,ntri
   results.append({k:v for k,v in r.items() if k not in ('simplices','face_dag')})
   if r['status']=='completed':
    for k,v in r['face_dag'].items():
     if k in faces and faces[k]!=v:raise ArithmeticError('conflicting global face DAG')
     faces[k]=v
    for s in r['simplices']:
     key=tuple(s['vertices'])
     if key in seen:raise ArithmeticError('repeated top simplex across cells')
     seen.add(key);ntri+=1
    totalvol+=Fraction(r['volume'])
    cellrec={k:v for k,v in r.items() if k!='face_dag'}
    with (D/'cells.jsonl').open('a') as f:f.write(json.dumps(cellrec)+'\n');f.flush()
   with (D/'terminal.jsonl').open('a') as f:f.write(json.dumps(results[-1])+'\n');f.flush();os.fsync(f.fileno())
  for cell in cells[:4]:
   r=cell_task(cell);record(r)
   if r['status']!='completed':raise RuntimeError('pulling pilot failed')
  failed=False
  with ProcessPoolExecutor(max_workers=15) as ex:
   fs={ex.submit(cell_task,cell):cell['id'] for cell in cells[4:]}
   for f in as_completed(fs):
    try:r=f.result()
    except BaseException as e:r=dict(id=fs[f],status='failed',error=repr(e))
    record(r)
    if r['status']!='completed':failed=True
  if failed:raise RuntimeError('one or more cell tasks failed')
  assert len(results)==len(cells) and totalvol==1
  (D/'face-dag.json').write_text(json.dumps(faces,sort_keys=True)+'\n');status='completed'
 except BaseException as e:print('ERROR',repr(e),file=sys.stderr)
 finally:stop.set();th.join()
 out=dict(status=status,start_utc=started,finish_utc=now(),wall_seconds=time.perf_counter()-t,stageA_parent_cpu_seconds=stagecpu,sum_cell_cpu_seconds=sum(x.get('cpu_seconds',0) for x in results),peak_cell_process_rss_bytes=max((x.get('peak_rss_bytes',0) for x in results),default=0),levels_completed=29 if cells else (len((D/'levels.jsonl').read_text().splitlines()) if (D/'levels.jsonl').exists() else 0),cells=len(cells),cell_terminals=len(results),unique_faces=len(faces),top_simplices=ntri,volume=str(totalvol),pilot_levels=4,pilot_cells=min(4,len(cells)),workers=15,python=sys.version,flint=flint.__version__)
 (D/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));sys.exit(0 if status=='completed' else 1)
if __name__=='__main__':main()
