"""Independent support re-summation and coordinate-compression union crosscheck."""
import json,math,time,os,importlib.util,argparse
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
S=Path(__file__).resolve().parent;P=S.parents[1];parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=S/'generated');D=parser.parse_args().output_dir.resolve();Z=json.loads((D/'decomposition.json').read_text());c=json.loads((P/'round002/verification/candidate002.json').read_text());dx=int(Z['Dx']);dm=int(Z['Dm']);dr=int(Z['Dr']);L=int(Z['L']);sx=L//dx;sr=L//dr;px,py,reach=map(F,(c['Px'],c['Py'],c['reach']));pxI,pyI,reachI=int(px*dx),int(py*dx),int(reach*dx); rmI=math.ceil(F(Z['rho_max'])*dx)+1;ext=reachI+rmI;mx=ext//pxI+2;my=ext//pyI+2
pts=[tuple(map(F,p)) for p in c['pts']];pts=[(int(x*dx),int(y*dx),int(m*dm)) for x,y,m in pts]
W=[(x+i*pxI,y+j*pyI,m,k) for i in range(-mx,mx+1) for j in range(-my,my+1) for k,(x,y,m) in enumerate(pts) if -rmI-reachI<=x+i*pxI<=pxI+rmI+reachI and -rmI-reachI<=y+j*pyI<=pyI+rmI+reachI];X0={x for x,y,m,k in W if k<50};Y0={y for x,y,m,k in W if k<50};records=Z['records']
spec=importlib.util.spec_from_file_location('lb',P/'checks/lower/lower_check.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def task(job):
 mode,arg=job; start=time.time()
 if mode=='area':
  boxes=[]
  for row in records:
   cl=row[6]
   if arg=='A0' and row[8] is not None:boxes.append(tuple(map(int,row[8])))
   if arg=='A1' and cl=='old_positive':boxes.append(tuple(map(int,row[7])))
   if arg=='A12' and cl in ('old_positive','old_empty_new_positive'):boxes.append(tuple(map(int,row[7])))
   if arg=='Afull':boxes.append(tuple(map(int,row[7])))
  value=F(mod.union_area(boxes),L**2);assert value==F(Z['areas'][arg]);return dict(mode=mode,name=arg,ok=True,pid=os.getpid(),seconds=time.time()-start,n=len(boxes))
 for row in records[arg::11]:
  xa,xb,ya,yb,m0,m1=map(int,row[:6]);M0=sum(m for x,y,m,k in W if k<50 and xa<=x<=xb and ya<=y<=yb);M1=sum(m for x,y,m,k in W if xa<=x<=xb and ya<=y<=yb);assert M0==m0 and M1==m1
  old=xa in X0 and xb in X0 and ya in Y0 and yb in Y0
  def rect(m):
   if m<=0:return None
   k=math.isqrt(m*dr**2//(4*dm))-1; r=(max(xb*sx-k*sr,0),min(xa*sx+k*sr,pxI*sx),max(yb*sx-k*sr,0),min(ya*sx+k*sr,pyI*sx))
   return r if k>0 and r[0]<r[1] and r[2]<r[3] else None
  ro=rect(M0) if old else None;rn=rect(M1);cl='old_positive' if old and ro else 'old_empty_new_positive' if old else 'new_boundary_positive';assert cl==row[6] and rn==tuple(map(int,row[7])) and ro==(tuple(map(int,row[8])) if row[8] is not None else None)
 return dict(mode=mode,name=arg,ok=True,pid=os.getpid(),seconds=time.time()-start,n=len(records[arg::11]))
if __name__=='__main__':
 start=time.time();jobs=[('area',n) for n in ('A0','A1','A12','Afull')]+[('support',i) for i in range(11)]
 with ProcessPoolExecutor(max_workers=15) as pool: res=list(pool.map(task,jobs))
 out=dict(jobs=res,wall=time.time()-start,actual_pids=len(set(x['pid'] for x in res)));(D/'independent_verify.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
