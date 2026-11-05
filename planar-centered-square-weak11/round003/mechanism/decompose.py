"""Exact receiver-union decomposition for one fixed periodic finite certificate.
Uses same reach and strict Dr-root construction as public check_exact2.py.
No uniqueness of area attribution is claimed: increments depend on class order.
"""
import sys,json,math,time,os,hashlib,importlib.util,argparse
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
P=Path(__file__).resolve().parents[2]; D=Path(__file__).resolve().parent
f0=P/'checks/lower/exactinputs/c2lower_1p70272.json'; f1=P/'round002/verification/candidate002.json'; core=P/'checks/lower/exact_check.py'
spec=importlib.util.spec_from_file_location('core',str(core));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def union(job):
 name,boxes=job;t=time.time();v=mod.union_area_int(boxes);return name,str(v),os.getpid(),time.time()-t

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,default=D/'generated'); args=parser.parse_args();O=args.output_dir.resolve();O.mkdir(parents=True,exist_ok=True)
 start=time.time();a=json.loads(f0.read_text());b=json.loads(f1.read_text());assert a['pts']==b['pts'][:-1] and (a['Px'],a['Py'],a['reach'])==(b['Px'],b['Py'],b['reach']); assert len(b['pts'])==51
 pts=[tuple(map(F,z)) for z in b['pts']];px,py,reach=map(F,(b['Px'],b['Py'],b['reach']));dr=2**40
 lcm=math.lcm;dx=math.lcm(*[v.denominator for v in [px,py,reach]+[z for p in pts for z in p[:2]]]);dm=math.lcm(*[p[2].denominator for p in pts]);L=lcm(dx,dr);sx=L//dx;sr=L//dr
 pxI,pyI,reachI=int(px*dx),int(py*dx),int(reach*dx);Pint=[(int(x*dx),int(y*dx),int(m*dm)) for x,y,m in pts];mass=sum(p[2] for p in Pint);nb=(int(reach/px)+2)*(int(reach/py)+2);rho_max=F(math.isqrt(nb*mass*dr**2//(4*dm))+1,dr);rmI=math.ceil(rho_max*dx)+1;ext=reachI+rmI;mx=ext//pxI+2;my=ext//pyI+2
 W=[(x+i*pxI,y+j*pyI,m,k,i,j) for i in range(-mx,mx+1) for j in range(-my,my+1) for k,(x,y,m) in enumerate(Pint) if -rmI-reachI<=x+i*pxI<=pxI+rmI+reachI and -rmI-reachI<=y+j*pyI<=pyI+rmI+reachI]; old=[z for z in W if z[3]<50];X=sorted({z[0] for z in W});Y=sorted({z[1] for z in W});oldx={z[0] for z in old};oldy={z[1] for z in old};xi={v:i for i,v in enumerate(X)};yi={v:i for i,v in enumerate(Y)};nx=len(X);ny=len(Y)
 print('PID',os.getpid(),'window old/new',len(old),len(W),'grid',nx,ny,'new coordinate globally existing?',Pint[-1][0] in oldx,Pint[-1][1] in oldy,flush=True)
 G=[[[0]*(ny+1) for _ in range(nx+1)] for _ in range(2)]
 for x,y,m,k,_,_ in W:
  G[1][xi[x]+1][yi[y]+1]+=m
  if k<50:G[0][xi[x]+1][yi[y]+1]+=m
 for g in G:
  for i in range(1,nx+1):
   acc=0
   for j in range(1,ny+1):acc+=g[i][j];g[i][j]=g[i-1][j]+acc
 def rectangle(xa,xb,ya,yb,M):
  if not M:return None
  k=math.isqrt(M*dr**2//(4*dm))-1
  if k<=0:return None
  rho=k*sr; r=(max(xb*sx-rho,0),min(xa*sx+rho,pxI*sx),max(yb*sx-rho,0),min(ya*sx+rho,pyI*sx))
  return r if r[0]<r[1] and r[2]<r[3] else None
 R0=[];R1=[];R2=[];R3=[];records=[];counts={'old_positive':0,'old_empty_new_positive':0,'new_boundary_positive':0}
 for ia,xa in enumerate(X):
  if xa< -rmI:continue
  for ib in range(ia,nx):
   xb=X[ib]
   if xb-xa>reachI or xb>pxI+rmI:break
   oldxx=xa in oldx and xb in oldx
   for ja,ya in enumerate(Y):
    if ya< -rmI:continue
    for jb in range(ja,ny):
     yb=Y[jb]
     if yb-ya>reachI or yb>pyI+rmI:break
     mass0=G[0][ib+1][jb+1]-G[0][ia][jb+1]-G[0][ib+1][ja]+G[0][ia][ja]
     mass1=G[1][ib+1][jb+1]-G[1][ia][jb+1]-G[1][ib+1][ja]+G[1][ia][ja]
     if mass1<=0:continue
     r1=rectangle(xa,xb,ya,yb,mass1)
     oldbbox=oldxx and ya in oldy and yb in oldy
     r0=rectangle(xa,xb,ya,yb,mass0) if oldbbox else None
     if oldbbox and r0 is not None:R0.append(r0)
     if r1 is None:continue
     if oldbbox and r0 is not None:cls='old_positive';R1.append(r1)
     elif oldbbox:cls='old_empty_new_positive';R2.append(r1)
     else:cls='new_boundary_positive';R3.append(r1)
     counts[cls]+=1;records.append([str(xa),str(xb),str(ya),str(yb),str(mass0),str(mass1),cls,[str(v) for v in r1],None if r0 is None else [str(v) for v in r0]])
 print('counts',counts,'old area rects',len(R0),'enumeration wall',time.time()-start,flush=True)
 jobs=[('A0',R0),('A1',R1),('A12',R1+R2),('Afull',R1+R2+R3)]
 with ProcessPoolExecutor(max_workers=4) as pool: out=list(pool.map(union,jobs))
 A={name:F(int(v),L*L) for name,v,_,_ in out};oldrec=json.loads((P/'checks/lower/exactinputs/c2lower_1p70272.check2.json').read_text());newrec=json.loads((P/'round002/verification/candidate002.check2.json').read_text());oldmass=sum(p[2] for p in pts[:-1]);newmass=sum(p[2] for p in pts)
 oldratio=F(int(oldrec['ratio_num']),int(oldrec['ratio_den']));newratio=F(int(newrec['ratio_num']),int(newrec['ratio_den']));print('exact old/new match',A['A0']==oldratio*oldmass,A['Afull']==newratio*newmass,flush=True);assert A['A0']==oldratio*oldmass and A['Afull']==newratio*newmass and A['A0']<=A['A1']<=A['A12']<=A['Afull']
 outdata={'Dx':str(dx),'Dm':str(dm),'Dr':str(dr),'L':str(L),'rho_max':str(rho_max),'new_x_in_old_periodic_grid':Pint[-1][0] in oldx,'new_y_in_old_periodic_grid':Pint[-1][1] in oldy,'counts':counts,'areas':{k:str(v) for k,v in A.items()},'increments':{'class1':str(A['A1']-A['A0']),'class2_given_class1':str(A['A12']-A['A1']),'class3_given_1_and_2':str(A['Afull']-A['A12'])},'old_new_match':[True,True],'union_workers':[{'name':n,'pid':pid,'elapsed':elapsed} for n,_,pid,elapsed in out],'records':records}
 (O/'decomposition.json').write_text(json.dumps(outdata,separators=(',',':'))); summary={k:v for k,v in outdata.items() if k!='records'}; summary['records_count']=len(records); summary['records_sha256']=hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest(); (O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('workers',outdata['union_workers']);print('areas',[(k,float(v)) for k,v in A.items()]);print('increments',[(k,float(F(v))) for k,v in outdata['increments'].items()]);print('wall',time.time()-start,flush=True)
if __name__=='__main__':main()
