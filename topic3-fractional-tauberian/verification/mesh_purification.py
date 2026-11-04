from fractions import Fraction as F
from pathlib import Path
import random,json,time,platform
from concurrent.futures import ProcessPoolExecutor

def task(arg):
 n,qs=arg;h=F(1,n);E=[]
 for i in range(n):
  for j in range(n):
   a=i*h;c=j*h;mass=((a+h)**2-a*a)*((c+h)**2-c*c)/4;E.append((a,a+mass/h,c,c+h))
 assert sum((b-a)*(d-c) for a,b,c,d in E)==F(1,4)
 errors=[]
 for xx,yy,kk in qs:
  x=F(xx,128);y=F(yy,128);ell=F(kk,128);val=sum(max(F(0),min(b,x+ell)-max(a,x))*max(F(0),min(d,y+ell)-max(c,y)) for a,b,c,d in E)
  exact=((x+ell)**2-x*x)*((y+ell)**2-y*y)/4;error=abs(val-exact);bound=(ell+2*h)**2-max(F(0),ell-2*h)**2;assert error<=bound<=8*h*(1+2*h);errors.append(str(error))
 return {'mesh_n':n,'cells':n*n,'mass':'1/4','max_exact_error':str(max(map(F,errors))),'errors':errors,'uniform_bound':str(8*h*(1+2*h))}
if __name__=='__main__':
 rng=random.Random(20261105);qs=[]
 for _ in range(64):
  k=rng.randint(1,96);qs.append([rng.randint(0,128-k),rng.randint(0,128-k),k])
 start=time.monotonic()
 with ProcessPoolExecutor(max_workers=4) as pool:rows=list(pool.map(task,[(n,qs) for n in [4,8,16,32]]))
 out={'environment':{'python':platform.python_version()},'density':'xy on [0,1]^2','queries_numerators_over_128':qs,'rows':rows,'seed':20261105,'workers':4,'seconds':time.monotonic()-start,'scope':'256 exact rational test integrals; global approximation bound proved analytically'}
 (Path(__file__).resolve().parents[1]/'results'/'mesh_purification.json').write_text(json.dumps(out,indent=2));print([{k:v for k,v in r.items() if k!='errors'} for r in rows],out['seconds'])
