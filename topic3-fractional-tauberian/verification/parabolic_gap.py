"""Exact rational union regressions and full-interval patch certificates."""
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
import random,json,time,platform
import mpmath as mp
from witness_capacity import atoms
POLYS=[list(map(F,["495/4096","2295808/262144","1016832/262144","986047/262144"])),
       list(map(F,["64512/262144","-2202624/262144","-824512/262144","-982081/262144"])),
       list(map(F,["490496/262144","-1505280/262144","1014784/262144","-229441/262144","949378/262144","-982081/262144"])),
       list(map(F,["1/2","21951/4096","1856/4096","15872/4096","15376/4096"]))]
def params(u):
 t=u*u;a=1-t;r=t/2+t*t/2+F(31,32)*t**3;q=a+F(3,16)*t*t;nu=t**3/1024;eta=u**3/128
 return t,a,r,q,nu,eta

def run(job):
 u,ns,ss=job;t,a,r,q,nu,eta=params(u)
 pure=[(r+ss[0]*eta,r-ss[0]*eta),(-r+ss[1]*eta,r+ss[1]*eta),(-r-ss[2]*eta,-r+ss[2]*eta),(r-ss[3]*eta,-r-ss[3]*eta)]
 normals=[(1,1),(-1,1),(-1,-1),(1,-1)]
 centers=[(x+sg[0]*n*nu,y+sg[1]*n*nu) for (x,y),sg,n in zip(pure,normals,ns)]
 def area(cs):return sum(atoms([(x-a/2,x+a/2,y-a/2,y+a/2) for x,y in cs]).values())
 pure_area=area(pure);full=area(centers);base=(a+2*r)**2
 assert a>=2*r+4*eta
 assert pure_area>=base-4*eta**2
 assert full>=base-4*eta**2-16*a*nu
 return (pure_area-base+4*eta**2)/t**3,(full-base+4*eta**2+16*a*nu)/t**3

if __name__=='__main__':
 start=time.monotonic();rng=random.Random(20261105);jobs=[]
 for u in [F(1,10),F(1,20),F(1,100),F(1,1000)]:
  for ns in product([-1,1],repeat=4):
   for ss in product([-1,1],repeat=4):jobs.append((u,ns,ss))
 for _ in range(2000):
  jobs.append((F(rng.randrange(1,1001),10000),[F(rng.randrange(-32,33),32) for _ in range(4)],[F(rng.randrange(-32,33),32) for _ in range(4)]))
 with ProcessPoolExecutor(max_workers=8) as pool: rows=list(pool.map(run,jobs,chunksize=32))
 lower=[cs[0]+sum(min(c,F(0))/100**j for j,c in enumerate(cs[1:],1)) for cs in POLYS];assert min(lower)>0
 mp.mp.dps=80;cap=[];moment=[]
 for ts in ['.01','.001','.0001','.000001']:
  t=mp.mpf(ts);a=1-t;c=2*a*a-1
  for k in [mp.mpf('0.001'),mp.mpf('.1'),mp.mpf('1')]:
   d=k*t**3;L=mp.sqrt(2*mp.sqrt(2*c)+2*d)
   normalized=mp.quad(lambda z:mp.sqrt(2)-z*z/(mp.sqrt(c+d*z*z)+mp.sqrt(c)),[0,L]);assert normalized>mp.mpf('0.125')
   cap.append({'t':ts,'d_over_t3':str(k),'cap_over_d32':mp.nstr(normalized,65)})
 for pv in [256,1024,4096,16384]:
  p=mp.mpf(pv);v=mp.quad(lambda u:mp.exp((p-1)*mp.log1p(-u/p))*u**mp.mpf('4.5'),[0,p/100])/2**24
  moment.append({'p':pv,'p_4p5_scaled_correction':mp.nstr(v,65)})
 out={'python':platform.python_version(),'mpmath':mp.__version__,'workers':8,'seed':20261105,'union_cases':len(rows),'endpoint_cases':1024,'random_rational_cases':2000,'min_pure_normalized_margin':str(min(r[0] for r in rows)),'min_full_normalized_margin':str(min(r[1] for r in rows)),'whole_interval_polynomial_lower':[str(v) for v in lower],'cap_precision':80,'cap_rows':cap,'moment_rows':moment,'moment_limit':mp.nstr(mp.gamma(mp.mpf('5.5'))/2**24,65),'seconds':time.monotonic()-start,'scope':'exact finite union regressions plus whole-interval polynomial signs; convex geometry transfer proved separately'}
 (Path(__file__).resolve().parents[1]/'results/parabolic_gap_python.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['cap_rows','moment_rows']},indent=2))
