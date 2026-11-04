import mpmath as m,numpy as np,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def run(arg):
 tt,n=arg;m.mp.dps=50;t=m.mpf(tt);nodes,weights=np.polynomial.legendre.leggauss(n);mass=1-4*t*t
 totalQ=m.mpf(0);totalG=m.mpf(0);max_point_res=0
 for node,weight in zip(nodes,weights):
  s=m.mpf(str((node+1)/4));w=m.mpf(str(weight/4));lo=m.mpf(-40);hi=m.mpf(40)
  for _ in range(150):
   q=(lo+hi)/2;r=1+2*s*t+q*t*t;v=[r*(1-s),r*s]
   candidates=[z for u in v for z in [-2*u,2*u,-2-2*u,-2+2*u]];val=[]
   for a in candidates:
    b=[max(0,-a,u-a/2) for u in v]
    U=[min(1,max(0,1+u+a/2))+min(1,max(0,1-u+a/2)) for u in v]
    val.append(((1-t*b[0])*(1-t*b[1])-t*t*U[0]*U[1])/((1+t*a)*m.sqrt(mass)))
   if max(val)>1-t:lo=q
   else:hi=q
  q=(lo+hi)/2;totalQ+=w*q;totalG+=w*(2*s+t*q)**2
  max_point_res=max(max_point_res,float(abs(q-(-2-4*s+8*s*s))))
 coeff=(8*totalQ+4*totalG+8+8*t)/mass
 return {'t':tt,'Gauss_nodes':n,'quartic_estimate':str(coeff),'error_from_minus_two_thirds':str(coeff+m.mpf(2)/3),'max_radial_second_error':max_point_res}
if __name__=='__main__':
 start=time.monotonic();jobs=[(t,n) for t in ['.02','.01','.005','.002','.001'] for n in [64,128]]
 with ProcessPoolExecutor(max_workers=8) as pool:rows=list(pool.map(run,jobs))
 out={'rows':rows,'mpmath_version':m.__version__,'precision':50,'bisection_steps':150,'workers':8,'elapsed_seconds':time.monotonic()-start,'scope':'finite high precision evidence; unsplit angular quadrature not a certificate'}
 (Path(__file__).resolve().parents[1]/'results'/'k1_fourth_python.json').write_text(json.dumps(out,indent=2));print(out)
