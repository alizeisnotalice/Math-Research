from fractions import Fraction as F
from witness_capacity import atoms,solve
from pathlib import Path
import mpmath as m,json,time,platform,numpy,scipy,sympy as sy

def layer(a,l):
 S=[l*l+4*a*l,l*l,l*l-4*a*a,(l-2*a)**2]
 if S[3]>=1:M=F(4)
 elif S[2]>=1:M=3+S[3]
 elif S[1]>=1:M=2+S[3]+S[2]
 elif S[0]>=1:M=1+S[3]+S[2]+S[1]
 else:M=4*l*l
 return S,M
if __name__=='__main__':
 start=time.monotonic();rows=[];alpha=F(19,20)
 for a in [F(1,100),F(1,20),F(1,10),F(1,8)]:
  for l in [F(1,2),F(9,10),F(99,100),F(1),F(1001,1000),F(101,100),F(21,20),F(6,5),F(3,2),F(2)]:
   C=[(a,0),(-a,0),(0,a),(0,-a)];R=[(x-l/2,x+l/2,y-l/2,y+l/2) for x,y in C];caps=atoms(R);ss,MM=layer(a,l)
   actual=[sum(v for mask,v in caps.items() if mask.bit_count()>=k) for k in range(1,5)];assert ss==actual
   budget=F(1);top=F(0)
   for k in range(4,0,-1):
    cap=sum(v for mask,v in caps.items() if mask.bit_count()==k);take=min(cap,budget);top+=k*take;budget-=take
   assert top==MM
   res=solve(caps,[alpha*l]*4);gap=abs(res['eta_float']-float(MM/4-alpha*l));assert gap<1e-10
   rows.append({'a':str(a),'ell':str(l),'layer_capacities':[str(v) for v in ss],'M':str(MM),'LP_error':gap})
 m.mp.dps=80;roots=[]
 for al in ['.91','.95','.99','.999','.9999']:
  al=m.mpf(al);t=1-al;lo=m.mpf(0);hi=m.mpf(1)/8
  for _ in range(280):
   mid=(lo+hi)/2;H=(1+2*mid*mid)/m.sqrt(1+4*mid*mid)-mid
   if H>al:lo=mid
   else:hi=mid
  a=(lo+hi)/2;poly=t+2*t**4-8*t**6+16*t**7+30*t**8-160*t**9+64*t**10
  roots.append({'alpha':str(al),'a':str(a),'fourth_ratio':str((a-t)/t**4),'sixth_ratio':str((a-t-2*t**4)/t**6),'order10_residual':str(a-poly),'rectangle_radius':str((1/al-al)/2)})
 u=sy.symbols('u');H=(1+2*u*u)/sy.sqrt(1+4*u*u)-u;inverse_kernel=sy.series(u/(1-H),u,0,10).removeO();coeffs=[sy.expand(inverse_kernel**j).coeff(u,j-1)/j for j in range(1,11)];assert coeffs==[1,0,0,2,0,-8,16,30,-160,64]
 out={'lagrange_coefficients':[str(x) for x in coeffs],'environment':{'python':platform.python_version(),'mpmath':m.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'sympy':sy.__version__},'layer_cases':rows,'roots':roots,'precision':80,'bisection_steps':280,'max_LP_error':max(r['LP_error'] for r in rows),'elapsed_seconds':time.monotonic()-start,'scope':'exact rational finite layer tests and high precision roots; general theorem proved analytically'}
 (Path(__file__).resolve().parents[1]/'results'/'common_scale_python.json').write_text(json.dumps(out,indent=2));print({'cases':len(rows),'max_LP_error':out['max_LP_error'],'roots':roots,'seconds':out['elapsed_seconds']})
