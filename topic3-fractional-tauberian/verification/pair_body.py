from fractions import Fraction as F
from witness_capacity import atoms,solve
from pathlib import Path
import mpmath as m,json,time,platform,sys,numpy,scipy

def enclosure(al,N):
 W=1/al-al;c=2*al*al-1;h=W/N;scale=10**40
 def f(u):return al-c/(al-u)
 lo=hi=0
 for j in range(N):
  a=j*h;b=(j+1)*h;trap=h*(f(a)+f(b))/2;mid=h*f((a+b)/2)
  lo+=(trap*scale).numerator//(trap*scale).denominator
  z=mid*scale;hi+=-((-z.numerator)//z.denominator)
 return F(lo,scale),F(hi,scale)
def replay(data):
 for r in data['pair_cases']:
  al=F(r['alpha']);u=F(r['u']);v=F(r['v']);G=lambda z:(z-u)*(z-v)-2*al*z+1
  assert G(al)==F(r['G_L']) and G(1/al)==F(r['G_U'])
  if r['kind']=='inside':
   ell=F(r['ell']);caps=atoms([(-ell/2,ell/2,-ell/2,ell/2),(u-ell/2,u+ell/2,v-ell/2,v+ell/2)]);mm={int(s):F(x) for s,x in r['primal_masses'].items()}
   assert sum(mm.values())<=1 and all(0<=x<=caps.get(s,F(0)) for s,x in mm.items())
   slack=min(sum(x for s,x in mm.items() if s>>i&1)-al*ell for i in range(2));assert slack==F(r['primal_strict_slack']) and slack>0
  else:assert G(al)<=0 and G(1/al)<=0
 for r in data['exact_area_enclosures']:
  lo,hi=enclosure(F(r['alpha']),r['subintervals']);assert lo==F(r['lower_rational']) and hi==F(r['upper_rational']) and lo<=hi
 print({'rational_pair_certificates_replayed':len(data['pair_cases']),'full_interval_area_enclosures_replayed':len(data['exact_area_enclosures'])})
if __name__=='__main__':
 if '--check-only' in sys.argv:
  replay(json.loads((Path(__file__).resolve().parents[1]/'results'/'pair_body_python.json').read_text()));raise SystemExit(0)
 start=time.monotonic();pairs=[];intervals=[]
 for al in [F(71,100),F(3,4),F(9,10),F(99,100)]:
  W=1/al-al;c=2*al*al-1
  for frac in [F(1,10),F(1,4),F(1,2),F(3,4),F(9,10)]:
   u=frac*W;vb=al-c/(al-u)
   for kind,shift in [('inside',-1),('boundary',0),('outside',1)]:
    v=vb+shift*(1-al)/100;assert 0<=u<al and 0<=v<al
    G=lambda z:(z-u)*(z-v)-2*al*z+1
    row={'alpha':str(al),'u':str(u),'v':str(v),'kind':kind,'G_L':str(G(al)),'G_U':str(G(1/al))}
    if kind=='inside':
     dz=(1/al-al)/2
     while G(al+dz)<=0:dz/=2
     ell=al+dz;Q=[(-ell/2,ell/2,-ell/2,ell/2),(u-ell/2,u+ell/2,v-ell/2,v+ell/2)];res=solve(atoms(Q),[al*ell]*2,denominator_limit=10**12,solver_options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10});assert res['status']=='exact-feasible';row.update({'ell':str(ell),**res})
    else:assert G(al)<=0 and G(1/al)<=0
    pairs.append(row)
  previous=None
  for N in [16,64,256,1024]:
   lo,hi=enclosure(al,N);assert lo<=hi
   if previous:assert previous[0]<=lo<=hi<=previous[1]
   previous=(lo,hi);intervals.append({'alpha':str(al),'subintervals':N,'lower_rational':str(lo),'upper_rational':str(hi),'gap':str(hi-lo),'coverage':'all u in [0,W] by concavity; exact outward rounding per strip'})
 m.mp.dps=80;values=[];tails=[]
 for aa in ['.71','.75','.9','.99','.999','.9999']:
  a=m.mpf(aa);W=1/a-a;c=2*a*a-1;U=1-a*a+c*m.log(c/(a*a));integral=m.quad(lambda y:a*W*W*(1-y)/(a-W*y),[0,1]);t=1-a;R=(1/a-a)**2/2
  values.append({'alpha':aa,'upper':str(U),'lower':str(R),'integral_residual':str(U-integral),'upper_cubic_ratio':str((U-2*t*t)/t**3),'lower_cubic_ratio':str((R-2*t*t)/t**3)})
  x=1-a*a
  for N in [2,3,4,6,8]:
   trunc=sum(m.mpf(2**j-j-1)/(j*(j-1))*x**j for j in range(2,N+1));rem=U-trunc;low=m.mpf(2**(N+1)-(N+1)-1)/(N*(N+1))*x**(N+1);up=(2*x)**(N+1)/(N*(N+1)*(1-2*x));assert low<=rem<=up;tails.append({'alpha':aa,'N':N,'remainder':str(rem),'lower_bound':str(low),'upper_bound':str(up)})
 for r in intervals:
  a=m.mpf(F(r['alpha']).numerator)/F(r['alpha']).denominator;c=2*a*a-1;U=1-a*a+c*m.log(c/(a*a));l=F(r['lower_rational']);h=F(r['upper_rational']);assert m.mpf(l.numerator)/l.denominator<=U<=m.mpf(h.numerator)/h.denominator
 out={'environment':{'python':platform.python_version(),'mpmath':m.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__},'pair_cases':pairs,'exact_area_enclosures':intervals,'values':values,'tail_checks':tails,'precision':80,'solver_options':{'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10,'rational_denominator_limit':10**12},'outward_rounding_denominator':str(10**40),'elapsed_seconds':time.monotonic()-start,'scope':'rational finite pair certificates and full-interval concavity enclosures; general theorems analytic'}
 (Path(__file__).resolve().parents[1]/'results'/'pair_body_python.json').write_text(json.dumps(out,indent=2));print({'pair_cases':len(pairs),'area_enclosures':len(intervals),'max_integral_residual':str(max(abs(m.mpf(r['integral_residual'])) for r in values)),'seconds':out['elapsed_seconds']})
