from fractions import Fraction as Q
from itertools import permutations,product
from collections import defaultdict
from math import prod,comb,isqrt
from pathlib import Path
from flint import fmpq_poly, fmpq
import sympy as s, mpmath as mp, hashlib, json,time,datetime,sys,resource
D=Path(__file__).resolve().parent;start=datetime.datetime.now(datetime.timezone.utc).isoformat();w=time.perf_counter();cpu=time.process_time()
sources=[Path('/Users/zhengzhihao/Desktop/T/cube8-source-package-20261025/cube-finite-dim-bounds/upper/u11_green.py'),Path('/Users/zhengzhihao/Desktop/T/cube8-source-package-20261025/cube-finite-dim-bounds/audit/second_review/n4/check.py')]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
try:
 V=tuple((1,)+v for v in product((-1,1),repeat=3));R=((1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(2,1,1,0),(3,1,1,1));base=[((0,1,4,5),1),((4,1,2,3),1),((4,1,3,5),2)]
 dot=lambda a,b:sum(x*y for x,y in zip(a,b)); rate=lambda r:sum(abs(dot(v,r)) for v in V)
 cones=[]
 for perm in permutations(range(4)):
  for signs in product((-1,1),repeat=4):
   for inds,det in base:cones.append((tuple(tuple(signs[j]*R[i][perm[j]] for j in range(4)) for i in inds),det))
 assert len(cones)==1152 and len({tuple(sorted(c)) for c,d in cones})==1152
 z,t,u=s.symbols('z t u', real=True);rat=lambda a:s.Rational(a.numerator,a.denominator)
 forms={};groupcounts={};tail={}
 for name,x in [('A',(1,1,1,0)),('B',(1,1,0,0))]:
  groups=defaultdict(Q)
  for cone,det in cones:
   q=[rate(r) for r in cone];bs=tuple(sorted(Q(dot(x,r),qq) for r,qq in zip(cone,q)));groups[bs]+=Q(det,prod(q))
  assert len(groups)<200000
  P=fmpq_poly([0]); H=fmpq_poly([1])
  for bs,c in groups.items():
   den=prod((fmpq_poly([-fmpq(str(b)),1]) for b in bs),start=fmpq_poly([1]))
   num=P*den+H*fmpq(str(c)); common=num.gcd(H*den); P=num//common; H=(H*den)//common
   assert len(P)+len(H)<200000
  f=s.cancel(sum(rat(Q(str(c)))*z**i for i,c in enumerate(P))/sum(rat(Q(str(c)))*z**i for i,c in enumerate(H)))
  assert len(s.Poly(s.numer(f),z).terms())+len(s.Poly(s.denom(f),z).terms())<200000
  p=s.cancel(f.subs(z,-s.I*t)/16); assert s.simplify(s.im(p))==0 and s.simplify(p-p.subs(t,-t))==0
  assert s.denom(p).subs(t,0)!=0
  forms[name]=p;groupcounts[name]=len(groups)
  # Triangle inequality on raw cone sum: |t-i*b|>=t, t>0.
  tail[name]=sum((Q(det,16*prod(rate(v) for v in cone)) for cone,det in cones),Q(0))
 rows=[]
 for name,x,y in [('AA','A','A'),('AB','A','B'),('BB','B','B')]:
  f=s.cancel(forms[x]*forms[y]); N=s.Poly(s.numer(f),t); De=s.Poly(s.denom(f),t); assert all(e[0]%2==0 for e,c in N.terms()+De.terms())
  nu=sum(c*u**(e[0]//2) for e,c in N.terms());du=sum(c*u**(e[0]//2) for e,c in De.terms());assert len(s.Poly(nu,u).terms())+len(s.Poly(du,u).terms())<200000
  assert s.degree(nu,u)<s.degree(du,u)
  pf=s.apart(nu/du,u);terms=s.Add.make_args(pf);assert len(terms)<200000
  coef=Q(0);termrows=[]
  for term in terms:
   n,d=map(lambda v:s.Poly(v,u),s.fraction(term)); fac=s.factor_list(d.as_expr(),u);assert len(fac[1])==1
   linear,k=fac[1][0];L=s.Poly(linear,u);assert L.degree()==1 and n.degree()==0
   a2=Q(L.nth(0)/L.nth(1)); sn=isqrt(a2.numerator);sd=isqrt(a2.denominator);assert sn*sn==a2.numerator and sd*sd==a2.denominator and a2>0
   a=Q(sn,sd);c=Q(n.nth(0))/Q(fac[0])/Q(L.nth(1))**k
   coef+=c*Q(comb(2*k-2,k-1),2*4**(k-1))/a**(2*k-1);termrows.append({'coefficient':str(c),'a':str(a),'k':k})
  assert coef>0
  vals=[]
  for digits in (80,160):
   with mp.workdps(digits):
    fn=s.lambdify(t,f,'mpmath');v=mp.quad(fn,[0,mp.mpf(1)/8,mp.mpf(1)/4,1,mp.inf])/mp.pi**8
    expect=(mp.mpf(coef.numerator)/coef.denominator)/mp.pi**7
    assert mp.almosteq(v,expect,rel_eps=mp.mpf(10)**(-digits+20));vals.append({'digits':digits,'quadrature':mp.nstr(v,digits-10),'closed_form':mp.nstr(expect,digits-10)})
  rows.append({'point':name,'integrand_pi_minus8':str(f),'partial_fractions':termrows,'G8_pi_minus7_coefficient':str(coef),'quadrature':vals,'small_t':'finite and continuous: complete simplified denominator nonzero at zero','t_ge_1_abs_bound_pi_minus8':str(tail[x]*tail[y])+' * t^-8'})
 out={'status':'completed','cone_count':len(cones),'group_counts':groupcounts,'p4_pi_minus4':{k:str(v) for k,v in forms.items()},'p4_t_ge_1_abs_bound_pi_minus4':{k:str(v)+' * t^-4' for k,v in tail.items()},'anchors':rows,'scope':'point values only; no global cube minimum or C8'}
 (D/'result.json').write_text(json.dumps(out,indent=2)+'\n');status='completed';err=None;exitcode=0
except BaseException as e:status='failed';err=repr(e);exitcode=1
receipt={'status':status,'error':err,'exit_code':exitcode,'start_utc':start,'finish_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wall_seconds':time.perf_counter()-w,'cpu_seconds':time.process_time()-cpu,'peak_rss':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'sympy':s.__version__,'mpmath':mp.__version__,'source_sha256':{str(p):sha(p) for p in sources+[D/'registration.md',D/'run.py']}}
(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));sys.exit(exitcode)
