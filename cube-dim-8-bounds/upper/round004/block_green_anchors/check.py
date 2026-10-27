from fractions import Fraction as Q
from itertools import permutations,product
from math import prod,comb
from pathlib import Path
import sympy as s,json,time,sys
D=Path(__file__).resolve().parent;start=time.perf_counter();r=json.loads((D/'result.json').read_text());t,u=s.symbols('t u',real=True)
# Independently evaluate the unsimplified, complete cone sum at exact rational t.
V=list(product((1,),(-1,1),(-1,1),(-1,1)));R=[(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1),(2,1,1,0),(3,1,1,1)];B=[((0,1,4,5),1),((4,1,2,3),1),((4,1,3,5),2)]
dot=lambda a,b:sum(x*y for x,y in zip(a,b));rate=lambda a:sum(abs(dot(v,a)) for v in V)
def mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
for name,x in [('A',(1,1,1,0)),('B',(1,1,0,0))]:
 f=s.sympify(r['p4_pi_minus4'][name],locals={'t':t})
 for q in [Q(1,10),Q(1,3),Q(1),Q(3)]:
  total=Q(0)
  for perm in permutations(range(4)):
   for signs in product((-1,1),repeat=4):
    for inds,det in B:
     den=(Q(1),Q(0))
     for j in inds:
      v=tuple(signs[k]*R[j][perm[k]] for k in range(4));den=mul(den,(q*rate(v),-dot(x,v)))
     total+=Q(det)*den[0]/(den[0]**2+den[1]**2)/16
  assert total==Q(f.subs(t,s.Rational(q.numerator,q.denominator)))
 # Verify known exact 4D anchor by partial fractions followed by integration.
 num,den=s.fraction(f);N=s.Poly(num,t);H=s.Poly(den,t)
 nu=sum(c*u**(i[0]//2) for i,c in N.terms());du=sum(c*u**(i[0]//2) for i,c in H.terms())
 def integrate(terms):
  ans=Q(0)
  for c,a,k in terms:ans+=c*Q(comb(2*k-2,k-1),2*4**(k-1))/a**(2*k-1)
  return ans
 terms=[]
 for term in s.Add.make_args(s.apart(nu/du,u)):
  n,d=s.fraction(term);factor,parts=s.factor_list(d,u);assert len(parts)==1;poly,k=parts[0];P=s.Poly(poly,u);a2=Q(P.nth(0)/P.nth(1));a=Q(s.sqrt(s.Rational(a2.numerator,a2.denominator)));c=Q(n)/Q(factor)/Q(P.nth(1))**k;terms.append((c,a,k))
 assert integrate(terms)=={'A':Q(17,640),'B':Q(1,24)}[name]
for row in r['anchors']:
 x,y=row['point'];f=s.sympify(r['p4_pi_minus4'][x],locals={'t':t})*s.sympify(r['p4_pi_minus4'][y],locals={'t':t})
 fractions=row['partial_fractions'];partial=sum((s.Rational(Q(z['coefficient']).numerator,Q(z['coefficient']).denominator)/(t*t+s.Rational(Q(z['a']).numerator,Q(z['a']).denominator)**2)**z['k'] for z in fractions),s.S.Zero)
 assert s.cancel(f-partial)==0
 ans=sum((Q(z['coefficient'])*Q(comb(2*z['k']-2,z['k']-1),2*4**(z['k']-1))/Q(z['a'])**(2*z['k']-1) for z in fractions),Q(0));assert ans==Q(row['G8_pi_minus7_coefficient'])
print(json.dumps({'status':'PASS','exact_cone_sum_sample_points_per_anchor':4,'known_n4_integrals':['17/640','1/24'],'partial_fraction_identities':3,'wall_seconds':time.perf_counter()-start},indent=2))
