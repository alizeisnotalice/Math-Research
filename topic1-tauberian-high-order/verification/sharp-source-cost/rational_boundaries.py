from fractions import Fraction as F
from pathlib import Path
import itertools,json
positive=lambda x:max(F(0),x)
vals=[F(0),F(1,4),F(1),F(2)]; zs=[F(-3),F(-2),F(-1),F(-1,4),F(0),F(1,4),F(1),F(2),F(3)];rows=[]
for x,y,z in itertools.product(vals,vals,zs):
 M=max(x,y);c=max(-M,min(M,z));Pp=positive(x+c);Pm=positive(x-c);Qp=positive(y+c);Qm=positive(y-c)
 excess=(Pp*Qp*(Pp+Qp)+Pm*Qm*(Pm+Qm))/2-x*y*(x+y)
 assert excess>=0
 if z<c:assert Pp==Qp==0
 if z>c:assert Pm==Qm==0
 rows.append({'x':str(x),'y':str(y),'z':str(z),'clip':str(c),'excess':str(excess)})

def integral_breaks(pairs):
 points=sorted(set([F(0)]+[q for c,p,q in pairs]));out=F(0)
 for l,h in zip(points[:-1],points[1:]):
  for c,p,q in pairs:
   if q>=h:out+=c*p*(q*(h-l)-(h*h-l*l)/2)
 return out
examples=[[(F(1),F(1),F(2)),(F(2),F(3),F(1))],[(F(1,3),F(0),F(2)),(F(4),F(1,2),F(0))],[(F(2),F(3),F(5)),(F(7),F(2),F(1))]]
for pairs in examples:
 a=integral_breaks(pairs);b=integral_breaks([(c,q,p) for c,p,q in pairs]);target=sum((c*p*q*(p+q)/2 for c,p,q in pairs),F(0));assert a+b==target
for n in range(2,17):
 v=F(1,8);r=F(3,2)
 for p,q in ((F(1,32),F(1,64)),(F(1,128),F(1,256))):
  H=min(q,p*q/v**2);lhs=(v/v**n)*(v*r)**(n-2)*(p*r)*(H*r);rhs=r**n*p*q*min(1/v,p/v**3)
  assert lhs==rhs and H<=q and v**2*H<=p*q
  X=H/3;assert p*q<=3*v*X
result={'clip_cases':len(rows),'clip_equality_cases':sum(F(r['excess'])==0 for r in rows),'strict_left_cases':sum(F(r['z'])<F(r['clip']) for r in rows),'strict_right_cases':sum(F(r['z'])>F(r['clip']) for r in rows),'targets_examples':len(examples),'large_angle_cases':30,'all_pass':True,'rows':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
