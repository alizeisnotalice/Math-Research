from fractions import Fraction as F
import json
R=F(1);a=F(1);L=F(1);rows=[]
for n in (2,3):
 B=F(n*2**(n-1)+1)
 for l,b,m in [(F(0),F(0),F(0)),(F(1,2),F(1,4),B/4),(F(1),F(0),B/4),(F(1,4),F(2),F(0))]:
  z=b+m/B; cap=lambda x:a*((R+x)**n-R**n)
  lhs=cap(l)-max(B*max(l-b,0)-m,0);rhs=cap(z)
  rows.append({'n':n,'l':str(l),'b':str(b),'m':str(m),'B':str(B),'z':str(z),'outside_L':z>L,'identity':max(B*max(l-b,0)-m,0)==B*max(l-z,0),'margin':str(rhs-lhs),'pass':lhs<=rhs})
print(json.dumps(rows,indent=2))
