"""Exact local checks for the two-atom proof; no numerical scan certifies geometry."""
import sympy as s
from fractions import Fraction as F

a,b=s.symbols("a b", real=True)
H=b*(a+b-1)/2
H0=(1-a)/2
ts=[0,(1-a)/2,(1-b)/2,(1+b)/2,(1+a)/2,1]
la=[a,a,(a+b)/2,(a-b)/2,0,0]
lb=[b,b,b,0,0,0]
expected=[
[0],
[-H0,-H0-b*b/4],
[-H,-((1-a)*(1-b)+2*b*b)/4,-H-b*b/4],
[H,H-(1-a)*(1-b)/4,b*b/4,H-b*b/4],
[H0,b*b/4,H0*(1+b)/2,H0*(1-b)/2,H0-b*b/4],
[0,0,0,0,0,0],
]
checks=[]
for i in range(6):
 for j in range(i+1):
  actual=(1-ts[i])*(1-ts[j])-la[i]*la[j]-lb[i]*lb[j]
  diff=s.rem(s.Poly(s.expand(actual-expected[i][j]),a),
             s.Poly(a*a+b*b-1,a)).as_expr()
  checks.append(diff==0)
print("all 21 exact vertex-table identities:", all(checks))
def circle_zero(expr):
 return s.rem(s.Poly(s.expand(expr),a),s.Poly(a*a+b*b-1,a)).as_expr()==0
print("H-H0 identity:",circle_zero(H-H0-(1-a)*(a-b)/2))
print("H-b^2/4 identity:",s.expand(H-b*b/4-b*(2*a+b-2)/4)==0)
print("a^2-(1-b/2)^2 identity:",circle_zero(a*a-(1-b/2)**2-b*(1-5*b/4)))
c=F(3,2)-s.sqrt(2)/4
print("N2 constant:",c, "decimal:",s.N(c,30))
area_single_union=F(39,5)
r=F(7,5)
inner_width=2*r-F(19,10)
extra=inner_width*2*(r-1)
print("rational inner witness area:",area_single_union+extra)
print("rational inner witness ratio:",(area_single_union+extra)/8)
print("strict inner average:",F(8,1)/(4*r*r))
area=(98-39*s.sqrt(2))/5
print("full irrational area:",area)
print("full ratio > rational certificate:",bool(area/8 > F(213,200)))
print("all checks passed:",all(checks) and extra==F(18,25))

