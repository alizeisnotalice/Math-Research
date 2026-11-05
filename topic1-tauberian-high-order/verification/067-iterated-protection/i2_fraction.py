from fractions import Fraction as F
import json

def add(*ps):
 out=[F(0)]*max(map(len,ps))
 for p in ps:
  for i,x in enumerate(p):out[i]+=x
 return out
def sc(c,p):return [c*x for x in p]
def mul(p,q):
 r=[F(0)]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):r[i+j]+=x*y
 return r
def diff(p):return [F(i)*p[i] for i in range(1,len(p))] or [F(0)]
def primitive(p):return [F(0)]+[x/F(i+1) for i,x in enumerate(p)]
def val(p,x):return sum(c*x**i for i,c in enumerate(p))
def integrate(p,l,h):
 P=primitive(p);return val(P,h)-val(P,l)
def calc(d):
 hl=[F(0),F(0),(1-d)/2]
 hh=[-d*d*(1-2*d)/(2*(1-d)),d-d*d/(1-d),d*d/(2*(1-d))]
 cl=primitive(hl);ch=primitive(hh);ch[0]+=val(cl,d)-val(ch,d)
 total=val(ch,F(1));T=[F(0),F(1)];one=[F(1)];results=[]
 for h,cum,L,R in ((hl,cl,F(0),d),(hh,ch,d,F(1))):
  for kind in ('inverse','cross'):
   if kind=='inverse':
    q=h;m=add(sc(F(1,2),mul(h,h)),mul(h,[F(1),F(-1)]),cum);jac=diff(h)
   else:
    q=add([F(1),F(-1)],h);m=add([2*total],sc(-1,cum),sc(F(1,2),mul(h,h)));jac=add(one,sc(-1,diff(h)))
   z=add(q,sc(-1,m));integrand=mul(add(sc(2*d,mul(add(one,q),z)),sc(-2*d*d,mul(z,z))),jac)
   results.append(integrate(integrand,L,R))
 return sum(results),results
coeff=[-50400,98280,-199584,228816,-216972,195792,-137916,73413,-21099,-3367,2617,186,-26,20]
def expected(d):return -d*sum(F(x)*d**j for j,x in enumerate(coeff))/30240
rows=[]
for d in (F(1,20),F(1,40),F(1,80)):
 got,parts=calc(d);rows.append({'delta':str(d),'I2':str(got),'four_branch_integrals':[str(x) for x in parts],'closed_polynomial_matches':got==expected(d)})
print(json.dumps({'method':'Fraction polynomial-in-threshold integration, independent implementation of branch substitution','rows':rows,'all_pass':all(x['closed_polynomial_matches'] for x in rows)},indent=2))
