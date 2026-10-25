"""Independent Acb640 jets and exact rational density Hermite checks; no kernel."""
import hashlib,json
from pathlib import Path
from fractions import Fraction as Q
from flint import arb,acb,fmpq,ctx
P=Path(__file__).resolve().parent
ctx.prec=640
def B(q):
 q=Q(q);return arb(fmpq(q.numerator,q.denominator))
def bounds(d):
 f=Q(2)**int(d['common_binary_exponent']);return int(d['lower_mantissa'])*f,int(d['upper_mantissa'])*f
def endpoints(v):
 l,e=v.lower().man_exp();u,f=v.upper().man_exp();e,f=int(e),int(f);g=min(e,f)
 return dict(lower_mantissa=str(int(l)*2**(e-g)),upper_mantissa=str(int(u)*2**(f-g)),common_binary_exponent=str(g))
def contained(n,r):
 nl,nu=bounds(n);rl,ru=bounds(r);return nl<=rl<=ru<=nu
def atanjet(x,u,k):
 if x==0:
  return B(0 if k==0 or k%2==0 else (-1)**((k-1)//2)/(Q(k)*u**k))
 if k==0:return (B(x)/B(u)).atan()
 i=acb(0,1);z=acb(B(x));v=acb(B(u));s=1 if k%2 else -1
 return (s*((i/(v+i*z))**k-(-i/(v-i*z))**k)/(2*k*i)).real
def hdd(nodes,jet):
 z=sorted(nodes);d=[jet(v,0) for v in z]
 for k in range(1,len(z)):
  d=[jet(z[j],k) if z[j]==z[j+k] else (d[j+1]-d[j])/(B(z[j+k]-z[j]) if isinstance(d[j],arb) else z[j+k]-z[j]) for j in range(len(z)-k)]
 return d[0]
def cmul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
def cpow(z,k):
 out=(Q(1),Q(0))
 for _ in range(k):out=cmul(out,z)
 return out
def pole(x,u):return (u/(u*u+x*x),x/(u*u+x*x))
def rationaljet(x,u,k):
 w=cmul(cpow((Q(0),Q(1)),k),cpow(pole(x,u),k+1));return w[1]
def product(nodes,u):
 out=(Q(1),Q(0))
 for x in nodes:out=cmul(out,pole(x,u))
 return out[0]
def main():
 reg=json.loads((P/'registration.json').read_text());native=json.loads((P/'anchors_native.json').read_text())
 assert native['parser_selftests_passed'] and len(reg['anchors'])==len(native['anchors'])==6
 a,T=Q(reg['a']),Q(reg['T']);rows=[]
 for j,(raw,n) in enumerate(zip(reg['anchors'],native['anchors'])):
  z=[Q(*v) for v in raw]
  def jet(x,k):return atanjet(x,a,k)-atanjet(x,T,k)
  v=endpoints(hdd(z,jet));ds=[endpoints(hdd(z+[z[i]],jet)) for i in range(6)]
  assert contained(n['value'],v),(j,'value')
  assert all(contained(nd,rd) for nd,rd in zip(n['partials'],ds)),(j,'partials')
  densities={}
  for label,u in [('p_a_cone',a),('p_T_cone',T)]:
   q=hdd(z,lambda x,k:rationaljet(x,u,k));prod=product(z,u);assert q==prod,(j,label,'DD/product')
   lo,hi=bounds(n[label]);assert lo<=q<=hi,(j,label,'containment')
   densities[label]=str(q)
  rows.append(dict(case=j,nodes=[str(v) for v in z],band_DD5_reference=v,partials_DD6_reference=ds,density_exact_rationals=densities,all_pass=True))
 out=dict(status='PASS',reference_precision_bits=640,method='Acb log-derived atan jets, sorted Hermite; density exact Fraction Hermite vs independent complex product',cases=rows,native_sha256=hashlib.sha256((P/'anchors_native.json').read_bytes()).hexdigest())
 (P/'anchors_validation.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS 6 cases,42 band quantities,12 exact density identities')
if __name__=='__main__':main()
