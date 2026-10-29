"""Actual finite sums: multi-basis, same-field, same-query anchor characteristic functions."""
from pathlib import Path
from fractions import Fraction as F
import math,json,hashlib,time
import numpy as np
HERE=Path(__file__).resolve().parent

def mul(a,b):
 out=0
 while b:
  if b&1:out^=a
  b>>=1;a<<=1
  if a>>7:a^=0x83
 return out

def power(a,n):
 out=1
 while n:
  if n&1:out=mul(out,a)
  a=mul(a,a);n>>=1
 return out

def trace(a):
 out=0
 for _ in range(7):out^=a;a=mul(a,a)
 assert out in [0,1]
 return out

def rank(rows):
 p={}
 for x in rows:
  while x:
   k=x.bit_length()-1
   if k in p:x^=p[k]
   else:p[k]=x;break
 return len(p)

def walshnorm(f):
 a=np.array(f,dtype=float).copy();n=len(a);h=1
 while h<n:
  for j in range(0,n,2*h):
   x=a[j:j+h].copy();y=a[j+h:j+2*h].copy();a[j:j+h]=x+y;a[j+h:j+2*h]=x-y
  h*=2
 return float(abs(a/n).sum())

def bound(aa,eta,V,anchor):
 U=8.;W=U+V;r=W*aa.max();s=abs(eta);s1=float(s@aa);s2=float(s@(aa*aa))
 return math.sqrt(128)/127*U*V*math.exp(W*s1/(1-r))*((1+aa[anchor]*W)*(s1*s1+s2)/(1-r)**2+2*aa[anchor]*s1/(1-r))

def pair(z):return [float(z.real),float(z.imag)]

def main():
 start=time.perf_counter();K=127;B=8;eps=F(1,128);bases=[[1,2,4],[3,8,16],[5,32,64]]
 assert all(rank(beta)==3 for beta in bases)
 states=np.array([power(2,j) for j in range(K)]);assert len(set(states))==127
 labels=np.array([[sum(trace(mul(beta,x))<<i for i,beta in enumerate(basis)) for x in range(128)] for basis in bases])
 scales=[F(10,9),F(7,6),F(11,9)];aa=np.tile(np.array([float(F(1,2)/(K*L)) for L in scales]),3)
 m=np.arange(1,5);bm=4*729/29*((np.sin(np.pi*m*10/27)/(np.pi*m))**4+3*(np.sin(np.pi*m/9)/(np.pi*m))**4);cm=np.sqrt(2*bm/3);V=float(cm.sum())
 anchors=[0,4,8];omega=np.array([.5,1/3,1/6]);cases=[];norms=[];ranks=[]
 for mode,phases in [('shared',[F(0)]*3),('shifted',[F(0),F(1,9),F(-1,7)])]:
  lograw=[];logdec=[]
  for theta in [F(1,7),F(4,9)]:
   H=[];BB=[];allrows=[]
   for bi,basis in enumerate(bases):
    for L in scales:
     def endpoint(z):
      off=z.numerator//z.denominator;r=z-off;g=power(2,off%K)
      shifted=np.array([labels[bi,mul(g,x)] for x in range(128)])
      rows=[mul(g,beta) for beta in basis];allrows.extend(rows)
      return np.clip((float(r)-(shifted+.5)/B+float(eps)/2)/float(eps),0,1)-float(r)
     hh=endpoint(theta+K*L/2)-endpoint(theta-K*L/2)
     H.append(hh[states]);norms.append(walshnorm(hh));assert norms[-1]<=8+1e-12
     xx=(np.arange(K)+float(theta))/K+float(phases[bi])
     BB.append(np.sum((cm*np.sin(np.pi*m*float(L)))[:,None]*np.cos(2*np.pi*m[:,None]*xx),axis=0))
   ranks.append(dict(mode=mode,theta=str(theta),all_endpoint_rank=rank(allrows),row_count=len(allrows)))
   H=np.array(H);BB=np.array(BB)
   lograw.append(np.log1p(aa[:,None]*(H+BB)))
   logdec.append(np.log1p(aa[:,None,None]*(H[:,:,None]+BB[:,None,:])))
  def moment(tau,dec=False):
   return sum(np.exp(np.tensordot(tau,L,axes=(0,0))).mean() for L in (logdec if dec else lograw))/2
  etas=[np.zeros(9),np.ones(9)/12,np.array([1,0,-1,0,2,0,-1,0,1])/4,np.arange(-4,5)/20]
  for eta in etas:
   fs=[];fts=[];caps=[];As=[]
   for anchor in anchors:
    e=np.zeros(9);e[anchor]=1;A=moment(e).real;At=moment(e,True).real
    assert abs(A-At)<2e-15
    f=moment(e+1j*eta)/A;ft=moment(e+1j*eta,True)/At;cap=bound(aa,eta,V,anchor)/A
    assert abs(f-ft)<=cap+1e-13
    fs.append(f);fts.append(ft);caps.append(cap);As.append(A)
    for n in [1,100,10000]:
     gap=abs(f**n-ft**n);bnd=min(2,n*cap);assert gap<=bnd+1e-11
     cases.append(dict(mode=mode,eta=eta.tolist(),kind='single_anchor',anchor=anchor,n=n,A=A,
                       one_coordinate_gap=abs(f-ft),CF_gap=gap,bound=bnd,eta_l1=float(abs(eta).sum())))
   for n in [1,100,10000]:
    # Fixed weights on normalized anchor laws: one global bank label, versus iid per-coordinate labels.
    global_gap=abs(sum(omega*np.array(fs)**n)-sum(omega*np.array(fts)**n))
    iid_gap=abs(sum(omega*np.array(fs))**n-sum(omega*np.array(fts))**n)
    bnd=min(2,n*sum(omega*np.array(caps)))
    assert max(global_gap,iid_gap)<=bnd+1e-11
    cases.append(dict(mode=mode,eta=eta.tolist(),kind='mixture_anchors',n=n,weights=omega.tolist(),
                       global_mixture_CF_gap=global_gap,iid_mixture_CF_gap=iid_gap,bound=bnd))
 out=dict(scope='Floating finite sums only; no Haar quadrature certificate, event bound, or HL improvement.',
          field=dict(d=7,modulus=0x83,generator=2,b=3,K=K,bases=bases),phases=['1/7','4/9'],
          scales=list(map(str,scales)),epsilon=str(eps),p='1/2',M=4,U_bound=8,V=V,r=float((8+V)*aa.max()),
          max_computed_Walsh_norm=max(norms),endpoint_ranks=ranks,cases=cases,
          elapsed_seconds=time.perf_counter()-start,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (HERE/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'cases':len(cases),'max_single_CF_gap':max([r['CF_gap'] for r in cases if r['kind']=='single_anchor']),
                   'max_Walsh_norm':max(norms),'endpoint_ranks':ranks,'seconds':out['elapsed_seconds']}))
if __name__=='__main__':main()
