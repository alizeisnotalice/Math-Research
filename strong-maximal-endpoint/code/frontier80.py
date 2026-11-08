import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import math,json,sys,time
from pathlib import Path
from fractions import Fraction as Q
from concurrent.futures import ProcessPoolExecutor
from sliver79 import logi
def add(a,b):
 c=[Q(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c
def mul(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def coeff(alpha,delta):
 v=mul([-alpha,Q(1)],[-alpha,Q(1)]);w=[Q(1),-2*alpha,4*alpha*alpha]
 neg=[-z for z in mul(v,w)]
 pos=[Q(0),4*delta*alpha**2,-8*delta*alpha**3,4*delta*alpha**2]
 return add(neg,pos)
def bernstein(poly,l,r):
 n=len(poly)-1
 power=[sum((poly[k]*math.comb(k,j)*l**(k-j)*(r-l)**j for k in range(j,n+1)),Q(0)) for j in range(n+1)]
 return [sum((power[j]*Q(math.comb(i,j),math.comb(n,j)) for j in range(i+1)),Q(0)) for i in range(n+1)]
def cell_cert(alpha,eps,l,r,depth=0):
 bs=bernstein(coeff(alpha,eps),l,r)
 if min(bs)>0:return min(bs),1,depth
 assert depth<10
 c=(l+r)/2;a=cell_cert(alpha,eps,l,c,depth+1);b=cell_cert(alpha,eps,c,r,depth+1)
 return min(a[0],b[0]),a[1]+b[1],max(a[2],b[2])
def chunk(job):
 m,start,end=job;h=Q(1,12*m);eps=h*h;mins=[];leaves=0;depth=0
 for k in range(start,end):
  a=Q(2,3)+k*h;l=max(Q(2,3),a-h/2);r=min(Q(3,4),a+h/2)
  v,c,d=cell_cert(a,eps,l,r);mins.append(v);leaves+=c;depth=max(depth,d)
 return dict(M=m,minimum=str(min(mins)),leaves=leaves,depth=depth,nodes=end-start)
def X(u):return (u+1/u)/2
def Y(u):return u+1/(4*u)
def grid_k(m,x):
 guess=x-math.sqrt(float((x-1)*(x+1-Q(1,m))));i=max(0,min(m-1,math.floor(m*guess)))
 def der(k):
  a=Q(k,m);return a*a-2*x*a+1+(x-1)/m
 while der(i)<0:i-=1
 while der(i+1)>0:i+=1
 assert 0<=i<m and der(i)>=0 and der(i+1)<=0
 def kval(k):return Q(m*(m-1)-k*(k-1),2*m*m)/(x-Q(k,m))
 a,b=kval(i),kval(i+1)
 return max(a,b),i if a>=b else i+1
def stepjob(m):
 rows=[];delta=Q(4,m)
 for k in range(33):
  u=Q(2,3)+Q(k,384);x=X(u);y=Y(u)-delta
  kx,ix=grid_k(m,x);ky,iy=grid_k(m,y);val=2*kx*ky
  assert val>=1+Q(68,15*m)
  rows.append(dict(u=str(u),x=str(x),y=str(y),maximizing_anchor_indices=[ix,iy],exact_maximal_value=str(val)))
 return dict(mesh=m,uniform_source_error_upper=str(Q(4,m)),halo_depth=str(delta),budget_lower=str(Q(math.isqrt(m),147456)),samples=rows)
if __name__=='__main__':
 out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True);t=time.monotonic();ms=[2**k for k in range(13)]
 jobs=[(m,j,min(j+32,m+1)) for m in ms for j in range(0,m+1,32)]
 with ProcessPoolExecutor(max_workers=12) as p:blocks=list(p.map(chunk,jobs));steps=list(p.map(stepjob,[4**k for k in range(5,21)]))
 ln2=logi(Q(2));jf=((1+5*ln2[0])/8,(1+5*ln2[1])/8);coverage=[]
 for z in steps:
  z['source_modular_interval']=[str(max(Q(0),jf[0]-Q(12,z['mesh']))),str(jf[1])]
  z['normalized_budget_lower']=str(Q(z['budget_lower'])/jf[1])
  z['actual_halo_ratio_upper']=str((1+8*ln2[1])/(jf[0]-Q(12,z['mesh'])))
 for m in ms:
  c=[z for z in blocks if z['M']==m];eps=Q(1,144*m*m)
  coverage.append(dict(M=m,epsilon=str(eps),anchor_count=m+1,certified_nodes=sum(z['nodes'] for z in c),minimum_Bernstein_coefficient=str(min(Q(z['minimum']) for z in c)),leaves=sum(z['leaves'] for z in c),maximum_subdivision_depth=max(z['depth'] for z in c),
    unnormalized_margin_lower=str(Q(229,768)*eps),budget_lower=str(Q(m,6144)),budget_upper_interval=[str((m+1)*v) for v in jf],discarded_layer_volume=str(eps/24)))
 ns=[2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576]
 high=[]
 for n in ns:
  lc=math.log(3072) if n==2 else math.log(1179648)+(n-2)*math.log(4)
  high.append(dict(N=n,depth_dyadic_exponent=8*n,log_section_constant=lc,log_truncated_cover_budget_lower=-math.log(24)-lc+4*n*math.log(2)))
 doc=dict(workers=12,seconds=time.monotonic()-t,source_modular_interval=[str(v) for v in jf],curvature_constants=dict(x_length='1/24',maximum_depth='1/96',source_column_charge_lower='1/768',section_budget_coefficient=3072),coverage=coverage,rational_step=steps,high_dimension=high)
 (out/'frontier80.json').write_text(json.dumps(doc,indent=2)+'\n')
 print('PASS',sum(z['certified_nodes'] for z in coverage),'rational anchor cells; maximum Bernstein subdivision',max(z['maximum_subdivision_depth'] for z in coverage),'and',len(steps)*33,'exact step-source halo tests;',doc['seconds'],'seconds')
