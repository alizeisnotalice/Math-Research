#!/usr/bin/env python3
"""Exact geometry and rational-log certificates for concave geometric capacity."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import json,math,time,sys,resource,hashlib
import numpy as np
import scipy
from scipy.optimize import linprog,minimize
HERE=Path(__file__).resolve().parent
SOURCE=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r716'/'wolfram.json'
TERMS=48

def atomize(boxes):
 k=len(boxes);d=len(boxes[0]);masses={}
 for s in range((1<<k)-1,0,-1):
  ids=[i for i in range(k) if s>>i&1];v=F(1)
  for j in range(d):v*=max(F(0),min(F(boxes[i][j][1]) for i in ids)-max(F(boxes[i][j][0]) for i in ids))
  masses[s]=v-sum(a for t,a in masses.items() if t&s==s)
 masks=sorted(s for s,a in masses.items() if a>0);mass=[masses[s] for s in masks]
 inc=[[int(s>>i&1) for s in masks] for i in range(k)];vol=[sum(a*x for a,x in zip(mass,row)) for row in inc]
 assert all(a>=0 for a in masses.values())
 return mass,inc,vol,masks

def hall_check(pi,mass,inc,masks):
 om=sum(mass);k=len(pi)
 return sum(pi)==1 and all(p>=0 for p in pi) and all(sum(pi[i] for i in range(k) if s>>i&1)<=sum(mass[j] for j,t in enumerate(masks) if s&t)/om for s in range(1,1<<k))

def lf_weights(mass,inc,vol,masks,p):
 k=len(vol);a=len(mass);om=sum(mass);vf=np.array([float(x) for x in vol]);mf=np.array([float(x) for x in mass]);ic=np.array(inc,dtype=float)
 aub=[];bub=[]
 for j in range(a):aub.append(list(ic[:,j])+[-int(t==j) for t in range(a)]);bub.append(1.)
 for s in range(1,1<<k):
  coverage=sum(mass[j] for j,t in enumerate(masks) if s&t)/om
  aub.append([float(vol[i]*(int(s>>i&1)-coverage)) for i in range(k)]+[0.]*a);bub.append(0.)
 res=linprog(np.r_[-float(p)*vf,mf],A_ub=np.array(aub),b_ub=np.array(bub),bounds=(0,None),method='highs')
 assert res.success,res.message
 lam=[max(F(0),F(str(x)).limit_denominator(10**7)) for x in res.x[:k]];scale=sum(v*l for v,l in zip(vol,lam));assert scale>0;pi=[v*l/scale for v,l in zip(vol,lam)]
 assert hall_check(pi,mass,inc,masks),'exact Hall reconstruction failed'
 return pi,res.message

def testing_price(pi,mass,inc,vol,p):
 density=[sum(pi[i]*inc[i][j]/vol[i] for i in range(len(pi))) for j in range(len(mass))]
 assert min(density)>0
 remain=p;price=F(0)
 for j in sorted(range(len(mass)),key=lambda j:density[j],reverse=True):
  hit=mass[j]*density[j]
  if remain<=hit:price+=remain/density[j];remain=0;break
  price+=mass[j];remain-=hit
 assert remain==0
 return price

def primal_source(mass,inc,vol,p,price):
 a=len(mass);k=len(vol);om=sum(mass);aub=[];bub=[]
 for i in range(k):
  for j in range(a):
   if inc[i][j]:aub.append([-float(inc[i][t]*mass[t]/vol[i]) for t in range(a)]+[float(t==j) for t in range(a)]);bub.append(0.)
 aub.append([0.]*a+[-float(x) for x in mass]);bub.append(-float(p*om))
 rr=linprog(np.r_[[float(x) for x in mass],np.zeros(a)],A_ub=np.array(aub),b_ub=np.array(bub),bounds=(0,1),method='highs');assert rr.success
 f=[F(str(x)).limit_denominator(10**7) for x in rr.x[:a]];ell=[F(str(x)).limit_denominator(10**7) for x in rr.x[a:]]
 av=[sum(mass[j]*inc[i][j]*f[j] for j in range(a))/vol[i] for i in range(k)]
 assert all(0<=x<=1 for x in f+ell) and all(ell[j]<=av[i] for i in range(k) for j in range(a) if inc[i][j])
 assert sum(mass[j]*ell[j] for j in range(a))>=p*om
 assert sum(mass[j]*f[j] for j in range(a))==price
 return f,ell

@lru_cache(None)
def logunit(y):
 assert 1<=y<=2
 t=(y-1)/(y+1);t2=t*t;term=t;lower=F(0)
 for j in range(TERMS):lower+=2*term/(2*j+1);term*=t2
 return lower,lower+2*term/((2*TERMS+1)*(1-t2))
@lru_cache(None)
def logbounds(x):
 assert x>0
 k=x.numerator.bit_length()-x.denominator.bit_length();y=x/(F(2)**k)
 while y<1:k-=1;y*=2
 while y>=2:k+=1;y/=2
 lo,hi=logunit(y);l2,h2=logunit(F(2))
 return lo+k*(l2 if k>=0 else h2),hi+k*(h2 if k>=0 else l2)

def linear_interval(terms):
 lo=F(0);hi=F(0)
 for coeff,x in terms:
  l,h=logbounds(x)
  if coeff>=0:lo+=coeff*l;hi+=coeff*h
  else:lo+=coeff*h;hi+=coeff*l
 return lo,hi

def decimal_bound(x,upper=False,places=30):
 scale=10**places;v=x.numerator*scale//x.denominator
 if upper and F(v,scale)<x:v+=1
 sg='-' if v<0 else '';v=abs(v)
 return sg+str(v//scale)+'.'+str(v%scale).zfill(places)

def run(item):
 idx,rec=item;start=time.monotonic();mass,inc,vol,masks=atomize(rec['boxes']);om=sum(mass);p=F(rec['p']);pi,msg=lf_weights(mass,inc,vol,masks,p);price=testing_price(pi,mass,inc,vol,p)
 assert price==F(rec['source_envelope']),(idx,price,rec['source_envelope'])
 source_f,source_ell=primal_source(mass,inc,vol,p,price)
 active=[i for i,v in enumerate(pi) if v>0];pir=[pi[i] for i in active];amat=[[mass[j]*inc[i][j]/vol[i] for j in range(len(mass))] for i in active];q=[a/om for a in mass]
 ar=np.array([[float(x) for x in row] for row in amat]);pr=np.array([float(x) for x in pir]);qr=np.array([float(x) for x in q]);k=len(pr)
 def fg(nu):
  pp=nu@ar;logs=np.log(pp/qr);ln=np.log(nu/pr);g=ar@logs-ln;value=float(np.dot(nu,g));return value,g
 if k==1:best=pr;solver='single active row';success=True
 else:
  starts=[pr,np.full(k,1/k)]+[.98*np.eye(k)[i]+.02/k for i in range(k)];out=[]
  for x0 in starts:
   rr=minimize(lambda z:-fg(z)[0],x0,jac=lambda z:-fg(z)[1],bounds=[(1e-12,1)]*k,constraints=[{'type':'eq','fun':lambda z:sum(z)-1,'jac':lambda z:np.ones(k)}],method='SLSQP',options={'ftol':1e-13,'maxiter':2000})
   if np.isfinite(rr.fun):out.append(rr)
  rr=min(out,key=lambda z:z.fun);best=rr.x;solver=str(rr.message);success=bool(rr.success)
 assert fg(best)[0]>=fg(pr)[0]-1e-8
 counts=[max(1,int(round(z*10**10))) for z in best];nu=[F(z,sum(counts)) for z in counts];pn=[sum(nu[i]*amat[i][j] for i in range(k)) for j in range(len(mass))]
 fterms=[(pn[j],pn[j]/q[j]) for j in range(len(q))]+[(-nu[i],nu[i]/pir[i]) for i in range(k)];flo,fhi=linear_interval(fterms)
 gis=[linear_interval([(amat[i][j],pn[j]/q[j]) for j in range(len(q))]+[(F(-1),nu[i]/pir[i])]) for i in range(k)];upper=max(h for l,h in gis);assert upper>=flo
 pl=[sum(pir[i]*amat[i][j] for i in range(k)) for j in range(len(q))];klo,khi=linear_interval([(pl[j],pl[j]/q[j]) for j in range(len(q))]);vlo,vhi=logbounds(om/price)
 result=dict(case=idx,dimension=rec['dimension'],target=str(p),boxes=rec['boxes'],atom_masks=masks,atom_masses=[str(v) for v in mass],row_volumes=[str(v) for v in vol],union_volume=str(om),pi=[str(v) for v in pi],active_rows=active,nu=[str(v) for v in nu],exact_Hall=True,exact_NP_price=str(price),price_matches_r716=True,source_f=[str(x) for x in source_f],source_ell=[str(x) for x in source_ell],exact_primal_dual_price=True,log_series_terms=TERMS,capacity_lower=decimal_bound(flo),capacity_upper=decimal_bound(upper,True),certified_capacity_gap_upper=decimal_bound(upper-flo,True),KL_pi_lower=decimal_bound(klo),KL_pi_upper=decimal_bound(khi,True),test_minus_capacity_lower=decimal_bound(vlo-upper),test_minus_capacity_upper=decimal_bound(vhi-flo,True),solver_status=solver,solver_success=success,elapsed_seconds=time.monotonic()-start,worker_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 return result
if __name__=='__main__':
 start=time.monotonic();records=json.loads(SOURCE.read_text())['results']
 with ProcessPoolExecutor(max_workers=12) as pool:rows=list(pool.map(run,enumerate(records)))
 out=dict(status='complete',scope='27 finite responsible-LF catalogues from r716; not near-critical dimensional asymptotics',source='../r716/wolfram.json',source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),workers=12,blas_threads=1,log_series_terms=TERMS,certificate_method='exact Fraction atanh-log series with geometric positive remainder; candidate nu rational; concavity upper max gradient',versions=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__),wall_seconds=time.monotonic()-start,results=rows)
 (HERE/'capacity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='complete',cases=len(rows),wall_seconds=out['wall_seconds'],max_capacity_interval_width=max(float(r['certified_capacity_gap_upper']) for r in rows),solver_non_success=sum(not r['solver_success'] for r in rows))))
