"""Independent direct one-dimensional moment check at two already-computed triples."""
import json,time,os,hashlib,resource
from pathlib import Path
from fractions import Fraction as Q
import mpmath as mp
mp.mp.dps=250;D=Path(__file__).resolve().parent;start=time.monotonic();cpu=time.process_time();checks=[]
for mm in (32,512):
 source=D/('case_m'+str(mm)+'_d256.json');doc=json.loads(source.read_text());row=next(r for r in doc['results'] if r['u_over_d']=='1/4' and r['theta1']=='-1' and r['theta2']=='1');n=mp.mpf(2)**mm;d=mp.mpf(1)/256;x=mm*mp.log(2)/2;t=n/x;h=mp.mpf(2)**-40/(x*d);r=1+d/4
 def src(theta):
  q=r-d+theta*h;rho=t*q;R=q+d-h;g=d-h+rho/n**3;a=R-g;ell=1-g/(2*R)
  def E(m):
   z=g-m
   return mp.mpf(0) if z<=0 else -2*R*mp.log1p(-z/(2*R))-z
  def J(w):
   m=max(mp.mpf(0),abs(w)-a);return 1+(E(m)+E(0)-m)/(2*rho)
  def S(w):return (1+int(abs(w)<a))/ell
  return J,S,a,R
 J1,S1,a1,R1=src(-1);J2,S2,a2,R2=src(1);cuts=sorted(set([mp.mpf(0),a1,a2,R1,R2,r]))
 def ave(f):return mp.fsum(mp.quad(f,[a,b]) for a,b in zip(cuts,cuts[1:]) if b>a)/r
 M=ave(lambda w:J1(w)*J2(w));A=ave(lambda w:S1(w)*J2(w));B=ave(lambda w:S2(w)*J1(w));C=ave(lambda w:S1(w)*S2(w));j1,j2,s1,s2=J1(r),J2(r),S1(r),S2(r);k=n-1
 value=(s1*s2*M**k+k*(s1*j2*B+s2*j1*A+j1*j2*C)*M**(k-1)+k*(k-1)*j1*j2*A*B*M**(k-2))/(4*n*n)
 err=abs(value-mp.mpf(row['S_cross_over_n2']));assert err<mp.mpf('1e-30')
 checks.append({'m':mm,'d':'1/256','u_over_d':'1/4','theta1':'-1','theta2':'1','direct_value':mp.nstr(value,40),'recorded':row['S_cross_over_n2'],'absolute_difference':mp.nstr(err,10),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()})
out={'status':'PASS','precision_dps':250,'checks':checks,'scope':'independent moment integration at two previously computed triples; numerical verification, not an interval proof','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'pid':os.getpid(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(D/'check_moments.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
