"""Directed rational Jensen upper bounds on the one-dimensional Lorenz moments."""
import argparse,functools,hashlib,json,math,multiprocessing as mp,sys,time
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from sliver79 import logi
D=10**24

def iroot(n,k):
 if n==0:return 0
 if k==2:return math.isqrt(n)
 x=1<<((n.bit_length()+k-1)//k)
 while True:
  y=((k-1)*x+n//pow(x,k-1))//k
  if y>=x:break
  x=y
 while pow(x+1,k)<=n:x+=1
 while pow(x,k)>n:x-=1
 return x

def power_upper(a,power):
 assert a>=0 and 0<power<1
 k,q=power.numerator,power.denominator
 t=pow(a.numerator,k)*pow(D,q)//pow(a.denominator,k)
 n=iroot(t,q)+1
 assert pow(n,q)*pow(a.denominator,k)>=pow(a.numerator,k)*pow(D,q)
 return Q(n,D)
@functools.lru_cache(None)
def primitive(v):
 if not v:return Q(0),Q(0)
 lo,hi=logi(v);base=v+v*v/2
 return base-2*v*hi,base-2*v*lo

def job(spec):
 r0,m,mode=spec;r=Q(r0);total=Q(0);prev=Q(0);rows=[]
 for j in range(1,m+1):
  v=Q(j,m)**(2 if mode=='squared' else 1);h=v-prev;fhi=primitive(v)[1];flo=primitive(prev)[0];avg=(fhi-flo)/h;pu=power_upper(avg,r-1);total+=h*pu
  if m==16 and r==Q(3,2) and mode=='uniform':rows.append(dict(a=str(prev),b=str(v),primitive_difference_upper=str(fhi-flo),mean_upper=str(avg),power_upper=str(pu)))
  prev=v
 upper=r*total;mu_lo=Q('1.2784645427610');mu_hi=Q('1.2784645427612');assert logi(mu_lo-1)[1]+mu_lo<0 and logi(mu_hi-1)[0]+mu_hi>0
 log_upper=logi(upper)[1];coef=((r-1)*mu_lo-log_upper)/r
 return dict(r=str(r),mesh=m,mode=mode,J_upper=str(upper),log_J_upper=str(log_upper),mu_lower=str(mu_lo),mu_upper=str(mu_hi),information_coefficient_lower=str(coef),coefficient_lower_float=float(coef),J_upper_float=float(upper),certificate_rows=rows,all_root_upper_checks=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('output');a=ap.parse_args();start=time.monotonic();out=Path(a.output);out.mkdir(parents=True,exist_ok=True);specs=[(str(r),n,mode) for r in (Q(3,2),Q(153,100),Q(77,50),Q(31,20),Q(39,25),Q(157,100),Q(8,5)) for n in (16,64,256,1024) for mode in ('uniform','squared')]
 with mp.get_context('spawn').Pool(12) as pool:rows=list(pool.imap_unordered(job,specs))
 rows.sort(key=lambda x:(Q(x['r']),x['mesh'],x['mode']));base=next(x for x in rows if x['r']=='3/2' and x['mesh']==16 and x['mode']=='uniform');assert Q(base['information_coefficient_lower'])>Q(3,25);best=max(rows,key=lambda x:Q(x['information_coefficient_lower']));assert Q(best['information_coefficient_lower'])>Q('0.1232')
 z=dict(pass_all=True,cases=len(rows),workers=12,base_certificate=base,best_certificate=best,rows=rows,seconds=time.monotonic()-start,scope='Directed rational Jensen certificates for the Lorenz envelope information lower bound; no unrestricted physical endpoint upper or finite critical-dimension certification.');(out/'information96.json').write_text(json.dumps(z,indent=2)+'\n');print(json.dumps({k:v for k,v in z.items() if k not in ['rows','base_certificate','best_certificate']}));print('BEST',best['r'],best['mesh'],best['mode'],best['coefficient_lower_float'])
if __name__=='__main__':main()
