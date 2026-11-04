from fractions import Fraction as F
from math import comb
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path
import argparse,json,time

def out(x,mode):
 with localcontext() as c:c.prec=30;c.rounding=mode;return str(Decimal(x.numerator)/Decimal(x.denominator))
def calc(n,j,C):
 u=F(1,2**j);g=1/(1+2*u);a=(2*u)**n;B=comb(n,2);bar=F(8*n,3);k=bar*g*g;e=k*u**3;lamstar=B*k**3/(8*n*n*g);beta=F(5*B,96*n**3)*bar**4;lam=lamstar-beta*u+C*u*u;w=u*g+lam*u**6;rho=1-2*w;K=a/(2*w)**(n-2);f=n*K*rho**(n-1);h=e*e/f;delta=h/e;zmax=1-2*delta;fmax=f*(1+h/rho)**(n-1);D=(2*w)**n-a*rho**n;H=(3**n)*h
 integral=zmax-zmax**3/6-17*delta*zmax
 Rec=4*B*(2*w-4*e)**(n-2)*e*e*h/fmax*integral
 B0=2*B*(2*w)**(n-2)*e**3/f
 if n>=3:B0+=8*n*comb(n-1,2)*H*(e+H)**2*(2*w)**(n-3)
 slack=D-B0+Rec;guards={'n_ge_2':n>=2,'C_u2_le_one':C*u*u<=1,'lambda_bounded':abs(lam)<=B*bar**3/(8*n*n)+beta+2,'w_positive_lt_quarter':0<w<F(1,4),'e_lt_w_over_4':0<e<w/4,'rho_positive':rho>2*w,'f_positive':f>0,'delta_lt_1_over_35':0<delta<F(1,35),'zmax_positive':zmax>0,'integral_positive':integral>0,'fmax_relative_increment_le_delta':fmax/f-1<=delta,'tangential_positive':2*w-4*e>0,'H_lt_w_minus_e':H<w-e,'recovery_nonnegative':Rec>=0,'budget_strict':slack>0}
 return {'n':n,'u':str(u),'j':j,'C':C,'lambda':str(lam),'guards':guards,'status':'PASS' if all(guards.values()) else 'FAIL','slack':[out(slack,ROUND_FLOOR),out(slack,ROUND_CEILING)],'recovery':[out(Rec,ROUND_FLOOR),out(Rec,ROUND_CEILING)],'note':'Exact arithmetic for supplied analytic lower/upper bounds; no independent geometric proof.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--j',type=int,default=30);ap.add_argument('--C',type=int,default=1000000);ar=ap.parse_args();t=time.perf_counter();rows=[calc(n,ar.j,ar.C) for n in (2,3)];data={'cases':rows,'all_pass':all(x['status']=='PASS' for x in rows),'seconds':time.perf_counter()-t};Path('corner_finite.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'all_pass':data['all_pass'],'statuses':[(x['n'],x['status']) for x in rows],'seconds':data['seconds']}))
