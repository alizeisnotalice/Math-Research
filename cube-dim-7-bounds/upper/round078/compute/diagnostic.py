#!/usr/bin/env python3
from fractions import Fraction as F
from decimal import Decimal as D, getcontext
from pathlib import Path
import hashlib,json,os,platform,time
getcontext().prec=80
here=Path(__file__).resolve().parent
conv=lambda x:D(x.numerator)/D(x.denominator)
x=F(7392); chi=F(63,1280); p=F(183,1000); eJ=F(1,2**20)
m=F(2211,50000)/(1-eJ); v=F(531,100000); u=(chi-m)*x-p; V=v*x; bound=p
E=u*u/(2*(V+bound*u/3))
term=total=F(1)
for i in range(1,65):
 term*=E/i; total+=term
assert u>0 and total>10**6
fail_bennett=(-(conv(V)/conv(bound)**2)*((D(1)+conv(bound*u/V))*(D(1)+conv(bound*u/V)).ln()-conv(bound*u/V))).exp()
eB=eN=eT=F(1,2**19); eR=eS=eD=F(1,2**20); mass=F(3043,10000)
P0=1-eB-F(6,768); Ulo=P0*(1-eN)/(1+eJ); Uhi=(1+eN)/(P0*(1-eJ)); M=2*(1+eR)*Uhi
price=F(257,2500); price_exact=(D(256)/D(231)).ln()
fails={'R073_fixed':F(283,10000),'cantelli_unrounded':V/(V+u*u),'bernstein_fixed':F(1,10**6),'bernstein_S64':1/total}
rows={}
for name,fail in fails.items():
 ret=1-fail; row=2*ret*(1-eR)*Ulo-eS-eT; gamma=mass*row/M; budget=(1-eD)*gamma
 rows[name]={k:str(conv(val)) for k,val in {'tail_upper':fail,'ret_lower':ret,'row_b':row,'M':M,'Gamma':gamma,'gamma_after_eD':budget,'price_upper_R073':price,'gamma_minus_price_upper':budget-price}.items()}
 rows[name]['gamma_minus_exact_price']=str(conv(budget)-price_exact)
 rows[name]['exact_gamma_after_eD']=str(budget)
result={'scope':'Conditional diagnostic; Fraction Bernstein arithmetic is exact, no Cn claim','pid':os.getpid(),'python':platform.python_version(),'input_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [here/'registration.json',here/'amendment.json',Path(__file__)]},'gap_exact':str(u),'gap_decimal':str(conv(u)),'variance_exact':str(V),'increment_upper_exact':str(bound),'bernstein_exponent_exact':str(E),'bernstein_exponent_decimal':str(conv(E)),'S64_lower_exp_exact':str(total),'S64_failure_upper':str(conv(1/total)),'bennett_failure_decimal_not_certified':str(fail_bennett),'rows':rows,'tail_monotonicity':'u/x=a-p/x increases, so Bernstein exponent x*(a-p/x)^2/(2*(v+b*(a-p/x)/3)) increases; Bennett exponent x*v/b²*h(b*(a-p/x)/v) increases.','stop':'fixed batch complete; no grid','cpu_seconds':time.process_time()}
(here/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'pid':os.getpid(),'gap':str(conv(u)),'E':str(conv(E)),'prefix_failure':str(conv(1/total)),'bennett_failure':str(fail_bennett),'rows':rows},indent=2))
