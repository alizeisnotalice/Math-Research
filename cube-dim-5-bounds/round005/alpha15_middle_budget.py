from pathlib import Path
import hashlib,json
from flint import arb,ctx
ctx.prec=256
alpha=arb(3)/2;q=arb(10)/7
A=(arb(2)/3).gamma()/(alpha*arb.pi())
K=arb(25)/6*(arb(10)/(3*arb.pi())).sqrt()
Cf=K*arb(16)**(arb(5)/2)*A**4/3
B=arb(5)/6*Cf**(arb(3)/10)*K**(arb(7)/10)
# det48 five-direction density upper bound; remainder is a probability convolution.
peak=A**5/48
I0=arb(11)/1000000;norm_limit=I0**(arb(7)/10)
rows=[]
for denominator in (4096,16384,65536):
 for T in (1,4):
  eps=arb(1)/denominator
  short=B*eps**(arb(6)/5)
  long=arb(3)/7*peak*arb(T)**(-arb(7)/3)
  remain=norm_limit-short-long
  rows.append({'epsilon':'1/'+str(denominator),'T':T,'short_Lq_upper_ball':str(short),'long_Lq_upper_ball':str(long),'remaining_middle_Lq_limit_ball':str(remain),'positive_remaining_budget':bool(remain>0),'middle_integral_limit_ball':str(remain**q) if remain>0 else None})
out={'status':'conditional_middle_integral_budget_only_not_C5_bound','g_target':'1/10000','I_total_target':'11/1000000','q':'10/7','peak_upper_formula':'A^5/48','long_time_bound':'(3/7)*(A^5/48)*T^(-7/3)','rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
for r in rows:print(r['epsilon'],r['T'],r['middle_integral_limit_ball'])
