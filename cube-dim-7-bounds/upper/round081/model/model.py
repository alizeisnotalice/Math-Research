#!/usr/bin/env python3
from fractions import Fraction as F
from decimal import Decimal as D, getcontext
from pathlib import Path
import hashlib,json,os,sys,time
sys.set_int_max_str_digits(0); getcontext().prec=80
p=Path(__file__).resolve().parent; reg=json.loads((p/'registration.json').read_text()); calls=0

def log_bounds(v):
 global calls
 calls+=1
 assert v>=1
 z=(v-1)/(v+1); acc=F(0); power=z
 for j in range(64):
  acc+=2*power/(2*j+1); power*=z*z
 return acc,acc+2*power/(129*(1-z*z))
def decimal(v):return D(v.numerator)/D(v.denominator)
rows=[]
for raw in reg['points']:
 s=F(raw); tlo,thi=log_bounds(2*s/(s+1))
 mu_lo=s-1-(s+1)*thi; mu_hi=s-1-(s+1)*tlo
 assert mu_lo>0
 lslo,lshi=log_bounds(s); lplo,_=log_bounds(1+2*mu_lo); _,lphi=log_bounds(1+2*mu_hi)
 clo=lslo-lphi; chi=lshi-lplo
 dlo=1-2*thi; dhi=1-2*tlo
 rows.append({'s':str(s),'mu_lower':str(mu_lo),'mu_upper':str(mu_hi),'c_lower':str(clo),'c_upper':str(chi),'c_display':str((decimal(clo)+decimal(chi))/2),'derivative_numerator_lower':str(dlo),'derivative_numerator_upper':str(dhi),'derivative_sign':('positive' if dlo>0 else 'negative' if dhi<0 else 'undetermined')})
assert calls<=200 and rows[0]['derivative_sign']=='positive' and rows[2]['derivative_sign']=='negative'
e=D(1).exp(); sstar=1/(2/e.sqrt()-1); mu=sstar-1+(sstar+1)*((sstar+1)/(2*sstar)).ln(); c=(sstar/(1+2*mu)).ln()
out={'status':'PASS exact rational series enclosures and Decimal stationary computation; model only','pid':os.getpid(),'python':sys.version,'precision_digits':getcontext().prec,'function_calls_log_bounds':calls,'registration_sha256':hashlib.sha256((p/'registration.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sstar_decimal':str(sstar),'mu_star_decimal':str(mu),'cstar_decimal':str(c),'rows':rows,'cpu_seconds':time.process_time(),'stopping':'5 points; 20 log-enclosure calls; no further optimization'}
(p/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'calls':calls,'sstar':str(sstar),'mu':str(mu),'cstar':str(c),'points':[{k:r[k] for k in ['s','c_display','derivative_sign']} for r in rows],'cpu_seconds':out['cpu_seconds']},indent=2))
