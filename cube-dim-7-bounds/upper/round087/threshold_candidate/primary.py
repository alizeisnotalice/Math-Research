from fractions import Fraction as F
from pathlib import Path
import json,hashlib,platform,os,sys
sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent; r=json.loads((p/'primary_registration.json').read_text()); v=lambda k:F(r[k]); source=Path(r['proof_path'])
assert hashlib.sha256(source.read_bytes()).hexdigest()==r['proof_sha256']
def logs(a):
 assert a>1
 z=(a-1)/(a+1); w=z; total=F(0)
 for j in range(96):total+=2*w/(2*j+1);w*=z*z
 return total,total+2*w/(193*(1-z*z))
s=v('s');chi=v('chi');sig=v('sigma');q=v('qmax');gamma=v('gamma');c=v('c')
lo,hi=logs(2*s/(s+1)); mu_lo=s-1-(s+1)*hi; mu_hi=s-1-(s+1)*lo
bandlo,bandhi=logs(s*q);pricelo,pricehi=logs(1/(1-sig)); ct=sig/2-chi*(1-sig)
g={'s_range':4<=s<=5,'chi_above_mu':chi>mu_hi,'sigma_domain':0<sig<1,'tau_rate_positive':ct>0,'qmax_domain':1/s<q<1,'gamma_positive':0<gamma,'gamma_below_log_sqmax':gamma<bandlo,'c_below_gamma_minus_price':c<gamma-pricehi}
b={'mu_lower':mu_lo,'mu_upper':mu_hi,'chi_minus_mu_lower':chi-mu_hi,'tau_rate':ct,'log_sqmax_lower':bandlo,'gamma_margin':bandlo-gamma,'price_upper':pricehi,'cost_slack_lower':gamma-pricehi-c}
out={'status':'PASS' if all(g.values()) else 'FAIL','pid':os.getpid(),'python':platform.python_version(),'registration_sha256':hashlib.sha256((p/'primary_registration.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gates':g,'failed':[k for k,val in g.items() if not val],'bounds_exact':{k:str(val) for k,val in b.items()},'display':{k:float(val) for k,val in b.items()},'scope':'primary parameters only, not Cn or threshold'}
(p/'primary_result.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({'status':out['status'],'pid':out['pid'],'gates':g,'display':out['display']},indent=2))
