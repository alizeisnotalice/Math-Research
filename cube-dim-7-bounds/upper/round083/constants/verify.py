#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,platform,sys,time
sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent; r=json.loads((p/'registration.json').read_text()); v=lambda k:F(r[k])
def log_bounds(y):
 assert y>1
 z=(y-1)/(y+1); power=z; total=F(0)
 for j in range(64): total+=2*power/(2*j+1); power*=z*z
 return total,total+2*power/(129*(1-z*z))
s=v('s'); qmin=v('qmin'); qmax=v('qmax'); sig=v('sigma'); chi=v('chi'); delta=v('delta'); e=v('e')
ret=v('Ret'); mass=v('mL'); gamma=v('gamma'); price=v('price'); c=v('c'); slack=v('scalar_slack')
target=v('target_first_profile'); transported=v('transported_first_profile'); face=v('inverse_face_rate'); pin=v('pin')
rlo=1+delta; rhi=s*qmax-delta; P0=1-2*e; lo=P0*(1-e)/(1+e); hi=(1+e)/(P0*(1-e)); M=2*(1+e)*hi; b=2*ret*(1-e)*lo-2*e; Gamma=mass*b/M
loglo,loghi=log_bounds(2*s/(s+1)); mu_lo=s-1-(s+1)*loghi; mu_hi=s-1-(s+1)*loglo
masslo,masshi=log_bounds(rhi/rlo); pricelo,pricehi=log_bounds(1/(1-sig)); qlo,qhi=log_bounds(2*s/(s+1))
eta=(1-qmax)/(2*qmax); nu=qmin-(1+qmin)*e; tau_rate=sig/2-chi*(1-sig)
gates={'fixed_domain':s>1 and qmin==1/s and 0<qmax<1 and qmin<qmax and 0<sig<1 and 0<e<1/10,
 'target_band':1<rlo<rhi<s,'mu_lower_positive':mu_lo>0,'mu_below_7914':mu_hi<F(7914,10000),
 'mu_below_event_target':mu_hi<target,'target_to_transport_gap':target<transported,'transport_to_face':transported/(1-e)<face,
 'pin_log_bound':qhi<pin,'face_to_chi_gap':face<chi,'nu_positive':nu>0,'source_eta_positive':eta>0,
 'tau_rate_positive':tau_rate>0,'price_log_upper':pricehi<price,'mass_lower':(1-e)*masslo>mass,
 'positive_degrees':0<P0<1 and M>2 and b>0,'complete_paid_gamma':(1-e)*Gamma>gamma,
 'scalar_slack_exact':gamma-price-c==slack and slack>0,'scalar_exponential_rates':sig>0 and sig/2>0}
bounds={'mu_lower':mu_lo,'mu_upper':mu_hi,'target':target,'target_margin':target-mu_hi,'transported':transported,'face_mean_after_inverse':transported/(1-e),'face_rate':face,'face_rate_margin':face-transported/(1-e),'pin_log_upper':qhi,'pin':pin,'positive_gap_rate':chi-face,'tau_rate':tau_rate,'eta':eta,'nu':nu,'rhi':rhi,'target_mass_log_lower':masslo,'target_mass_after_error':(1-e)*masslo,'mass':mass,'P0':P0,'Ulo':lo,'Uhi':hi,'M':M,'b':b,'Gamma':Gamma,'paid_Gamma':(1-e)*Gamma,'gamma':gamma,'gamma_margin':(1-e)*Gamma-gamma,'price_log_upper':pricehi,'price':price,'price_margin':price-pricehi,'c':c,'slack':slack}
out={'status':'PASS' if all(gates.values()) else 'FAIL','scope':'Fixed constants only, no finite threshold and no Cn theorem. Profile, geometry, probability, finite history and operator transfer remain to be proved.', 'pid':os.getpid(),'python':platform.python_version(),'registration_sha256':hashlib.sha256((p/'registration.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gates':gates,'failed':[k for k,v in gates.items() if not v],'bounds_exact':{k:str(v) for k,v in bounds.items()},'bounds_display':{k:float(v) for k,v in bounds.items()},'cpu_seconds':time.process_time()}
(p/'result.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps({'status':out['status'],'pid':os.getpid(),'gates':len(gates),'failed':out['failed'],'bounds':out['bounds_display']},indent=2))
