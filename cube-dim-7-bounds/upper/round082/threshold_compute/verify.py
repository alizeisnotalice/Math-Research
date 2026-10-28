#!/usr/bin/env python3
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,os,platform,sys,time
sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent; reg=json.loads((p/'registration.json').read_text()); src=Path(reg['proof_path'])
assert hashlib.sha256(src.read_bytes()).hexdigest()==reg['proof_sha256']
def logs(v):
 assert v>1
 z=(v-1)/(v+1); power=z; total=F(0)
 for j in range(64): total+=2*power/(2*j+1); power*=z*z
 return total,total+2*power/(129*(1-z*z))
def S(y):
 assert y>=0
 term=total=F(1)
 for j in range(1,17): term*=y/j; total+=term
 return total
L=F(reg['L0']); x=2*L/5; n=L**16/factorial(16); h=1/L**2; e=F(1,10000); sig=F(3,5); chi=F(7,10)
qmin=F(1,4); qmax=F(99,100); delta=F(1,100000); rlo=1+delta; rhi=F(99,25)-delta
P0=1-2*e; Ulo=P0*(1-e)/(1+e); Uhi=(1+e)/(P0*(1-e)); M=2*(1+e)*Uhi
ret=F(999,1000); b=2*ret*(1-e)*Ulo-2*e; mL=F(137,100); Gamma=mL*b/M; gamma=F(13,10); c=F(3,8); price=F(23,25)
nu=F(1,4)-F(5,4)*e; eta=F(1,198); y=(nu-1/n)*eta*x/(1+e); gap=F(29,1000)*x-F(1,2)
rho=n/(4*x); mean=576*x*x/n+48*x/n; var=2304*x*x/n
cap=M/S(y); tau=rhi*L**2*4096/S(L/50); root=1/n+1/S(sig*L)+8192/S(sig*L/2)
masslo,masshi=logs(rhi/rlo); pricelo,pricehi=logs(F(5,2)); qlo,qhi=logs(F(8,5))
gates={
 'registration_L0_and_sqrt':L==2**24 and F(4096)**2==L,
 'n_lower_ge_300':n>=300,'x_ge_2':x>=2,'beta_lt_2':x/n<1,'kappa_gt_1':x*n*n>1,
 'matching_trim':2*h<delta,'h_below_001':h<=F(1,1000),'active_R_above_one':rlo-2*h>1,
 'full_R_positive':4*qmin-h>0,'g_positive':3*qmin-h>0,'eps_positive_small':1/(x*n*n)<1,
 'eps_q_below_h_over_4':qmax/(x*n*n)<=h/4,'t_above_target':n/x>rhi,'source_above_4rhi':rho>4*rhi,
 'source_above_12':rho>=12,'target_source_separation':4*rhi>5,
 'relative_gap_log8_5':qhi<F(1,2),'logprod_payment':80/L<=e,
 'individual_J_error':6/rho<=e,'individual_J_local_guard':6/rho<=F(1,2),
 'baseline_mean':mean<=e/2,'baseline_chebyshev':4*var/e**2<=e/8,
 'count_minus_probability':4/(e*e*n)<=e/8,'count_plus_probability':4/(e*e*n)<=e/8,
 'target_cap_fraction':1+x/(2*n)-1/(2*rhi)<=F(9,10),
 'target_cap_probability':F(9,10)**299<=e/8,
 'profile_mu_upper':qlo>F(47,100),'first_profile_pin':4/n<=F(1,200),
 'first_profile_probability':F(640000)/n<=e/8,'count_whole_matching':2*h+1/n<=3*e/4,
 'profile_transport':80*h<=F(1,100),'profile_mean_rounding':F(67,100)/(1-e)<F(671,1000),
 'profile_variance_rounding':F(4,3)/(1-e)<2,'chebyshev_positive_gap':gap>0,
 'retention_chebyshev':2*x/(gap*gap)<=F(1,1000),
 'weighted_cap_nu':nu>0 and nu-1/n>0,'weighted_cap_positive_domain':y>=4,'weighted_cap_payment':cap<=e,
 'radial_width':1/(1-(h/rlo)**2)<=1+e,'tau_positive_rate':sig/2-chi*(1-sig)==F(1,50),
 'tau_prefix_tail_domain':L/50>=4,'tau_payment':tau<=e,
 'logt_lower':S(L/2)>=x,'B_minus_gamma_positive':L/4>gamma,
 'target_mass':(1-e)*masslo>mL,'ambient_degree':P0>0 and M>2 and b>0,
 'gamma_full_deficit':(1-e)*Gamma>gamma,'price_log_upper':pricehi<price,
 'root_prefix_domains':sig*L/2>=2,'root_full_payment':root<F(1,200),
 'scalar_slack_exact':gamma-price-c==F(1,200),
 'finite_history_dilation':n>=2 and x/n<1 and n/(4*x)>5,
 'finite_history_displacement':F(1,1)/n**3-F(4,1)/(8*n**4)>0,
 'finite_history_norm_guard':rhi<4 and 4*rhi>5,
}
bounds={'N_lower':n,'x':x,'h':h,'rho_lower':rho,'baseline_mean_upper':mean,'baseline_variance_upper':var,'baseline_failure_upper':4*var/e**2,'count_failure_upper':4/(e*e*n),'first_failure_upper':F(640000)/n,'chebyshev_gap':gap,'source_face_bad_upper':2*x/(gap*gap),'weighted_cap_exponent_lower':y,'weighted_cap_payment_upper':cap,'tau_payment_upper':tau,'root_payment_upper':root,'mass_log_lower':masslo,'mL':mL,'Ulo':Ulo,'Uhi':Uhi,'M':M,'b':b,'Gamma':Gamma,'gamma_budget':(1-e)*Gamma,'gamma':gamma,'price_upper_enclosure':pricehi,'price_budget':price,'final_c':c,'scalar_slack':F(1,200)}
out={'status':'PASS' if all(gates.values()) else 'FAIL','scope':'Exact base-envelope arithmetic only, conditional on frozen proof; does not audit operator chain or establish all-tail by itself.','pid':os.getpid(),'python':platform.python_version(),'L0':str(L),'registration_sha256':hashlib.sha256((p/'registration.json').read_bytes()).hexdigest(),'proof_sha256':reg['proof_sha256'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gates':gates,'failed':[k for k,v in gates.items() if not v],'rational_bounds':{k:str(v) for k,v in bounds.items()},'display':{k:float(v) for k,v in bounds.items()},'cpu_seconds':time.process_time()}
(p/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'gate_count':len(gates),'failed':out['failed'],'display':out['display'],'cpu_seconds':out['cpu_seconds']},indent=2))
