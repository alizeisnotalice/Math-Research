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
L=F(reg['L0']); x=F(387,1000)*L; n=L**16/factorial(16); h=1/L**2; e=F(1,100000); sig=F(613,1000); chi=F(1583,2000)
s=F(9,2); qmin=F(2,9); qmax=F(9999,10000); delta=F(1,1000000); rlo=1+delta; rhi=F(4499549,1000000)
P0=1-2*e; Ulo=P0*(1-e)/(1+e); Uhi=(1+e)/(P0*(1-e)); M=2*(1+e)*Uhi
ret=F(99999,100000); b=2*ret*(1-e)*Ulo-2*e; mL=F(1503,1000); Gamma=mL*b/M; gamma=F(1501,1000); c=F(11,20); price=F(19,20)
nu=F(22221,100000); eta=F(1,19998); y=(nu-1/n)*eta*x/(1+e); gap=F(3,100000)*x-F(3,5)
rho=n/(s*x); mean=49*s*s*x*x/n+14*s*x/n; var=196*s*s*x*x/n
cap=M/S(y); sqrtL=F(reg['sqrt_L']); tau=rhi*L**2*sqrtL/S(F(379,2000000)*L); root=1/n+1/S(sig*L)+2*sqrtL/S(sig*L/2)
masslo,masshi=logs(rhi/rlo); pricelo,pricehi=logs(F(1000,387)); qlo,qhi=logs(F(18,11)); muhi=F(7,2)-F(11,2)*qlo
expected_count=F(1,2)+(F(1,2)*4-F(5,2))*4**0 # unused
# exact polynomial lower bound for y S16'-(5/2)S16: negative degrees 0,1,2, positive degree 3,4.
y0=F(4); tau_poly=-F(5,2)-F(3,2)*y0-F(1,2)*y0**2+F(1,2)*y0**3/F(factorial(3))+F(3,2)*y0**4/F(factorial(4))
gates={
 'registration_L0_and_sqrt':L==2**50 and sqrtL**2==L,
 'n_lower_ge_506250':n>=506250,'x_ge_2':x>=2,'beta_lt_2':x/n<1,'kappa_gt_1':x*n*n>1,
 'matching_trim':2*h<delta,'active_R_above_one':rlo-2*h>1,'full_R_positive':s*qmin-h>0,'g_positive':(s-1)*qmin-h>0,
 'eps_positive_small':1/(x*n*n)<1,'eps_q_below_h_over_s':qmax/(x*n*n)<=h/s,
 't_at_least_50':n/x>=50,'t_above_target':n/x>rhi,'source_above_4rhi':rho>4*rhi,
 'individual_J_local_guard':7/rho<=F(1,2),'individual_J_error':7/rho<=e,
 'source_target_separation':4*rhi>6 and rhi<5,
 'relative_gap_log18_11':qhi<F(1,2),'logprod_payment':180*x*h<=e,
 'baseline_mean':mean<=e/2,'baseline_chebyshev':4*var/e**2<=e/8,
 'count_minus_probability':4/(e*e*n)<=e/8,'count_plus_probability':4/(e*e*n)<=e/8,
 'target_cap_fraction':1+x/(2*n)-1/(2*rhi)<=F(9,10), 'target_cap_probability':F(9,10)**299<=e/8,
 'first_profile_mu':muhi<F(7914,10000),'first_profile_pin':F(81,16)/n<=F(1,100000),
 'first_profile_probability':F(64072265625)/n<=e/8,'count_whole_matching':2*h+1/n<=3*e/4,
 'profile_transport':80*h<=F(1,100000),'profile_mean_rounding':F(79144,100000)/(1-e)<F(79147,100000),
 'profile_variance_rounding':F(27,16)/(1-e)<2,'chebyshev_positive_gap':gap>0,'retention_chebyshev':2*x/(gap*gap)<=e,
 'weighted_cap_nu':nu==1/s-(1+1/s)*e and nu-1/n>0,'weighted_cap_eta':eta==(1-qmax)/(2*qmax),
 'weighted_cap_positive_domain':y>=4,'weighted_cap_payment':cap<=e,
 'radial_width':1/(1-(h/rlo)**2)<=1+e,'tau_positive_rate':sig/2-chi*(1-sig)==F(379,2000000),
 'tau_prefix_tail_domain':F(379,2000000)*L>=4,'tau_four_term_poly_positive_at_4':tau_poly>0,'tau_payment':tau<=e,
 'logt_lower':S(L/2)>=x,'B_minus_gamma_positive':L/4>gamma,
 'target_mass':(1-e)*masslo>mL,'ambient_degree':P0>0 and M>2 and b>0,
 'gamma_full_deficit':(1-e)*Gamma>gamma,'price_log_upper':pricehi<price,
 'root_prefix_domains':sig*L/2>=2,'root_full_payment':root<F(1,1000),
 'scalar_slack_exact':gamma-price-c==F(1,1000),
 'finite_history_dilation':n>=2 and x/n<1 and n/(s*x)>6,
 'finite_history_displacement':1/n**3-F(4,1)/(8*n**4)>0,'finite_history_norm_guard':1+rhi<6 and 4*rhi>6,
}
bounds={'N_lower':n,'x':x,'h':h,'rho_lower':rho,'baseline_mean_upper':mean,'baseline_variance_upper':var,'baseline_failure_upper':4*var/e**2,'count_failure_upper':4/(e*e*n),'first_failure_upper':F(64072265625)/n,'chebyshev_gap':gap,'source_face_bad_upper':2*x/(gap*gap),'weighted_cap_exponent_lower':y,'weighted_cap_payment_upper':cap,'tau_payment_upper':tau,'root_payment_upper':root,'mass_log_lower':masslo,'mL':mL,'M':M,'b':b,'Gamma':Gamma,'gamma_budget':(1-e)*Gamma,'gamma':gamma,'price_upper_enclosure':pricehi,'price_budget':price,'mu_upper':muhi,'final_c':c,'scalar_slack':F(1,1000),'tau_four_term_poly_at_4':tau_poly}
out={'status':'PASS' if all(gates.values()) else 'FAIL','scope':'Exact base arithmetic conditional on frozen R083 proof; no Cn theorem or full proof audit.','pid':os.getpid(),'python':platform.python_version(),'L0':str(L),'registration_sha256':hashlib.sha256((p/'registration.json').read_bytes()).hexdigest(),'proof_sha256':reg['proof_sha256'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gates':gates,'failed':[k for k,v in gates.items() if not v],'rational_bounds':{k:str(v) for k,v in bounds.items()},'display':{k:float(v) for k,v in bounds.items()},'cpu_seconds':time.process_time()}
(p/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'gate_count':len(gates),'failed':out['failed'],'display':out['display'],'cpu_seconds':out['cpu_seconds']},indent=2))
