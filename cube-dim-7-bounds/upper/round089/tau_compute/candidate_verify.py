#!/usr/bin/env python3
from fractions import Fraction as F
from decimal import Decimal as D, getcontext
from math import factorial
from pathlib import Path
import hashlib,json,os,platform,sys,time
sys.set_int_max_str_digits(0);getcontext().prec=30
p=Path(__file__).resolve().parent; r=json.loads((p/'candidate_registration.json').read_text())
v=lambda k:F(r[k]); w=v

def logs(y):
 assert y>1
 z=(y-1)/(y+1); power=z; total=F(0)
 for j in range(96):total+=2*power/(2*j+1);power*=z*z
 return total,total+2*power/(193*(1-z*z))
def S(y):
 assert y>=0
 term=total=F(1)
 for j in range(1,17):term*=y/j;total+=term
 return total
def disp(f): return str(D(f.numerator)/D(f.denominator))
L=w('L0'); sqrtL=w('sqrt_L'); N=L**16/factorial(16); s=v('s'); chi=v('chi'); sig=v('sigma'); qmax=v('qmax'); gamma=v('gamma'); c=v('c')
delta=w('trim');e=w('e');ret=w('ret'); mass=w('mass'); face=w('m_face');pin=F(3,5);price=w('price');slack=w('scalar_allowance')
x=(1-sig)*L; h=1/L**2; qmin=1/s; rlo=1+delta;rhi=s*qmax-delta;rho=N/(s*x)
qlo,qhi=logs(2*s/(s+1));mulo=s-1-(s+1)*qhi;muhi=s-1-(s+1)*qlo; dlo=(chi-muhi)/4;dhi=(chi-mulo)/4
masslo,masshi=logs(rhi/rlo);pricelo,pricehi=logs(1/(1-sig));P0=1-2*e; Ulo=P0*(1-e)/(1+e); Uhi=(1+e)/(P0*(1-e));M=2*(1+e)*Uhi;b=2*ret*(1-e)*Ulo-2*e;Gamma=mass*b/M
nu=1/s-(1+1/s)*e;eta=(1-qmax)/(2*qmax);weighted=(nu-1/N)*eta*x/(1+e);facegap=(chi-face)*x-pin; dt=w('delta_tau')
mean=81*s*s*x*x/N+18*s*x/N;var=324*s*s*x*x/N;cap=M/S(weighted);ct=dt-chi*(1-sig);tau=rhi*L**2/S(ct*L);root=1/N+1/S(dt*L)+1/S(sig*L)+L/S((sig-dt)*L)
# coefficients for lower polynomial y S16' - 5/2 S16; deliberately lower quadratic.
y=F(4); tau_poly=-F(5,2)-F(3,2)*y-F(1,2)*y*y+y**3/12+y**4/16
gates={'sqrt_L_exact':sqrtL**2==L,'N_at_least_400':N>=400,'x_at_least_2':x>=2,'beta_lt_2':x/N<1,'kappa_gt_1':x*N*N>1,
 'trim_band':2*h<delta and rlo<rhi,'active_R_above_1':rlo-2*h>1,'full_R_positive':s*qmin-h>0,'g_positive':(s-1)*qmin-h>0,
 'epsilon_small':1/(x*N*N)<1,'relative_gap':qmax/(x*N*N)<=h/s,'t_at_least_2rhi':N/x>=2*rhi,
 'source_radius_gt20':rho>20,'source_radius_gt_4rhi':rho>4*rhi,'rhi_gt1_lt5':1<rhi<5,
 'local_pin_bound':qhi<pin,'individual_J_local':9/rho<=F(1,2),'individual_J_error':9/rho<=e,'logprod_bridge':240/L<=e,
 'baseline_mean':mean<=e/2,'baseline_failure':4*var/e**2<=e/8,'count_fail_minus':4/(e*e*N)<=e/8,'count_fail_plus':4/(e*e*N)<=e/8,
 'target_cap_fraction':1+x/(2*N)-1/(2*rhi)<=F(19,20),'target_cap_probability':F(19,20)**399<=e/8,
 'mu_gap_positive':dlo>0,'first_profile_pin':10/N<=dlo/2,'first_profile_bad':400/(dlo*dlo*N)<=e/8,
 'whole_matching_count':2*h+1/N<=3*e/4,'profile_transport':120*h<=dlo,
 'face_upper_strict':(muhi+2*dhi)/(1-e)<face<chi,'source_variance':s/(1-e)<10,
 'face_gap_positive':facegap>0,'face_retention':10*x/(facegap*facegap)<=1-ret,
 'nu_positive':nu>0 and nu-1/N>0,'eta_positive':eta>0,'weighted_exp_domain':weighted>=4,'weighted_cap_payment':cap<=e,
 'radial_width':1/(1-(h/rlo)**2)<=1+e,'tau_rate_positive':ct>0 and sig-dt>0 and dt>0,'tau_prefix_domain':ct*L>=4,'tau_four_term_lower_at_4':tau_poly>0,'tau_payment':tau<=e,
 'logt_lower':S(L/2)>=x,'B_minus_gamma_positive':L/4>gamma,
 'mass_trim_band':(1-e)*masslo>mass>gamma,'retained_mass':ret*mass>gamma,'ambient_degree':0<P0<1 and M>2 and b>0,'gamma_full_deficit':(1-e)*Gamma>gamma,
 'price_log_upper':pricehi<price,'root_prefix_domain':min(sig,dt,sig-dt)*L>=2,'root_full_payment':root<slack,'scalar_slack':gamma-price-c==slack,
 'finite_history_dilation':N>=2 and x/N<1 and rho>20,'finite_history_displacement':1/N**3-F(4,1)/(8*N**4)>0,'finite_history_radius':1+rhi<6 and rho>6}
bounds={'N':N,'x':x,'rlo':rlo,'rhi':rhi,'dlo':dlo,'dhi':dhi,'mu_upper':muhi,'target_mu_plus_d_upper':muhi+dhi,'transport_mu_plus_2d_upper':muhi+2*dhi,'face':face,'gap_rate':chi-face,'mean_upper':mean,'baseline_failure':4*var/e**2,'count_failure':4/(e*e*N),'target_cap_failure':F(19,20)**399,'profile_failure':400/(dlo*dlo*N),'face_gap':facegap,'face_failure':10*x/(facegap*facegap),'nu':nu,'eta':eta,'weighted_exponent':weighted,'cap_upper':cap,'tau_rate':ct,'tau_upper':tau,'root_upper':root,'mass_log_lower':masslo,'mass':mass,'M':M,'b':b,'gamma_budget':(1-e)*Gamma,'gamma':gamma,'price_upper':pricehi,'price_budget':price,'c':c,'scalar_slack':slack,'tau_poly_at_4':tau_poly}
out={'status':'PASS' if all(gates.values()) else 'FAIL','pid':os.getpid(),'python':platform.python_version(),'scope':'fixed new-tau conditional candidate; proof and tail transfer require independent review','registration_sha256':hashlib.sha256((p/'candidate_registration.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gates':gates,'failed':[k for k,val in gates.items() if not val],'bounds_exact':{k:str(val) for k,val in bounds.items()},'display':{k:disp(val) for k,val in bounds.items()},'cpu_seconds':time.process_time()}
(p/'candidate_result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'pid':out['pid'],'gates':len(gates),'failed':out['failed'],'display':out['display'],'cpu_seconds':out['cpu_seconds']},indent=2))
