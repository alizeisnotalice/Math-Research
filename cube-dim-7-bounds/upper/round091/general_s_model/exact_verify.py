#!/usr/bin/env python3
from fractions import Fraction as F
from decimal import Decimal as D,getcontext
from math import factorial
from pathlib import Path
import json,hashlib,os,sys,time
sys.set_int_max_str_digits(0);getcontext().prec=32
p=Path(__file__).resolve().parent;r=json.loads((p/'exact_registration.json').read_text());cfg=r['parameters'];src=Path(r['proof_path']);assert hashlib.sha256(src.read_bytes()).hexdigest()==r['proof_sha256']
assert hashlib.sha256((p/'summary.json').read_bytes()).hexdigest()==r['source_summary_sha256']
L=F(sys.argv[1]);assert str(L) in r['thresholds'];v=lambda k:F(cfg[k])
def logs(y):
 assert y>1
 z=(y-1)/(y+1);power=z;total=F(0)
 for j in range(96):total+=2*power/(2*j+1);power*=z*z
 return total,total+2*power/(193*(1-z*z))
def S(y):
 assert y>=0
 term=total=F(1)
 for j in range(1,25):term*=y/j;total+=term
 return total
def disp(f):return str(D(f.numerator)/D(f.denominator))
s=v('s');chi=v('chi');sigma=v('sigma');dt=v('delta');qmax=v('qmax');gamma=v('gamma');c=v('c');trim=v('trim');e=v('e');ret=v('ret');m=v('m');face=v('m_face');price=v('price');slack=gamma-price-c
N=L**16/factorial(16);x=(1-sigma)*L;h=1/L**2;qmin=1/s;rlo=1+trim;rhi=s*qmax-trim;rho=N/(s*x)
qlo,qhi=logs(2*s/(s+1));mu_lo=s-1-(s+1)*qhi;mu_hi=s-1-(s+1)*qlo;dlo=(chi-mu_hi)/4;dhi=(chi-mu_lo)/4
masslo,masshi=logs(rhi/rlo);rawlo,rawhi=logs(s*qmax);pricelo,pricehi=logs(1/(1-sigma));log2lo,log2hi=logs(F(2))
P0=1-2*e;ulo=P0*(1-e)/(1+e);uhi=1/ulo;M=2*(1+e)*uhi;b=2*ret*(1-e)*ulo-2*e;Gamma=m*b/M
nu=1/s-(1+1/s)*e;eta=(1-qmax)/(2*qmax);wexp=(nu-1/N)*eta*x/(1+e);cap=M/S(wexp);ct=dt-chi*(1-sigma);tau=rhi*L**2/S(ct*L);root=1/N+1/S(dt*L)+1/S(sigma*L)+L/S((sigma-dt)*L)
mean=4*s**4*x*x/N+4*s*s*x/N;var=16*s**4*x*x/N;Ds=4*s*s+6*s+1;gap=(chi-face)*x-1;targetcap=1/S(F(10000)/(4*s))
gates={'registered_L0':str(L) in r['thresholds'],'fixed_domain':s>=4 and 0<sigma<1 and 1/s<qmax<1 and 0<gamma<rawlo and c<gamma-pricehi,
 'primary_gap':chi>mu_hi and mu_lo>0,'tau_domain':chi*(1-sigma)<dt<sigma,'error_domain':0<e<min(F(1,12),1/(2*(s+1))),
 'trim_band':2*h<trim and rlo<rhi,'active_R_above_1':rlo-2*h>1,'h_local':0<h<=F(1,1000),'full_R_positive':s*qmin-h>0,'g_positive':(s-1)*qmin-h>0,
 'epsilon_small':1/(x*N*N)<min(1,h),'relative_gap':qmax/(x*N*N)<=h/s,'x_ge2_beta':x>=2 and x/N<1,'kappa_gt_1':x*N*N>1,
 't_2rhi':N/x>=2*rhi,'source_rho_4_s_plus_1':rho>4*(s+1),'target_radius':1+rhi<s+1,
 'Q_bound':log2hi<1,'individual_J_local':2*s/rho<F(1,2),'individual_J_error':2*s/rho<=e,'logprod_bridge':2*Ds*x*h<=e,
 'baseline_mean':mean<=e/2,'baseline_failure':4*var/e**2<=e/8,'count_fail_minus':4/(e*e*N)<=e/8,'count_fail_plus':4/(e*e*N)<=e/8,
 'target_cap_fraction':1+x/(2*N)-1/(2*rhi)<=1-1/(4*s),'target_cap_N':N>=10001,'target_cap_probability':targetcap<=e/8,
 'd_mu_positive':dlo>0,'first_profile_pin':s*s/N<=dlo/2,'first_profile_bad':4*s**4/(dlo*dlo*N)<=e/8,
 'whole_matching_count':2*h+1/N<=3*e/4,'profile_transport':10*s*h<=dlo,
 'face_upper_strict':(mu_hi+2*dhi)/(1-e)<face<chi,'source_variance':s*s/(1-e)<2*s*s,
 'face_gap_positive':gap>0,'face_retention':2*s*s*x/(gap*gap)<=1-ret,
 'nu_lower':nu>1/(2*s) and nu-1/N>0,'eta_positive':eta>0,'weighted_exp_domain':wexp>=4,'weighted_cap_payment':cap<=e,
 'radial_width':1/(1-(h/rlo)**2)<=1+e,'tau_rate_positive':ct>0,'tau_prefix_domain':ct*L>=4,'tau_payment':tau<=e,
 'logt_lower':S(L/2)>=x,'B_minus_gamma_positive':L/4>gamma,'target_mass':(1-e)*masslo>m>gamma,'retained_mass':ret*m>gamma,
 'ambient_degree':0<P0<1 and M>2 and b>0,'gamma_full_deficit':(1-e)*Gamma>gamma,
 'price_log_upper':pricehi<price,'root_prefix_domain':min(sigma,dt,sigma-dt)*L>=2,'root_full_payment':root<slack,'scalar_slack':slack>0,
 'finite_history_dilation':N>=2 and x/N<1 and rho>4*(s+1),'finite_history_displacement':1/N**3-F(4)/(8*N**4)>0,'finite_history_radius':1+rhi<s+1 and rho>s+1}
bounds={'N':N,'x':x,'rhi':rhi,'mu_lower':mu_lo,'mu_upper':mu_hi,'dmu_lower':dlo,'dmu_upper':dhi,'mface':face,'face_gap_rate':chi-face,'source_face_failure':2*s*s*x/gap**2,'target_cap_failure':targetcap,'baseline_failure':4*var/e**2,'count_failure':4/(e*e*N),'profile_failure':4*s**4/(dlo*dlo*N),'weighted_exponent':wexp,'weighted_cap_payment':cap,'tau_rate':ct,'tau_payment':tau,'root_payment':root,'mass_log_lower':masslo,'m':m,'gamma_budget':(1-e)*Gamma,'gamma':gamma,'price_log_upper':pricehi,'price':price,'c':c,'slack':slack}
out={'status':'PASS' if all(gates.values()) else 'FAIL','pid':os.getpid(),'threshold':str(L),'scope':'conditional exact base gates, R090 whole-chain not independently certified','registration_sha256':hashlib.sha256((p/'exact_registration.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'gates':gates,'failed':[k for k,val in gates.items() if not val],'bounds_exact':{k:str(val) for k,val in bounds.items()},'display':{k:disp(val) for k,val in bounds.items()},'cpu_seconds':time.process_time()}
file=p/('exact_result_L'+str(L)+'.json');file.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'pid':out['pid'],'threshold':str(L),'gates':len(gates),'failed':out['failed'],'display':out['display'],'cpu_seconds':out['cpu_seconds']},indent=2))
