#!/usr/bin/env python3
"""One newly registered R073 tuple: exact finite gates, one process only."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import platform

HERE = Path(__file__).resolve().parent
REG = json.loads((HERE/'registration.json').read_text())
AMD = json.loads((HERE/'registration_amendment.json').read_text())
P = REG['chosen_parameters']
sources = {
    'relative_J':(Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md'),'c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1'),
    'first_profile':(Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round066/first_profile/LEMMA.md'),'67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030'),
    'second_profile':(Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile/PROOF.md'),'4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8'),
}
for path,digest in sources.values():
    assert hashlib.sha256(path.read_bytes()).hexdigest()==digest


def logs(v,terms=64):
    assert v>1
    z=(v-1)/(v+1)
    power,total=z,F(0)
    for j in range(terms):
        total+=2*power/(2*j+1)
        power*=z*z
    return total,total+2*power/((2*terms+1)*(1-z*z))


def prefix(y):
    assert y>0
    term=total=F(1)
    for j in range(1,65):
        term*=y/j
        total+=term
    return total


L=F(P['L0']); sigma=F(P['sigma']); chi=F(P['chi'])
s=F(P['s']); qmin=F(P['qmin']); qmax=F(P['qmax'])
D=F(P['hden']); trim=F(P['delta'])
eB,eN,eT,eJ,eR,eS,eD,ea=(F(P[k]) for k in ('eB','eN','eT','eJ','eRad','eS','eD','epsA'))
x=(1-sigma)*L; N=L**64/F(factorial(64)); h=1/(D*L)
rlo,rhi=1+trim,s*qmax-trim
Ret=F(AMD['Ret']); gamma=F(AMD['gamma']); cfinal=F(AMD['final_c'])
mL=F(AMD['target_mass_lower']); price=F(AMD['price_upper'])
P0=1-eB-6/D
Ulo=P0*(1-eN)/(1+eJ)
Uhi=(1+eN)/(P0*(1-eJ))
M=2*(1+eR)*Uhi
b=2*Ret*(1-eR)*Ulo-eS-eT
Gamma=mL*b/M

masslo,masshi=logs(rhi/rlo)
Qlo,Qhi=logs(F(6,5))
mu1hi=F(1,2)-F(5,2)*Qlo
mu2lo=1-5*Qhi-F(5,2)*Qhi*Qhi
mu2hi=1-5*Qlo-F(5,2)*Qlo*Qlo
mean2lo=(1-1/N)*(mu2lo-(mu2hi+F(1,25))*h)
mean2hi=mu2hi+1/(250*N*qmin)
source_log_lo,source_log_hi=logs(qmax/qmin)
price_lo,price_hi=logs(1/(1-sigma))

ML=F(3,4)*L*L/(qmin*qmin*N)+F(1,2)*L/(qmin*N)
VL=F(3,2)*L*L/(qmin*qmin*N)
meanface=F(2211,50000)/(1-eJ)
pin=F(183,1000)
vface=F(531,100000)
gap0=(chi-meanface)*x-pin
Ret0=gap0*gap0/(vface*x+gap0*gap0)
eta=(1-qmax)/(2*qmax)
cS=(eta-chi)*(1-sigma)
cT=sigma/2-chi*(1-sigma)
cn=F(1025,1024)
cap=cn*source_log_hi*rhi*D*L/prefix(cS*L)
tau=rhi*D*L*91/prefix(cT*L)
scalar=1/N+1/prefix(sigma*L)+182/prefix(sigma*L/2)

gates={
 'registered_L0':L==8192,
 'parameter_domain':s==F(3,2) and F(0)<sigma<F(1,2) and F(0)<chi,
 'bands':qmin==F(2,3) and qmin<qmax<1 and 1<rlo<rhi<F(3,2),
 'trim_full_interval':2*h<trim,
 'local_width_H2':h<=F(1,2**20),
 'relative_epsilon':qmax/(x*N*N)<=h/s,
 'epsilon_small':1/(x*N*N)<F(1,2),
 'source_radius':N*qmin/x>4*rhi,
 'target_inside_t':N/x>rhi,
 'beta_lt_two':x/N<1,
 'kappa_gt_one':x*N*N>1,
 'x_range':2<=x<=L,
 'n_for_H2':N>=2000,
 'x_over_n_for_H2':x/N<=F(1,1000),
 'individual_J':x/(3*qmin*N)<=eJ<=F(1,1024),
 'baseline_log_mean':ML<=eB/2,
 'baseline_log_probability':4*VL/(eB*eB)<=ea/8,
 'both_count_probabilities':4/(eN*eN*N)<=ea/8,
 'whole_interval_count':2*h+1/N<=3*eN/4,
 'target_cap_fraction':1+x/(2*N)-1/(2*rhi)<=F(3,4),
 'target_cap_probability':N>=75 and F(3,4)**74<=ea/8,
 'first_profile_log_direction':mu1hi<F(221,5000),
 'first_profile_pin':N>=9000,
 'first_profile_probability':F(324000000)/N<=ea/8,
 'first_profile_transport':h<=F(1,100000),
 'second_profile_constant_positive':0<mu2lo<mu2hi,
 'second_profile_mean_lower':mean2lo>F(528,100000),
 'second_profile_mean_upper':mean2hi<F(5295,1000000),
 'second_profile_probability':F(5760000)/N<=ea/8,
 'second_profile_transport':h/5<F(1,1000000),
 'true_face_variance_upper':F(5301,1000000)/(1-eJ)<vface,
 'small_pin':Qhi<pin,
 'positive_face_gap_rate':chi>meanface,
 'positive_face_gap_at_base':gap0>0,
 'fixed_Cantelli_retention':F(0)<Ret<Ret0<1,
 'radial_width':1/(1-(h/rlo)**2)<=1+eR,
 'product_P0':0<P0<1,
 'ambient_degree':M>2 and b>0,
 'new_target_mass':(1-ea)*masslo>mL,
 'new_deficit':0<gamma<(1-eD)*Gamma<Gamma<1,
 'cap_exponent':cS>0,
 'tau_exponent':cT>0,
 'cap_finite_n_correction':eta*x/N<=F(1,1025),
 'complete_cap_prefix':cap<=eS,
 'complete_tau_prefix':tau<=eT,
 'cap_prefix_tail_domain':cS*L>=4,
 'tau_prefix_tail_domain':cT*L>=4,
 'logt_lower':(L/2)**2/2>=L,
 'tau_positive':L/4>gamma,
 'exact_log_price':price_hi<price,
 'scalar_prefix_tail_domain':sigma*L/2>=4,
 'complete_root_remainder':scalar<F(AMD['complete_scalar_slack']),
 'complete_cost_slack':gamma-price-cfinal==F(AMD['complete_scalar_slack'])>0,
}
assert all(gates.values()), [k for k,v in gates.items() if not v]
key_bounds={
 'N_lower_at_base':N,'x0':x,'h0':h,'P0':P0,'Ulo':Ulo,'Uhi':Uhi,'M':M,'b':b,
 'Gamma':Gamma,'gamma_budget':(1-eD)*Gamma,'gamma_chosen':gamma,
 'Cantelli_gap0':gap0,'Cantelli_Ret0':Ret0,'Ret_chosen':Ret,
 'Cantelli_mean_rate':meanface,'Cantelli_variance_rate':vface,
 'source_log_upper':source_log_hi,'source_eta':eta,'cap_rate':cS,'tau_rate':cT,
 'source_finite_n_factor':cn,'cap_ratio_upper':cap,'tau_ratio_upper':tau,
 'target_mass_log_lower':masslo,'pin_log_upper':Qhi,'profile2_mu_lower':mu2lo,
 'profile2_mu_upper':mu2hi,'profile2_mean_lower':mean2lo,'profile2_mean_upper':mean2hi,
 'price_upper':price_hi,'price_budget':price,'full_root_remainder':scalar,
 'gamma_budget_margin':(1-eD)*Gamma-gamma,'final_c':cfinal,
 'full_cost_slack':gamma-price-cfinal,
}
out={
 'status':'PASS',
 'scope':'One new registered finite tuple and its exact base gates. The all-tail/operator claims are proved analytically in PROOF.md, pending independent review.',
 'python':platform.python_version(),'arithmetic':'fractions.Fraction, one process',
 'L0':str(L),'gamma':str(gamma),'final_c':str(cfinal),
 'log_method':'64 positive atanh terms and exact geometric tail upper bound',
 'exponential_method':'full positive S64 prefix, not a single Taylor term',
 'registration_sha256':hashlib.sha256((HERE/'registration.json').read_bytes()).hexdigest(),
 'amendment_sha256':hashlib.sha256((HERE/'registration_amendment.json').read_bytes()).hexdigest(),
 'source_sha256':{k:d for k,(p,d) in sources.items()},
 'gates':gates,'key_rational_bounds':{k:str(v) for k,v in key_bounds.items()},
 'no_production_search':True,
}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','gates':len(gates),'L0':str(L),'gamma':str(gamma),
                  'c':str(cfinal),'Gamma_budget_display':float((1-eD)*Gamma),
                  'cap_upper_display':float(cap),'tau_upper_display':float(tau),
                  'source_locks':out['source_sha256']},indent=2))
