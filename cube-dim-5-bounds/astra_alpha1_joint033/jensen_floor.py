#!/usr/bin/env python3
"""A strict lower bound on a method budget, not on the actual I."""
import json,os,time
from fractions import Fraction as F
from flint import arb,ctx
from levy_dirichlet import P,S,RHO,digest,ar,dump

START=time.monotonic()

def main():
    reg=P/'jensen_registration.json'
    assert json.loads(reg.read_text())['status']=='registered_before_computation'
    data=json.loads((P/'results.json').read_text())
    assert data['parameters']=={'beta':'1/2','b':'1/10','rho':'247/256','s':'247/2560'}
    ctx.prec=256
    a=F(1,2)+S;s=ar(S);rho=ar(RHO)
    rows=[];prod=arb(1)
    for model in data['models']:
        m=model['m'];degree=m-2
        coeffsum=sum((F(t['coefficient']) for t in model['normalized_R_terms']),F(0))
        denom=F(1)
        for j in range(degree):denom*=m*a+j
        mean=coeffsum*a**degree/denom
        pref=(ar(a).gamma()/(arb(1)/2).gamma())**m*ar(F(m,2)+2*S).gamma()/ar(m*a).gamma()
        lower=pref*ar(mean)**(-s)
        prod*=lower**(ar(model['theta'])/rho)
        rows.append({'m':m,'Gram_det':model['Gram_det'],'coefficient_sum':str(coeffsum),
                     'Dirichlet_polynomial_mean_exact':str(mean),'K_Jensen_lower':str(lower)})
    moment=(1+2*s).gamma()/(1+s).gamma()
    J=4*(arb(3)/2).gamma()*moment**(3/rho)*prod
    c1=1/(32*arb.pi()**3)
    I=c1**(arb(5)/4)*arb.pi()**2/60*J**(arb(5)/4)
    check68=bool(I>ar(F(68,1000000)))
    check689=bool(I>ar(F(689,10000000)))
    assert check68
    out={'status':'certified_default_joint_Holder_budget_floor','models':rows,
         'J_budget_floor':str(J),'I_budget_floor':str(I),
         'budget_floor_above_68e_minus6':check68,'budget_floor_above_689e_minus7':check689,
         'not_a_lower_bound_on_actual_J_or_I':True,'new_C5_upper_bound':False,
         'script':digest(__file__),'registration':digest(reg),'source_exact_coefficients':digest(P/'results.json'),
         'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-START,'Arb_bits':ctx.prec}}
    dump('jensen_floor.json',out);print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
