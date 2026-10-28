"""One R071 tuple: exact finite gates and full Taylor prefixes, no search."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
def logs(v, terms=64):
    y=(v-1)/(v+1); power=y; total=F(0)
    for j in range(terms):
        total+=2*power/(2*j+1); power*=y*y
    return total,total+2*power/((2*terms+1)*(1-y*y))
def S64(y):
    term,total=F(1),F(1)
    for j in range(1,65):
        term*=y/j;total+=term
    return total
L=F(8192); x=F(91,100)*L; N=L**64/factorial(64)
s=F(3,2); qmin=F(2,3);qmax=F(91,100)
trim=F(1,2**19); rlo=1+trim;rhi=s*qmax-trim
h=1/(256*L); eB=F(1,2**19);eN=eB;eT=eB
eJ=F(1,2**20);eR=eJ;eS=eJ;eD=eJ
ea=F(1,65536); eret=F(14,25); retain=1-eret
sigma=F(9,100);chi=F(9,200);eta=(1-qmax)/(2*qmax)
P0=1-eB-F(6,256)
Ulo=P0*(1-eN)/(1+eJ);Uhi=(1+eN)/(P0*(1-eJ))
M=2*(1+eR)*Uhi;b=2*retain*(1-eR)*Ulo-eS-eT
mL=F(311,1000);gamma=F(13,100);Gamma=mL*b/M
masslo,_=logs(rhi/rlo);qlo,qhi=logs(F(6,5))
mu0hi=F(1,2)-F(5,2)*qlo
c2lo=3*(2-F(5,6)*(qhi*qhi+2*qhi+2))
c2hi=3*(2-F(5,6)*(qlo*qlo+2*qlo+2))
ES2lo=(1-1/N)*(c2lo-(c2hi+F(1,25))*h)
ES2hi=c2hi+1/(N*qmin)
c=F(81,20000)
prefix=S64(c*L)
Csrc=(qmax-qmin)/qmin
capratio=2*Csrc*rhi*256*L/prefix
tauratio=rhi*256*L*91/prefix # sqrt(8192)<91, only this base bound
meanrate=F(2211,50000)/(1-eJ)
gface=chi-meanrate-F(183,1000)/x
vface=F(531,100000)
ML=6*F(1,8)*L**2/(qmin*qmin*N)+F(1,2)*L/(qmin*N)
VL=12*F(1,8)*L**2/(qmin*qmin*N)
scalar=1/N+1/S64(sigma*L)+182/S64(sigma*L/2)
price_lo,price_hi=logs(1/(1-sigma))
gates={
 'L0':L==8192,
 'active_trim':2*h<trim,
 'relative_domain':1<rlo<rhi<F(3,2) and h<F(1,1000),
 'relative_epsilon':qmax/(x*N*N)<=2*h/3,
 'source_radius':qmin*N/x>4*rhi,
 'target_inside_t':N/x>rhi,
 'source_inside_t':qmax<1,
 'positive_kappa':1/(x*N*N)<F(1,2),
 'individual_J':x/(3*qmin*N)<=eJ,
 'baseline_mean':ML<=eB/2,
 'baseline_probability':4*VL/(eB*eB)<=ea/8,
 'count_probability':4/(eN*eN*N)<=ea/8,
 'count_calibration':2*h+1/N<=3*eN/4,
 'target_cap_fraction':1+x/(2*N)-1/(2*rhi)<=F(3,4),
 'target_cap_probability':F(3,4)**74<=ea/8 and N>=75,
 'profile1_mean':mu0hi<F(221,5000),
 'profile1_pin':N>=9000,
 'profile1_probability':324000000/N<=ea/8,
 'profile1_transport':h<=F(1,100000),
 'profile2_mean_lower':ES2lo>=F(5275,1000000),
 'profile2_mean_upper':ES2hi<=F(5295,1000000),
 'profile2_probability':1/(N*qmin*qmin*F(1,200000)**2)<=ea/8,
 'profile2_transport':2*h<=F(1,1000000),
 'profile2_facevariance':F(5301,1000000)/(1-eJ)<vface,
 'small_pin':qhi<F(183,1000),
 'Cantelli_gap':gface>F(1,1330),
 'Cantelli_ratio':vface*1330**2/x<=F(14,11),
 'radial_width':1/(1-(h/rlo)**2)<=1+eR,
 'product_P':0<P0<1,
 'degree':M>2 and b>0,
 'mass':(1-ea)*masslo>mL,
 'deficit':(1-eD)*Gamma>gamma,
 'cap_exponent':(eta-chi)*(1-sigma)==c,
 'tau_exponent':sigma/2-chi*(1-sigma)==c,
 'cap_finite_n':eta*x/N<=F(1,2),
 'cap_prefix':capratio<=eS,
 'tau_prefix':tauratio<=eT,
 'prefix_y':c*L>=4,
 'logt_lower':(L/2)**2/2>=L,
 'tau_positive':L/4>gamma,
 'price':price_hi<F(189,2000),
 'root_scalar':scalar<F(1,2000),
 'complete_cost':gamma-F(189,2000)-F(7,200)==F(1,2000),
}
assert all(gates.values()),{k:v for k,v in gates.items() if not v}
out={'status':'PASS','scope':'one fixed finite tuple; H2 local law and full operator transfer are analytic responsibilities',
 'gates':gates,'L0':'8192','gamma':'13/100','final_c':'7/200',
 'Gamma':str(Gamma),'Gamma_margin':str((1-eD)*Gamma-gamma),
 'cap_ratio_bound':str(capratio),'tau_ratio_bound':str(tauratio),
 'face_gap_rate':str(gface),'Cantelli_ratio_upper':str(vface*1330**2/x),
 'profile2_constant_interval':[str(c2lo),str(c2hi)],
 'full_scalar_remainder':str(scalar),'price_upper':str(price_hi),
 'no_production_search':True}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','gates':len(gates),'L0':8192,'gamma':'13/100','final_c':'7/200'}))
