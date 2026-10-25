from pathlib import Path
from flint import arb,ctx
from fractions import Fraction as F
import hashlib,json,time
p=Path(__file__).resolve().parents[1];start=time.monotonic();ctx.prec=256
sources=['astra_upper002/report.md','audit_g011/report.md','round006/orbit_minkowski_proof.md','audit_upper019/report.md']
def ar(x):
 x=F(x);return arb(x.numerator)/x.denominator
q=arb(10)/7;c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2))
J=2*(arb(7)/4).gamma()*((arb(17)/16).gamma()/(arb(67)/64).gamma())**16
I=c**q*arb.pi()**2/60*J**q
g3=c*ar(F(648,1463))**(arb(5)/2)/(ar(F(126404864,474619145))**(arb(5)/8)*ar(F(1517312,149017))**(arb(7)/8))
U=(arb(10)/3)**q*I/g3**q
U_coarse=(arb(10)/3)**q*ar(F(157929,10**10))/ar(F(1,27000))**q
oldQ=arb(16)**(arb(5)/14)*(arb(7)/4).gamma()*(arb(4)/7).gamma()/(arb(19)/28).gamma()*((arb(61)/56).gamma()/(arb(239)/224).gamma())**16
oldI=c**q*arb.pi()**2/60*oldQ
oldU=(arb(10)/3)**q*oldI/ar(F(1,65000))**q
assert U>32 and U_coarse>32 and oldU>U_coarse
out={'status':'complete_auxiliary_same_model_bound_dominated_by_32','I_exact_formula':'c^(10/7)*pi^2/60*(2*Gamma(7/4)*(Gamma(17/16)/Gamma(67/64))^16)^(10/7)','g_exact_formula':'c*(648/1463)^(5/2)/((126404864/474619145)^(5/8)*(1517312/149017)^(7/8))','I_ball':I.str(65,more=True),'g_ball':g3.str(65,more=True),'auxiliary_cost_ball':U.str(65,more=True),'coarse_auxiliary_cost_ball':U_coarse.str(65,more=True),'preceding_auxiliary_cost_ball':oldU.str(65,more=True),'current_C5_upper':'32','new_C5_upper':False,'source_hashes':{s:hashlib.sha256((p/s).read_bytes()).hexdigest() for s in sources},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
(p/'round006/joint_cost.json').write_text(json.dumps(out,indent=2)+'\n')
print('complete model cost:',U);print('coarse rational inputs:',U_coarse);print('previous complete auxiliary:',oldU)
