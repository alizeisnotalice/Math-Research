"""New compressed independent R056 certificate. No author helpers or programs."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,os,sys

ROOT=Path('/Users/zhengzhihao/Desktop/T');OUT=Path(__file__).parent
L=F(2)**20;sqrtL=F(2)**10;d=F(3,10);sig=F(7,100);chi=F(7,200)
e=F(1,1024);eA=F(1,64);eps=F(1,2**24);x=(1-sig)*L
N=L**8/factorial(8);h=eps/(L*d);qm=F(31,40);qM=F(37,40)
rl=F(43,40)+F(1,2**20);rh=F(49,40)-F(1,2**20)
P=1-F(1,256)-1280*L*d*h
ulo=P*(1-e)/(1+e);uhi=(1+e)/(P*(1-e));degree=F(21,8)
z=(rh-rl)/rl;mass=(1-eA)*(z-z*z/2+z**3/3-z**4/4)
bracket=(1-e)*2*ulo/qM-2*e;Gamma=mass*bracket/degree;gamma=(1-e)*Gamma
ML=(6*d**3*L**2/qm**2+d*L/qm)/N;VL=12*d**3*L**2/(qm**2*N)
bad_log=VL/(F(1,256)-ML)**2
bad_total=bad_log+8/(e*e*N)+F(3,4)**74
u=F(500,3083);Q=u-u*u/2+u**3/3
meanrate=F(43,31)/(1-e)*Q**2;varrate=F(2,3)*F(43,31)/(1-e)*Q**3
gap=chi-meanrate-Q/x;inverse_bad=varrate/(x*gap*gap)
cap8=factorial(8)*d*d/(eps*F(93,20000)**8*L**7)
tau8=factorial(8)*d/(eps*F(49,20000)**8*L**6*sqrtL)
scalar8=factorial(8)/(sig*L)**8+2*sqrtL*factorial(8)/(chi*L)**8+F(factorial(8),L**8)
checks={
 'trim_and_local_domain':h<min(F(1,2**20),d/8,F(1,10000)) and rl-2*h>F(10749,10000) and rl-d-h>=qm and rh-d+h<=qM,
 'source_and_kappa':N>75 and N/L>4 and qm*N/L>4*rh and 1/(x*N*N)<h/2,
 'target_cap_uniform_face':1/rh-F(1,4)>F(1,2) and F(3,4)**74<eA/4,
 'target_probability_direct_union':ML<F(1,512) and bad_total<eA,
 'individual_coordinate_J':2*d*L/(qm*N)<e,
 'endpoint_count_uniform_lower':2*h+1/N+e/4<e,
 'endpoint_count_uniform_upper':2*h*d+e/4<e,
 'product':P==F(65275,65536),
 'mass_log_bound_domain':0<z<1 and mass>0,
 'pin_ratio_exact':2*F(10749,10000)/(2*F(10749,10000)-d)==F(3583,3083),
 'perface_moment_rates':Q<F(1,6) and meanrate<F(4,125) and varrate<F(1,300),
 'perface_direct_Cheb':gap>0 and inverse_bad<e,
 'coarse_perface_tail':chi-F(4,125)-F(1,2000)==F(1,400) and 1600/(3*x)<e and Q/x<F(1,2000),
 'exact_degree_budget':2*uhi/qm<degree and 2<degree and bracket>0,
 'new_cap_exponent':F(74,75)*(1-qM)/(2*qM)==F(1,25) and F(1,25)-chi==F(1,200),
 'cap_Taylor8_margin':cap8<e,
 'tau_exponent_and_Taylor8':sig/2-chi*(1-sig)==F(49,20000) and tau8<e,
 'safe_budget_gamma':0<gamma<Gamma<1,
 'complete_scalar_margin':scalar8<gamma-F(7,93)-F(1,40),
 # Finite-a centre guard at the worst permitted integer n=75; it improves afterwards.
 'finite_dilation_guard':F(1,75**3)-F(4,8*75**4)>0,
}
assert all(checks.values()),checks
paths={
 'final_source':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md',
 'local_lemma':ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md',
 'parameter_proof':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round056/parameter_gates/PROOF.md',
}
hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
assert hashes['final_source']=='6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a'
assert hashes['local_lemma']=='97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,'sources':{k:str(p) for k,p in paths.items()},
 'exact_values':{k:str(v) for k,v in {'mass':mass,'ulo':ulo,'uhi':uhi,'Gamma':Gamma,'gamma':gamma,'net_coefficient':gamma-F(7,93),
 'final_slack':gamma-F(7,93)-F(1,40),'target_bad_union':bad_total,'inverse_bad_direct':inverse_bad,'cap8':cap8,'tau8':tau8,'scalar8':scalar8}.items()},
 'scope':'New independent compression: aggregate target union and direct exact-rate per-face Chebyshev, plus degree, Taylor8 and full scalar margins. No author 29/49 checks, local baseline program, production diagnostics, grid, MATLAB or numerical exponentials.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'new compressed exact certificates; pid',os.getpid())
