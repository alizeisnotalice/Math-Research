"""One standalone simple tuple, independent exact budget and threshold checks."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,json,os,sys

OUT=Path(__file__).parent;ROOT=Path('/Users/zhengzhihao/Desktop/T')
a=F(9);s=F(19);d=F(1,1024);e0=F(1,2**20)
eJ=eD=eS=eT=F(1,1024);eA=F(1,64);P=F(3579,4096)
qmin=1-3*d/4;ulo=P/(2*(1+eJ));uhi=1/(P*(1-d/2)*(1-eJ))
degree=max(F(4),2*uhi/qmin);mass=(1-eA)*d/(2*(1+3*d/4))
bracket=(1-3/a)*ulo-eS-eT;Gamma=mass*bracket/degree;gamma=(1-eD)*Gamma
sigma=s*d*d;raw=gamma-sigma/(1-sigma);c=F(1,2**16)
L=F(2)**27;N2=L**2/2;h=e0/(L*d)
M=factorial(8)*(24*d**3/L**6+4*d/L**7);V=48*factorial(8)*d**3/L**6
eta=d/32-a*d*d;nu=(s/2-a)*d*d
cap=24*2**20*d*d/(eta**4*L**3)
tau=F(2)**51*factorial(64)/F(2)**384
scalar=F(2)**-53+F(factorial(16),2**176)+F(2**15*factorial(16),2**160)
# Raw R053 budget independently reconstructed, not read from a producer helper.
ar=F(73,8);sr=F(293,16);dr=F(1023,2**20)
ur=P/(2*(1+eJ));vr=1/(P*(1-dr/2)*(1-eJ))
mr=(1-eA)*dr/(2*(1+3*dr/4));br=(1-3/ar)*ur-eS-eT
gr=(1-eD)*mr*br/max(F(4),2*vr/(1-3*dr/4));sigr=sr*dr*dr
rawr=gr-sigr/(1-sigr)
Rr=F(2)**-63+F(factorial(16),2**256)+F(2**17*factorial(16),2**240)
checks={
 'domain':a>3 and s>2*a and 0<d<=F(1,16) and d<1/(32*a) and sigma<F(1,2) and bracket>0,
 'sharp_degree':degree==4,
 'mass_and_bracket':mass==F(63,131168) and bracket==F(151679,524800),
 'deficit':Gamma==F(9555777,275347865600) and gamma==F(9775559871,281956214374400),
 'net_budget':raw==F(4893063658542547,295647162275777740800) and raw>c and 0<gamma<1,
 'extended_geometry':h<d/8 and 1+d/4-2*h>1 and 1-3*d/4-h>F(1,2) and 1-d/4+h<1,
 'rho_and_epsilon':N2/L*(1-3*d/4-h)>8 and 2*d/(e0*N2**2)<F(1,2),
 'n_nd_Ld':N2>33 and N2*d>2 and L*d>1,
 'mean_and_direct_Cheb':M<F(1,32) and V/(F(1,8)-M)**2<eA/3,
 'Kw_probability':8/(9*N2*d)<eA/3,
 'target_corner':1/(1+3*d/4)-2/L>F(1,2) and F(3,4)**32<eA/3,
 'product_log_budget':F(1,8)+1280*L*d*h==F(517,4096),
 'individual_J':4*d*factorial(8)/L**7<eJ,
 'distinct_rates':eta==F(23,2**20) and nu==F(1,2**21),
 'cap_threshold':eta*L==2944 and cap==F(24*2**27,2944**4) and cap<eS,
 'tau_threshold':nu*L==64 and tau<eT,
 # The tau base prefactor 2^(101/2) < 2^51 is checked by squaring.
 'tau_prefactor':(2**20*d)**2*L**3<F(2)**102,
 'scalar_exponents':sigma*L==2432 and sigma*L>2**11 and sigma*L/2>2**10,
 'full_scalar_slack':scalar<raw-c,
 'full_raw_vs_simple_cost':Rr+4/F(2)**32<rawr-(gamma-sigma),
}
assert all(checks.values()),checks
paths={
 'threshold':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round055/simple_tuple_threshold/PROOF.md',
 'family':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md',
 'R053':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round053/explicit_candidate_threshold/THRESHOLD_PROOF.md'}
hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
assert hashes['threshold']=='5d2c61adf0e341b2a11daabec61daedd6e009b65a0139bb468711175f0288e9f'
assert hashes['family']=='aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6'
assert hashes['R053']=='26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,'sources':{k:str(p) for k,p in paths.items()},
 'budget':{k:str(v) for k,v in {'ulo':ulo,'uhi':uhi,'M':degree,'mass':mass,'bracket':bracket,'Gamma':Gamma,'gamma':gamma,'sigma':sigma,'raw':raw}.items()},
 'base_bounds':{k:str(v) for k,v in {'Mlog':M,'Vlog':V,'cap':cap,'tau':tau,'scalar':scalar,'comparison_slack':rawr-(gamma-sigma)-Rr-4/F(2)**32}.items()},
 'scope':'Standalone tuple and exact base margins. Infinite-tail and finite-history arguments are analytic in review.md. No author 35 checks, search/list traversal or production helper.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'independent exact certificates; pid',os.getpid())
