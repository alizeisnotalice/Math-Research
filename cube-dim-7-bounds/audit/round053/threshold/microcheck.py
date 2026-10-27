"""Independent R053 reconstruction and compressed exact base certificates."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,json,os,sys

OUT=Path(__file__).parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
a=F(73,8);s=F(293,16);d=F(1023,2**20)
eps0=F(1,2**20);eJ=eD=eS=eT=F(1,1024);eA=F(1,64)
P=F(3579,4096);qmin=1-3*d/4
ulo=P/(2*(1+eJ));uhi=1/(P*(1-d/2)*(1-eJ))
degree=max(F(4),2*uhi/qmin)
mass=(1-eA)*d/(2*(1+3*d/4))
bracket=(1-3/a)*ulo-eS-eT
Gamma=mass*bracket/degree;gamma=(1-eD)*Gamma
sigma=s*d*d;raw=gamma-sigma/(1-sigma);c=F(1,2**16)
expected={'ulo':F(3579,8200),'uhi':F(8796093022208,7674592741893),
 'mass':F(64449,134315936),'bracket':F(11148919,38310400),
 'Gamma':F(718536680631,20582788938137600),
 'gamma':F(735063024285513,21076775872652902400),
 'sigma':F(306632997,17592186044416),
 'raw':F(6468305130094850733787980147,370780099533215320770594052505600)}
actual={'ulo':ulo,'uhi':uhi,'mass':mass,'bracket':bracket,'Gamma':Gamma,'gamma':gamma,'sigma':sigma,'raw':raw}
L=F(2)**32;h=eps0/(L*d);N2=L**2/2
M=factorial(8)*(24*d**3/L**6+4*d/L**7)
V=48*factorial(8)*d**3/L**6
eta=d/32-a*d*d;nu=(s/2-a)*d*d
cap=24*2**20*d*d/(eta**4*L**3)
tau=F(2)**58*factorial(64)/F(2)**384
scalar=F(2)**-63+F(factorial(16),2**256)+F(2**17*factorial(16),2**240)
checks={
 'exact_budget_reconstruction':actual==expected,
 'full_domain':a>3 and s>2*a and 0<d<=F(1,16) and d<1/(32*a) and sigma<F(1,2) and bracket>0,
 'error_domain':all(0<e<1 for e in (eJ,eA,eD)) and eS>0 and eT>0,
 'degree_branch':degree==4,
 'strict_net_budget':raw>c and Gamma>gamma>0 and gamma<1,
 'extended_geometry':h<d/8 and 1+d/4-2*h>1 and 1-3*d/4-h>F(1,2) and 1-d/4+h<1,
 'epsilon_guard':2*d/(eps0*N2**2)<F(1,2),
 'rho_guard':N2/L*(1-3*d/4-h)>8,
 'nd_and_Ld':N2*d>2 and L*d>1 and N2>33,
 'baseline_mean':M<F(1,32),
 'baseline_direct_Cheb':V/(F(1,8)-M)**2<eA/3,
 'Kw_probability':8/(9*N2*d)<eA/3,
 'target_cap_coordinate':1/(1+3*d/4)-2/L>F(1,2),
 'target_cap_face':F(3,4)**32<eA/3,
 'product_exact_log_budget':F(1,8)+1280*L*d*h==F(517,4096) and P==1-F(517,4096),
 'individual_J':4*d*factorial(8)/L**7<eJ,
 'source_cap_error':eta>0 and cap<eS,
 'tau_base_exponent':nu*L>64 and 2**20*d*2**48<2**58,
 'tau_error':tau<eT,
 'scalar_base_exponent':sigma*L>2**16,
 'scalar_complete_slack':scalar<raw-c,
 'cost_comparison_49_42':5*F(2)**-40<c,
 'cost_comparison_older':F(2)**-38+F(14,2**20)<c,
}
assert all(checks.values()),checks
paths={
 'threshold':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round053/explicit_candidate_threshold/THRESHOLD_PROOF.md',
 'family':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md',
 'family_audit':ROOT/'cube-r030-soft/cube-dim-7-bounds/audit/round051/flexible_budget/proof_review.md',
}
hashes={key:hashlib.sha256(path.read_bytes()).hexdigest() for key,path in paths.items()}
assert hashes['threshold']=='26adb4e31f0026637ea17d8f7943bc2933dce97b8655e121f2604b6ddeea2775'
assert hashes['family']=='aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6'
assert hashes['family_audit']=='66355a1d58d3b6a832e844645b9e61ca9e3c5a2f27bc0958947054aa9144f960'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'sources':{k:str(v) for k,v in paths.items()},'hashes':hashes,
 'reconstructed_budget':{k:str(v) for k,v in actual.items()},
 'independent_bounds':{k:str(v) for k,v in {'M':M,'V':V,'cap':cap,'tau':tau,'scalar':scalar,'net_slack':raw-c}.items()},
 'scope':'One fixed tuple, exact base constants only; analytic all-tail/operator audit is review.md. No author 38 checks, helper imports, R052 enumeration or production run.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'independent exact certificates; pid',os.getpid())
