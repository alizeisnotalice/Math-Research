"""Independent compressed R049 base-point certificate; only exact rationals."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib, json, os, sys

OUT=Path(__file__).parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
L=Q(2)**52; d=Q(2)**-24; eps0=Q(2)**-20
h=eps0/(L*d); sigma=32*d*d; Delta=d/65536
N2=L*L/2; N8=L**8/factorial(8)
qlo=1-3*d/4-h
M=factorial(8)*(24*d**3/L**6+4*d/L**7)
V=48*factorial(8)*d**3/L**6
cap=24*Q(2)**44/(L**3*d**2)
tau=Q(2)**74*factorial(64)/Q(2)**384
scalar=Q(2)**-103+Q(factorial(16),2**144)+Q(2**27*factorial(16),2**128)
checks={
 'endpoint_parameters':(sigma,Delta)==(Q(2)**-43,Q(2)**-40),
 'bridge_domain':h<d/8 and 1+d/4-2*h>1 and qlo>Q(1,2) and 1-d/4+h<1,
 'rho_uniform':N2/L*qlo>8,
 'epsilon_over_h':2*d/(eps0*N2**2)<Q(1,2),
 'n_and_nd':N2>17 and N2*d>8 and L*d>1,
 'baseline_mean':M<Q(1,32),
 # Direct Chebyshev, rather than copying the author's separate V gate.
 'baseline_bad_probability':V/(Q(1,8)-M)**2<Q(1,8),
 'Kw_worst_z8':Q(16,15**2)<Q(1,8),
 'target_cap_coordinate':1/(1+3*d/4)-2/L>Q(1,2),
 'target_cap_face':Q(3,4)**16<Q(1,8),
 'product_bridge':1280*L*d*h<Q(1,8),
 'source_cap_stronger':cap<Q(1,64),
 'tau_stronger':tau<Q(1,128),
 'complete_scalar_stronger':scalar<Q(2)**-43,
 'endpoint_scalar_margin':2*sigma+scalar<Delta/2,
}
assert all(checks.values()),checks
source_paths={
 'threshold':ROOT/'cube-r030-parameters/cube-dim-7-bounds/upper/round049/explicit_threshold/THRESHOLD_PROOF.md',
 'family':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md',
 'local':ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round045/one_coordinate/LEMMA.md',
}
hashes={key:hashlib.sha256(path.read_bytes()).hexdigest() for key,path in source_paths.items()}
assert hashes['threshold']=='b8c385a90cefbe3b234f5269be749ad13984243e01d9109658bb38e27304e33c'
assert hashes['family']=='add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,
 'sources':{key:str(path) for key,path in source_paths.items()},
 'independent_bounds':{key:str(value) for key,value in {'M':M,'V':V,'cap':cap,'tau':tau,'scalar':scalar}.items()},
 'scope':'Base-point constants only. Infinite tails and all-q/history/operator bridges are analytic in review.md. No author verifier or production code was run.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'independent exact gates; pid',os.getpid())
