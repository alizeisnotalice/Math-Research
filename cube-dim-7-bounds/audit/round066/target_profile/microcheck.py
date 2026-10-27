"""R066 new compressed profile/budget/price certificates; no author program."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,sys

ROOT=Path('/Users/zhengzhihao/Desktop/T');OUT=Path(__file__).parent
s=F(3,2);qmin=F(2,3);qmax=F(917,1000);sig=F(83,1000);chi=F(9,200)
delta=F(1,2**30);rl=1+delta;rh=F(2751,2000)-delta
e=F(1,2**20);eret=F(1,4096);ea=F(1,65536);mL=F(637,2000);gamma=F(317,1000)
P=1-e-F(1,2**60);ulo=P*(1-e)/(1+e);uhi=1/ulo
M=2*(1+e)*uhi;b=2*(1-eret)*(1-e)*ulo-2*e;Gamma=mL*b/M
# Six even alternating terms are an exact strict LOWER for log(6/5).
S6=sum(((-1)**(j+1)*F(1,5)**j/j for j in range(1,7)),F(0))
Q5=S6+F(1,5)**6/6
# A four-term positive atanh lower is enough for this target mass.
u=(rh-rl)/(rh+rl);mass_log_lower=2*sum((u**(2*j+1)/(2*j+1) for j in range(4)),F(0))
# Different exact-price recipe from the separate odd-seven log1p certificate.
v=sig;price_upper=sum((v**j/j for j in range(1,6)),F(0))+v**6/(6*(1-v))
checks={
 'new_relative_geometry':s*qmin==1 and s*qmax==F(2751,2000) and 1<rl<rh<F(3,2),
 'profile_mean_strict_lower_log':S6==F(2279,12500) and F(1,2)-F(5,2)*S6==F(221,5000),
 'pin_log_upper':Q5==F(34187,187500) and Q5<F(1,5),
 'mass_lower':(1-ea)*mass_log_lower>mL,
 'whole_profile_bridge':F(8)/qmin*F(2)/s==16 and F(256)/qmin*F(2)/s==512 and F(512,2**80)<F(1,2**60),
 'source_mean_interface':F(2211,50000)/(1-e)<F(223,5000),
 'source_variance_interface':F(2)*s*Q5**3/(3*(1-e))<F(1,160)<F(1,150),
 'source_tail_gap':chi-F(223,5000)-F(1,5000)==F(1,5000),
 'degree_and_deficit':M>2 and b>0 and 0<gamma<(1-e)*Gamma<Gamma<1,
 'two_strict_exponents':(1-qmax)/(2*qmax)-chi==F(47,183400) and sig/2-chi*(1-sig)==F(47,200000),
 'author_full_scalar_slack':gamma-sig/(1-sig)-F(9,40)==F(341,229250)>0,
 'exact_price_corollary':price_upper<F(4333,50000) and gamma-F(4333,50000)-F(23,100)==F(17,50000)>0,
}
assert all(checks.values()),checks
paths={
 'full_source':ROOT/'cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md',
 'first_profile':ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round066/first_profile/LEMMA.md',
 'relative_J':ROOT/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md',
 'R064':ROOT/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md',
 'R064_audit':ROOT/'cube-r030-soft/cube-dim-7-bounds/audit/round064/radial_calibration/REVIEW.md',
 'exact_scalar':ROOT/'cube-r030-parameters/cube-dim-7-bounds/audit/round066/exact_scalar/SCALAR_PROOF.md'}
hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
assert hashes['full_source']=='4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1'
assert hashes['first_profile']=='67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030'
assert hashes['relative_J']=='c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1'
receipt={'pid':os.getpid(),'python':sys.version,'checks':checks,'hashes':hashes,'sources':{k:str(p) for k,p in paths.items()},
 'exact_values':{k:str(x) for k,x in {'S6':S6,'mass_log_lower':mass_log_lower,'Gamma':Gamma,'Gamma_gap':(1-e)*Gamma-gamma,'price_upper':price_upper,'price_corollary_gap':gamma-price_upper-F(23,100)}.items()},
 'scope':'New compressed exact constants, even-six target-log direction, four-term mass and positive-five price tail. No author 16 checks, profile/local verifier, old baseline, grid, production or MATLAB.'}
(OUT/'microcheck.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS',len(checks),'new compressed certificates; pid',os.getpid())
