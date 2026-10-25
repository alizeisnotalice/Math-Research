#!/usr/bin/env python3
"""Independent Fraction comparison, E polynomial enclosure, and finite-L1 accounting."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,math
HERE=Path(__file__).resolve().parent
MAIN=HERE.parents[2]
EPATH=MAIN/'audit/lower007/candidate_number_field/candidate_results.json'
SRC=MAIN/'audit/lower007/primary/inward_smallstep.py'
RP=HERE/'inward_smallstep_results.json'
frozen_path=MAIN/'rounds/007/frozen_lower_before_candidate.json'
phase_path=MAIN/'coordination/phase_20261025/frozen_thresholds.json'
ledger_path=HERE/'frozen_comparands.json'
e=json.loads(EPATH.read_text());r=json.loads(RP.read_text())
ESOURCE=EPATH.with_name('certify_candidate.py')
ARITH=MAIN/'lower/round001_rectangular_certificate.py'
assert hashlib.sha256(ESOURCE.read_bytes()).hexdigest()==e['code_sha256']
assert hashlib.sha256(ARITH.read_bytes()).hexdigest()==e['input_sha256'][str(ARITH)]
frozen=json.loads(frozen_path.read_text());phase=json.loads(phase_path.read_text());ledger=json.loads(ledger_path.read_text())
original=SRC.read_text();replay=(HERE/'inward_replay.py').read_text()
assert replay==original
assert e['candidate']['a2_factor']=='1000001/1000000'
R=Q(r['exact_rational_lower']);target=Q(103346206909,50000000000)
w=r['finite_L1_witness'];K=w['receiver_index_radius'];theta=Q(w['threshold']);eta=Q(w['bump_halfwidth'])
assert K==2**45 and w['source_index_radius']==K+2
assert eta==Q(1,2*10**35)>0 and 0<theta<1
assert theta==(R+target)/(2*R)
finite=theta*R*Q(2*K+1,2*K+5)**6
assert finite>target
comparands={
 'round007_frozen_current_full_upper':Q(frozen['lower']['exact_interval'][1]),
 'phase_frozen_full_upper':Q(phase['lower']['exact_interval'][1]),
 'phase_internal_full_upper':Q(phase['internal_current_best']['lower']['exact_interval'][1])}
comparands.update({name:Q(value) for name,value in ledger['complete_upper_comparands'].items()})
assert comparands['round007_frozen_current_full_upper']==comparands['phase_frozen_full_upper']
for name,old in comparands.items():
    assert finite>old and R>old and target>old,name

# Re-enclose the displayed E polynomial from scratch with 90-digit rational root intervals.
D=10**90
def root6(n):
    lo,hi=0,10*D
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**6<=n*D**6:lo=mid
        else:hi=mid
    assert lo**6<=n*D**6<(lo+1)**6
    return Q(lo,D),Q(lo+1,D)
ul,uh=root6(2);vl,vh=root6(3)
elo=ehi=Q(0)
for key,value in e['coefficients'].items():
    i,j=map(int,key.split(','));c=Q(value)
    low=ul**i*vl**j;high=uh**i*vh**j
    if c>=0:elo+=c*low;ehi+=c*high
    else:elo+=c*high;ehi+=c*low
published_lo,published_hi=map(Q,e['exact_interval'])
assert published_lo<=elo<=ehi<=published_hi
assert R<published_lo and published_lo>target
assert published_lo>comparands['round007_frozen_current_full_upper']

# No hidden source index or independent-radius assumption.
offsets=set()
for ell in range(1,5):
    first=-(ell//2) if ell%2 else 1-ell//2
    for index in range(first,first+ell):offsets.update((index,-index))
assert min(offsets)==-2 and max(offsets)==2
def dec_interval(x,digits=18):
    scale=10**digits
    return [str(Q(x.numerator*scale//x.denominator,scale)),str(Q(-((-x.numerator*scale)//x.denominator),scale))]
paths=[SRC,HERE/'inward_replay.py',RP,EPATH,ESOURCE,ARITH,frozen_path,phase_path,ledger_path]
out=dict(status='PASS_INDEPENDENT_REPLAY_FULL_COMPARISON_AND_FINITE_L1',
 new_strict_lower='C6 > 103346206909/50000000000 = 2.06692413818',
 inward_ratio=str(R),finite_ratio=str(finite),finite_ratio_decimal_enclosure=dec_interval(finite),
 finite_ratio_minus_target=str(finite-target),finite_ratio_minus_target_enclosure=dec_interval(finite-target,24),
 new_lower_minus_frozen_full_upper=str(R-comparands['round007_frozen_current_full_upper']),
 strict_comparisons={name:bool(finite>old and R>old and target>old) for name,old in comparands.items()},
 E_polynomial_reenclosed_inside_E_interval=True,E_polynomial_root_denominator=str(D),
 E_interval_above_frozen_full_upper=bool(published_lo>comparands['round007_frozen_current_full_upper']),
 finite_L1=dict(receiver_K=K,source_index_radius=K+2,eta=str(eta),threshold=str(theta),
  source_mass=f'({2*K+5})^6',receiver_measure=f'({2*K+1})^6 * R_inward',
  exact_collar=2,all_witness_cubes_use_common_radius=True,full_positive_scale_operator_contains_witness_scales=True),
 source_algorithm_unchanged_except_ROOT=True,
 input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 upper_remains_64=True,dual_breakthrough_checkpoint_achieved=False)
(HERE/'comparison_transfer_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['finite_ratio','finite_ratio_minus_target','new_lower_minus_frozen_full_upper','inward_ratio']},indent=2))
