"""Exact finite case check of the 24 main edge-half envelopes.

Each interval endpoint is slope*ell_i/h+constant. All independent ell_i/h>=256.
This checks envelope separation only, not greedy selection or strict witnesses.
"""
from fractions import Fraction as F
from itertools import product,combinations
import json

K=8
w=K+1
W=w+K
D=2*W
tau_min=F(256)
interval_near=((F(0),F(0)),(F(0),F(W)))
interval_far=((F(1),F(-W)),(F(1),F(0)))
interval_left=((F(0),F(D)),(F(1,2),F(-K)))
interval_right=((F(1,2),F(K)),(F(1),F(-D)))
regions=[]
for free in range(3):
    fixed=[i for i in range(3) if i!=free]
    for faces in product((0,1),repeat=2):
        for half in (0,1):
            axes=[None]*3
            axes[free]=interval_left if half==0 else interval_right
            for axis,face in zip(fixed,faces):axes[axis]=interval_near if face==0 else interval_far
            regions.append({'free':free,'faces':faces,'half':half,'axes':axes})

def positive_gap(upper,lower):
    slope=lower[0]-upper[0]
    intercept=lower[1]-upper[1]
    if slope<0:return None
    minimum=slope*tau_min+intercept
    return minimum if minimum>0 else None

minimum_gap=None
checked=0
for left,right in combinations(regions,2):
    gaps=[]
    for axis in range(3):
        a,b=left['axes'][axis],right['axes'][axis]
        for gap in (positive_gap(a[1],b[0]),positive_gap(b[1],a[0])):
            if gap is not None:gaps.append(gap)
    assert gaps
    sep=max(gaps)
    minimum_gap=sep if minimum_gap is None else min(minimum_gap,sep)
    checked+=1
assert checked==276
large_tube_fraction=F(12,K*K)
remaining_D=1-F(1,4)-large_tube_fraction
assert large_tube_fraction<=F(1,4) and remaining_D>F(1,2)
print(json.dumps({'K':K,'w_over_h':w,'W_over_h':W,'endpoint_collar_over_h':D,
                  'free_mid_source_collar_halfwidth_over_h':2*K,
                  'independent_side_lengths_lower_bound_over_h':str(tau_min),
                  'main_regions':len(regions),'universally_separated_region_pairs':checked,
                  'minimum_certified_envelope_gap_over_h':str(minimum_gap),
                  'large_P_tube_fraction_bound':str(large_tube_fraction),
                  'large_P_D_residual_fraction_a_le_1_over_4':str(remaining_D),
                  'scope':'Finite affine-envelope cases, all independent side lengths >=256h',
                  'strict_witness_or_cover_certified':False,'all_assertions_passed':True},indent=2))
