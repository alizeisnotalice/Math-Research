"""Exact spatial version of the existing second-moment Holder certificate."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb,floor
import importlib.util,hashlib,json
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent
source=ROOT.parent/'astra_upper001'/'verify.py'
spec=importlib.util.spec_from_file_location('upper001_bridge',source)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def ar(q):
    q=F(q);return arb(q.numerator)/q.denominator


def main():
    ctx.prec=256
    sets=[{} for _ in range(5)]
    total=0
    for indices in combinations(range(16),4):
        w=old.bareiss_det([old.V[i] for i in indices]+[(1,0,0,0,0)])**2
        if not w:continue
        total+=w
        for j in range(5):
            for subset in combinations(indices,j):sets[j][subset]=sets[j].get(subset,0)+w
    assert total==16**4
    sums=[sum(w*w for w in group.values()) for group in sets]
    overlaps=[sum((-1)**(j-k)*comb(j,k)*sums[j] for j in range(k,5)) for k in range(5)]
    assert min(overlaps)>=0 and sum(overlaps)==total**2
    A=sum((w*old.laplace_moment([2]*k+[1]*(8-2*k)) for k,w in enumerate(overlaps)),F(0))
    corner=F(1517312,149017)
    B=(corner-5*A)/10
    assert A>0 and B>0
    m1,m2=F(648,1463),F(126404864,474619145)
    c=3/(2**(arb(15)/2)*arb.pi()**(arb(5)/2))
    numerator=c*ar(m1)**(arb(5)/2)/ar(m2)**(arb(5)/8)
    thresholds=[]
    for g in [F(1,10000),F(3,20000)]:
        exact=(numerator/ar(g))**(arb(8)/7)
        lower=F(floor(float(exact)*10**6)-1,10**6)
        assert ar(lower)<exact
        thresholds.append({'g_target':str(g),'Q4_threshold_Arb':str(exact),'certified_Q4_threshold':str(lower)})
    assert A-B/2<0
    balls=[]
    for radius_squared,row in zip([F(14,5),F(11,5)],thresholds):
        upper=B/2*radius_squared**2
        assert upper<F(row['certified_Q4_threshold'])
        balls.append({'norm_squared_at_most':str(radius_squared),'g_at_least':row['g_target'],'Q4_upper':str(upper)})
    out={'status':'certified_spatial_Green_lower_interface_no_full_domain_target_claim',
         'Q4':'E b(x)^2 = A sum_i x_i^4 + B sum_{i<j} x_i^2 x_j^2',
         'A':str(A),'B':str(B),'basis_subsets':1820,'overlap_weights':overlaps,
         'g3_formula':'G(x) >= numerator / Q4(x)^(7/8)',
         'numerator_Arb':str(numerator),'thresholds':thresholds,
         'certified_simple_ball_regions':balls,
         'axis_g3_Arb':str(numerator/ar(A)**(arb(7)/8)),
         'point_1_1_1_1_0_g3_Arb':str(numerator/ar(4*A+6*B)**(arb(7)/8)),
         'remaining_weak_region':'Q4(x) exceeds the relevant certified threshold; no ridge coverage has been proved here',
         'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,ROOT/'registration.md',Path(__file__)]}}
    (ROOT/'spatial_bridge.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
