"""Six fixed zero-kernel cost configurations; exact original endpoints only."""
import hashlib
import json
import signal
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

signal.alarm(170)
sys.set_int_max_str_digits(20000)
START = time.monotonic()
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
OUT = HERE/'results.json'
assert not OUT.exists(), 'No evidence overwrite'
read = lambda p: json.loads(p.read_text())


def dyadic(record):
    return [Q(int(record[s+'_mantissa']))*Q(2)**record[s+'_binary_exponent']
            for s in ('lower', 'upper')]


ROOTS = []
D = 10**100


def fifth_root(value, label):
    assert value >= 0
    if not value:
        return Q(0), Q(0)
    target = value.numerator*D**5
    den = value.denominator
    lo, hi = 0, D
    while hi**5*den <= target:
        hi *= 2
    steps = 0
    while hi-lo > 1:
        mid = (hi+lo)//2
        if mid**5*den <= target:
            lo = mid
        else:
            hi = mid
        steps += 1
        assert steps < 2048
    assert lo**5*den <= target < hi**5*den
    ROOTS.append(dict(label=label, lower=str(Q(lo,D)), upper=str(Q(hi,D)),
                      denominator=str(D), steps=steps))
    return Q(lo,D), Q(hi,D)


def main():
    paths = dict(peak=PROJECT/'upper/round036/exact_peak_budget/results.json',
                 roots=PROJECT/'audit/round058/a2_next_event_review/review_results.json',
                 point=PROJECT/'upper/round032/green_target_diagnostic/results.json',
                 point_review=PROJECT/'upper/round032/green_target_diagnostic/checked_results.json',
                 reference=PROJECT/'audit/round061/cutoff_gradient_pilot_review/manifest.json',
                 guidance=PROJECT/'coordination/guidance_20261026/next_steps_original.md')
    for name in ('low-primary','S3-primary'):
        base=PROJECT/'upper/round061/cutoff_gradient_pilot/production'/name
        paths[name+'_box']=base/'box_2.json'
        paths[name+'_jet']=base/'jet.json'
    pins={k:dict(path=str(p.relative_to(PROJECT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
          for k,p in paths.items()}
    expected=read(HERE/'input_manifest.json')['inputs']
    assert pins == expected
    roots=read(paths['roots'])
    ul,uh=map(Q,roots['root_enclosures']['u'])
    vl,vh=map(Q,roots['root_enclosures']['v'])
    assert ul**6<2<uh**6 and vl**6<3<vh**6
    # C1: analytic uniform same-time fixed-basis floor, no basis/kernel rerun.
    assert Q(32)+52+Q(677,32)==Q(3365,32)
    assert Q(16)+52+Q(677,16)==Q(1765,16)>Q(3365,32)
    gL=Q(7**6,2**14*38830**6)
    assert gL>0
    # C2: new all-time consequence of the already proved036 bound.
    peak=read(paths['peak'])['configurations'][0]
    Al,Au=dyadic(peak['A']);pl,pu=dyadic(peak['pi'])
    assert 0<Al<Au and 3<pl<pu<Q(22,7)
    BU=Au*pu/(32*ul**2*(vl**3+1))
    _,IU=fifth_root(BU**6,'complete_I_upper_from_alltime036')
    assert IU**5>=BU**6
    # C3: two disjoint REAL whole boxes, with fixed-jet gap closed by065.
    reference=read(paths['reference'])
    assert reference['full_gradient_gap_closed_for_registered_cases']
    assert not reference['universal_all_points_windows_gradient_crosscheck_claim']
    boxes=[];IL=Q(0)
    for name in ('low-primary','S3-primary'):
        box=read(paths[name+'_box']);jet=read(paths[name+'_jet'])
        assert jet['model']=='body32_unit_cauchy_v1'
        assert box['halfwidth_exact']=='1/1000'
        intervals=[list(map(Q,row))for row in box['face_box']]
        assert len(intervals)==5 and all(0<=a<b<=1 for a,b in intervals)
        volume=Q(1)
        for a,b in intervals:volume*=b-a
        lower=Q(box['certificate']['tangent_bound_lower'])
        assert lower>0
        powered,_=fifth_root(lower**6,name+'_lower_power')
        IL+=volume*powered
        boxes.append(dict(name=name,box=[[str(a),str(b)]for a,b in intervals],
                          volume=str(volume),band_lower=str(lower),
                          I_contribution_lower=str(volume*powered)))
    b0=[[Q(a),Q(b)]for a,b in boxes[0]['box']]
    b1=[[Q(a),Q(b)]for a,b in boxes[1]['box']]
    assert any(b<=c or d<=a for (a,b),(c,d) in zip(b0,b1))
    # C4: original FULL point endpoint, never rational display reconstruction.
    point=read(paths['point'])['points'][0]
    gU=Q(point['check']['G_interval'][1])
    reviewed=read(paths['point_review'])
    assert reviewed['full_old_G_chain_verified'] and reviewed['refutes_global_g_at_least_1_over_60000']
    assert 0<gL<gU<Q(1,60000)
    structural,_=fifth_root(Q(6**6),'structural_model_cost_floor')
    box_cost,_=fifth_root(Q(6**6)*IL**5/gU**6,'model_cost_lower_from_boxes')
    _,cost_upper=fifth_root(Q(6**6)*IU**5/gL**6,'model_cost_upper_coarse_pair')
    # C5: actual S interval vs unproved-target plan, and fixed-IU limitation.
    glp,_=fifth_root((gL/6)**6,'gL_div6_power_lower')
    _,gup=fifth_root((gU/6)**6,'gU_div6_power_upper')
    g0=Q(1,65000)
    g0l,g0u=fifth_root((g0/6)**6,'unproved_g0_div6_power')
    planl,planu=fifth_root(Q(6**6)*IU**5/g0**6,'conditional_plan_U')
    fixedIU_Supper=64*gup-IU
    assert fixedIU_Supper<0 and 64*g0u-IU<0
    assert max(structural,box_cost)<64<cost_upper
    assert 64*glp-IU<0<64*gup-IL
    assert len(ROOTS)<=12
    # C6: immutable endpoints and all output classifications.
    assert pins=={k:dict(path=str(p.relative_to(PROJECT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest())for k,p in paths.items()}
    result=dict(status='PASS_SIX_FIXED_CERTIFICATE_CONFIGURATIONS',model='body32_unit_cauchy_v1',q='6/5',
                g_L_global_coarse=str(gL),g_U_full_original_point=str(gU),
                B_U_complete_alltime_reenvelope=str(BU),I_U_complete=str(IU),I_L_disjoint_boxes=str(IL),
                model_U_interval=[str(max(structural,box_cost)),str(cost_upper)],
                actual_model_S_interval=[str(64*glp-IU),str(64*gup-IL)],
                conditional_g0_unproved=str(g0),conditional_plan_U_interval=[str(planl),str(planu)],
                conditional_plan_S_interval=[str(64*g0l-IU),str(64*g0u-IU)],
                fixed_IU_even_gU_plan_S_upper=str(fixedIU_Supper),
                true_model_excluded=False,fixed_coarse_IU_plan_excluded=True,
                g0_global_proved=False,new_C6_bound=False,public_internal_C6_upper='64',
                boxes=boxes,fifth_root_certificates=ROOTS,input_sha256=pins,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                proof_sha256=hashlib.sha256((HERE/'theorem.md').read_bytes()).hexdigest(),
                new_kernel_calls=0,old_producer_executions=0,old_point_recomputations=0,
                configurations=6,elapsed_seconds=time.monotonic()-START)
    with OUT.open('x')as f:f.write(json.dumps(result,indent=2)+'\n')
    print('PASS six zero-kernel cost configurations; same-model complete endpoints; no C6 improvement')


if __name__=='__main__':main()
