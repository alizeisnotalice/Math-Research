"""Exact rational piecewise-polynomial and independent cone checks, 2026-10-07.

All production outputs are confined to this directory. No old runner is executed.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parent
PAPER=ROOT.parent/'cross_shell_pair_geometry_20261007.md'
OLD_CODE=ROOT.parent/'strict_winner_shell_guard_pressure_20261007.py'
OLD_DATA=ROOT.parent/'strict_winner_shell_guard_pressure_20261007_results.json'
DIMENSIONS=(8,32,128)
MC_SAMPLES=65536


def pack(x):
    return {'exact':str(x),'float':float(x)}


def mul(p,q):
    ans=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): ans[i+j]+=a*b
    return ans


def integral(p,lo,hi):
    return sum(c*(hi**(j+1)-lo**(j+1))/F(j+1) for j,c in enumerate(p))


def p_formula(s,Q):
    a=abs(s)
    if a<=Q-1: return F(1)
    if a<=Q+1: return (Q+1-a)/2
    return F(0)


def p_intersection(s,Q):
    # Independent definition: length([-1,1] intersect [s-Q,s+Q])/2.
    return max(F(0), min(F(1),s+Q)-max(F(-1),s-Q))/2


def b_formula(s,Q):
    return (F(abs(s-1)<=Q)+F(abs(s+1)<=Q))/2


def interval_p(s,Q):
    if abs(s)<Q-1: return [F(1)]
    return [(Q+1)/2,F(-1,2) if s>0 else F(1,2)]


def one_dimensional_integrals(Q):
    assert Q>=1
    breaks=sorted(set([-Q-1,-Q+1,F(0),Q-1,Q+1]))
    sums={k:F(0) for k in ['p2','b2','bp','Ip_minus','Ip_plus','I_minus_minus','I_minus_plus','I_plus_minus','I_plus_plus']}
    intervals=[]
    for lo,hi in zip(breaks,breaks[1:]):
        mid=(lo+hi)/2
        p=interval_p(mid,Q);b=[b_formula(mid,Q)]
        ints={k:[F(abs(mid-eps)<=Q)] for k,eps in [('minus',-1),('plus',1)]}
        polynomials={'p2':mul(p,p),'b2':mul(b,b),'bp':mul(b,p)}
        for name,I in ints.items():polynomials['Ip_'+name]=mul(I,p)
        for a,I in ints.items():
            for bname,J in ints.items():polynomials['I_'+a+'_'+bname]=mul(I,J)
        vals={name:integral(poly,lo,hi) for name,poly in polynomials.items()}
        for name,val in vals.items():sums[name]+=val
        intervals.append({'lo':str(lo),'hi':str(hi),'p_coefficients':list(map(str,p)),'b_value':str(b[0]),'integrals':{name:str(val) for name,val in vals.items()}})
    A=2*Q-F(2,3);C=2*Q-1
    assert sums['p2']==A and sums['b2']==sums['bp']==C
    assert sums['Ip_minus']==sums['Ip_plus']==C
    assert sums['I_minus_minus']==sums['I_plus_plus']==2*Q
    assert sums['I_minus_plus']==sums['I_plus_minus']==2*(Q-1)
    return sums,intervals


def cone_faces(delta,Q):
    # Explicit 2n face/sign sum using interval-intersection probabilities.
    n=len(delta);p=[p_intersection(s,Q) for s in delta]
    terms=[]
    for j in range(n):
        for eps in [-1,1]:
            val=F(abs(delta[j]-eps)<=Q)
            for k in range(n):
                if k!=j:val*=p[k]
            terms.append(val)
    return terms


def point_check(delta,Q):
    n=len(delta)
    assert all(p_formula(s,Q)==p_intersection(s,Q) for s in delta)
    p=[p_formula(s,Q) for s in delta]
    L=F(0)
    for j in range(n):
        term=b_formula(delta[j],Q)
        for k in range(n):
            if k!=j:term*=p[k]
        L+=term/n
    left=cone_faces(delta,Q);right=cone_faces([-s for s in delta],Q)
    L_enum=sum(left)/F(2*n)
    # Full face-pair expansion, preserving all same/different-axis terms.
    K_enum=sum(a*b for a in left for b in right)/F((2*n)**2)
    assert L==L_enum and K_enum==L*L
    return {'L':pack(L),'K':pack(K_enum),'face_terms':2*n,'face_pair_terms':(2*n)**2}


def kernel_case(n,Q):
    sums,pieces=one_dimensional_integrals(Q)
    A=sums['p2'];C=sums['b2'];D=sums['bp']
    target=C*A**(n-1)/n+(1-F(1,n))*C*C*A**(n-2)
    # Independent integration by individual fixed faces and fixed signs.
    same=F(0);different=F(0)
    for j in range(n):
        for k in range(n):
            for a in ['minus','plus']:
                for b in ['minus','plus']:
                    if j==k:same+=sums['I_'+a+'_'+b]*A**(n-1)
                    else:different+=sums['Ip_'+a]*sums['Ip_'+b]*A**(n-2)
    direct=(same+different)/F((2*n)**2)
    assert direct==target
    if Q==1:
        simplified=F(4,3)**(n-2)*(1+F(1,3*n))
        scaled=F(9,16)*(1+F(1,3*n))*F(2,3)**n
        assert target==simplified and target/F(2)**n==scaled
    tail=Q-1+F(1,2*n)
    cases={'origin':[F(0)]*n,'near_central_tail':[tail]*n,
           'alternating_tail':[tail if j%2==0 else -tail for j in range(n)],
           'central_closed_boundary':[Q-1]*n,
           'one_outer_closed_boundary':[Q+1]+[F(0)]*(n-1),
           'one_outside':[Q+1+F(1,n)]+[F(0)]*(n-1)}
    points={label:point_check(delta,Q) for label,delta in cases.items()}
    assert points['origin']['L']['exact']=='1'
    assert points['central_closed_boundary']['L']['exact']=='1'
    assert points['one_outer_closed_boundary']['L']['exact']==str(F(1,2*n))
    assert points['one_outside']['L']['exact']=='0'
    return {'n':n,'Q':str(Q),'A':pack(A),'C':pack(C),'D':pack(D),'piecewise_integrals':pieces,
            'Lebesgue_kernel_volume':pack(direct),'scaled_displacement_volume_at_r1':pack(direct/F(2)**n),
            'same_axis_volume':pack(same/F((2*n)**2)),'different_axis_volume':pack(different/F((2*n)**2)),
            'point_checks':points,'scope':'All Q>=1; Q=1+1/n is an algebra check and lies below exp(1/n), so does not upper-bound actual unit-short-shell mutual capture.'}


def cone_mc(n,Q,index):
    # A fixed moderate-probability displacement avoids exponential rejection.
    delta=Q-1+F(1,2*n)
    exact=point_check([delta]*n,Q)
    rng_a=np.random.default_rng(709118+n*100+index*10)
    rng_b=np.random.default_rng(709119+n*100+index*10)
    ca=cb=ck=0
    def draw_capture(rng,neg=False):
        omega=rng.uniform(-1,1,size=(2048,n))
        face=rng.integers(0,n,size=2048);sign=rng.choice([-1.,1.],size=2048)
        omega[np.arange(2048),face]=sign
        return np.all(np.abs(float(delta)+(omega if neg else -omega))<=float(Q),axis=-1)
    for _ in range(MC_SAMPLES//2048):
        a=draw_capture(rng_a);b=draw_capture(rng_b,True)
        ca+=int(a.sum());cb+=int(b.sum());ck+=int((a&b).sum())
    error=math.sqrt(math.log(40)/(2*MC_SAMPLES))
    estimates={'L_left':ca/MC_SAMPLES,'L_right':cb/MC_SAMPLES,'K_joint':ck/MC_SAMPLES}
    targets={'L_left':exact['L']['float'],'L_right':exact['L']['float'],'K_joint':exact['K']['float']}
    intervals={key:[max(0.,v-error),min(1.,v+error)] for key,v in estimates.items()}
    return {'n':n,'Q':str(Q),'delta_each_coordinate':str(delta),'samples':MC_SAMPLES,
            'seeds':[709118+n*100+index*10,709119+n*100+index*10],'estimates':estimates,'targets':targets,
            'absolute_errors':{key:abs(v-targets[key]) for key,v in estimates.items()},
            'marginal95_Hoeffding_intervals':intervals,
            'targets_inside_intervals':{key:lo<=targets[key]<=hi for key,(lo,hi) in intervals.items()},
            'scope':'Independent A/B cone draws, fixed deterministic rational displacement. MC is a floating algorithm cross-check, not an exact certificate or a source-averaged energy computation.'}


def source_box_case(n):
    # Preserve the old b=2 normalization; M=n^2 b, a=1, c=1.
    a=F(1);b=F(2);M=n*n*b;ell=M-b
    def mass(h):
        I0=integral([F(1)],-h,h)
        I2=integral([F(0),F(0),F(1)],-h,h)
        return I0**n+F(1,M*M)*I2*I0**(n-1)
    W=mass(M);core=mass(ell);q=core/W
    formula=(ell/M)**n*(3+(ell/M)**2)/4
    assert q==formula and W==F(4,3)*(2*M)**n
    assert q>=F(3,4)*(1-F(1,n))
    near=F(3,2)*(b/M)**n
    lower=q*q-near
    assert lower>0
    # All x in the receiver box: this is a global range calculation.
    min_response=1+a*a/(12*M*M)
    max_response=1+(M-b/2)**2/(M*M)+b*b/(12*M*M)
    assert 1<min_response<=max_response<2
    strict_growth=(b*b-a*a)/(12*M*M)
    assert strict_growth>0 and ell+b/2==M-b/2
    # B=n log2 >= n/2 >=1: exact lower bound log2>=1/2 from integral_1^2 1/x dx.
    assert F(n,2)>=1
    return {'n':n,'a':str(a),'b':str(b),'M':str(M),'source_box_halfwidth':str(ell),'q0':pack(q),
            'q0_squared':pack(q*q),'near_pair_upper_bound':pack(near),'non_mutual_energy_lower_bound_23':pack(lower),
            'weak_q0_guard':pack(F(3,4)*(1-F(1,n))),'receiver_response_min':str(min_response),'receiver_response_max':str(max_response),
            'response_growth_from_a_to_b':str(strict_growth),'checks':['source mass from exact one-coordinate moments','q0 matches old source-box formula','hardband c<u<2c over whole receiver box','unique b winner by strictly positive quadratic R growth','inner source box outputs stay inside receiver box','B>=1','far pair exclusion and subtraction direction for (23)'],
            'scope':'Exact algebraic validation of the specified continuous L1 source guard; no old numerical runner executed and no full general energy bound claimed.'}


def old_data_formula_check():
    old=json.loads(OLD_DATA.read_text());checks=[]
    for row in old['rows']:
        n=row['n'];q=(1-F(1,n*n))**n*(3+(1-F(1,n*n))**2)/4
        assert q==F(row['inner_source_mass_fraction_exact'])
        assert q*q==F(row['Gamma_energy_lower_per_v_exact'])
        checks.append({'n':n,'equal_coordinate_x':row['equal_coordinate_x'],'saved_q0_matches':True,'saved_q0_squared_matches':True})
    return checks


def main():
    for name in ['exact_results.json','cone_mc_results.json','receipt.json']:
        assert not (ROOT/name).exists(), 'Refuse overwrite of new validation output.'
    inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PAPER,OLD_CODE,OLD_DATA]}
    kernels=[];boxes=[];mcs=[]
    for stage,n in enumerate(DIMENSIONS,1):
        qs=[F(1),1+F(1,n),1+F(2,n)]
        stage_k=[kernel_case(n,Q) for Q in qs]
        stage_mc=[cone_mc(n,Q,index) for index,Q in enumerate(qs)]
        box=source_box_case(n)
        kernels.extend(stage_k);mcs.extend(stage_mc);boxes.append(box)
        (ROOT/f'round{stage}_n{n}.json').write_text(json.dumps({'stage':stage,'kernel_checks':stage_k,'cone_mc':stage_mc,'source_box_nonmutual':box},indent=2))
        print(json.dumps({'stage':stage,'n':n,'exact_kernel_cases':3,'nonmutual_lower':float(F(box['non_mutual_energy_lower_bound_23']['exact'])),'mc_all_targets_in_marginal_intervals':all(all(row['targets_inside_intervals'].values()) for row in stage_mc)}),flush=True)
    old_checks=old_data_formula_check()
    exact={'status':'PASS_EXACT_RATIONAL_PAIR_KERNEL_AND_NONMUTUAL_GUARD','rounds':3,'kernel_cases':kernels,'source_box_checks':boxes,'old_saved_formula_checks':old_checks}
    mc={'status':'COMPLETE_FIXED_MODERATE_PROBABILITY_CONE_MC','samples_per_case':MC_SAMPLES,'cases':mcs,'confidence_scope':'Each interval is marginal 95%; not simultaneous across 27 estimates. No exponential small-probability blind sampling.'}
    (ROOT/'exact_results.json').write_text(json.dumps(exact,indent=2))
    (ROOT/'cone_mc_results.json').write_text(json.dumps(mc,indent=2))
    assert inputs=={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PAPER,OLD_CODE,OLD_DATA]}, 'Input files changed during validation.'
    receipt={'status':exact['status'],'read_scope':'cross_shell_pair_geometry_20261007.md sections 4-6; prior strictwinner source-box code/results read only','input_sha256_unchanged':inputs,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rounds':3,'dimensions':DIMENSIONS,'Q_rules':['1','1+1/n','1+2/n'],'rational_kernel_cases':9,'exact_point_cases':54,'nonmutual_source_cases':3,'old_saved_formula_cases_checked':len(old_checks),'old_data_rerun':False,'MC_cases':9,'MC_targets_all_in_marginal_intervals':all(all(row['targets_inside_intervals'].values()) for row in mcs),'limitations':['Free mutual kernel removes original hardband/winner eligibility.','Not actual FIRST/softFIRST geom certification.','(17) is Lebesgue volume, not a source-pair density bound.','(23) only certifies this nonmutual source-box guard, not super-polylog growth.','Q=1+1/n is not an upper relaxation of exp(1/n).','MC intervals are marginal and floating cross-check only.']}
    (ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
