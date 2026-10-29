"""Reproduce exact-local-phase and independent direct-source tests (no large-n MC)."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse
import hashlib
import json
import os
import sys
import time
import numpy as np
from source_oracle import SourceOracle, NPZLabels, GFLabels, gf_pow, gf_mul, gf_trace, xor_rank, mp, mpf, LOCAL_NUMERICS

OUT = Path(__file__).resolve().parent


def serial(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, mp.mpf): return mp.nstr(x, 60)
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x).__name__)


def exact_mean_H(oracle, L):
    cuts = {F(0), F(1)}
    for j in range(oracle.K):
        a = F(j, oracle.K)+F(2*oracle.labels.label(j)+1, 2*oracle.B*oracle.K)
        for sign in [-1, 1]:
            for smooth_sign in [-1, 1]:
                cuts.add((a+sign*L/2+smooth_sign*oracle.delta/2) % 1)
    cuts = sorted(cuts)
    # H is affine between these exact knots; midpoint rule is exact.
    integral = F(0)
    for a, b in zip(cuts, cuts[1:]):
        x = (a+b)/2
        integral += (b-a)*(oracle.primitive(x+L/2)-oracle.primitive(x-L/2))
    return integral


def run_small(source):
    source = Path(source).resolve()
    labels = NPZLabels(source, 3)
    K, B = labels.K, labels.B
    oracle = SourceOracle(labels, F(1, 64*K*1024), M=32)
    a = lambda j: F(j,K)+F(2*labels.label(j)+1,2*B*K)
    # All lengths at one query use the same x; x is kept exact through endpoints.
    queries = [
        ('common_three_scales', oracle.point(112,F(27,64)), [F(1,3),F(7,6),F(29,11)]),
        ('cross_grid_boundary', oracle.point(K-1,F(1)-oracle.h/3),[F(1,K),F(7,6)]),
        ('half_peak_from_center_to_far', a(318)+F(1,2*K),[F(1,K)]),
        ('peak_edge_to_peak_edge', a(441),[oracle.delta]),
        ('peak_center_half_box', a(441)+oracle.delta/4,[oracle.delta/2]),
        ('negative_period_boundary', -oracle.delta/3,[F(2,K)]),
        ('many_complete_periods', oracle.point(500,F(3,7)),[F(10**20)+F(7,6)]),
    ]
    rows=[]
    max_diff=mp.mpf(0)
    for name,x,lengths in queries:
        values=oracle.query(x,lengths)
        for val in values:
            direct=oracle.direct_small(x,val['L'])
            discrepancy=abs(val['Y']-direct['Y'])
            assert discrepancy < mp.mpf('1e-60'), (name,discrepancy)
            max_diff=max(max_diff,discrepancy)
            # Compare peak mass independently and exactly before any Fourier arithmetic.
            recovered=val['L']+val['H']/K
            assert recovered==direct['mass']
            rows.append(dict(case=name,x=x,L=val['L'],H=val['H'],B=val['B'],
                             Y=val['Y'],direct_Y=direct['Y'],exact_peak_mass=direct['mass'],
                             absolute_difference=discrepancy))
    # Cross-grid primitive near 0 and huge period translations remain identical exactly.
    for x in [F(0), oracle.delta/3, F(1)-oracle.delta/3, a(0),a(K-1)]:
        assert oracle.primitive(x)==oracle.primitive(x+10**30)
    clip_checks=0
    for j in range(K):
        u=F(2*labels.label(j)+1,2*B)
        for phase in [F(0),oracle.h/3,u-oracle.h,u-oracle.h/3,u,
                      u+oracle.h/3,u+oracle.h,F(1)-oracle.h/3]:
            x=oracle.point(j,phase)
            assert oracle.primitive(x)==oracle.primitive_one_cell(x)
            clip_checks+=1
    Ls=[F(1,3),F(7,6),F(29,11)]
    mean_checks=[dict(L=L,exact_integral_H=exact_mean_H(oracle,L),
                      exact_integral_B=0,exact_integral_d=0) for L in Ls]
    assert all(r['exact_integral_H']==0 for r in mean_checks)
    # Retain the full H/B cross term. The identity follows by Fourier coefficient
    # integration; a second independent quadrature below validates its sign/value.
    L=F(76385,65472)
    atoms=[a(j) for j in range(K)]
    cross=mp.mpf(0)
    for m,c in enumerate(oracle.coeff,1):
        atom_cos=sum((mp.cos(2*mp.pi*mpf((m*aa)%1)) for aa in atoms),mp.mpf(0))
        sinc=mp.sin(mp.pi*mpf(m*oracle.delta))/(mp.pi*mpf(m*oracle.delta))
        cross+=2*c*mp.sin(mp.pi*mpf((m*L)%2))**2*sinc*atom_cos/(mp.pi*m)
    # Every ramp breakpoint is included; longdouble arithmetic is adequate for
    # this bounded d10 validation. No float conversion is used by the actual oracle.
    ld=np.longdouble
    av=np.array([ld(v.numerator)/ld(v.denominator) for v in atoms])
    width=ld(oracle.delta.numerator)/ld(oracle.delta.denominator)
    ll=ld(L.numerator)/ld(L.denominator)
    cuts=np.unique(np.r_[ld(0),ld(1),((av[:,None]+np.array([-1,1],dtype=ld)*width/2)[:,None,:]
                         +np.array([-1,1],dtype=ld)[None,:,None]*ll/2).ravel()%1])
    gx,gw=np.polynomial.legendre.leggauss(8)
    nodes=(cuts[:-1,None]+cuts[1:,None])/2+(cuts[1:]-cuts[:-1])[:,None]*gx/2
    weights=(cuts[1:]-cuts[:-1])[:,None]*gw/2
    def pv(x):
        x=x%1
        left,right=av-width/2,av+width/2
        beg=np.searchsorted(left,x,side='right');end=np.searchsorted(right,x,side='right')
        partial=np.where(beg>end,(x-left[np.minimum(end,K-1)])/width,0)
        return end+partial-K*x
    hh=pv(nodes+ll/2)-pv(nodes-ll/2)
    pi=np.arccos(ld(-1));bb=np.zeros_like(nodes)
    for m,c in enumerate(oracle.coeff,1):
        bb+=ld(str(c))*np.sin(pi*m*ll)*np.cos(2*pi*m*nodes)
    integrals={
        'HH':float(np.sum(weights*hh*hh)),
        'BB':float(np.sum(weights*bb*bb)),
        'cross_2HB':float(np.sum(weights*2*hh*bb)),
        'full_HplusB_squared':float(np.sum(weights*(hh+bb)**2)),
        'analytic_cross_2HB':cross,
    }
    assert abs(integrals['cross_2HB']-float(cross)) < 1e-12
    assert abs(integrals['HH']+integrals['BB']+integrals['cross_2HB']-integrals['full_HplusB_squared'])<1e-12
    return dict(K=K,B=B,delta=oracle.delta,p=oracle.p,M=oracle.M,
                source_npz_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                test_rows=rows,max_absolute_difference=max_diff,zero_mean_checks=mean_checks,
                clipped_one_cell_exact_checks=clip_checks,
                full_covariance_at_L=L,covariance_components=integrals,
                note='Historical d10,B8 is a finite implementation diagnostic, not the growing-family qualification.')


def run_gf_small():
    labels=GFLabels(7,2,0x83,2,[1,2]);K,B=labels.K,labels.B
    states=[labels.state(j) for j in range(K)]
    assert len(set(states))==K and set(states)==set(range(1,128))
    W=[0,1,2,3]
    ratios={}
    for w in W[1:]:
        for z in W[1:]:
            if w!=z:
                ratio=gf_mul(w,gf_pow(z,126,labels.modulus,7),labels.modulus,7)
                assert ratio not in ratios
                ratios[ratio]=(w,z)
    defects=[k for k in range(1,K) if labels.state(k) in ratios]
    pair_checks=[]
    for k in range(K):
        rank=labels.endpoint_rank([0,k]);hist=Counter((labels.label(j),labels.label(j+k)) for j in range(K))
        relation=None
        if k in defects:
            w,z=ratios[labels.state(k)]
            relation=(w,z) # W basis [1,2], so field bit pattern equals character index.
        for y in range(B):
            for z in range(B):
                if k==0: expected=(128//B if y==z else 0)-int(y==z==0)
                elif relation:
                    aa,cc=relation
                    sign=(-1)**(((aa&y).bit_count()+(cc&z).bit_count())&1)
                    expected=(128//(B*B))*(1+sign)-int(y==z==0)
                else: expected=128//(B*B)-int(y==z==0)
                assert hist[y,z]==expected
        pair_checks.append(dict(k=k,rank=rank,kind='zero' if k==0 else ('defect' if relation else 'ordinary')))
    # Several actual endpoint offsets at the same query phase. Rank bounds the
    # support of their joint law; this does not make Fourier phase independent.
    r=F(3,7);lengths=[F(1,3),F(7,6),F(29,11)]
    shifts=[(r+sign*K*L/2).__floor__() for L in lengths for sign in [-1,1]]
    rank=labels.endpoint_rank(shifts)
    hist=Counter(tuple(labels.label(j+s) for s in shifts) for j in range(K))
    allzero=tuple(0 for _ in shifts)
    image_support=set(hist)|{allzero}
    for output in image_support:
        assert hist[output]==(1 << (7-rank))-int(output==allzero)
    assert len(image_support)==(1 << rank)
    # A state character (trace label bit) correlated with the discrete-log circle.
    correlations=[]
    for bit in range(2):
        cs=sum(((-1)**((labels.label(j)>>bit)&1)*mp.cos(2*mp.pi*j/K) for j in range(K)),mp.mpf(0))/K
        sn=sum(((-1)**((labels.label(j)>>bit)&1)*mp.sin(2*mp.pi*j/K) for j in range(K)),mp.mpf(0))/K
        correlations.append(dict(bit=bit,mean_character=F(-1,K),cos_correlation=cs,
                                 sin_correlation=sn,magnitude=mp.sqrt(cs*cs+sn*sn)))
    assert all(abs(row['magnitude']-mp.sqrt(128)/127)<mp.mpf('1e-60') for row in correlations)
    return dict(d=7,b=2,modulus=0x83,generator=2,basis=[1,2],K=K,B=B,
                all_pair_shifts_checked=K,defect_shifts=defects,
                pair_rank_histogram=dict(Counter(row['rank'] for row in pair_checks)),
                common_query_phase=r,endpoint_lengths=lengths,endpoint_shifts=shifts,
                endpoint_rank=rank,endpoint_full_field_support=len(image_support),
                nonzero_state_joint_law_exact=True,fourier_correlations=correlations,
                conclusion='Trace labels are functions of g^j; Fourier phase uses j/K. Their product cannot be replaced by independent marginals.')


def run_manifest(path):
    start=time.time();labels=GFLabels.from_manifest(path)
    parameters=json.loads(Path(path).read_text())['parameters']
    # Actual finite width at a large period, with no table: local geometry tests.
    n=int(parameters['n']);p=F(parameters['p']);w=int(parameters['w']);M=int(parameters['M'])
    v=9*p*p*n/(29*labels.K**2);delta=1/(w*labels.K*v)
    oracle=SourceOracle(labels,delta,p=p,M=M)
    direct_trace_checks=0
    for jtest in [0,1,2,3,17,123,23456,labels.K//3,labels.K//2,labels.K-1]:
        state=labels.state(jtest)
        direct=sum(gf_trace(gf_mul(beta,state,labels.modulus,labels.d),labels.modulus,labels.d)<<i
                   for i,beta in enumerate(labels.basis))
        assert direct==labels.label(jtest)
        direct_trace_checks+=1
    j=labels.K//3;u=F(2*labels.label(j)+1,2*labels.B);a=oracle.point(j,u)
    local=[a-delta/2,a-delta/4,a,a+delta/4,a+delta/2]
    primitive=[oracle.primitive(x) for x in local]
    expected=[-u+oracle.h, F(1,4)-u+oracle.h/2,F(1,2)-u,
              F(3,4)-u-oracle.h/2,1-u-oracle.h]
    assert primitive==expected
    assert all(oracle.primitive_one_cell(x)==oracle.primitive(x) for x in local)
    binary64_local_values=[float(x) for x in local]
    assert len(set(binary64_local_values))==1
    Ls=[F(7,6),F(7,6)+F(1,labels.K*8192),F(7,6)+F(1,labels.K*2048)]
    calls_before=labels.calls
    rows=oracle.query(a,Ls)
    calls_for_query=labels.calls-calls_before
    # More than 10^8 periods: complete periods cancel from H exactly.
    for x in local:
        assert oracle.primitive(x)==oracle.primitive(x+10**30)
    next_v=9*p*p*n/(29*(2*labels.K+1)**2)
    assert mpf(v)>=16*mp.log(n) and mpf(next_v)<16*mp.log(n)
    assert labels.B>=v and labels.B/2<v and M*M>=v and (M-1)*(M-1)<v
    assert labels.K+1>3*labels.B**4
    endpoint_shifts=[(u+sign*labels.K*L/2).__floor__() for L in Ls for sign in [-1,1]]
    endpoint_rank=labels.endpoint_rank(endpoint_shifts)
    alternate_endpoint_tests=[]
    for name,alt_lengths in [('cell_separated',[F(7,6)+F(k,labels.K) for k in [0,1,2]]),
                            ('widely_separated',[F(1,3),F(7,6),F(29,11)])]:
        shifts=[(u+sign*labels.K*L/2).__floor__() for L in alt_lengths for sign in [-1,1]]
        alternate_endpoint_tests.append(dict(name=name,lengths=alt_lengths,shifts=shifts,
                                             trace_rank=labels.endpoint_rank(shifts)))
    return dict(manifest=str(Path(path).resolve()),manifest_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                d=labels.d,b=labels.b,K=labels.K,B=labels.B,delta=delta,p=oracle.p,M=M,
                v=v,n=n,w=w,dimension_and_resolution_gates_passed=True,
                basis_and_trace_masks_integers=2*labels.b,no_K_table=True,
                direct_trace_vs_mask_checks=direct_trace_checks,
                local_peak_primitive_values=primitive,local_peak_exact_checks_passed=True,
                binary64_distinct_coordinates_for_five_peak_points=len(set(binary64_local_values)),
                common_query_x=a,query=rows,label_calls_for_three_scale_query=calls_for_query,
                fixed_phase_for_endpoint_rank=u,endpoint_shifts=endpoint_shifts,
                endpoint_rank=endpoint_rank,endpoint_rank_full_field_support=2**endpoint_rank,
                alternate_endpoint_tests=alternate_endpoint_tests,
                elapsed_seconds=time.time()-start,scope='Finite number of exact-local-phase oracle queries; not coordinate Monte Carlo in n dimensions.')


def main():
    parser=argparse.ArgumentParser()
    inputs=parser.add_mutually_exclusive_group()
    inputs.add_argument('--manifest',help='GF manifest; alternatively CUBE_GF_MANIFEST')
    inputs.add_argument('--labels',help='Historical B8 source_d10.npz; alternatively CUBE_SOURCE_LABELS')
    args=parser.parse_args()
    if args.manifest is None and args.labels is None:
        args.manifest=os.environ.get('CUBE_GF_MANIFEST')
        if args.manifest is None:
            args.labels=os.environ.get('CUBE_SOURCE_LABELS')
    if args.manifest is None and args.labels is None:
        if LOCAL_NUMERICS is None:
            parser.error('Provide --labels/CUBE_SOURCE_LABELS or --manifest/CUBE_GF_MANIFEST.')
        args.labels=LOCAL_NUMERICS/'inputs/source_d10.npz'
    start=time.time()
    if args.manifest:
        result=run_manifest(args.manifest);dest=OUT/'large_manifest_results.json'
    else:
        result=dict(scope='Actual exact-local-phase source-average oracle and finite validation; no n~10^32-coordinate MC.',
                    small_npz=run_small(args.labels),small_gf=run_gf_small());dest=OUT/'validation_results.json'
    result['elapsed_seconds']=time.time()-start
    result['runtime']=dict(python=sys.version,executable=sys.executable,numpy=np.__version__,
                           mpmath=mp.__version__,mpmath_decimal_digits=mp.mp.dps,
                           mpmath_interval_arithmetic=False)
    result['script_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'source_oracle.py',Path(__file__)]}
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2,default=serial)+'\n')
    print(json.dumps(dict(output=str(dest),elapsed_seconds=result['elapsed_seconds'],all_assertions_passed=True)))

if __name__=='__main__': main()
