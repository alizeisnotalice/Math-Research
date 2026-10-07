#!/usr/bin/env python3
"""Exact guards for diagnostic-mask support and inside-weight marginal bounds.

Source boxes certify local ORIGINAL kernel/source support, not actual residual
FIRST or CP/GP qualification. Weight arrays check a proved interval lemma only.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib,datetime,time

HERE=Path(__file__).resolve().parent
OUT=HERE/'actual_face_mask_exchange_20261007_results.json'


def elementary(weights,degree):
    c=[F(1)]+[F(0)]*degree
    for w in weights:
        for j in range(degree,0,-1):c[j]+=w*c[j-1]
    return c


def marginal(weights,l,j):
    if l==0:return F(0)
    denominator=elementary(weights,l)[l]
    numerator=weights[j]*elementary(weights[:j]+weights[j+1:],l-1)[l-1]
    return numerator/denominator


def source_support(n):
    eps=F(1,100*n);L=R=F(1);sigma=F(2,n)
    cy=[F(1,2*n)]*n
    x=cy.copy();x[0]+=1
    cz=cy.copy();cz[0]+=F(3,5)
    # Three positive-volume axis boxes of halfside eps: original soft source
    # y in Cy, hard source z in Cz, and receiver x in X. A common complete
    # source is 1/2 times the uniform density on each source box; W=1.
    assert cy[0]-eps>0 and cy[0]+eps<F(1,n)
    soft_first_lower=F(1)-2*eps
    soft_other_upper=2*eps
    hard_first_lower=F(2,5)-2*eps
    hard_first_upper=F(2,5)+2*eps
    hard_other_upper=2*eps
    assert soft_first_lower>L/2 and soft_other_upper<L/2
    assert hard_first_lower>hard_other_upper and hard_first_upper<R/2
    # h_L(x-y) forces A to contain coordinate 1. Conditional K=1 only
    # A={1} can have positive original P^A; phi is positive everywhere.
    allowed_size_one_masks=[1]
    return dict(n=n,L=L,R=R,sigma=sigma,epsilon=eps,
        full_source_weights=[F(1,2),F(1,2)],full_source_mass=F(1),
        soft_source_box_center=cy,hard_source_box_center=cz,
        receiver_box_center=x,all_three_box_halfsides=eps,
        original_soft_outside_set_for_all_points=[1],
        original_hard_face_for_all_points=1,
        soft_outside_margin_lower=soft_first_lower-L/2,
        hard_cube_interior_margin_lower=R/2-hard_first_upper,
        hard_face_uniqueness_margin_lower=hard_first_lower-hard_other_upper,
        positive_size_one_mask_support=allowed_size_one_masks,
        posterior_intersection_conditional_K_one=F(1),
        uniform_prior_intersection_conditional_K_one=F(1,n),
        kernel_support_counterexample=True,
        actual_residual_counterexample=False,
        unverified_global_fields=['FIRST and fullfuture on complete source',
            'actual hard maximal winner and common threshold band',
            'strict CP/GP','original far and forest/coin gates',
            'GOOD/score/owntrace/birth and all actual histories',
            'lowS and K dual cutoff plus adaptive-core complement'])


def main():
    start=time.perf_counter();rounds=[];checks=0
    for r,n in enumerate([512,4096,32768],1):
        tests=[]
        for o in [0,1,3]:
            d=n-o
            # Values in [3/8,1] are exact interval-lemma tests, NOT claims
            # that phi at arbitrary geometric positions equals these values.
            w=[F([3,4,6,8][j%4],8) for j in range(d)]
            for l in [0,1,4]:
                k=o+l
                bound=F(8*l,3*d+5*l)
                global_bound=F(8*k,3*n+5*k)
                assert bound<=global_bound<=F(8*k,3*n)
                for j in [0,min(3,d-1),d-1]:
                    p=marginal(w,l,j)
                    assert p<=bound
                    # Also validate marginal by the exact partition identity
                    # e_l(all)=e_l(omit j)+w_j e_(l-1)(omit j).
                    allc=elementary(w,l);omit=elementary(w[:j]+w[j+1:],l)
                    assert allc[l]==omit[l]+(w[j]*omit[l-1] if l else 0)
                    tests.append(dict(o=o,d=d,k=k,l=l,j=j+1,weight=w[j],
                        posterior_inclusion=p,sharp_interval_upper=bound,
                        global_count_upper=global_bound,coarse_count_upper=F(8*k,3*n),
                        exact_verified=True));checks+=2
        source=source_support(n);checks+=5
        rounds.append(dict(round=r,n=n,source_support_guard=source,
            interval_lemma_checks=tests,weight_period=[F(3,8),F(1,2),F(3,4),F(1)],
            source_sha256=hashlib.sha256(json.dumps(source,default=str,sort_keys=True).encode()).hexdigest()))
        print('ROUND',r,'n',n,'interval_records',len(tests),'forced_hit',1,'prior',F(1,n),flush=True)
    data=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        rounds=rounds,exact_checks=checks,all_passed=True,
        arithmetic='Fraction source boxes, support margins, polynomial coefficients, elementary-symmetric marginals and bounds; no floating kernel integration, Monte Carlo, fitting or radius grid.',
        scope='Original P^A support on a common complete positive L1 two-box source is certified. It rejects unconditional K/n transfer, not a bound on fully qualified actual R_dagger. Inside coefficient arrays only test the universal interval lemma phi in [3/8,1]. No source is normalized after restriction.',
        runtime_seconds=time.perf_counter()-start,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with OUT.open('x') as f:json.dump(data,f,default=str,ensure_ascii=False,indent=2)
    print('OUTPUT',OUT,'CHECKS',checks,'RUNTIME',data['runtime_seconds'],flush=True)


if __name__=='__main__':main()
