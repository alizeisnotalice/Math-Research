"""Exact finite block-moment/proposal audit, not a source experiment."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json

Y = [
    [[F(1,2),F(1),F(3,2),F(1)], [F(3,2),F(1),F(1,2),F(1)]],
    [[F(1,4),F(5,4),F(7,4),F(3,4)], [F(3,4),F(1,4),F(5,4),F(7,4)]],
]
GS=[1,2,4]; J=2; R=2; n=4
assert all(sum(y)==4 and min(y)>0 for row in Y for y in row)
def prod(z):
    a=F(1)
    for v in z:a*=v
    return a
def blocks(g):return [list(range(k*n//g,(k+1)*n//g)) for k in range(g)]
means=[[[sum(y[2*q:2*q+2])/2 for q in range(2)] for y in row] for row in Y]
secs=[[[sum(x*x for x in y[2*q:2*q+2])/2 for q in range(2)] for y in row] for row in Y]
def value(x,j,g):return prod(sum(prod(Y[j][r][x[i]] for i in b) for r in range(R))/R for b in blocks(g))
def conditional(c,j,g):
    ms=[]; rs=[]
    for b in blocks(g):
        a=[prod(means[j][r][c[i]] for i in b) for r in range(R)]
        m=sum(a)/R
        # Triangle bound below is checked with squares, avoiding irrational roots.
        rv=[prod(secs[j][r][c[i]]/means[j][r][c[i]]**2 for i in b)-1 for r in range(R)]
        # Cauchy gives a rational, weaker upper bound on the stated Minkowski bound.
        bound=sum((a[r]/(R*m))*rv[r] for r in range(R))
        ms.append(m);rs.append(bound)
    return prod(ms),prod(1+r for r in rs)-1

rows=[]; density_mass=F(0); max_mean_gap=F(0); exact_checks=0
for c in product(range(2),repeat=n):
    hidden=[tuple(2*c[i]+u[i] for i in range(n)) for u in product(range(2),repeat=n)]
    for j in range(J):
        for g in GS:
            m,rv=conditional(c,j,g)
            z=[value(x,j,g) for x in hidden]
            true_m=sum(z)/len(z); true_var=sum((a-true_m)**2 for a in z)/len(z)
            assert m==true_m and true_var<=rv*m*m
            # Verify the tighter square-root Minkowski separately numerically.
            import math
            rbs=[]
            for b in blocks(g):
                aa=[prod(means[j][r][c[i]] for i in b) for r in range(R)]
                mb=sum(aa)/R
                rr=[prod(secs[j][r][c[i]]/means[j][r][c[i]]**2 for i in b)-1 for r in range(R)]
                rbs.append(sum(float(aa[r]/(R*mb))*math.sqrt(float(rr[r])) for r in range(R))**2)
            tight=math.prod(1+t for t in rbs)-1
            assert float(true_var/m**2)<=tight+1e-12
            exact_checks+=1
        max_mean_gap=max(max_mean_gap,abs(conditional(c,j,1)[0]-conditional(c,j,4)[0]))
    ds=[]
    for g in GS:
        for j in range(J):
            formula=conditional(c,j,g)[0]
            latent=sum(prod(prod(means[j][rr[k]][c[i]] for i in b) for k,b in enumerate(blocks(g)))
                       for rr in product(range(R),repeat=g))/R**g
            assert formula==latent
            ds.append(formula)
    d=F(1,2)+F(1,2)*sum(ds)/len(ds)
    assert d>=F(1,2)
    density_mass+=d/2**n
assert density_mass==1
tail=[]
for lam in [F(1),F(5,4),F(3,2),F(2)]:
    co={str(g):lam*sum(max(value(x,j,g) for j in range(J))>lam for x in product(range(4),repeat=n))/4**n for g in GS}
    tail.append(dict(threshold=str(lam),coefficients={g:str(x) for g,x in co.items()}))
loss=[]
for g in GS:
    ratios=[]
    for x in product(range(4),repeat=n):
        block_values=[[sum(prod(Y[j][r][x[i]] for i in b) for r in range(R))/R
                       for j in range(J)] for b in blocks(g)]
        upper=prod(max(row) for row in block_values)
        common=max(prod(row[j] for row in block_values) for j in range(J))
        ratios.append(upper/common)
    loss.append(dict(G=g,probability_no_common_scale_loss=str(F(sum(t==1 for t in ratios),len(ratios))),
                     mean_log_blockwise_over_common=sum(math.log(float(t)) for t in ratios)/len(ratios),
                     maximum_ratio=str(max(ratios))))
out=dict(scope='Four-coordinate rational toy; arbitrary positive Y, not original cube source; no MC.',
         exact_conditional_checks=exact_checks,proposal_total_mass=str(density_mass),
         max_G1_vs_coordinate_mixture_conditional_mean_gap=str(max_mean_gap),
         W_upper=2,common_query_true_tail_coefficients=tail,
         blockwise_vs_common_query_diagnostic=loss,
         variance_check='Exact rational weaker Cauchy bound and numerical tighter Minkowski bound both pass; coarse blocks independent conditional on their fine counts.')
Path(__file__).with_name('block_exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
