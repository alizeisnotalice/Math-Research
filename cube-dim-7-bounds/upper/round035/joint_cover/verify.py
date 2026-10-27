#!/usr/bin/env python3
"""Three exact finite interface diagnostics; no parameter/dimension scan."""
import hashlib
import itertools
import json
import platform
from fractions import Fraction as F
from pathlib import Path

import scipy
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def record(x):
    return {'exact':str(x),'decimal':float(x)}

def rat(x):
    return F(float(x)).limit_denominator(10**8)

def box(center,r):
    return tuple((c-r,c+r) for c in center)

def volume(v):
    z = F(1)
    for lo,hi in v:
        z *= hi-lo
    return z

def intersects(a,b):
    return all(max(x[0],y[0])<min(x[1],y[1]) for x,y in zip(a,b))

def small_box(v,a):
    return tuple(((lo+hi)/2-(hi-lo)/(2*a),(lo+hi)/2+(hi-lo)/(2*a)) for lo,hi in v)

def union_and_increments(boxes,order):
    n = len(boxes[0])
    cuts = [sorted({edge for v in boxes for edge in v[k]}) for k in range(n)]
    increments = [F(0)]*len(boxes)
    for edges in itertools.product(*(list(zip(c[:-1],c[1:])) for c in cuts)):
        point = [(lo+hi)/2 for lo,hi in edges]
        members = [i for i in order if all(lo<x<hi for x,(lo,hi) in zip(point,boxes[i]))]
        if members:
            increments[members[0]] += volume(edges)
    return sum(increments,F(0)),increments

def cliques(vs):
    nn = len(vs)
    edges = {(i,j) for i in range(nn) for j in range(i+1,nn) if intersects(vs[i],vs[j])}
    def adjacent(i,j):
        return i==j or (min(i,j),max(i,j)) in edges
    all_cs = []
    for size in range(1,nn+1):
        for c in itertools.combinations(range(nn),size):
            if all(adjacent(i,j) for i,j in itertools.combinations(c,2)):
                all_cs.append(c)
    maximal = [c for c in all_cs if not any(all(adjacent(i,j) for i in c) for j in range(nn) if j not in c)]
    return adjacent,all_cs,maximal

def peak(q,all_cs,allowed=None):
    if allowed is None:
        allowed = set(range(len(q)))
    return max([F(0)]+[sum((q[i] for i in c),F(0)) for c in all_cs if set(c)<=allowed])

def exact_mass_lp(h,maximal):
    primal = linprog([1.0]*len(maximal),
        A_ub=[[-int(i in c) for c in maximal] for i in range(len(h))],
        b_ub=[-float(x) for x in h],bounds=(0,None),method='highs')
    assert primal.success
    m = [rat(x) for x in primal.x]
    z = [rat(-x) for x in primal.ineqlin.marginals]
    assert all(x>=0 for x in m+z)
    assert all(sum((x for x,c in zip(m,maximal) if i in c),F(0))>=h[i] for i in range(len(h)))
    assert all(sum((z[i] for i in c),F(0))<=1 for c in maximal)
    pp = sum(m,F(0))
    dd = sum((x*y for x,y in zip(h,z)),F(0))
    assert pp==dd
    return pp,{'maximal_cliques':[list(c) for c in maximal],
               'mass':record(pp),'primal_dual_equal_exact':True,
               'dual_z':[str(x) for x in z]}

def diagnose(vs,a,cc,name):
    qs = [small_box(v,a) for v in vs]
    h = [volume(v) for v in vs]
    order = sorted(range(len(vs)),key=lambda i:(-(vs[i][0][1]-vs[i][0][0]),i))
    u,inc = union_and_increments(qs,order)
    q = [x/y for x,y in zip(inc,h)]
    assert sum((x*y for x,y in zip(q,h)),F(0))==u
    adjacent,all_cs,maximal = cliques(vs)
    c = set(cc)
    assert all(adjacent(i,j) for i,j in itertools.combinations(c,2))
    anchor = max(c,key=lambda i:(h[i],-i))
    outside = set(range(len(vs)))-c
    good = {j for j in outside if all(adjacent(j,i) for i in c)}
    boundary = {j for j in outside if adjacent(j,anchor) and j not in good}
    is_module = all(all(adjacent(j,i) for i in c) or not any(adjacent(j,i) for i in c) for j in outside)
    qnew = q.copy()
    for i in c:
        qnew[i] = F(0)
    qnew[anchor] = sum((h[i]*q[i] for i in c),F(0))/h[anchor]
    assert sum((x*y for x,y in zip(qnew,h)),F(0))==u
    old,new = peak(q,all_cs),peak(qnew,all_cs)
    savings = sum((q[i] for i in c),F(0))-qnew[anchor]
    delta = peak(q,all_cs,boundary)
    rest = peak(q,all_cs,outside)
    lgood = peak(q,all_cs,good)
    upper = max(rest,old-savings+delta)
    assert new<=upper
    if is_module:
        assert old==max(rest,sum((q[i] for i in c),F(0))+lgood)
        assert new==max(rest,qnew[anchor]+lgood)
        assert new<=old
    mass,mass_receipt = exact_mass_lp(h,maximal)
    assert u/mass<=old and u/mass<=new
    unionv,_ = union_and_increments(vs,list(range(len(vs))))
    eta = F(11,10)
    actual_ratio = u/(eta*unionv)
    assert actual_ratio<=old/eta and actual_ratio<=new/eta
    return {'name':name,'n':len(vs[0]),'a':str(a),'cube_count':len(vs),
            'U_volume':record(u),'budgets':[str(x) for x in h],
            'initial_q':[str(x) for x in q],'reassigned_q':[str(x) for x in qnew],
            'compression_clique':sorted(c),'anchor':anchor,'clique_module':is_module,
            'boundary_vertices':sorted(boundary),'savings_s':record(savings),'boundary_Delta':record(delta),
            'old_true_peak':record(old),'new_true_peak':record(new),'certified_boundary_upper':record(upper),
            'payment_conservation_exact':True,'boundary_bound_exactly_attained':new==upper,
            'actual_common_f_ratio_eta_11_over_10':record(actual_ratio),
            'mass_LP_ratio':record(u/mass),'mass_LP':mass_receipt,
            'interpretation':'Finite witness-payment certificate diagnostic, not an operator lower bound or uniform upper bound.'}

a = F(3)
n = 2
nn = 6
nested = [box((2*F(1,4)**i,F(0)),a*F(1,4)**i) for i in range(nn)]
r0 = diagnose(nested,a,list(range(nn)),'nested_own_payment_obstruction')
closed_form = (1-F(1,4)**(n*nn))/(a**n*(1-F(1,4)**n))
assert F(r0['old_true_peak']['exact'])==F(nn,a**n)
assert F(r0['new_true_peak']['exact'])==closed_form
r0['all_n_formula_checked_at_n2_N6'] = True

good = [((F(-1),F(5)),),((F(-1,4),F(5,4)),),
        ((F(1),F(6)),),((F(11,2),F(8)),)]
bad = good[:-1]+[((F(9,2),F(7)),)]
r1 = diagnose(good,a,[0,1],'nontrivial_clique_module')
r2 = diagnose(bad,a,[0,1],'sharp_boundary_counterexample')
assert F(r1['old_true_peak']['exact'])==F(14,15)
assert F(r1['new_true_peak']['exact'])==F(41,60)
assert F(r2['old_true_peak']['exact'])==F(14,15)
assert F(r2['new_true_peak']['exact'])==F(61,60)
assert F(r2['savings_s']['exact'])==F(1,4)
assert F(r2['boundary_Delta']['exact'])==F(1,3)
sources = [HERE/'proof.md',HERE/'verify.py',ROOT/'upper/round032/joint_cover_theory/proof.md',ROOT/'paper/main.tex']
receipt = {'status':'passed','python':platform.python_version(),'scipy':scipy.__version__,
           'family_count':3,'dimensions':[1,2],'examples':[r0,r1,r2],
           'source_checkout':'f65a5bf92911ebed4cf803071bd7693978e5e582',
           'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
           'limitations':'No n7 specialization, no high-dimensional parameter scan, no operator lower bound, no new uniform Cn constant, no MATLAB, no Git writes.'}
(HERE/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'passed','examples':[{k:r[k] for k in ('name','old_true_peak','new_true_peak','savings_s','boundary_Delta','clique_module')} for r in (r0,r1,r2)]},indent=2))
