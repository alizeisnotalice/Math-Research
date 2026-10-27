#!/usr/bin/env python3
"""Two rational finite examples; no dimensional/parameter search."""
import hashlib
import itertools
import json
import platform
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

import scipy
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def rat(x):
    return F(float(x)).limit_denominator(10**8)

def fs(x):
    x = F(x)
    return {'exact':str(x), 'decimal':float(x)}

def cube(center, radius, dilation=F(1)):
    return tuple((x-radius*dilation, x+radius*dilation) for x in center)

def volume(box):
    out = F(1)
    for lo, hi in box:
        out *= hi-lo
    return out

def inside(point, box):
    return all(lo<x<hi for x,(lo,hi) in zip(point,box))

def cells(groups):
    """Exact coordinate partition; midpoint only labels an entire open cell."""
    all_boxes = [box for group in groups for box in group]
    n = len(all_boxes[0])
    cuts = [sorted({edge for box in all_boxes for edge in box[k]}) for k in range(n)]
    intervals = [list(zip(c[:-1],c[1:])) for c in cuts]
    for edges in itertools.product(*intervals):
        point = tuple((lo+hi)/2 for lo,hi in edges)
        vol = volume(edges)
        incidences = tuple(tuple(i for i,box in enumerate(group) if inside(point,box)) for group in groups)
        yield vol, incidences

def aggregate(group):
    result = defaultdict(F)
    for vol,(inc,) in cells([group]):
        if inc:
            result[inc] += vol
    return dict(result)

def select(data, a, b, tau, eta):
    boxes = [cube(center,r) for center,r in data]
    q_atoms = aggregate(boxes)
    remaining = set(range(len(data)))
    selected, increments, weights = [], [], []
    while remaining:
        eligible = []
        for i in remaining:
            center_i,_ = data[i]
            load = sum((p for j,p in zip(selected,weights)
                        if inside(center_i,cube(data[j][0],data[j][1],b))),F(0))
            if load <= tau:
                eligible.append(i)
        if not eligible:
            break
        i = min(eligible,key=lambda j:(-data[j][1],j))
        inc = sum((vol for pattern,vol in q_atoms.items()
                   if i in pattern and not any(j in pattern for j in selected)),F(0))
        h = eta*volume(cube(data[i][0],data[i][1],a))
        selected.append(i)
        increments.append(inc)
        weights.append(inc/h)
        remaining.remove(i)
    return selected,increments,weights

def witness_lp(v_all, atoms):
    patterns = list(atoms)
    vols = [volume(v) for v in v_all]
    primal = linprog([1.0]*len(patterns),
                     A_ub=[[-int(i in pat) for pat in patterns] for i in range(len(v_all))],
                     b_ub=[-float(v) for v in vols],bounds=(0,None),method='highs')
    assert primal.success
    m = [rat(x) for x in primal.x]
    z = [rat(-x) for x in primal.ineqlin.marginals]
    assert all(x>=0 for x in m+z)
    assert all(sum((mass for mass,pat in zip(m,patterns) if i in pat),F(0))>=vols[i] for i in range(len(v_all)))
    assert all(sum((z[i] for i in pat),F(0))<=1 for pat in patterns)
    primal_obj = sum(m,F(0))
    dual_obj = sum((zi*vi for zi,vi in zip(z,vols)),F(0))
    assert primal_obj==dual_obj
    return primal_obj, {'positive_cells':len(patterns),'primal_dual_equal_exact':True,
                        'mass_optimum':fs(primal_obj),'dual_z':[str(x) for x in z]}

def flow_lp(v, h, p, tau, cover_atoms, witness_atoms):
    patterns = list(cover_atoms)
    rows = list(witness_atoms)
    variables = [(si,j) for si,pat in enumerate(patterns) for j in pat]
    lcol = len(variables)
    objective = [0.0]*lcol+[1.0]
    aeq,beq = [],[]
    for si,pat in enumerate(patterns):
        assert sum((p[j] for j in pat),F(0))>tau
        aeq.append([float(k==si) for k,j in variables]+[0.0])
        beq.append(float(cover_atoms[pat]))
    aub,bub = [],[]
    for rho in rows:
        aub.append([float(F(1)/h[j]) if j in rho else 0.0 for si,j in variables]+[-1.0])
        bub.append(-float(sum((v[j]/h[j] for j in rho),F(0))))
    bounds = [(0,float(cover_atoms[patterns[si]]*p[j]/tau)) for si,j in variables]+[(0,None)]
    primal = linprog(objective,A_ub=aub,b_ub=bub,
                     A_eq=aeq or None,b_eq=beq or None,bounds=bounds,method='highs')
    assert primal.success
    zz = [rat(x) for x in primal.x[:-1]]
    ll = rat(primal.x[-1])
    yy = [rat(-x) for x in primal.ineqlin.marginals]
    assert all(y>=0 for y in yy) and sum(yy,F(0))==1
    c = [F(0)]*len(v)
    for value,(si,j) in zip(zz,variables):
        assert 0<=value<=cover_atoms[patterns[si]]*p[j]/tau
        c[j] += value
    for si,pat in enumerate(patterns):
        assert sum((value for value,(ki,j) in zip(zz,variables) if ki==si),F(0))==cover_atoms[pat]
    assert all(sum(((v[j]+c[j])/h[j] for j in rho),F(0))<=ll for rho in rows)
    kappa = [sum((y for y,rho in zip(yy,rows) if j in rho),F(0))/h[j] for j in range(len(v))]
    dual = sum((kappa[j]*v[j] for j in range(len(v))),F(0))
    for pat,d in cover_atoms.items():
        remaining,price = F(1),F(0)
        for j in sorted(pat,key=lambda j:(kappa[j],j)):
            fill = min(remaining,p[j]/tau)
            price += fill*kappa[j]
            remaining -= fill
        assert remaining==0
        dual += d*price
    assert dual==ll
    return ll,c,{'cover_cells':len(patterns),'witness_cells':len(rows),
                 'flow_variables':len(variables),'load_optimum':fs(ll),
                 'flow_primal_and_knapsack_dual_equal_exact':True,
                 'c':[str(x) for x in c]}

def example(n):
    a,b,tau0,eta = F(3),F(27,5),F(4,25),F(11,10)
    if n==1:
        data = [((F(0),),F(1)),((F(1),),F(9,10)),((F(3),),F(4,5)),((F(7),),F(3,5))]
    else:
        data = [((F(0),F(0)),F(1)),((F(1),F(0)),F(9,10)),
                ((F(3),F(0)),F(4,5)),((F(0),F(3)),F(7,10)),((F(7),F(0)),F(3,5))]
    q_all = [cube(center,r) for center,r in data]
    v_all = [cube(center,r,a) for center,r in data]
    sel,v,p = select(data,a,b,tau0/eta,eta)
    sel0,v0,p0 = select(data,a,b,tau0,F(1))
    assert sel==sel0 and v==v0 and all(x==y/eta for x,y in zip(p,p0))
    v_selected = [v_all[i] for i in sel]
    e_selected = [cube(data[i][0],data[i][1],b+1) for i in sel]
    cover_atoms,witness_atoms,all_atoms = defaultdict(F),defaultdict(F),defaultdict(F)
    uvol,union_vvol = F(0),F(0)
    for vol,(qi,vi,rho,ei) in cells([q_all,v_all,v_selected,e_selected]):
        if qi:
            uvol += vol
        if vi:
            union_vvol += vol
            all_atoms[vi] += vol
        if rho:
            witness_atoms[rho] += vol
        if qi and not any(i in qi for i in sel):
            assert ei
            cover_atoms[ei] += vol
    optimum,mass_receipt = witness_lp(v_all,all_atoms)
    h = [eta*volume(box) for box in v_selected]
    ll,c,flow_receipt = flow_lp(v,h,p,tau0/eta,cover_atoms,witness_atoms)
    ll0,c0,flow0 = flow_lp(v,[volume(box) for box in v_selected],p0,tau0,cover_atoms,witness_atoms)
    assert ll*eta==ll0
    assert sum(v,F(0))+sum(c,F(0))==uvol
    # Actual f=eta*1_unionV maps all masses and makes the joint identity exact.
    joint_integral = sum(((v[j]+c[j])/h[j]*eta*volume(v_selected[j]) for j in range(len(v))),F(0))
    assert joint_integral==uvol
    actual_mass = eta*union_vvol
    actual_ratio = uvol/actual_mass
    mass_ratio = uvol/optimum
    assert actual_ratio<=mass_ratio and actual_ratio<=ll
    max_p_load = max(sum((p[j] for j in rho),F(0)) for rho in witness_atoms)
    old_recovery = (1+((b+1)/a)**n/(tau0/eta))*max_p_load
    assert ll<=old_recovery
    return {'n':n,'input_cubes':len(data),'selected_indices':sel,
            'a':str(a),'b':str(b),'geometric_tau':str(tau0),'eta':str(eta),
            'U_volume':fs(uvol),'actual_f_mass':fs(actual_mass),'actual_ratio':fs(actual_ratio),
            'exact_mass_LP_ratio_certificate':fs(mass_ratio),'witness_mass_LP':mass_receipt,
            'margin_flow':flow_receipt,'geometry_flow':flow0,
            'flow_scaling_identity_exact':True,'joint_integral_identity_exact':True,
            'old_separated_recovery_same_instance':fs(old_recovery),
            'interpretation':'Finite interface verification only; no uniform-dimensional upper bound follows.'}

results = [example(1),example(2)]
sources = [HERE/'proof.md',HERE/'verify.py',ROOT/'paper/main.tex',ROOT/'revision_v3/asymptotic/proof.md']
receipt = {'status':'passed','python':platform.python_version(),'scipy':scipy.__version__,
           'batches':1,'family_count':2,'dimensions':[1,2],
           'source_checkout':'f65a5bf92911ebed4cf803071bd7693978e5e582',
           'examples':results,
           'sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
           'limitations':'No high-dimensional scan, no n=7 optimization, no operator lower bound, no new uniform constant, no MATLAB.'}
(HERE/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':'passed','examples':results},indent=2))
