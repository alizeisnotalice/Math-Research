#!/usr/bin/env python3
"""Actual n=1 common-beta two-source diagnostics; finite tests only.

All geometry, W, lost masses, theta and costs are exact Fraction values.
Entropy expressions use 80-digit Decimal and are numerical diagnostics.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
from collections import defaultdict
import argparse, importlib.util, json, sys, time, hashlib
sys.dont_write_bytecode = True
LABELS = {'minimal_four_atoms','four_atom_perturbation_s471007',
          'four_atom_perturbation_s471016','geometric_contact_clouds_s491101',
          'geometric_contact_clouds_s491102','geometric_contact_clouds_s491108',
          'chain_L4_alpha1','chain_L8_alpha1','chain_L8_alpha2'}


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def decimal(q):
    return D(q.numerator)/D(q.denominator)


def entropy(q, binary=False):
    if q in (0, 1):
        return D(0)
    z = decimal(q)
    return -z*z.ln() - ((1-z)*(1-z).ln() if binary else 0)


def band(q):
    if q <= F(1, 4): return '0<q<=1/4'
    if q <= F(1, 2): return '1/4<q<=1/2'
    if q <= F(3, 4): return '1/2<q<=3/4'
    if q < 1: return '3/4<q<1'
    return 'q=1'


def evaluate(flow, record, globalmod, old, batch, label, atoms, alpha):
    started = time.monotonic()
    saved, radii, cells, _ = flow.context(old, atoms, alpha)
    lo, hi = alpha/4, alpha/2
    I = hi-lo
    X = alpha*F(saved['eligible_volume'])
    data = []
    # context splits at ALL source/radius endpoints, so even z's coarse-j
    # members (which need not be record labels at z) are fixed on a cell.
    for a, b, J, trace in cells:
        mid = (a+b)/2
        rec = record.intervals(trace, lo, hi)
        members = {k: frozenset(i for i,(y,w) in enumerate(atoms)
                                if abs(y-mid) < radii[k])
                   for k in range(1, saved['D']+1)}
        signature = (rec, members)
        if data and data[-1][1] == a and data[-1][2:] == signature:
            data[-1] = (data[-1][0], b, rec, members)
        else:
            data.append((a,b,rec,members))
        # Independent membership check at two interior points.
        for pt in ((2*a+b)/3, (a+2*b)/3):
            for k in range(1, saved['D']+1):
                assert members[k] == frozenset(i for i,(y,w) in enumerate(atoms)
                                              if abs(y-pt) < radii[k])
    O = H = square_theta = square_p = F(0)
    theta_bins, p_bins, gap_cost = defaultdict(F), defaultdict(F), defaultdict(F)
    pair_bins = defaultdict(F)
    pair_by_theta = defaultdict(lambda: defaultdict(F))
    gap_by_theta = defaultdict(lambda: defaultdict(F))
    source_depth = {rho: F(0) for rho in (F(1,2),F(3,4),F(7,8),F(15,16))}
    theta_depth = {b: {rho:F(0) for rho in source_depth} for b in
                   ('0<q<=1/4','1/4<q<=1/2','1/2<q<=3/4','3/4<q<1','q=1')}
    likelihood_theta = binary_theta = likelihood_p = binary_p = D(0)
    lost_entropy_on_O = D(0)
    small_O = large_O = large_square = F(0)
    small_entropy = D(0)
    rows = []
    source_pairs = 0
    # ln(2)=2*atanh(1/3): this partial sum is an exact strict lower bound.
    # Acceptance checks below prove a stronger finite-instance inequality.
    ln2_lower = 2*sum((F(1,3**(2*r+1)*(2*r+1)) for r in range(80)),F(0))
    acceptance_max = F(0)
    acceptance_checks = 0
    for xc, (xa,xb,ix,mx) in enumerate(data):
        for zc, (za,zb,iz,mz) in enumerate(data):
            acceptance = F(0)
            # Independent of positive geometry: include every record pair.
            for j,(ja,jb) in ix.items():
                gj = sum((atoms[i][1] for i in mx[j]),F(0))/(2*radii[j])
                assert gj > 0
                for k,(ka,kb) in iz.items():
                    if j < k:
                        width = max(F(0),min(jb,kb)-max(ja,ka))
                        acceptance += width/gj
            assert acceptance <= ln2_lower
            acceptance_max = max(acceptance_max,acceptance)
            acceptance_checks += 1
            rectangle = (xb-xa)*(zb-za)
            for j,(ja,jb) in ix.items():
                R, Vj = radii[j], 2*radii[j]
                outside = rectangle-globalmod.strip_area(xa,xb,za,zb,R)
                if outside == 0: continue
                for k,(ka,kb) in iz.items():
                    if k <= j: continue
                    W = min(jb,kb)-max(ja,ka)
                    if W <= 0: continue
                    r, Vk = radii[k], 2*radii[k]
                    common = mx[j]&mz[k]
                    if not common: continue
                    lost_set = mx[j]-mz[j]
                    assert common.isdisjoint(lost_set)
                    lost = sum((atoms[i][1] for i in lost_set),F(0))
                    coarse_mass = sum((atoms[i][1] for i in mx[j]),F(0))
                    assert lost > 0 and W*Vj <= lost
                    theta = W*Vj/lost
                    p = lost/coarse_mass
                    common_mass = sum((atoms[i][1] for i in common),F(0))
                    h = outside*common_mass*lost/(I*Vj*Vj*Vk)
                    cost = W*outside*common_mass/(I*Vj*Vk)
                    assert cost == theta*h and cost > 0
                    O += cost; H += h
                    square_theta += theta*cost; square_p += p*cost
                    likelihood_theta += entropy(theta)*decimal(h)
                    binary_theta += entropy(theta,True)*decimal(h)
                    likelihood_p += entropy(p)*decimal(h)
                    binary_p += entropy(p,True)*decimal(h)
                    lost_entropy_on_O += -decimal(p).ln()*decimal(cost)
                    if theta <= F(1,4):
                        small_O += cost
                        small_entropy += entropy(theta)*decimal(h)
                    else:
                        large_O += cost
                        large_square += theta*cost
                    tb = band(theta)
                    theta_bins[tb] += cost; p_bins[band(p)] += cost
                    gap_cost[k-j] += cost
                    gap_by_theta[tb][k-j] += cost
                    depth = {rho:F(0) for rho in source_depth}
                    pair_rows = []
                    subtotal = F(0)
                    for yi in sorted(common):
                        y, wy = atoms[yi]
                        per_y_cost = W*outside*wy/(I*Vj*Vk)
                        for rho in source_depth:
                            left, right = max(xa,y-rho*R), min(xb,y+rho*R)
                            if left < right:
                                a = (right-left)*(zb-za)-globalmod.strip_area(left,right,za,zb,R)
                                dcost = W*a*wy/(I*Vj*Vk)
                                depth[rho] += dcost
                                source_depth[rho] += dcost
                                theta_depth[tb][rho] += dcost
                        # Cellwise depth range is exact; minimum can be zero.
                        min_depth = F(0) if xa <= y <= xb else min(abs(xa-y),abs(xb-y))/R
                        max_depth = max(abs(xa-y),abs(xb-y))/R
                        assert 0 <= min_depth <= max_depth <= 1
                        y_pair_total = F(0)
                        for vi in sorted(lost_set):
                            v,wv = atoms[vi]
                            eta = abs(v-y)/R
                            assert 1-r/R < eta < 2
                            pcost = theta*wy*wv*outside/(I*Vj*Vj*Vk)
                            subtotal += pcost; y_pair_total += pcost
                            pair_rows.append(dict(y_index=yi,v_index=vi,source_separation=str(abs(v-y)),
                                                  normalized_separation=str(eta),
                                                  source_boundary_depth_range=[str(min_depth),str(max_depth)],
                                                  cost=str(pcost)))
                            eb = ('eta<=1' if eta <= 1 else '1<eta<=3/2' if eta <= F(3,2)
                                  else '3/2<eta<2')
                            pair_bins[eb] += pcost
                            pair_by_theta[tb][eb] += pcost
                            source_pairs += 1
                        assert y_pair_total == per_y_cost
                    assert subtotal == cost
                    rows.append(dict(xcell=xc,zcell=zc,x_interval=[str(xa),str(xb)],
                                     z_interval=[str(za),str(zb)],j=j,k=k,gap=k-j,
                                     common_beta_interval=[str(max(ja,ka)),str(min(jb,kb))],
                                     W=str(W),common_beta_probability=str(W/I),
                                     win_given_x_record_probability=str(W/(jb-ja)),
                                     coarse_mass=str(coarse_mass),lost=str(lost),
                                     lost_fraction=str(p),theta=str(theta),
                                     outside_area=str(outside),cost=str(cost),unrecorded_two_source_cost=str(h),
                                     source_depth_cost={str(rho):str(value) for rho,value in depth.items()},
                                     source_pairs=pair_rows))
    assert sum(theta_bins.values(),F(0)) == O
    assert sum(pair_bins.values(),F(0)) == O
    assert sum(gap_cost.values(),F(0)) == O
    assert source_depth[F(1,2)] == 0
    assert list(source_depth.values()) == sorted(source_depth.values())
    assert square_theta <= O and square_p <= O
    assert 4*large_square >= large_O
    assert small_entropy >= decimal(small_O)
    assert decimal(square_theta)+likelihood_theta >= decimal(O)
    return dict(batch=batch,label=label,alpha=str(alpha),D=saved['D'],X=str(X),O=str(O),
                O_over_X=str(O/X),merged_cell_count=len(data),positive_rectangles=len(rows),
                source_pair_count=source_pairs,theta_min=str(min(F(r['theta']) for r in rows)) if rows else None,
                acceptance_probability_checks=acceptance_checks,
                max_sum_W_over_gj=str(acceptance_max),
                ln2_rational_lower=str(ln2_lower),
                acceptance_max_le_ln2_lower_exact=True,
                theta_max=str(max(F(r['theta']) for r in rows)) if rows else None,
                lost_fraction_min=str(min(F(r['lost_fraction']) for r in rows)) if rows else None,
                lost_fraction_max=str(max(F(r['lost_fraction']) for r in rows)) if rows else None,
                theta_cost={key:str(v) for key,v in theta_bins.items()},
                lost_fraction_cost={key:str(v) for key,v in p_bins.items()},
                gap_cost={str(key):str(v) for key,v in sorted(gap_cost.items())},
                source_separation_cost={key:str(v) for key,v in pair_bins.items()},
                source_separation_by_theta={key:{eb:str(v) for eb,v in value.items()}
                                           for key,value in pair_by_theta.items()},
                gap_by_theta={key:{str(gap):str(v) for gap,v in value.items()}
                              for key,value in gap_by_theta.items()},
                source_depth_cost={str(key):str(v) for key,v in source_depth.items()},
                source_depth_by_theta={key:{str(rho):str(v) for rho,v in value.items()}
                                       for key,value in theta_depth.items()},
                candidate_diagnostics=dict(unrecorded_H=str(H),square_theta_H=str(square_theta),
                  O_over_square_theta=str(O/square_theta) if square_theta else None,
                  square_lost_fraction_on_O=str(square_p),
                  O_over_square_lost_fraction=str(O/square_p) if square_p else None,
                  likelihood_entropy_theta_H=str(likelihood_theta),
                  binary_entropy_theta_H=str(binary_theta),
                  likelihood_entropy_lost_fraction_H=str(likelihood_p),
                  binary_entropy_lost_fraction_H=str(binary_p),
                  lost_fraction_likelihood_entropy_on_O=str(lost_entropy_on_O),
                  small_theta_O=str(small_O),small_theta_entropy_H=str(small_entropy),
                  large_theta_O=str(large_O),large_theta_square_H=str(large_square),
                  small_entropy_ge_small_O_numerical=True,
                  four_large_square_ge_large_O_exact=True,
                  square_plus_likelihood_entropy_ge_O_numerical=True),
                atoms=[dict(location=str(y),mass=str(w)) for y,w in atoms],
                positive_cell_records=rows,elapsed_seconds=time.monotonic()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    flow = module(args.repo_root/'rounds/047/flow_probe.py','flow53')
    rec = module(args.repo_root/'rounds/046/threshold_probe.py','record53')
    pre = module(args.repo_root/'rounds/046/prefix_probe.py','prefix53')
    glob = module(args.repo_root/'rounds/049/global_probe.py','global53')
    old = rec.load(args.repo_root,args.output)
    refs = {c['label']:c for c in json.loads((args.repo_root/'rounds/049/verification.json').read_text())['global_records']['cases']}
    started = time.monotonic()
    cases = []
    with localcontext() as ctx:
        ctx.prec = 80
        for batch,label,atoms,alpha,params in flow.cases(old,pre):
            if label not in LABELS: continue
            c = evaluate(flow,rec,glob,old,batch,label,atoms,alpha)
            assert F(c['O_over_X']) == F(refs[label]['O_over_X'])
            c['independent_round49_O_check'] = 'passed'
            cases.append(c)
            print('batch',batch,label,'positive_cells',c['positive_rectangles'],
                  'theta',c['theta_min'],c['theta_max'],'gaps',','.join(c['gap_cost']),flush=True)
    assert len(cases) == 9
    out = dict(status='passed',dimension=1,cases=cases,
               batch_counts=[sum(c['batch']==b for c in cases) for b in range(3)],
               scope='actual full eligible E, original atomic P, one common uniform beta; no arbitrary mask substitution',
               exact_arithmetic='Fraction',entropy_precision_digits=80,
               total_positive_cells=sum(c['positive_rectangles'] for c in cases),
               total_source_pairs=sum(c['source_pair_count'] for c in cases),
               total_acceptance_probability_checks=sum(c['acceptance_probability_checks'] for c in cases),
               general_entropy_or_square_budget_proved=False,
               finite_tests_prove_uniform_bound=False,
               elapsed_seconds=time.monotonic()-started,
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',len(cases),'cases',out['total_positive_cells'],'positive cells',
          out['total_source_pairs'],'source pairs',round(out['elapsed_seconds'],2),'seconds')


if __name__ == '__main__': main()
