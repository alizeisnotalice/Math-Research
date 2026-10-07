#!/usr/bin/env python3
"""Finite exact actual n1 common-record lag relaxations. No uniform theorem.

S_support means S_fine with the common-source-exists indicator.
All costs, area integrals, cellwise essential suprema and row ratios use Fraction.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import argparse, importlib.util, json, hashlib, sys, time
sys.dont_write_bytecode = True
ORIGINAL = {'minimal_four_atoms','four_atom_perturbation_s471007',
            'four_atom_perturbation_s471016','geometric_contact_clouds_s491101',
            'geometric_contact_clouds_s491102','geometric_contact_clouds_s491108',
            'chain_L4_alpha1','chain_L8_alpha1','chain_L8_alpha2'}
EXTRA = {'four_atom_perturbation_s471008','four_atom_perturbation_s471010',
         'four_atom_perturbation_s471013','chain_L16_alpha2','chain_L32_alpha2',
         'chain_L8_alpha3/2'}
KINDS = ('coarse','fine','support')


def module(path,name):
    spec = importlib.util.spec_from_file_location(name,path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def interval_length(x,za,zb,h):
    return max(F(0),min(zb,x+h)-max(za,x-h))


def lag_length(x,za,zb,R,U):
    value = interval_length(x,za,zb,U)-interval_length(x,za,zb,R)
    assert value >= 0
    return value


def evaluate(flow,rec,glob,old,batch,label,atoms,alpha):
    started = time.monotonic()
    saved,radii,cells,_ = flow.context(old,atoms,alpha)
    lo,hi = alpha/4,alpha/2; I = hi-lo
    Evolume = F(saved['eligible_volume']); X = alpha*Evolume
    data = []
    for a,b,J,trace in cells:
        mid = (a+b)/2
        records = rec.intervals(trace,lo,hi)
        members = {k:frozenset(i for i,(y,w) in enumerate(atoms) if abs(y-mid)<radii[k])
                   for k in range(1,saved['D']+1)}
        signature = (records,members)
        if data and data[-1][1] == a and data[-1][2:] == signature:
            data[-1] = (data[-1][0],b,records,members)
        else: data.append((a,b,records,members))
    assert sum((b-a for a,b,ix,mx in data),F(0)) == Evolume
    S = {kind:F(0) for kind in KINDS}
    gaps = {kind:defaultdict(F) for kind in KINDS}
    O = O_fine = O_two_source = F(0)
    O_gap = defaultdict(F)
    point_max = {kind:F(0) for kind in KINDS}
    row_max = {kind:F(0) for kind in KINDS}
    row_witness = {kind:None for kind in KINDS}
    row_records = []
    active_pairs = positive_support_pairs = F(0)
    breakpoint_checks = 0
    for xc,(xa,xb,ix,mx) in enumerate(data):
        entries = []
        cell_rows = {j:{kind:F(0) for kind in KINDS} for j in ix}
        for zc,(za,zb,iz,mz) in enumerate(data):
            rectangle = (xb-xa)*(zb-za)
            for j,(ja,jb) in ix.items():
                R = radii[j]; Vj = 2*R; Uc = R+radii[j+1]
                out = rectangle-glob.strip_area(xa,xb,za,zb,R)
                Ac = glob.strip_area(xa,xb,za,zb,Uc)-glob.strip_area(xa,xb,za,zb,R)
                for k,(ka,kb) in iz.items():
                    if j >= k: continue
                    W = min(jb,kb)-max(ja,ka)
                    if W <= 0: continue
                    r = radii[k]; Vk = 2*r; Uf = R+r
                    Af = glob.strip_area(xa,xb,za,zb,Uf)-glob.strip_area(xa,xb,za,zb,R)
                    assert 0 <= Af <= Ac <= out
                    common = mx[j]&mz[k]
                    mass = sum((atoms[i][1] for i in common),F(0))
                    if mass:
                        # Fixed common source forces every interior point into fine lag.
                        assert Af == out and Ac == out
                        assert mass/Vk <= alpha
                        O += W*mass*out/(I*Vj*Vk)
                        O_fine += W*mass*Af/(I*Vj*Vk)
                        lost = sum((atoms[i][1] for i in mx[j]-mz[j]),F(0))
                        assert lost > 0 and W*Vj <= lost
                        theta = W*Vj/lost
                        O_two_source += theta*mass*lost*out/(I*Vj*Vj*Vk)
                        O_gap[k-j] += W*mass*out/(I*Vj*Vk)
                    if Ac == 0: continue
                    active_pairs += 1
                    positive_support_pairs += bool(common)
                    factor = W/(I*Vj)
                    values = dict(coarse=factor*Ac,fine=factor*Af,
                                  support=factor*Af if common else F(0))
                    for kind,value in values.items():
                        S[kind] += value; gaps[kind][k-j] += value
                        cell_rows[j][kind] += value
                    entries.append((j,za,zb,R,Uc,Uf,factor,bool(common)))
        # All length terms are continuous piecewise affine on this x-cell.
        # At its endpoints these are one-sided limits, hence essential suprema.
        knots = {xa,xb}
        for j,za,zb,R,Uc,Uf,factor,support in entries:
            for h in {R,Uc,Uf}:
                for z in (za,zb):
                    for sign in (-1,1):
                        p = z+sign*h
                        if xa < p < xb: knots.add(p)
        maxima = {j:{kind:F(0) for kind in KINDS} for j in ix}
        length_integral = {kind:F(0) for kind in KINDS}
        previous_x = None; previous_total = None
        for x in sorted(knots):
            total = {kind:F(0) for kind in KINDS}
            point_rows = {j:{kind:F(0) for kind in KINDS} for j in ix}
            for j,za,zb,R,Uc,Uf,factor,support in entries:
                values = dict(coarse=factor*lag_length(x,za,zb,R,Uc),
                              fine=factor*lag_length(x,za,zb,R,Uf))
                values['support'] = values['fine'] if support else F(0)
                for kind,value in values.items():
                    total[kind] += value
                    # Undo uniform beta averaging and condition on this x-record.
                    ratio = value*I/(ix[j][1]-ix[j][0])
                    point_rows[j][kind] += ratio
            for kind in KINDS:
                point_max[kind] = max(point_max[kind],total[kind])
                if previous_x is not None:
                    length_integral[kind] += (x-previous_x)*(previous_total[kind]+total[kind])/2
            for j in ix:
                for kind in KINDS:
                    value = point_rows[j][kind]
                    maxima[j][kind] = max(maxima[j][kind],value)
                    if value > row_max[kind]:
                        row_max[kind] = value
                        row_witness[kind] = dict(xcell=xc,j=j,x_one_sided_limit=str(x),
                                                 record_interval=[str(v) for v in ix[j]],
                                                 value=str(value))
            breakpoint_checks += 1
            previous_x = x; previous_total = total
        for kind in KINDS:
            assert length_integral[kind] == sum((cell_rows[j][kind] for j in ix),F(0))
        for j in ix:
            width = ix[j][1]-ix[j][0]
            average = {kind:cell_rows[j][kind]*I/(width*(xb-xa)) for kind in KINDS}
            assert all(average[kind] <= maxima[j][kind] for kind in KINDS)
            row_records.append(dict(xcell=xc,x_interval=[str(xa),str(xb)],j=j,
                                    record_interval=[str(v) for v in ix[j]],
                                    conditional_row_average={kind:str(v) for kind,v in average.items()},
                                    conditional_row_max={kind:str(v) for kind,v in maxima[j].items()}))
    assert O == O_fine == O_two_source
    assert O <= alpha*S['support'] <= alpha*S['fine'] <= alpha*S['coarse']
    assert all(sum(gaps[kind].values(),F(0)) == S[kind] for kind in KINDS)
    assert sum(O_gap.values(),F(0)) == O
    assert all(S[kind]/Evolume <= point_max[kind] for kind in KINDS)
    return dict(batch=batch,label=label,params=dict(N=len(atoms),alpha=str(alpha),D=saved['D']),
                Evolume=str(Evolume),X=str(X),O=str(O),O_over_X=str(O/X),
                S={kind:str(v) for kind,v in S.items()},
                S_over_Evolume={kind:str(v/Evolume) for kind,v in S.items()},
                pointwise_max_Sx={kind:str(v) for kind,v in point_max.items()},
                record_row_max={kind:str(v) for kind,v in row_max.items()},
                record_row_witness=row_witness,record_rows=row_records,
                false_coarse_to_fine_fraction=str((S['coarse']-S['fine'])/S['coarse']) if S['coarse'] else None,
                false_fine_without_source_fraction=str((S['fine']-S['support'])/S['fine']) if S['fine'] else None,
                false_coarse_without_source_fraction=str((S['coarse']-S['support'])/S['coarse']) if S['coarse'] else None,
                S_gap={kind:{str(k):str(v) for k,v in sorted(gap.items())} for kind,gap in gaps.items()},
                O_gap={str(k):str(v) for k,v in sorted(O_gap.items())},
                merged_cells=len(data),geometrically_active_record_pairs=int(active_pairs),
                common_source_active_record_pairs=int(positive_support_pairs),
                breakpoint_checks=breakpoint_checks,
                actual_O_independent_fine_lag_and_two_source_checks='passed',
                strip_area_equals_length_trapezoid_exact=True,
                O_le_alpha_Ssupport_le_alpha_Sfine_le_alpha_Scoarse=True,
                elapsed_seconds=time.monotonic()-started)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[2])
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    flow = module(args.repo_root/'rounds/047/flow_probe.py','flow54')
    rec = module(args.repo_root/'rounds/046/threshold_probe.py','record54')
    prefix = module(args.repo_root/'rounds/046/prefix_probe.py','prefix54')
    glob = module(args.repo_root/'rounds/049/global_probe.py','global54')
    old = rec.load(args.repo_root,args.output)
    refs = {c['label']:c for c in json.loads((args.repo_root/'rounds/053/verification.json').read_text())['record_tests']['cases']}
    started = time.monotonic(); cases = []
    for batch,label,atoms,alpha,params in flow.cases(old,prefix):
        if label not in ORIGINAL|EXTRA: continue
        c = evaluate(flow,rec,glob,old,batch,label,atoms,alpha)
        if label in ORIGINAL:
            assert F(c['O']) == F(refs[label]['O'])
            c['round53_actual_O_reference_check'] = 'passed'
        cases.append(c)
        print(batch,label,'cells',c['merged_cells'],'O/X',float(F(c['O_over_X'])),
              'Scoarse/E',float(F(c['S_over_Evolume']['coarse'])),
              'max row',float(F(c['record_row_max']['coarse'])),flush=True)
    assert len(cases) == len(ORIGINAL|EXTRA)
    out = dict(status='passed',dimension=1,arithmetic='Fraction exact',cases=cases,
               batch_counts=[sum(c['batch']==b for c in cases) for b in range(3)],
               original_round53_cases_checked=len(ORIGINAL),
               scope='full eligible E, actual first-exit records, original atomic P, common uniform beta',
               general_constant_proved=False,finite_tests_prove_uniform_bound=False,
               elapsed_seconds=time.monotonic()-started,
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',len(cases),out['batch_counts'],'seconds',round(out['elapsed_seconds'],2),flush=True)


if __name__ == '__main__': main()
