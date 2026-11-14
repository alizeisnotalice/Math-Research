#!/usr/bin/env python3
"""Optimizer-free exact replay of all finite-catalogue LF capacities.

Usage: python3 replay_lf.py [path/to/r732]
The default source is sibling ../r732, so this also works in public r736.
Only Python's standard library is used. Reconstructs geometry from box cells,
reconstructs the LF constraints and independently clips them using every crossing
vertex pair, retaining only full-active-rank points. Reuses r732/replay_exact.py
for independent 64-term rational logarithm enclosures. Numerical solver status
is never read. The stored Wolfram report is an additional consistency gate;
completeness is also recomputed here, not inferred from its success flags.
"""
import hashlib
import importlib.util
import itertools
import json
import math
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))


def rref(rows, columns):
    a = [list(map(Q, row)) for row in rows]
    pivots = []
    for j in range(columns):
        pivot = next((i for i in range(len(pivots), len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        n = len(pivots)
        a[n], a[pivot] = a[pivot], a[n]
        c = a[n][j]
        a[n] = [x/c for x in a[n]]
        for i in range(len(a)):
            if i != n and a[i][j]:
                c = a[i][j]
                a[i] = [x-c*y for x, y in zip(a[i], a[n])]
        pivots.append(j)
    return a, pivots


def rank(rows, columns):
    return len(rref(rows, columns)[1])


def normalize(row):
    row = list(map(Q, row))
    scale = math.lcm(*(x.denominator for x in row))
    ints = [int(x*scale) for x in row]
    if not any(ints[:-1]):
        assert ints[-1] >= 0, 'infeasible constant constraint'
        return None
    divisor = math.gcd(*ints)
    return tuple(x//divisor for x in ints)


def atomize(boxes):
    boxes = [[[Q(str(v)) for v in side] for side in box] for box in boxes]
    dimension = len(boxes[0])
    assert all(len(box) == dimension and all(a < b for a, b in box) for box in boxes)
    axes = [sorted({v for box in boxes for v in box[j]}) for j in range(dimension)]
    intervals = [list(zip(axis, axis[1:])) for axis in axes]
    masses = {}
    for cell in itertools.product(*intervals):
        midpoint = [(a+b)/2 for a, b in cell]
        mask = sum(1 << i for i, box in enumerate(boxes)
                   if all(a < x < b for (a, b), x in zip(box, midpoint)))
        if mask:
            masses[mask] = masses.get(mask, Q(0)) + math.prod(b-a for a, b in cell)
    masks = sorted(masses)
    return masks, [masses[s] for s in masks], [math.prod(b-a for a, b in box) for box in boxes]


def verify_face(base, rec, wolfram):
    masks, mass, vol = atomize(base['boxes'])
    assert masks == base['atom_masks']
    assert mass == list(map(Q, base['atom_masses']))
    assert vol == list(map(Q, base['row_volumes']))
    assert len(base['boxes'][0]) == base['dimension']
    assert sum(mass) == Q(base['union_volume'])
    k, om = len(vol), sum(mass)
    f, ell = list(map(Q, base['source_f'])), list(map(Q, base['source_ell']))
    target, price = Q(base['target']), Q(base['exact_NP_price'])
    assert len(f) == len(ell) == len(mass) and 0 < target < 1
    av = [sum(m*f[j] for j, (m, s) in enumerate(zip(mass, masks)) if s & (1 << i))/vol[i] for i in range(k)]
    assert all(0 <= x <= 1 for x in f+ell)
    assert all(ell[j] <= av[i] for j, s in enumerate(masks) for i in range(k) if s & (1 << i))
    assert dot(mass, ell) >= target*om and dot(mass, f) == price
    density = [[Q(int(bool(s & (1 << i))))/vol[i] for i in range(k)] for s in masks]
    equations = [[Q(1)]*k+[Q(0), Q(1)], av+[Q(0), target]]
    equations += [density[j]+[Q(-1), Q(0)] for j in range(len(f)) if 0 < f[j] < 1]
    nv = k+1
    reduced, pivots = rref(equations, nv)
    assert all(any(row[:nv]) or row[-1] == 0 for row in reduced)
    free = [i for i in range(nv) if i not in pivots]
    d = len(free)
    origin = [Q(0)]*nv
    basis = [[Q(0)]*d for _ in range(nv)]
    for j, i in enumerate(free):
        basis[i][j] = Q(1)
    for row, i in zip(reduced, pivots):
        origin[i] = row[-1]
        basis[i] = [-row[j] for j in free]
    assert rec['case'] == base['case'] and Q(rec['target']) == target
    assert rec['free_dimension'] == d and rec['free_variables'] == free
    assert list(map(Q, rec['affine_origin'])) == origin
    assert [list(map(Q, row)) for row in rec['affine_basis']] == basis
    dmax = max(1/v for v in vol)
    constraints = []
    for i in range(k):
        c = [Q(0)]*nv
        c[i] = Q(-1)
        constraints.append((c, Q(0)))
    constraints += [([Q(0)]*k+[Q(-1)], Q(0)), ([Q(0)]*k+[Q(1)], dmax)]
    for j, value in enumerate(f):
        if value in (0, 1):
            c = density[j]+[Q(-1)]
            constraints.append((c if value == 0 else [-v for v in c], Q(0)))
    for s in range(1, 1 << k):
        constraints.append(([Q(int(bool(s & (1 << i)))) for i in range(k)]+[Q(0)],
                            sum(m for m, t in zip(mass, masks) if s & t)/om))
    limits = [dmax if i == k else Q(1) for i in free]
    rows = []
    for j, bound in enumerate(limits):
        for sign, rhs in [(-1, Q(0)), (1, bound)]:
            c = [Q(0)]*d
            c[j] = Q(sign)
            rows.append(normalize(c+[rhs]))
    seen = set(rows)
    for c, rhs in constraints:
        transformed = normalize([dot(c, [basis[i][j] for i in range(nv)]) for j in range(d)] + [rhs-dot(c, origin)])
        if transformed is not None and transformed not in seen:
            rows.append(transformed)
            seen.add(transformed)
    assert rows == [tuple(row) for row in rec['constraint_normals']]
    # Independent clipping: all crossing pairs are considered. Rank filtering
    # removes interior/diagonal intersections. This does not use the producer's
    # edge test or any recorded clipping trace.
    vertices = set(itertools.product(*[(Q(0), bound) for bound in limits]))
    done = rows[:2*d]
    for row in rows[2*d:]:
        residuals = {v: dot(row[:-1], v)-row[-1] for v in vertices}
        inside = [(v, t) for v, t in residuals.items() if t < 0]
        outside = [(v, t) for v, t in residuals.items() if t > 0]
        candidates = {v for v, t in residuals.items() if t <= 0}
        for v, a in inside:
            for w, b in outside:
                t = -a/(b-a)
                candidates.add(tuple(x+t*(y-x) for x, y in zip(v, w)))
        done.append(row)
        vertices = {v for v in candidates if rank([r[:-1] for r in done if dot(r[:-1], v) == r[-1]], d) == d}
        assert vertices, 'empty LF face'
    listed = {tuple(map(Q, v['free_coordinates'])) for v in rec['vertices']}
    assert len(listed) == len(rec['vertices']) == rec['extended_vertex_count']
    assert vertices == listed, 'incomplete or incorrect extended vertex list'
    projected = set()
    for v in rec['vertices']:
        z = tuple(map(Q, v['free_coordinates']))
        x = [origin[i]+dot(basis[i], z) for i in range(nv)]
        assert x == list(map(Q, v['pi']))+[Q(v['threshold'])]
        assert all(dot(c, x) <= rhs for c, rhs in constraints)
        assert all(dot(eq[:-1], x) == eq[-1] for eq in equations)
        projected.add(tuple(x[:k]))
    unique = [tuple(map(Q, pi)) for pi in rec['unique_pi']]
    assert unique == sorted(projected) and len(unique) == rec['unique_pi_count']
    affine_dimension = rank([[x-y for x, y in zip(v, unique[0])] for v in unique[1:]], k)
    assert wolfram['case'] == base['case'] and wolfram['unique_pi_count'] == len(unique)
    assert wolfram['projected_affine_dimension'] == affine_dimension
    assert wolfram['listed_vertices_exact_LF'] is True and wolfram['all_LF_pi_in_listed_convex_hull'] is True
    assert wolfram['affine_hull_tests'] == 2*(k-affine_dimension)
    assert len(wolfram['exact_linear_maxima']) == wolfram['affine_hull_tests']+wolfram['complete_hull_facets']
    assert all(Q(value) == 0 for value in wolfram['exact_linear_maxima'])
    return dict(case=base['case'], exact_box_geometry=True, exact_shared_primal=True,
                exact_LF_constraints=True, independent_complete_vertex_replay=True,
                projected_vertices=len(unique), extended_vertices=len(vertices),
                projected_affine_dimension=affine_dimension, wolfram_report_consistent=True)


def main():
    if not __debug__:
        raise RuntimeError('Run without -O or PYTHONOPTIMIZE; assertions are required.')
    started = time.monotonic()
    source_dir = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent/'r732'
    paths = {'source': source_dir/'capacity.json', 'inner_replay': source_dir/'replay_exact.py',
             'vertices': HERE/'lf_vertices.json', 'capacities': HERE/'lf_capacities.json',
             'wolfram': HERE/'verify_face_wolfram.json', 'wolfram_code': HERE/'verify_face.wl'}
    raw = {key: path.read_bytes() for key, path in paths.items()}
    hashes = {key: hashlib.sha256(value).hexdigest() for key, value in raw.items()}
    source, enum, cap, wl = [json.loads(raw[key]) for key in ['source', 'vertices', 'capacities', 'wolfram']]
    assert enum['status'] == cap['status'] == wl['status'] == 'complete'
    assert hashes['source'] == enum['source_sha256'] == cap['source_sha256'] == wl['source_sha256']
    assert hashes['vertices'] == cap['vertices_sha256']
    def indexed(rows):
        out = {r['case']: r for r in rows}
        assert len(out) == len(rows)
        return out
    bases, faces, reports, summaries = [indexed(rows) for rows in [source['results'], enum['results'], wl['results'], cap['LF_capacity_summaries']]]
    assert set(bases) == set(faces) == set(reports) == set(summaries)
    assert len(bases) == 27
    geometry = [verify_face(bases[i], faces[i], reports[i]) for i in sorted(bases)]
    spec = importlib.util.spec_from_file_location('r732_exact_replay', paths['inner_replay'])
    inner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inner)
    assert inner.TERMS == 64
    expected = {(i, j) for i, face in faces.items() for j in range(len(face['unique_pi']))}
    seen, results = set(), []
    same_fields = ['dimension', 'target', 'boxes', 'atom_masks', 'atom_masses', 'row_volumes',
                   'union_volume', 'exact_NP_price', 'source_f', 'source_ell']
    for row in cap['vertex_certificates']:
        i, j = row['source_case'], row['vertex']
        assert (i, j) in expected and (i, j) not in seen
        seen.add((i, j))
        assert all(row[key] == bases[i][key] for key in same_fields)
        assert list(map(Q, row['pi'])) == list(map(Q, faces[i]['unique_pi'][j]))
        result = inner.check(row)
        result.update(source_case=i, vertex=j, geometry_and_shared_primal_match=True)
        results.append(result)
    assert seen == expected and len(results) == 35
    outer = []
    for i in sorted(bases):
        rows = [r for r in cap['vertex_certificates'] if r['source_case'] == i]
        lo, hi = max(Q(r['capacity_lower']) for r in rows), max(Q(r['capacity_upper']) for r in rows)
        summary = summaries[i]
        assert Q(summary['target']) == Q(bases[i]['target'])
        assert summary['LF_vertex_count'] == len(rows)
        assert Q(summary['LF_capacity_lower']) == lo and Q(summary['LF_capacity_upper']) == hi
        assert Q(summary['LF_capacity_interval_width_upper']) >= hi-lo
        # Older wider lower bounds remain valid after boundary refinement.
        assert Q(summary['improvement_over_original_fixed_pi_lower']) <= lo-Q(bases[i]['capacity_upper'])
        outer.append(dict(case=i, all_LF_vertices_accounted=True, max_vertex_endpoints_verified=True,
                          LF_capacity_lower=str(lo), LF_capacity_upper=str(hi)))
    assert Q(bases[11]['target']) == Q(3, 4) and len(faces[11]['unique_pi']) == 1
    out = dict(status='complete', all_pass=True, faces=len(geometry), vertex_certificates=len(results),
               dependencies='Python standard library only', log_series_terms=64,
               optimizer_called=False, numerical_solver_success_used=False,
               independent_geometry_and_LF_completeness=True,
               source_paths={key: str(path.relative_to(source_dir)) if path.parent == source_dir else path.name for key, path in paths.items()},
               sha256=hashes, python=sys.version.split()[0],
               scope='Complete LF-face capacity intervals for this 27-case finite catalogue only; no central asymptotic estimate.',
               method='Exact cell geometry; shared-primal inequalities; exact affine reduction; independent all-crossing-pairs clipping with active-rank filtering; r732 64-term Fraction logarithm replay; max vertex endpoints.',
               case11_q_three_quarters_LF_singleton=True,
               test_minus_capacity_rule='Stored intervals may be wider than the refined capacity-derived interval; enclosure is required, endpoint equality is not.',
               elapsed_seconds=time.monotonic()-started, geometry_results=geometry,
               vertex_results=results, outer_results=outer)
    (HERE/'replay_lf.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({key: out[key] for key in ['status', 'all_pass', 'faces', 'vertex_certificates', 'log_series_terms', 'independent_geometry_and_LF_completeness', 'elapsed_seconds']}))


if __name__ == '__main__':
    main()
