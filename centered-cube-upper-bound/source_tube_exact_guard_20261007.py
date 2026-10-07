"""Exact rational finite tube components; no experiment or endpoint claim.

Single segments: rational determinant sums for the finite zonotope parallelotope
partition, with an independent planar convex-hull check. Finite tree tubes:
exact planar polygon inclusion-exclusion, extruded in the other n-2 coordinates.
Orthogonal trees also use an independent coordinate-cell union computation.
No high-dimensional bounding box is substituted for a tube volume.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json


def rat(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def cross(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def hull(points):
    pts = sorted(set(points))
    if len(pts) <= 1:
        return pts
    low, up = [], []
    for p in pts:
        while len(low) >= 2 and cross(low[-2], low[-1], p) <= 0:
            low.pop()
        low.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return low[:-1]+up[:-1]


def area(poly):
    if len(poly) < 3:
        return F(0)
    return abs(sum(poly[i][0]*poly[(i+1) % len(poly)][1]
                   -poly[i][1]*poly[(i+1) % len(poly)][0]
                   for i in range(len(poly))))/2


def clip(subject, other):
    """Convex polygon intersection, exact half-plane clipping."""
    out = list(subject)
    for j, a in enumerate(other):
        b = other[(j+1) % len(other)]
        if not out:
            break
        old, out = out, []
        for i, p in enumerate(old):
            q = old[(i+1) % len(old)]
            cp, cq = cross(a, b, p), cross(a, b, q)
            if cp >= 0:
                out.append(p)
            if (cp >= 0) != (cq >= 0):
                t = cp/(cp-cq)
                out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
        if out:
            out = hull(out)
    return out


def segment_polygon(p, q, r):
    h = r/2
    return hull([(v[0]+dx, v[1]+dy) for v in (p, q)
                 for dx in (-h, h) for dy in (-h, h)])


def polygon_union_ie(polygons):
    total, terms = F(0), []
    for k in range(1, len(polygons)+1):
        subtotal = F(0)
        for ids in combinations(range(len(polygons)), k):
            poly = polygons[ids[0]]
            for j in ids[1:]:
                poly = clip(poly, polygons[j])
                if not poly:
                    break
            subtotal += area(poly)
        total += subtotal if k % 2 else -subtotal
        terms.append(rat(subtotal))
    return total, terms


def rectangle_union_cells(rects):
    """Independent exact axis-grid integration of a rectangle union."""
    xs = sorted({v for box in rects for v in (box[0], box[1])})
    ys = sorted({v for box in rects for v in (box[2], box[3])})
    total = F(0)
    for a, b in zip(xs, xs[1:]):
        for c, d in zip(ys, ys[1:]):
            x, y = (a+b)/2, (c+d)/2
            if any(lo <= x <= hi and bot <= y <= top for lo, hi, bot, top in rects):
                total += (b-a)*(d-c)
    return total


def sparse_det(matrix):
    """Exact Laplace elimination on sparse columns; fallback Gaussian."""
    if not matrix:
        return F(1)
    n = len(matrix)
    for j in range(n):
        ids = [i for i in range(n) if matrix[i][j]]
        if not ids:
            return F(0)
        if len(ids) == 1:
            i = ids[0]
            minor = [row[:j]+row[j+1:] for k, row in enumerate(matrix) if k != i]
            return (-1 if (i+j) % 2 else 1)*matrix[i][j]*sparse_det(minor)
    m = [row[:] for row in matrix]
    ans = F(1)
    for j in range(n):
        i = next((i for i in range(j, n) if m[i][j]), None)
        if i is None:
            return F(0)
        if i != j:
            m[i], m[j] = m[j], m[i]
            ans = -ans
        pivot = m[j][j]
        ans *= pivot
        for i in range(j+1, n):
            c = m[i][j]/pivot
            for k in range(j+1, n):
                m[i][k] -= c*m[j][k]
            m[i][j] = F(0)
    return ans


def zonotope_segment_volume(v, r):
    """Finite rational parallelotope determinant certificate for Q+[0,v]."""
    n = len(v)
    columns = [[r if i == j else F(0) for i in range(n)] for j in range(n)]+[v]
    dets = []
    for ids in combinations(range(n+1), n):
        mat = [[columns[j][i] for j in ids] for i in range(n)]
        dets.append(abs(sparse_det(mat)))
    return sum(dets), dets


def l1(p, q):
    return sum(abs(a-b) for a, b in zip(p, q))


def tree_record(n, name, edges, r, orthogonal):
    # Each edge begins at a previously reached vertex: exact finite graph guard.
    reached = {edges[0][0]}
    for p, q in edges:
        assert p in reached
        reached.add(q)
    length = sum(l1(p, q) for p, q in edges)
    polygons = [segment_polygon(p, q, r) for p, q in edges]
    volume2, terms = polygon_union_ie(polygons)
    exact = volume2*r**(n-2)
    bound = r**n+length*r**(n-1)
    assert exact <= bound
    independent = None
    if orthogonal:
        rects = []
        for p, q in edges:
            assert p[0] == q[0] or p[1] == q[1]
            rects.append((min(p[0], q[0])-r/2, max(p[0], q[0])+r/2,
                          min(p[1], q[1])-r/2, max(p[1], q[1])+r/2))
        independent = rectangle_union_cells(rects)
        assert independent == volume2
    # Projection product is exact: all segment vertices have coordinates 3..n=0.
    return dict(n=n, case=name, method="exact_2D_polygon_IE_times_r_power_n_minus_2",
                edges=[[[rat(x) for x in p], [rat(x) for x in q]] for p, q in edges],
                full_side=rat(r), length_l1=rat(length), volume2=rat(volume2),
                exact_volume_n=rat(exact), network_bound=rat(bound), slack=rat(bound-exact),
                inclusion_exclusion_unsigned_level_sums=terms,
                independent_axis_cell_union_area=rat(independent) if independent is not None else None,
                source_scope="finite_connected_segment_network_only", pass_=True)


def main():
    records = []
    O = (F(0), F(0))
    orth = [(O, (F(3,2), F(0))), ((F(3,2), F(0)), (F(3,2), F(1,3))),
            ((F(3,2), F(0)), (F(3,2), F(-2,3))), (O, (F(-1,2), F(0))),
            ((F(-1,2), F(0)), (F(-1,2), F(5,4)))]
    oblique = [(O, (F(3,4), F(2,3))), ((F(3,4), F(2,3)), (F(5,4), F(-1,4))),
               ((F(3,4), F(2,3)), (F(-1,3), F(7,6))), (O, (F(-2,3), F(-1,2)))]
    for n in (4, 16, 64):
        r = F(n+1, 2*n)
        vectors = {"axis": [F(3,2)]+[F(0)]*(n-1),
                   "all_coordinate_diagonal": [F((-1)**j*(j+1), n*(n+1)) for j in range(n)]}
        for name, v in vectors.items():
            exact, dets = zonotope_segment_volume(v, r)
            length = sum(abs(x) for x in v)
            expected = r**n+r**(n-1)*length
            assert exact == expected
            records.append(dict(n=n, case=name, method="rational_zonotope_parallelotope_determinants",
                                full_side=rat(r), vector=[rat(x) for x in v],
                                determinant_absolute_terms=[rat(x) for x in dets],
                                exact_volume_n=rat(exact), segment_identity_rhs=rat(expected), pass_=True))
        # Planar segment component independently reconstructed by convex hull/shoelace.
        v2 = (F(-3,7), F(5,6))
        poly = segment_polygon(O, v2, r)
        exact2 = area(poly)
        det2, _ = zonotope_segment_volume(list(v2), r)
        assert exact2 == det2 == r*r+r*l1(O, v2)
        records.append(dict(n=n, case="independent_planar_diagonal_slice",
                            method="convex_hull_shoelace_and_independent_determinants",
                            polygon=[[rat(x) for x in p] for p in poly],
                            volume2=rat(exact2), exact_volume_n=rat(exact2*r**(n-2)), pass_=True))
        records.append(tree_record(n, "orthogonal_branched_tree", orth, r, True))
        records.append(tree_record(n, "oblique_branched_tree", oblique, r, False))
    # Sanity guards for the independent polygon engine on elementary rationals.
    a = segment_polygon(O, O, F(1))
    b = segment_polygon((F(1,2), F(0)), (F(1,2), F(0)), F(1))
    assert area(clip(a, b)) == F(1,2)
    assert polygon_union_ie([a, b])[0] == F(3,2)
    result = dict(status="PASS_EXACT_FINITE_TUBE_COMPONENTS", rounds=[4,16,64],
                  records=records, count=len(records),
                  scope="Pure rational finite component certificates only. Single full-dimensional segments use a zonotope determinant partition; planar networks use exact 2D polygon union times a common cube slice. No bounding-box substitution, arbitrary rectifiable/fractal extension, actual FIRST certification, or new asymptotic theorem.",
                  numerical_precision="fractions.Fraction throughout; no float geometry",
                  output="source_tube_exact_guard_20261007_results.json")
    path = Path(__file__).with_name(result['output'])
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(dict(status=result['status'], rounds=result['rounds'], count=len(records), output=str(path))))


if __name__ == "__main__":
    main()
