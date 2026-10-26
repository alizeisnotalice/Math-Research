"""Independent exact precheck helpers. No producer imports or numeric kernel."""
import sys
sys.dont_write_bytecode = True
import hashlib, itertools, json, math, time
from collections import Counter
from fractions import Fraction as Q

class Budget:
    def __init__(self):
        self.started = time.monotonic()
        self.matrix_operations = 0
    def guard(self):
        if time.monotonic() - self.started > 900:
            raise RuntimeError("registered actual-time cap")
    def determinant(self, rows):
        self.matrix_operations += 1
        assert self.matrix_operations <= 30000
        self.guard()
        a = [list(map(int, row)) for row in rows]
        n = len(a)
        assert all(len(row) == n for row in a)
        sign, previous = 1, 1
        for k in range(n - 1):
            pivot = next((i for i in range(k, n) if a[i][k]), None)
            if pivot is None:
                return 0
            if pivot != k:
                a[k], a[pivot] = a[pivot], a[k]
                sign = -sign
            value = a[k][k]
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    numerator = value * a[i][j] - a[i][k] * a[k][j]
                    quotient, remainder = divmod(numerator, previous)
                    assert remainder == 0
                    a[i][j] = quotient
                a[i][k] = 0
            previous = value
        return sign * a[-1][-1]

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def canonical_line(v):
    s = next(x for x in v if x)
    return tuple(v) if s > 0 else tuple(-x for x in v)

def generate_frame():
    base = [(1,) + s for s in itertools.product((-1, 1), repeat=5)]
    sparse = []
    for support in itertools.combinations(range(6), 4):
        for signs in itertools.product((-1, 1), repeat=3):
            row = [0] * 6
            for i, s in zip(support, (1,) + signs):
                row[i] = s
            sparse.append(tuple(row))
    unscaled = base + sparse
    columns80 = [tuple(80 * x for x in row) for row in base] + sparse
    assert len(base) == 32 and len(sparse) == 120
    assert len(set(columns80)) == 152
    return unscaled, columns80

def check_frame(saved, budget):
    directions, columns = generate_frame()
    assert saved["common_denominator"] == 80
    assert saved["column_numerators"] == [list(v) for v in columns]
    assert saved["beta"] == "1/80"
    gram = [[Q(sum(v[i] * v[j] for v in columns), 6400)
             for j in range(6)] for i in range(6)]
    assert all(gram[i][j] == (Q(2561, 80) if i == j else 0)
               for i in range(6) for j in range(6))
    marginal = [Q(sum(abs(v[i]) for v in columns), 80) for i in range(6)]
    assert marginal == [Q(33)] * 6
    colset = set(columns)
    # Adjacent swaps plus one sign change generate the entire B6 group.
    for generator in range(6):
        images = []
        for v in columns:
            w = list(v)
            if generator < 5:
                w[generator], w[generator + 1] = w[generator + 1], w[generator]
            else:
                w[0] = -w[0]
            images.append(canonical_line(w))
        assert set(images) == colset and len(set(images)) == 152
    budget.guard()
    return directions, columns, {
        "frame_columns": 152, "weighted_Gram_diagonal": str(gram[0][0]),
        "coordinate_Cauchy_rates": list(map(str, marginal)),
        "group_generators_checked": 6, "weighted_column_images_checked": 912,
        "full_signed_permutation_invariance": True,
    }

def check_pairs(directions, point, denominator, producer, budget):
    point = tuple(point)
    n = len(directions)
    all_pairs = []
    digest = hashlib.sha256()
    line_distances = [
        Q(dot(point, point) * dot(v, v) - dot(point, v) ** 2,
          denominator ** 2 * dot(v, v))
        for v in directions
    ]
    for i, j in itertools.combinations(range(n), 2):
        vectors = [directions[i], directions[j], point]
        gram3 = [[dot(u, v) for v in vectors] for u in vectors]
        det2 = budget.determinant([row[:2] for row in gram3[:2]])
        det3 = budget.determinant(gram3)
        assert det2 > 0 and det3 >= 0
        distance = Q(det3, denominator ** 2 * det2)
        digest.update(f"{i},{j}:{distance}\n".encode())
        all_pairs.append((i, j, distance))
    assert len(all_pairs) == 11476
    smallest = min(d for i, j, d in all_pairs)
    witnesses = [[i, j] for i, j, d in all_pairs if d == smallest]
    assert smallest > 0, "HARD REJECT: true S2 point"
    assert producer["status"] == "PASS_STRICT_OFF_TRUE_S2"
    assert producer["point_numerators"] == list(point)
    assert producer["denominator"] == denominator
    assert producer["pair_count"] == len(all_pairs)
    assert Q(producer["pair_min_squared"]) == smallest
    assert Q(producer["line_min_squared"]) == min(line_distances)
    assert producer["all_pair_distance_digest"] == digest.hexdigest()
    assert [w["indices"] for w in producer["pair_min_witnesses"]] == witnesses
    return {
        "method": "integer Bareiss detGram3/detGram2, not normal-equation expansion",
        "pair_count": len(all_pairs), "matrix_operations": budget.matrix_operations,
        "minimum_squared_pair_distance": str(smallest),
        "minimum_squared_line_distance": str(min(line_distances)),
        "minimizing_pairs": witnesses,
        "all_pair_distance_digest": digest.hexdigest(),
        "true_S2_rejection_gate_passed": True,
        "all_pair_exact_values": [[i, j, str(d)] for i, j, d in all_pairs],
    }

def check_group(point, budget):
    counter = Counter()
    for perm in itertools.permutations(range(6)):
        x = tuple(point[i] for i in perm)
        for signs in itertools.product((-1, 1), repeat=6):
            counter[tuple(a * b for a, b in zip(x, signs))] += 1
    assert sum(counter.values()) == 46080
    assert len(counter) == 3840 and set(counter.values()) == {12}
    reps = Counter()
    for x, multiplicity in counter.items():
        negative = tuple(-a for a in x)
        assert x != negative and counter[negative] == multiplicity
        reps[canonical_line(x)] += multiplicity
    assert len(reps) == 1920 and set(reps.values()) == {24}
    zero_positions = [sum(x[i] == 0 for x in counter) for i in range(6)]
    assert zero_positions == [640] * 6
    zero_first_reps = sum(x[0] == 0 for x in reps)
    assert zero_first_reps == 320
    representatives = sorted(reps)
    budget.guard()
    return representatives, {
        "group_actions": 46080, "distinct_point_images": 3840,
        "image_multiplicity": 12, "negative_representatives": 1920,
        "representative_multiplicity": 24, "zero_positions_image_counts": zero_positions,
        "first_coordinate_zero_representatives": zero_first_reps,
        "representative_rule": "first nonzero coordinate positive, not coordinate0 positive",
        "representatives_sha256": hashlib.sha256(
            json.dumps(representatives, separators=(",", ":")).encode()).hexdigest(),
    }

def check_scaled_fan(fan, columns80, budget):
    assert fan["ray_count"] == len(fan["rays"]) == 478
    assert fan["cone_count"] == len(fan["cones"]) == 13229
    assert fan["parent_count"] == 21
    rays, scales = [], []
    for index, row in enumerate(fan["rays"]):
        assert row["id"] == index
        primitive = row["primitive"]
        assert len(primitive) == 6
        true_rate = Q(sum(abs(dot(v, primitive)) for v in columns80), 80)
        assert true_rate > 0 and Q(row["q_true_152"]) == true_rate
        d, q = true_rate.denominator, true_rate.numerator
        r = [d * x for x in primitive]
        assert all(-(1 << 63) < x < (1 << 63) for x in r)
        assert 0 < q < (1 << 63)
        rays.append({"R": r, "Q": q})
        scales.append(d)
    cones = []
    weights_hash = hashlib.sha256()
    for row in fan["cones"]:
        ids = row["ray_ids"]
        assert len(ids) == len(set(ids)) == 6
        assert all(0 <= i < 478 for i in ids)
        primitive_det = row["abs_primitive_determinant"]
        assert primitive_det > 0
        scaled_det = primitive_det * math.prod(scales[i] for i in ids)
        assert 0 < scaled_det < (1 << 63)
        old_weight = Q(primitive_det)
        new_weight = Q(scaled_det)
        for i in ids:
            old_weight /= Q(fan["rays"][i]["q_true_152"])
            new_weight /= rays[i]["Q"]
        assert old_weight == new_weight
        cones.append({"ray_ids": ids, "scaled_determinant": scaled_det})
        weights_hash.update((str(new_weight) + "\n").encode())
    budget.guard()
    return rays, cones, {
        "all_rays_checked": 478, "all_simplex_weights_checked": 13229,
        "primitive_determinants_inherited_from_F050": 13229,
        "new_simplex_matrix_determinants": 0,
        "simultaneous_column_scale_and_rate_identity": True,
        "max_scaled_coordinate": max(abs(x) for row in rays for x in row["R"]),
        "max_scaled_rate": max(row["Q"] for row in rays),
        "max_scaled_determinant": max(row["scaled_determinant"] for row in cones),
        "exact_weight_stream_sha256": weights_hash.hexdigest(),
    }

if __name__ == "__main__":
    raise SystemExit("Preparation helper only: use separately pinned review entry after stable producer input.")
