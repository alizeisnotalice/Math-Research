"""Exact new finite obstacle checks; no continuum/actual cube claim."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

BASE = Path(__file__).parent
REG = BASE / "fixed_jump_obstacle_interface_registration_20261007.json"
COUNT = 0

def check(ok, label):
    global COUNT
    COUNT += 1
    if not ok:
        raise AssertionError(label)

def mv(a, x):
    return [sum((v * w for v, w in zip(row, x)), F()) for row in a]

def solve(a, b):
    z = [row[:] + [v] for row, v in zip(a, b)]
    for k in range(len(z)):
        pivot = next(i for i in range(k, len(z)) if z[i][k])
        z[k], z[pivot] = z[pivot], z[k]
        scale = z[k][k]
        z[k] = [v / scale for v in z[k]]
        for i in range(len(z)):
            if i != k:
                scale = z[i][k]
                z[i] = [v - scale * w for v, w in zip(z[i], z[k])]
    return [row[-1] for row in z]

def components(j):
    unseen = set(range(len(j)))
    out = []
    while unseen:
        todo = [min(unseen)]
        found = set(todo)
        while todo:
            i = todo.pop()
            for k, v in enumerate(j[i]):
                if v and k not in found:
                    found.add(k)
                    todo.append(k)
        unseen -= found
        out.append(sorted(found))
    return out

def obstacle(j, nu):
    m = len(j)
    degree = list(map(sum, j))
    lap = [[degree[i] if i == k else -j[i][k] for k in range(m)] for i in range(m)]
    comps = components(j)
    if any(sum((nu[i] for i in c), F()) > len(c) for c in comps):
        return None
    u = [F()] * m
    trace = []
    for c in comps:
        mass = sum((nu[i] for i in c), F())
        if mass == len(c):
            a = c[:-1]
            values = solve([[lap[i][k] for k in a] for i in a], [nu[i] - 1 for i in a])
            values += [F()]
            shift = min(values)
            for i, v in zip(c, values):
                u[i] = v - shift
            trace.append({"regime": "critical", "component": len(c)})
            continue
        active = {i for i in c if nu[i] > 1}
        sizes = []
        while active:
            a = sorted(active)
            check(len(a) < len(c), "strict proper active set")
            values = solve([[lap[i][k] for k in a] for i in a], [nu[i] - 1 for i in a])
            previous = u[:]
            for i, v in zip(a, values):
                u[i] = v
            check(all(u[i] > 0 for i in active), "positive active solution")
            check(all(v >= w for v, w in zip(u, previous)), "monotone active expansion")
            sizes.append(len(a))
            lu = mv(lap, u)
            more = {i for i in c if i not in active and nu[i] - lu[i] > 1}
            if not more:
                break
            active |= more
        trace.append({"regime": "strict", "component": len(c), "active_sizes": sizes})
    return lap, u, trace

def digest(values):
    return hashlib.sha256("|".join(str(v) for v in values).encode()).hexdigest()

def run(m):
    j = [[F(1 + (i + k) % 3, 8 * m) if i != k else F() for k in range(m)] for i in range(m)]
    one = lambda mass: [F(mass)] + [F()] * (m - 1)
    profiles = [
        ("strict_single_3", one(3)),
        ("strict_single_half_capacity", one(F(m, 2))),
        ("strict_mixed", [F(2)] + [F(1, 4)] * (m - 1)),
        ("strict_alternating", [F(3, 2) if i % 3 == 0 else F(1, 8) for i in range(m)]),
        ("strict_arithmetic", [F((7 * i + 3) % 17, 12) for i in range(m)]),
        ("strict_uniform", [F(1, 4)] * m),
        ("critical_single", one(m)),
        ("critical_uniform", [F(1)] * m),
        ("over_capacity", one(m + 1)),
    ]
    receipts = []
    exterior_seen = False
    before = COUNT
    for name, nu in profiles:
        answer = obstacle(j, nu)
        if name == "over_capacity":
            check(answer is None, "global capacity rejection")
            receipts.append({"profile": name, "status": "REJECT_AS_REQUIRED"})
            continue
        check(answer is not None, "admissible source solved")
        lap, u, trace = answer
        lu = mv(lap, u)
        mu = [v - w for v, w in zip(nu, lu)]
        omega = [i for i, v in enumerate(u) if v > 0]
        outside = [i for i, v in enumerate(u) if v == 0]
        check(all(v >= 0 for v in u), "nonnegative odometer")
        check(all(0 <= v <= 1 for v in mu), "whole cap")
        check(all(mu[i] == 1 for i in omega), "saturation")
        check(sum(mu) == sum(nu), "whole mass conservation")
        check(len(omega) <= sum(nu), "total volume")
        good = [nu[i] if i in outside else F() for i in range(m)]
        bad = [nu[i] - good[i] for i in range(m)]
        mubad = [v - w for v, w in zip(bad, lu)]
        check(all(0 <= v <= 1 for v in good), "good density cap")
        check(all(0 <= v <= w for v, w in zip(mubad, mu)), "bad positive cap")
        check(sum(mubad) == sum(bad), "bad mass conservation")
        check(len(omega) <= sum(bad), "bad pays volume")
        check(sum(bad) + sum(good) == sum(nu), "source once")
        check(all(mubad[i] == sum((j[i][k] * u[k] for k in omega), F()) for i in outside), "full exterior absorption")
        flux = sum((mubad[i] for i in outside), F())
        exterior_seen |= flux > 0
        receipts.append({"profile": name, "status": "PASS", "mass": str(sum(nu)),
                         "bad_mass": str(sum(bad)), "good_mass": str(sum(good)),
                         "omega_states": len(omega), "exterior_bad_absorption": str(flux),
                         "u_sha256": digest(u), "mu_sha256": digest(mu), "trace": trace})
    check(exterior_seen, "actual nonzero exterior absorption tested")
    split = m // 2
    disconnected = [[v if (i < split) == (k < split) else F() for k, v in enumerate(row)] for i, row in enumerate(j)]
    nu = one(split + 1)
    check(sum(nu) <= m, "global capacity alone would allow source")
    check(obstacle(disconnected, nu) is None, "closed component capacity rejection")
    receipts.append({"profile": "disconnected_over_component_capacity", "status": "REJECT_AS_REQUIRED"})
    return {"states": m, "status": "PASS", "checks": COUNT - before, "profiles": receipts}

registration = json.loads(REG.read_text())
rounds = [run(m) for m in registration["rounds"]]
result = {"status": "PASS_EXACT", "predicates": COUNT, "rounds": rounds,
          "scope": registration["scope"],
          "registration_sha256": hashlib.sha256(REG.read_bytes()).hexdigest(),
          "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out = BASE / "fixed_jump_obstacle_interface_guard_results_20261007.json"
out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"status": result["status"], "predicates": COUNT, "output": str(out)}))
