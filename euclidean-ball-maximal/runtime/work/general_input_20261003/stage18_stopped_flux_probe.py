#!/usr/bin/env python3
"""Exact 1D stopped-mask flux calculations; no source/observer sampling.

Complete original atom superlevels/bands and canonical [J,K] are decomposed
in Fraction cells. An independent source-pair/U integration verifies the
opening/displacement decomposition, signed linear-quadratic cancellation,
entry/terminal clocks, and positive parts before/after the U mean.
All sums have the finite K<=L mask; every source-pair/layer mean is saved.
"""
from fractions import Fraction as F
from bisect import bisect_right
from pathlib import Path
import json


def merge(intervals):
    out = []
    for lo, hi in sorted(intervals):
        if lo >= hi:
            continue
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(hi, out[-1][1]))
        else:
            out.append((lo, hi))
    return out


def full_E(atoms, alpha):
    intervals = []
    for i in range(len(atoms)):
        m = F(0)
        for k in range(i, len(atoms)):
            m += atoms[k][1]
            r = m/(2*alpha)
            intervals.append((atoms[k][0]-r, atoms[i][0]+r))
    return merge(intervals)


def difference(A, B):
    out = []
    for lo, hi in A:
        current = lo
        for bl, bh in B:
            if bh <= current or bl >= hi:
                continue
            if bl > current:
                out.append((current, min(bl, hi)))
            current = max(current, bh)
            if current >= hi:
                break
        if current < hi:
            out.append((current, hi))
    return merge(out)


def flat(intervals):
    return [e for pair in intervals for e in pair]


def inside(x, endpoints):
    return bisect_right(endpoints, x) % 2 == 1


def run_scope(name, atoms, alpha, scope):
    beta = alpha/2
    E = full_E(atoms, alpha)
    if scope == "band":
        E = difference(E, full_E(atoms, 2*alpha))
    endpoints_E = flat(E)
    L = 1
    while F(1, 2**L) > min(w for y, w in atoms)/4:
        L += 1
    p = [F(1, 2**j) for j in range(L+1)]
    R = [pp/beta for pp in p]  # n=1
    points = set(endpoints_E)
    for r in R:
        for y, w in atoms:
            points.update((y-r, y+r))
    points = sorted(points)
    masks = [[] for _ in range(L+2)]
    block_integral, terminal_integral, entry_integral = F(0), F(0), F(0)
    eligible_volume, first_layer_volume = F(0), F(0)
    rows = []
    for lo, hi in zip(points, points[1:]):
        x = (lo+hi)/2
        if not inside(x, endpoints_E):
            continue
        m = [sum((w for y, w in atoms if abs(x-y) < r), F(0)) for r in R]
        J = next((j for j in range(1, L+1) if m[j] > p[j]), None)
        K = next((j for j in range(1, L+1) if m[j] > 2*p[j]), None)
        assert J is not None and K is not None and J <= K
        if J == 1:
            first_layer_volume += hi-lo
            continue
        bK = m[K]*(1-m[K])/p[K]
        bprev = m[J-1]*(1-m[J-1])/p[J-1]
        block = bK-bprev
        assert F(1, 2) < block <= 4
        block_integral += (hi-lo)*block
        terminal_integral += (hi-lo)*bK
        entry_integral += (hi-lo)*bprev
        eligible_volume += hi-lo
        for j in range(J, K+1):
            masks[j].append((lo, hi))
        rows.append(dict(lo=str(lo), hi=str(hi), J=J, K=K, block=str(block)))
    masks = [merge(v) for v in masks]
    fm = [flat(v) for v in masks]
    H, opening, linear, quadratic, terminal, entry = (F(0) for _ in range(6))
    after_U_positive, after_U_negative, before_U_positive, before_U_negative = (F(0) for _ in range(4))
    pair_layer_fluxes = []
    coupled_cells = 0
    for y, wy in atoms:
        base = {F(-1), F(1)}
        for j in range(2, L+1):
            for endpoint in fm[j]:
                for r in (R[j], R[j-1]):
                    u = (endpoint-y)/r
                    if -1 < u < 1:
                        base.add(u)
        for s, ws in atoms:
            layer_h = [F(0) for _ in range(L+1)]
            layer_before_positive = [F(0) for _ in range(L+1)]
            layer_before_negative = [F(0) for _ in range(L+1)]
            breaks = set(base)
            for j in range(1, L+1):
                for u in ((s-y)/R[j]-1, (s-y)/R[j]+1):
                    if -1 < u < 1:
                        breaks.add(u)
            breaks = sorted(breaks)
            for ul, uh in zip(breaks, breaks[1:]):
                u = (ul+uh)/2
                U_weight = (uh-ul)/2  # uniform U on (-1,1)
                weight = wy*ws*U_weight
                pathH, pathD, pathlin, pathquad, pathT, pathA = (0 for _ in range(6))
                for j in range(2, L+1):
                    xj, xp = y+R[j]*u, y+R[j-1]*u
                    fj, fp = int(inside(xj, fm[j])), int(inside(xp, fm[j]))
                    ep, ej = int(abs(xp-s) >= R[j-1]), int(abs(xj-s) >= R[j])
                    assert ej >= ep
                    displacement = ep*(fj-fp)
                    pathH += displacement
                    layer_h[j] += U_weight*displacement
                    layer_before_positive[j] += U_weight*max(0, displacement)
                    layer_before_negative[j] += U_weight*max(0, -displacement)
                    pathD += fj*(ej-ep)
                    pathlin += fj-fp
                    pathquad += fj*(1-ej)-fp*(1-ep)
                    # Fixed-center mask-clock jumps at the prior/current scale.
                    pathA += ep*max(0, fp-int(inside(xp, fm[j-1])))
                    pathT += ej*max(0, fj-int(inside(xj, fm[j+1])))
                assert 0 <= pathD <= 1
                assert pathH+pathD == pathlin-pathquad == pathT-pathA
                H += weight*pathH
                opening += weight*pathD
                linear += weight*pathlin
                quadratic += weight*pathquad
                terminal += weight*pathT
                entry += weight*pathA
                coupled_cells += 1
            assert all(layer_h[j] == layer_before_positive[j]-layer_before_negative[j]
                       for j in range(2, L+1))
            after_U_positive += wy*ws*sum((max(F(0), h) for h in layer_h[2:]), F(0))
            after_U_negative += wy*ws*sum((max(F(0), -h) for h in layer_h[2:]), F(0))
            before_U_positive += wy*ws*sum(layer_before_positive[2:], F(0))
            before_U_negative += wy*ws*sum(layer_before_negative[2:], F(0))
            pair_layer_fluxes.append(dict(y=str(y), s=str(s), y_mass=str(wy), s_mass=str(ws),
                layers=[dict(j=j, U_mean_signed_flux=str(layer_h[j]),
                             U_mean_positive_part=str(layer_before_positive[j]),
                             U_mean_negative_part=str(layer_before_negative[j]))
                        for j in range(2, L+1)]))
    assert 0 <= opening <= 1
    assert block_integral == (2/beta)*(H+opening)
    assert terminal_integral == (2/beta)*terminal
    assert entry_integral == (2/beta)*entry
    assert block_integral == terminal_integral-entry_integral
    assert after_U_positive-after_U_negative == H == before_U_positive-before_U_negative
    assert after_U_positive <= before_U_positive and after_U_negative <= before_U_negative
    if name == "separated_equal_atoms":
        assert opening == 0 and H > 0
    return dict(input=name, scope=scope, alpha=str(alpha), beta=str(beta), L=L,
                original_volume=str(sum((hi-lo for lo, hi in E), F(0))),
                J1_volume=str(first_layer_volume), eligible_volume=str(eligible_volume),
                stopped_signed_block_integral=str(block_integral),
                normalized_mask_displacement_flux=str(H),
                normalized_once_only_opening=str(opening),
                normalized_linear_flux=str(linear), normalized_quadratic_flux=str(quadratic),
                normalized_terminal_clock=str(terminal), normalized_entry_clock=str(entry),
                normalized_positive_after_U_mean=str(after_U_positive),
                normalized_negative_after_U_mean=str(after_U_negative),
                normalized_positive_before_U_mean=str(before_U_positive),
                normalized_negative_before_U_mean=str(before_U_negative),
                source_pair_layer_direction_means=pair_layer_fluxes,
                opening_space_budget=str(2/beta),
                coupling_U_cells=coupled_cells, observer_cells=rows,
                masks=[[list(map(str, interval)) for interval in v] for v in masks])


def main():
    p = F(1, 2**16)
    cases = [
        ("separated_equal_atoms", [(F(0), F(1, 2)), (F(10), F(1, 2))], F(2)),
        ("irregular_nine_atoms", [(y, F(1, 9)) for y in
            (F(-1), F(-39, 50), F(-8, 25), F(1, 10), F(3, 25),
             F(13, 25), F(1), F(26, 25), F(3, 2))], F(7, 10)),
        ("deep_near_threshold_cloud", [(F(0), F(2, 5)*p),
            (F(4, 5)*p, F(3, 5)*p+p**3),
            (F(8, 5)*p, p-p**3), (F(10), 1-2*p)], F(1)),
    ]
    rounds = []
    for name, atoms, alpha in cases:
        atoms = sorted(atoms)
        assert sum((w for y, w in atoms), F(0)) == 1
        rounds.append(dict(name=name, atoms=[dict(location=str(y), mass=str(w)) for y, w in atoms],
            scopes=[run_scope(name, atoms, alpha, scope) for scope in ("full", "band")]))
    result = dict(scope="Three actual canonical input calculations, full original superlevel and band for each; exact Fraction spatial/source-pair/U integrals with finite K<=L masks, including every source-pair/layer U mean and positive parts before/after U averaging; no Monte Carlo, entropy, or uniform weak certification",
                  rounds=rounds)
    root = Path(__file__).resolve().parents[2]
    out = root/"output/general_input_20261003/stage18_stopped_flux_probe.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(output=str(out), rounds=[dict(name=r["name"], scopes=[
        {key:s[key] for key in ("scope", "L", "stopped_signed_block_integral", "normalized_once_only_opening",
                                "normalized_mask_displacement_flux", "normalized_positive_after_U_mean",
                                "normalized_negative_after_U_mean", "normalized_positive_before_U_mean",
                                "normalized_negative_before_U_mean", "coupling_U_cells")} for s in r["scopes"]])
        for r in rounds]), indent=2))


if __name__ == "__main__":
    main()
