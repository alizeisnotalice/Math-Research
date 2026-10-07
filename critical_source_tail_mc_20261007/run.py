#!/usr/bin/env python3
"""Fixed-budget original-source, continuous-winner occupation Monte Carlo.
Exact dyadic response oracle; ideal-law MC bands are not rounding certificates.
"""
import gzip
import hashlib
import importlib.util
import json
import math
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
sys.set_int_max_str_digits(0)
sys.dont_write_bytecode = True  # Read-only import must not create old-folder caches.


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def ceildiv(a, b):
    return -((-a) // b)


class IntegerArrival:
    def __init__(self, inp, tau):
        self.n, self.A = inp['n'], inp['A']
        self.ks, self.qs = inp['ks'], inp['qs']
        self.weights = list(map(F, inp['weights']))
        assert self.weights == [F(1, 2), F(1, 6), F(1, 6), F(1, 6)]
        qp = [q**self.n for q in self.qs]
        self.den = 6 * math.prod(qp)
        self.coef = [(3 if l == 0 else 1) * math.prod(qp[j] for j in range(4) if j != l) for l in range(4)]
        self.tau = F(tau)
        self.old_data = dict(n=self.n, A=self.A, ks=self.ks, qs=self.qs, weights=self.weights)

    def evaluate(self, source_codes, source_k, offsets):
        ratios = [float(v).as_integer_ratio() for v in offsets]
        D = max(8, *(d for _, d in ratios))
        dx = [a * (D // d) for a, d in ratios]
        x = [int(y) * (D // source_k) + v for y, v in zip(source_codes, dx)]
        counts, events = [], defaultdict(Counter)
        for l, k in enumerate(self.ks):
            row = []
            for i, xi in enumerate(x):
                low = max(-self.A * k, ceildiv((xi - D // 2) * k, D))
                high = min(self.A * k, ((xi + D // 2) * k) // D)
                row.append(max(0, high - low + 1))
                low2 = max(-self.A * k, ceildiv((xi - D) * k, D))
                high2 = min(self.A * k, ((xi + D) * k) // D)
                for atom in range(low2, high2 + 1):
                    key = 2 * abs(xi - atom * (D // k))
                    if D < key <= 2 * D:
                        events[key][(l, i)] += 1
            counts.append(row)
        nonzero = [math.prod(c for c in row if c) for row in counts]
        zeros = [sum(c == 0 for c in row) for row in counts]
        bestmass, bestkey, bestpow, bestcounts = None, None, None, None
        log_values = []
        tie_count = 0
        for key in sorted({D, 2 * D, *events}):
            for (l, i), increment in events.get(key, {}).items():
                prev = counts[l][i]
                if prev:
                    nonzero[l] //= prev
                else:
                    zeros[l] -= 1
                counts[l][i] = prev + increment
                nonzero[l] *= counts[l][i]
            mass = sum(coef * (0 if zeros[l] else nonzero[l]) for l, coef in enumerate(self.coef))
            power = key**self.n
            if bestmass is None or mass * bestpow > bestmass * power:
                bestmass, bestkey, bestpow = mass, key, power
                bestcounts = [row[:] for row in counts]
                tie_count = 1
            elif mass * bestpow == bestmass * power:
                tie_count += 1
            log_values.append((key, math.log(mass) - math.log(self.den) - self.n * math.log(key / D) if mass else -math.inf))
        Dpow = D**self.n
        lhs = bestmass * Dpow * self.tau.denominator
        rhs = self.tau.numerator * self.den * bestpow
        E = lhs > rhs
        g = rhs / lhs if E else 0.0
        assert 0 <= g <= 1
        capture = 2 * max(abs(v) for v in dx) <= bestkey
        logM = math.log(bestmass) - math.log(self.den) - self.n * math.log(bestkey / D)
        loggap = logM - math.log(self.tau.numerator) + math.log(self.tau.denominator)
        other = [v for key, v in log_values if key != bestkey]
        runner_gap = max(0.0, logM - max(other)) if other else math.inf
        assert sum(c * math.prod(row) for c, row in zip(self.coef, bestcounts)) == bestmass
        return dict(R_key=bestkey, dyadic_denominator=D, R_float=bestkey / D,
                    full_mass_numerator=bestmass, full_mass_denominator=self.den,
                    full_capture_counts=bestcounts, E=E, source_captured=capture, g=g,
                    logM=logM, log_gap_tau=loggap, near_threshold=abs(loggap) <= 1e-10,
                    runnerup_log_gap=runner_gap, near_max=runner_gap <= 1e-10,
                    tied_max_candidates=tie_count, candidate_count=len(log_values),
                    simultaneous_arrival_groups=sum(sum(v.values()) > 1 for v in events.values()))


def partition(old, y, width, A):
    segments, p = old.high_sets(y, width, A)
    local = [(left - y, right - y, high) for left, right, high in segments]
    assert sum(right - left for left, right, _ in local) == width
    assert sum(right - left for left, right, high in local if high) == width * p
    return local, p


def sampler_guards(old, oracle, source_codes, source_k, T):
    normal_y = [F(int(v), source_k) for v in source_codes]
    plans = [(normal_y, F(1), None, 0),
             ([F(oracle.A)] * oracle.n, F(1), None, 0),
             ([F(-oracle.A)] * oracle.n, F(5, 4), 0, -1),
             (normal_y, F(7, 4), oracle.n - 1, 1)]
    records = []
    for y, width, face, sign in plans:
        free = [i for i in range(oracle.n) if i != face]
        ps = []
        for i in free:
            _, p = partition(old, y[i], width, oracle.A)
            ps.append(p)
        Z = [1 - p + p * T for p in ps]
        prior = old.status_law(ps)
        tilted = old.status_law([p * T / z for p, z in zip(ps, Z)])
        proposal = [F(1, 4) * a + F(3, 4) * b for a, b in zip(prior, tilted)]
        normalizer = math.prod(Z)
        RN = [1 / (F(1, 4) + F(3, 4) * T**h / normalizer) for h in range(len(prior))]
        assert sum(q * w for q, w in zip(proposal, RN)) == 1
        assert sum(h * q * w for h, (q, w) in enumerate(zip(proposal, RN))) == sum(ps)
        second = sum(q * w * w for q, w in zip(proposal, RN))
        assert second <= 4 and all(0 < w <= 4 for w in RN)
        records.append(dict(y=[str(v) for v in y], width=str(width), face=face, sign=sign,
                            p_i=[str(p) for p in ps], normalizer=str(normalizer),
                            RN_mass='1', RN_first_high_moment=str(sum(ps)),
                            RN_second_moment=str(second), all_exact_checks_passed=True))
    return records


def receiver(old, oracle, codes, source_k, rng, T, Cn):
    shell = rng.random() >= 1 / Cn
    width = math.exp(math.log(2) * rng.random()) if shell else 1.0
    facecode = int(rng.integers(0, 2 * oracle.n)) if shell else None
    face = facecode // 2 if shell else None
    sign = (1 if facecode % 2 else -1) if shell else 0
    tilted_branch = rng.random() >= .25
    y = [F(int(v), source_k) for v in codes]
    wf = F.from_float(width)
    offsets = [0.0] * oracle.n
    free, pvals, hs = [], [], []
    near_boundary = 0
    logZ = 0.0
    for i in range(oracle.n):
        if i == face:
            offsets[i] = sign * width / 2
            continue
        segments, pf = partition(old, y[i], wf, oracle.A)
        p = float(pf)
        Zi = 1 - p + p * T
        if tilted_branch:
            high = rng.random() < p * T / Zi
            choices = [(a, b) for a, b, flag in segments if bool(flag) == bool(high)]
            length = sum(b - a for a, b in choices)
            assert length > 0
            target = rng.random() * float(length)
            chosen = choices[-1]
            for a, b in choices:
                span = float(b - a)
                if target < span:
                    chosen = (a, b)
                    break
                target -= span
            a, b = chosen
            offset = float(a) + target
        else:
            offset = width * (rng.random() - .5)
            high = None
        offsetf = F.from_float(offset)
        actual_high = old.axis_count(y[i] + offsetf, F(3, 2), 1, oracle.A) == 2
        if high is not None:
            assert actual_high == high, 'sampler endpoint/status discrepancy'
        assert -wf / 2 <= offsetf <= wf / 2
        near_boundary += int(any(abs(float(offsetf - bound)) <= 1e-12 for a, b, _ in segments for bound in (a, b)))
        offsets[i] = offset
        free.append(i)
        pvals.append(p)
        hs.append(int(actual_high))
        logZ += math.log(Zi)
    logL = sum(hs) * math.log(T) - logZ
    assert max(abs(F.from_float(v)) for v in offsets) == wf / 2 if shell else max(abs(F.from_float(v)) for v in offsets) <= F(1, 2)
    RN = 1 / (.25 + .75 * math.exp(logL))
    assert 0 < RN <= 4
    face_high = None if face is None else int(old.axis_count(y[face] + F.from_float(offsets[face]), F(3, 2), 1, oracle.A) == 2)
    out = oracle.evaluate(codes, source_k, offsets)
    geometric = (width / out['R_float'])**oracle.n if shell else out['R_float']**(-oracle.n)
    Zscore = out['g'] * Cn * geometric * RN if out['E'] and out['source_captured'] else 0.0
    assert 0 <= Zscore <= 4 * Cn * (1 + 1e-12)
    out.update(core_or_shell='shell' if shell else 'core', r=width, r_hex=width.hex(),
               face=face, sign=sign, proposal_branch='tilted_product' if tilted_branch else 'base_product',
               free_coordinate_indices=free, free_p_i=pvals, free_high_status=hs,
               face_reference_high=face_high, log_normalizer=logZ, logL=logL, RN=RN,
               offsets_hex=[v.hex() for v in offsets], near_partition_boundary_count=near_boundary,
               score=Zscore, reference_high_count=sum(hs) + (face_high or 0))
    return out


def ess(values):
    a = np.asarray(values, dtype=float)
    return float(a.sum()**2 / np.square(a).sum()) if np.square(a).sum() else 0.0


def hoeffding(mean, B, count, union, delta, target_max):
    radius = B * math.sqrt(math.log(2 * union / delta) / (2 * count))
    return dict(mean=float(mean), radius=radius, lower=max(0.0, mean - radius),
                upper=min(target_max, mean + radius), bound_of_observation=B,
                target_max=target_max, count=count, union_count=union, delta=delta,
                scope='ideal exact continuous model; no floating-law interval certification')


def main():
    regpath = BASE / 'registration.json'
    reg = json.loads(regpath.read_text())
    assert sha(reg['oracle_module_path']) == reg['oracle_module_sha256']
    for inp in reg['source_inputs'].values():
        assert sha(inp['path']) == inp['sha256']
    spec = importlib.util.spec_from_file_location('read_only_original_arrival', reg['oracle_module_path'])
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)  # Pure definitions only, never old main/probes.
    started = time.perf_counter()
    rounds = []
    for plan in reg['rounds']:
        n, Ns, Nr = plan['n'], plan['source_samples'], plan['receivers_per_replica']
        meta = reg['source_inputs'][str(n)]
        inp = meta['input']
        oracle = IntegerArrival(inp, meta['tau'])
        T = float(F(plan['tilt_odds_T']))
        Cn = 1 + n * math.log(2)
        srng = np.random.default_rng(plan['seed_source'])
        rrng = [np.random.default_rng(plan[f'seed_replica_{j}']) for j in (1, 2)]
        codes_list, labels = [], []
        for _ in range(Ns):
            labelcoin = int(srng.integers(0, 6))
            label = 0 if labelcoin < 3 else labelcoin - 2
            k = inp['ks'][label]
            labels.append(label)
            codes_list.append([int(v) for v in srng.integers(-inp['A'] * k, inp['A'] * k + 1, size=n)])
        guards = sampler_guards(old, oracle, codes_list[0], inp['ks'][labels[0]], F(plan['tilt_odds_T']))
        guardpath = BASE / f'n{n}_sampler_guards.json.gz'
        with gzip.open(guardpath, 'xt', encoding='utf8') as f:
            json.dump(guards, f, separators=(',', ':'))
        profilepath = BASE / f'n{n}_source_profiles.jsonl.gz'
        S1, S2, allw, allscores, validation = [], [], [], [], []
        diagnostic = Counter()
        exact_ties = 0
        round_start = time.perf_counter()
        with gzip.open(profilepath, 'xt', encoding='utf8') as stream:
            for s, (codes, label) in enumerate(zip(codes_list, labels)):
                k = inp['ks'][label]
                samples = [[], []]
                for rep in (0, 1):
                    for draw in range(Nr):
                        out = receiver(old, oracle, codes, k, rrng[rep], T, Cn)
                        if (s, rep, draw) in {(0, 0, 0), (0, 1, 0), (1, 0, 0)}:
                            x = [F(c, k) + F.from_float(float.fromhex(v)) for c, v in zip(codes, out['offsets_hex'])]
                            cross = old.arrival_oracle(oracle.old_data, x)
                            Rex = F(out['R_key'], out['dyadic_denominator'])
                            Mex = F(out['full_mass_numerator'], out['full_mass_denominator']) / Rex**n
                            assert F(cross['winner_R']) == Rex and F(cross['M']) == Mex
                            validation.append(dict(source_index=s, replica=rep, draw=draw,
                                                   receiver_x=[str(v) for v in x], R=str(Rex), M=str(Mex),
                                                   original_oracle_candidate_count=cross['candidate_count'], passed=True))
                        samples[rep].append(out)
                        allw.append(out['RN'])
                        allscores.append(out['score'])
                        diagnostic['receivers'] += 1
                        diagnostic['E'] += out['E']
                        diagnostic['source_captured'] += out['source_captured']
                        diagnostic['E_and_source_capture'] += out['E'] and out['source_captured']
                        diagnostic['positive_score'] += out['score'] > 0
                        diagnostic['near_threshold'] += out['near_threshold']
                        diagnostic['near_max'] += out['near_max']
                        diagnostic['near_partition_boundary'] += out['near_partition_boundary_count']
                        diagnostic[out['core_or_shell']] += 1
                        diagnostic[out['proposal_branch']] += 1
                        exact_ties += out['tied_max_candidates'] > 1
                a, b = [sum(z['score'] for z in rep) / Nr for rep in samples]
                S1.append(a)
                S2.append(b)
                mean = (a + b) / 2
                conditional = hoeffding(mean, 4 * Cn, 2 * Nr, 3 * Ns, .025, Cn)
                memberships = [all((c * kk) % k == 0 for c in codes) for kk in inp['ks']]
                row = dict(source_index=s, source_component=label, source_k=k, source_integer_codes=codes,
                           source_y=[str(F(c, k)) for c in codes], source_grid_membership=memberships,
                           source_on_support_edge=any(abs(c) == inp['A'] * k for c in codes),
                           receiver_replicas=samples, S_replica_1=a, S_replica_2=b, S_estimate=mean,
                           conditional_S_ideal_hoeffding=conditional, I2_cross_score=a * b)
                stream.write(json.dumps(row, separators=(',', ':')) + '\n')
                if (s + 1) % 32 == 0:
                    stream.flush()
                    print(json.dumps(dict(n=n, completed_sources=s + 1, receiver_calls=diagnostic['receivers'],
                                          positive_scores=diagnostic['positive_score'], elapsed=time.perf_counter() - round_start)), flush=True)
        a, b = np.asarray(S1), np.asarray(S2)
        v1, v2 = (a + b) / 2, a * b
        bound1, bound2 = 4 * Cn, 16 * Cn * Cn
        summary = dict(n=n, source_samples=Ns, receivers_per_replica=Nr, tau=meta['tau'], Cn=Cn,
                       seeds={key: plan[key] for key in plan if key.startswith('seed')},
                       I1_over_W=hoeffding(float(v1.mean()), bound1, Ns, 3 * 2 * Ns, .025, Cn),
                       I2_over_W=hoeffding(float(v2.mean()), bound2, Ns, 3 * 2 * Ns, .025, Cn * Cn),
                       I2_over_I1_point_estimate=float(v2.mean() / v1.mean()) if v1.mean() else None,
                       replica_I1_means=[float(a.mean()), float(b.mean())],
                       source_S_estimate_min=float(v1.min()), source_S_estimate_max=float(v1.max()),
                       source_S_estimate_variance=float(v1.var(ddof=1)),
                       receiver_RN_ESS=ess(allw), receiver_score_ESS=ess(allscores),
                       outer_I1_score_ESS=ess(v1), outer_I2_cross_score_ESS=ess(v2),
                       RN_sample_mean=float(np.mean(allw)), diagnostics=dict(diagnostic),
                       exact_multiwinner_receivers=exact_ties,
                       source_label_counts={str(l): labels.count(l) for l in range(4)},
                       source_support_edge_count=sum(any(abs(c) == inp['A'] * inp['ks'][l] for c in codes) for codes, l in zip(codes_list, labels)),
                       plugin_source_tail=[dict(K=K, value=float(np.maximum(v1 - K, 0).mean()),
                                                biased_plugin_only=True) for K in (math.sqrt(n), 2 * math.sqrt(n))],
                       new_original_oracle_crosschecks=validation,
                       sampler_guards_path=str(guardpath), sampler_guards_sha256=sha(guardpath),
                       source_profiles_path=str(profilepath), source_profiles_sha256=sha(profilepath),
                       fixed_N_complete=True, ideal_model_unbiased_I2_cross=True,
                       floating_proposal_law_not_interval_certified=True,
                       elapsed_seconds=time.perf_counter() - round_start)
        roundpath = BASE / f'n{n}_results.json'
        with roundpath.open('x') as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
            f.write('\n')
        rounds.append(summary)
        print(json.dumps(dict(n=n, status='complete', I1=float(v1.mean()), I2=float(v2.mean()),
                              E=diagnostic['E'], positive_scores=diagnostic['positive_score'],
                              elapsed_seconds=summary['elapsed_seconds'])), flush=True)
    result = dict(status='complete', registration_sha256=sha(regpath), script_sha256=sha(Path(__file__)),
                  rounds=rounds, elapsed_seconds=time.perf_counter() - started,
                  old_oracle_main_or_old_probes_rerun=False, no_parameter_adaptation=True,
                  scope='New fixed-budget complete-source and continuous-winner occupation estimates, ideal-law statistical CI only. No asymptotic order or endpoint proof.')
    with (BASE / 'results.json').open('x') as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps(dict(status='complete', elapsed_seconds=result['elapsed_seconds'], results_sha256=sha(BASE / 'results.json'))), flush=True)


if __name__ == '__main__':
    main()
