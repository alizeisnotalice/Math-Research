#!/usr/bin/env python3
"""Bounded hard-band radial occupation diagnostic, g=1 only.

This script does not evaluate softFIRST or actual geom. It writes a new,
timestamped JSON and never overwrites earlier results. Requires NumPy.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, math, time
import numpy as np

SEED = 202610070941
CASES = [(1, 8, 8), (2, 32, 16), (3, 128, 32),
         ("N-control", 32, 8), ("N-control", 32, 32)]
FAMILIES = ["unequal_diagonal", "correlated_two_clusters", "random_sign_outlier"]
LAMBDA_FACTORS = [0.5, 1.0, 2.0]  # Fixed before observing any events.


def make_input(n, N, family, seed):
    rng = np.random.default_rng(seed)
    if family == "unequal_diagonal":
        atoms = np.repeat(np.linspace(-0.7, 0.7, N)[:, None], n, axis=1)
        weights = np.geomspace(1.0, 16.0, N)
    elif family == "correlated_two_clusters":
        centers = np.where(np.arange(N) < N // 2, -0.48, 0.48)
        shared = rng.normal(0, 0.065, (N, 1))
        atoms = centers[:, None] + shared + rng.normal(0, 0.025, (N, n))
        weights = np.where(np.arange(N) < N // 2, 1.0, 2.0)
    else:
        atoms = rng.choice([-0.42, 0.42], (N, n)) + rng.normal(0, 0.015, (N, n))
        atoms[-1] = 1.5
        weights = np.ones(N)
        weights[-1] = 0.2
    weights = weights / weights.sum()
    digest = hashlib.sha256(atoms.astype('<f8').tobytes() + weights.astype('<f8').tobytes()).hexdigest()
    return atoms, weights, digest


def maximal_from_q(q, weights, n, t):
    """q_j = full side needed for atom j / source radial side.

    Evaluate only the END of each exact equal-distance group. Choosing the
    smallest radius in an exact response tie is deterministic.
    """
    order = np.argsort(q, axis=1, kind="stable")
    qs = np.take_along_axis(q, order, axis=1)
    mass = np.cumsum(weights[order], axis=1)
    last = np.concatenate([qs[:, :-1] != qs[:, 1:], np.ones((len(q), 1), dtype=bool)], axis=1)
    with np.errstate(divide="ignore"):
        logq = np.log(qs)
        response = np.log(mass) - t - n * logq
    response = np.where(last, response, -np.inf)
    winner = np.argmax(response, axis=1)
    rows = np.arange(len(q))
    return qs, logq, mass, last, response, winner, response[rows, winner]


def scalar_direct(x, atoms, weights):
    ds = [float(np.max(np.abs(x - a))) for a in atoms]
    if 0.0 in ds:
        return float("inf"), 0.0
    options = [(sum(float(w) for d, w in zip(ds, weights) if d <= D)/(2*D)**len(x), 2*D)
               for D in sorted(set(ds))]
    return max(options, key=lambda p: (p[0], -p[1]))


def window_maximal(q, weights, n, t):
    """True continuous-window maximal function, full sides in [1,2]."""
    qs, logq, mass, last, _, _, _ = maximal_from_q(q, weights, n, t)
    r = math.exp(t/n)
    qa, qb = 1/r, 2/r
    inside = (qs >= qa) & (qs <= qb)
    endpoints = np.broadcast_to([qa, qb], (len(q), 2))
    endpoint_mass = np.stack([np.sum((q <= qa)*weights, axis=1),
                              np.sum((q <= qb)*weights, axis=1)], axis=1)
    qs = np.concatenate([qs, endpoints], axis=1)
    mass = np.concatenate([mass, endpoint_mass], axis=1)
    valid = np.concatenate([last & inside, np.ones((len(q), 2), dtype=bool)], axis=1)
    perm = np.argsort(qs, axis=1, kind='stable')
    qs = np.take_along_axis(qs, perm, axis=1)
    mass = np.take_along_axis(mass, perm, axis=1)
    valid = np.take_along_axis(valid, perm, axis=1)
    with np.errstate(divide='ignore'):
        logq = np.log(qs)
        resp = np.log(mass)-t-n*logq
    resp = np.where(valid, resp, -np.inf)
    ix = np.argmax(resp, axis=1)
    return qs, logq, mass, valid, resp, ix, resp[np.arange(len(q)), ix]


def validate():
    tests = [
        (np.array([0., 0.]), np.array([[0.5, 0.], [-0.5, 0.], [0., 1.]]), np.array([0.2, 0.3, 0.5])),
        (np.array([0.2, -0.1]), np.array([[0., 0.], [0.4, 0.2], [0.9, -0.1]]), np.array([0.2, 0.7, 0.1])),
        (np.array([0., 0.]), np.array([[0., 0.], [0.5, 0.]]), np.array([0.4, 0.6])),
    ]
    result = []
    for x, atoms, weights in tests:
        q = 2*np.max(np.abs(x[None, :] - atoms), axis=1)[None, :]
        qs, _, mass, last, response, winner, logu = maximal_from_q(q, weights, len(x), 0.)
        u, R = scalar_direct(x, atoms, weights)
        idx = winner[0]
        R2 = float(qs[0, idx])
        u2 = math.exp(float(logu[0])) if math.isfinite(logu[0]) else float("inf")
        ok = (u == u2 or math.isclose(u, u2, rel_tol=1e-12)) and math.isclose(R, R2, abs_tol=1e-12)
        assert ok, (u, R, u2, R2)
        result.append(dict(pass_=ok, scalar_u="infinity" if math.isinf(u) else u,
                           winning_R=R2, winning_mass=float(mass[0, idx]),
                           tie_group_candidate_count=int(last.sum())))
    # Exact response tie: mass 1/4 at side 1, mass 1 at side 2 in n=2.
    _, _, _, _, response, winner, _ = maximal_from_q(np.array([[1., 2.]]), np.array([0.25, 0.75]), 2, 0.)
    assert winner[0] == 0 and response[0, 0] == response[0, 1]
    result.append(dict(pass_=True, validation="exact response tie chooses smaller full side"))
    for x, atoms, weights in tests:
        q = 2*np.max(np.abs(x[None, :]-atoms), axis=1)[None, :]
        qs, _, _, _, _, ix, logu = window_maximal(q, weights, len(x), 0.)
        sides = sorted(set([1., 2.]+[float(v) for v in q[0] if 1. <= v <= 2.]))
        choices = [(sum(float(w) for d,w in zip(q[0], weights) if d <= side)/side**len(x),side)
                   for side in sides]
        u,R = max(choices, key=lambda p:(p[0],-p[1]))
        assert math.isclose(math.exp(logu[0]),u,rel_tol=1e-12) and qs[0,ix[0]] == R
    result.append(dict(pass_=True, validation='continuous window [1,2] candidate enumeration matches scalar endpoint-inclusive implementation'))
    # Delta_0 is a normalization unit test, not a new research family.
    rng = np.random.default_rng(91)
    n, Nsamples = 8, 12
    omega = rng.uniform(-0.5, 0.5, (Nsamples, n))
    face = rng.integers(0, 2*n, Nsamples)
    omega[np.arange(Nsamples), face//2] = np.where(face % 2, 0.5, -0.5)
    single_nodes = []
    for t in [0., math.log(2)/2, math.log(2), n*math.log(2)]:
        r = math.exp(t/n)
        q = 2*np.max(np.abs(omega), axis=1)[:, None]
        qs, _, _, _, _, ix, logu = maximal_from_q(q, np.array([1.]), n, t)
        assert np.all(qs == 1.) and np.all(logu == -t)
        u = math.exp(-t)
        # lambda=1/4 gives a lower closed t endpoint at 0 and upper open at log2.
        observed = int(0.5 < u <= 1.)
        assert observed == int(0. <= t < math.log(2))
        single_nodes.append(dict(t=t, R=r, u=u, Y=observed, all_samples_identical=True))
    lam = 0.25
    B = n*math.log(2)
    left, right = max(0., -math.log(4*lam)), min(B, -math.log(2*lam))
    length = max(0., right-left)
    assert length <= math.log(2) and length == math.log(2)
    assert math.isclose(math.expm1(B), 2.**n-1., rel_tol=1e-12)
    result.append(dict(pass_=True, validation="single atom normalization only",
                       n=n, lambda_=lam, t_interval_closed_open=[left, right],
                       analytic_L1=length, analytic_L2_energy=length, samples=single_nodes,
                       jacobian="dx=r^n dt dnu; normalized boundary [-1/2,1/2]^n cone measure, cube volume r^n; shell integral exp(B)-1=2^n-1 verified"))
    return result


def run_case(round_, n, N, family, caseid, M, nodes, sampling_seed_base):
    seed = SEED + caseid*10007
    atoms, weights, digest = make_input(n, N, family, seed)
    B = n*math.log(2)
    ts = np.linspace(0., B, nodes)
    logs = [-n*math.log(2)/2-math.log(N)+math.log(f) for f in LAMBDA_FACTORS]
    # Independent replicates; within each replicate the same source and cone
    # direction are reused at every node and for every lambda (CRN).
    lowers = np.zeros((2, 3, nodes))
    uppers = np.zeros_like(lowers)
    raw = np.zeros_like(lowers)
    fullscale_control_phi = np.zeros_like(lowers)
    ambiguity = np.zeros((3, nodes), dtype=int)
    winner_source = np.zeros(nodes)
    winner_mass = np.zeros(nodes)
    winner_count = np.zeros(nodes)
    exact_distance_tie_samples = 0
    close_winner_samples = 0
    for rep in range(2):
        rng = np.random.default_rng(sampling_seed_base+caseid*10007+rep+1000003)
        source = rng.choice(N, M, p=weights)
        omega = rng.uniform(-0.5, 0.5, (M, n))
        face = rng.integers(0, 2*n, M)
        omega[np.arange(M), face//2] = np.where(face % 2, 0.5, -0.5)
        differences = atoms[source, None, :] - atoms[None, :, :]
        for k, t in enumerate(ts):
            r = math.exp(t/n)
            # This algebraically identical calculation preserves the exact
            # source q=1 and avoids cancellation from forming x then x-z.
            q = 2*np.max(np.abs(differences/r + omega[:, None, :]), axis=2)
            assert np.all(q[np.arange(M), source] == 1.)
            fs_q, fs_logq, _, _, _, fs_ix, fs_logu = maximal_from_q(q, weights, n, t)
            qs, logq, mass, last, resp, ix, logu = window_maximal(q, weights, n, t)
            rows = np.arange(M)
            tol = 1e-10*(1+n)  # Numerical screening, not a rigorous interval certificate.
            near = last & (resp >= logu[:, None]-tol)
            close_winner_samples += int(np.count_nonzero(near.sum(axis=1)>1))
            exact_distance_tie_samples += int(np.count_nonzero(np.any(fs_q[:, :-1] == fs_q[:, 1:], axis=1)))
            Q = qs[rows, ix]
            L = t/n+logq
            # Window membership was imposed in candidate construction, so an
            # endpoint winner is legitimate even at an exact radial boundary.
            radial_lower = np.ones_like(last)
            radial_upper = np.ones_like(last)
            source_lower = qs >= 1.
            source_upper = qs >= 1.-tol
            value = np.exp(-n*np.maximum(logq, 0.))
            winner_source[k] += float(np.mean(Q >= 1.))/2
            winner_mass[k] += float(np.mean(mass[rows, ix]))/2
            counts = np.sum(q <= Q[:, None], axis=1)
            winner_count[k] += float(np.mean(counts))/2
            for li, loglambda in enumerate(logs):
                low, high = loglambda+math.log(2), loglambda+math.log(4)
                band_lower = (logu-tol > low) & (logu+tol <= high)
                band_upper = (logu+tol > low) & (logu-tol <= high)
                yl = value*(radial_lower & source_lower & band_lower[:, None])
                yu = value*(radial_upper & source_upper & band_upper[:, None])
                lower = np.min(np.where(near, yl, np.inf), axis=1)
                upper = np.max(np.where(near, yu, -np.inf), axis=1)
                lowers[rep, li, k] = np.mean(lower)
                uppers[rep, li, k] = np.mean(upper)
                ambiguity[li, k] += int(np.count_nonzero(upper-lower > 1e-14))
                raw[rep, li, k] = np.mean((logu > low) & (logu <= high) &
                                           (Q >= 1.))
                fsQ = fs_q[rows, fs_ix]
                fsL = t/n+fs_logq[rows, fs_ix]
                fs_event = (fs_logu > low) & (fs_logu <= high) & (fsL >= 0.) & (fsL <= math.log(2)) & (fsQ >= 1.)
                fullscale_control_phi[rep, li, k] = np.mean(fs_event*np.exp(-n*np.maximum(fs_logq[rows,fs_ix],0.)))
    records = []
    # Union bound over all 45 lambda-case records and all finite fine nodes.
    # Two-sided Hoeffding for EACH of lower/upper screened bounded RVs:
    # 2 tails * 2 RVs * all records * all nodes.
    epsilon = math.sqrt(math.log(4*len(CASES)*len(FAMILIES)*3*nodes/0.05)/(2*(2*M)))
    for li, (factor, loglambda) in enumerate(zip(LAMBDA_FACTORS, logs)):
        phi = lowers[:, li].mean(axis=0)
        phi_upper = uppers[:, li].mean(axis=0)
        cross = lowers[0, li]*lowers[1, li]
        fine = float(np.trapezoid(cross, ts))
        coarse = float(np.trapezoid(cross[::2], ts[::2]))
        sq = float(np.trapezoid(phi**2, ts))
        confidence_lo = np.maximum(0., phi-epsilon)
        confidence_hi = np.minimum(1., phi_upper+epsilon)
        records.append(dict(round=round_, n=n, N=N, family=family,
            seed=seed, sampling_seed_base=sampling_seed_base,
            input_sha256=digest, lambda_factor=factor, log_lambda=loglambda,
            energy_cross_replicates_fine=fine, energy_cross_replicates_coarse=coarse,
            coarse_fine_abs_difference=abs(fine-coarse),
            energy_empirical_square_fine=sq,
            fullscale_max_then_filter_control_energy_cross=float(np.trapezoid(fullscale_control_phi[0,li]*fullscale_control_phi[1,li],ts)),
            energy_upper_numerical_envelope_cross=float(np.trapezoid(uppers[0, li]*uppers[1, li], ts)),
            energy_finite_node_hoeffding_lower=float(np.trapezoid(confidence_lo**2, ts)),
            energy_finite_node_hoeffding_upper=float(np.trapezoid(confidence_hi**2, ts)),
            max_phi=float(phi.max()), max_raw_event_fraction=float(raw[:, li].mean(axis=0).max()),
            raw_event_observations=int(round(raw[:, li].sum()*M)),
            ambiguous_sample_node_observations=int(ambiguity[li].sum()),
            exact_distance_tie_sample_node_observations=exact_distance_tie_samples,
            close_winner_sample_node_observations=close_winner_samples,
            finite_node_band_epsilon=epsilon,
            t=ts.tolist(), phi_rep_A=lowers[0, li].tolist(), phi_rep_B=lowers[1, li].tolist(),
            phi_numerical_upper=phi_upper.tolist(), raw_event_fraction=raw[:, li].mean(axis=0).tolist(),
            fullscale_max_then_filter_control_phi=fullscale_control_phi[:,li].mean(axis=0).tolist(),
            winner_source_membership_fraction=winner_source.tolist(),
            winner_mean_mass=winner_mass.tolist(), winner_mean_atom_count=winner_count.tolist(),
            atoms=atoms.tolist(), weights=weights.tolist()))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples-per-replicate', type=int, default=None,
                        help='Override dimension-dependent bounded defaults: n8=256,n32=128,n128=64')
    parser.add_argument('--nodes', type=int, default=129)
    parser.add_argument('--sampling-seed-base', type=int, default=SEED,
                        help='Input atoms remain fixed; use a distinct sampling base to change Monte Carlo streams')
    parser.add_argument('--validate-only', action='store_true')
    args = parser.parse_args()
    assert args.nodes >= 3 and args.nodes % 2 == 1
    start = time.perf_counter()
    validations = validate()
    if args.validate_only:
        print(json.dumps(validations, indent=2, allow_nan=False))
        return
    records = []
    for i, (round_, n, N) in enumerate(CASES):
        for j, family in enumerate(FAMILIES):
            M = args.samples_per_replicate or {8:256, 32:128, 128:64}[n]
            current = run_case(round_, n, N, family, i*len(FAMILIES)+j, M, args.nodes,
                               args.sampling_seed_base)
            for record in current:
                record['samples_per_replicate'] = M
            records.extend(current)
            print(f'round={round_} n={n} N={N} family={family} energies='+
                  ','.join(f'{r["energy_cross_replicates_fine"]:.6g}' for r in current), flush=True)
    payload = dict(schema='radial-hard-band-pressure-v1', created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        seed=SEED, sampling_seed_base=args.sampling_seed_base,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime_seconds=time.perf_counter()-start,
        samples_per_replicate_by_n={str(n):(args.samples_per_replicate or {8:256, 32:128, 128:64}[n]) for n in [8,32,128]},
        independent_replicates=2, fine_nodes=args.nodes, coarse_nodes=(args.nodes+1)//2,
        a=1, b=2, g=1, source_total_mass=1, lambda_selection='Fixed factors 0.5,1,2 times (1/N)*2^(-n/2); no post-observation adaptation',
        quantity='PRIMARY continuous-window u=max_{1<=R<=2}mu(Q(x,R))/R^n; Y=1{2lambda<u<=4lambda}*(r/R)^n*1{r<=R}; phi=E_source,cone[Y], energy=integral phi^2 dt. Fullscale-max-then-filter is a separately labeled control.',
        full_maximal_response='Continuous window candidates: sides a=1, b=2, and every atom arrival side in [1,2]; exact floating equal-distance groups include all tied atoms; maximize cumulative mass / full_side^n; no radius grid. Sorted candidate enumeration is exact for computed floating distances, not certified exact real arithmetic.',
        sampling='Source label z sampled with probability its mass; uniformly choose one of 2n cone faces and uniform free coordinates; same draws across t and lambda',
        jacobian='For each fixed source, normalized cone measure nu on boundary [-1/2,1/2]^n: cube volume=r^n, dx=n*r^(n-1) dr dnu=r^n dt dnu. These source-weighted cone expectations count labels; they are not unweighted spatial/Lebesgue integrals. Overlapping source parameterizations retain distinct labels.',
        confidence='Cross-replicate phi_A*phi_B is unbiased for the squared screened expectation at each finite node, despite CRN correlation across nodes. Empirical mean square is biased upward. Hoeffding bounds use iid source-cone draws, Y in [0,1], global 95% union bound across all 45 records and finite fine nodes, with two-sided bounds on BOTH lower and upper screened random variables: epsilon=sqrt(log(4*45*K/.05)/(2*2M)). Their trapezoidal energy bounds only bracket the weighted finite-node target, not the continuum integral.',
        numerical_ambiguity='Tolerance (n+1)*1e-10 screens near band/source-membership boundaries and near-equal winning responses; lower=min eligible values over near winners, upper=max possible values. Exact source q=1 is retained. Exact a,b endpoints are explicit window candidates. This conservative numerical screening is not directed-rounding interval arithmetic.',
        limitations=['g=1 relaxation only; neither original softFIRST nor actual geom',
                      'Nested trapezoid difference is a diagnostic, not a certified quadrature error',
                      'Continuous-window pressure does not certify arbitrary finite radius sets J',
                      'Finite small N and fixed dimensions only; no asymptotic-order or endpoint-proof claim',
                      'Hoeffding confidence is conditional on correct numerical screening; rounding errors are not rigorously enclosed'],
        validation=validations, records=records)
    outdir = Path(__file__).resolve().parent
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    outfile = outdir/f'radial_hard_band_results_{stamp}.json'
    with outfile.open('x') as f:
        json.dump(payload, f, ensure_ascii=False, indent=2, allow_nan=False)
    print(f'OUTPUT={outfile}\nRUNTIME={payload["runtime_seconds"]:.3f}', flush=True)


if __name__ == '__main__':
    main()
