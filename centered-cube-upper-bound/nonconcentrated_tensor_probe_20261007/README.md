# Nonconcentrated implicit tensor probe, 2026-10-07

Applied `math-l02-adversarial-input-search`, including provenance, method and cube-interface references. Read B41 from the source table identified in `preregistration.json`, and the older tensor pressure README, code and repair receipt. Existing outputs were neither rerun nor modified.

This experiment provides floating numerical pressure on a B41 shared-resolution mechanism. It neither proves a special-family bound nor validates actualgeom or a general cube weak(1,1) bound. It does not reproduce a smoothed lower-family theorem. No FIRST, continuation, CP, GP or other geom gates appear.

## Complete source, not an empirical approximation

For each production dimension n=16,32,64, the original input is the finite probability measure

    mu = sum_{j=0}^2 p_j [ (64 q_j)^(-1)
          sum_{k=0}^{64 q_j-1} delta_{k/q_j} ]^tensor_n,
    q=(1,2,4), p=(1/4,1/2,1/4), W=1.

Every component is a complete tensor grid. Resolution is selected once per entire n-dimensional source draw, followed by independent coordinates conditional on that shared resolution. This is a positive mixture after tensorization; independently selecting resolutions for coordinates would be a different measure. The three grids overlap, and all coincident atoms have summed component mass. The union has exactly 256^n atoms, finite support in [0,64)^n and total mass one. Source sampling retains boundary sources and samples the complete mu/W; only the outer source/cone expectation is Monte Carlo. Maximal averages use the complete implicit measure, with no sparse replacement.

This equals a dilation by 64 of the dyadic grids of resolutions 64,128,256, so it realizes B41's finite shared-resolution mixture mechanism. No smoothing, positive-width L1 transfer, entropy inequality, weak quotient, or asymptotic lower-bound reproduction is inferred.

## Receiver, winner and microbox certificate

Take omega uniformly from the 2n faces of [-1,1]^n, with independent uniform free coordinates. Set x=z+r omega/2, r=exp(t/n), t in [0,n log2]. The receiver cube has full side R in [1,2]. Compute

    m(x,R)=sum_j p_j product_i #(component-j grid in [x_i-R/2,x_i+R/2])/(64q_j)^n,
    u(x)=max_{1<=R<=2} m(x,R)/R^n.

The candidate set contains both endpoints and every component coordinate lattice arrival in (1,2]. Count products are constant between arrivals, hence dividing by R^n decreases there. All equal floating arrival groups are completed before scoring; coincident arrivals from different components update every component before mixture evaluation. This evaluates the full continuous window without an R grid. The smallest candidate wins exactly equal computed scores. General mathematical ties and near equal scores remain floating uncertainties.

All source coordinates use a common integer quarter-unit Z. Component lattice offsets are formed relative to that local source. The source-tagged face arrival is exactly r in floating arithmetic, avoiding the older addition/subtraction cancellation bug. `validation.json` compares 48 n<=4 cases with an independent full all-atom construction and maximal scan: maximum log-score error 1.33e-15, maximum winner error 2.22e-16. A 100,000-draw prior check verifies the shared-resolution sampling and records empirical frequencies. `microbox_validation.json` independently scans 40 small-dimensional complete source/cube interval-placement cases and checks the safe bound; it also checks exact local source-face arrivals at n=16,32,64. These checks are floating implementation evidence, not interval proofs.

The main residual threshold is eta_n=(49/65536)/sqrt(n), with fixed eta0=49/65536 retained only as a diagnostic. Microboxes have side h=2/n. In each component and coordinate, an arbitrary closed interval of length h captures at most floor(h q_j)+1 points, and captures zero if that component's coordinate interval in Q is empty. Thus

    sup_B mu(B intersect Q) <= sum_j p_j product_i cap_ji/(64q_j)^n.

For every production n>=16, h q_j<1, so cap_ji=min(c_ji,1), where c_ji is the component's actual captured coordinate count. Dividing this bound by the actual winner m gives a safe upper bound on the microbox capture share, conditional on the computed winner/counts. Define Eexists as the exact existence event that some side-h box captures more than eta_n times m; its complement is the exact nonconcentrated residual used here. Strictly passing eta_n with numerical slack identifies a screened certified subset of that exact complement, conditional on the computed winner and counts. No particular finite greedy-cover complement E\U is used or claimed. The finite cover may overpay some exact nonconcentrated queries; it is only a fee construction in the separate analytic audit. Failing this sufficient test means **unknown**; it does not prove concentration and is not forcibly removed from the upper profile. The per-component interval bound itself is analytic and exact; the winner and log-ratio computations are floating and not formally outward rounded.

## Frozen lambda and fields

Each n has one common threshold fixed before production, log(lambda)=-n log64. There is no lambda or geometric parameter optimization. Let d=n log(Rwinner/r), vstar=log(65536/49). The tested fields are

    Gamma_v1:       1_{2lambda<u<=4lambda} 1_{0<=d<=1} 1_res,
    Theta_hard:     1_{2lambda<u<=4lambda} exp(-d) 1_{d>=0} 1_res,
    Theta_hard_v1:  1_{2lambda<u<=4lambda} exp(-d) 1_{0<=d<=1} 1_res,
    Gamma_vstar:    1_{2lambda<u<=4lambda} 1_{0<=d<=vstar} 1_res.

The normalized cone Jacobian is dx=r^n dt dcone; there is no additional 2^-n. Profiles are expectations over the original source and cone. Reported energy is the time integral of the squared profile.

Each NPZ has t, profile_A, profile_B and fixed_eta0_A/B. Profile shape is (4 quantities,3 axes,K). Quantities follow the above order. The axes are `screened_lower_cert_res`, `raw_cert_res`, `floating_upper_including_unknown`. A numerical tolerance (n+1)*1e-10 screens hardband boundaries, delay endpoints and near maximizing candidate scores. Lower profiles minimize eligible values over near maximizing candidates; upper profiles maximize values and retain all unknown residual queries. The upper profile omits the microbox gate and is a full hard-relaxation upper envelope; it must not be read as a computed cooperative residual energy. No tolerance widens the raw hardband, whose lower endpoint is strict and upper endpoint closed. General score/band ambiguity is retained, without symbolic or interval certification. Fixed eta0 diagnostics are raw comparisons and do not use a separate lambda.

## Three completed rounds and limits

| round | samples per side | fine nodes | nested coarse nodes | independent A/B pairs |
|---|---:|---:|---:|---:|
| 1 | 32 | 65 | 33 | 4 |
| 2 | 64 | 129 | 65 | 4 |
| 3 | 128 | 257 | 129 | 4 |

All three dimensions run in each round. Seed namespaces for different rounds, dimensions, pairs and sides are distinct and saved. Within a side, the same sampled source/cone paths are retained across t. Each pair's energy is `trapezoid(profile_A*profile_B,t)`, an unbiased cross-product estimator of the corresponding finite-node profile-square statistic when the two sides are independent. Each displayed mean averages four independent pair energies; SD is the sample SD across those four values, and SE=SD/2. The stored approximate 95% interval uses t3=3.182446305. Four pairs give weak distributional evidence; this approximation need not provide valid coverage. The bounded Hoeffding intervals use the pair-energy range [0,n log2] and are 95% marginal probabilistic bounds for the finite-node statistic. They are generally wide and are not simultaneous across all fields/dimensions/rounds.

Neither kind of Monte Carlo interval controls floating error, continuous-time quadrature, unknown residual certification, or actualgeom gates. Nested coarse/fine differences use the same sampled paths and are resolution diagnostics. Discontinuous band/winner/short-shell fields prevent claiming a deterministic quadrature error from those differences. Numerical lower/upper profiles screen the stated floating comparisons, not all possible machine error. No finite energy disproves a bound with an unspecified constant or exponent; no n-order is fitted from three dimensions.

`summary.md` gives the final numerical table; `completion_receipt.json` gives hashes and query-screen fractions. Immutable `round*.json` files save seeds, all per-node branch screens, per-pair energies, coarse energies and Monte Carlo intervals. Full profiles are in 36 NPZ files. `preregistration.json` was written before production.

Runtime:

    /Users/zhengzhihao/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3

Run `probe.py --validate`, `microbox_check.py`, then `probe.py --round 1`, `--round 2`, `--round 3` separately in a fresh copied output directory. Existing round tags refuse overwrite. `summarize.py` reads completed rounds only.

## Independent dimension-128/512 controls

An additional registered control uses n128 with M32/K513 and n512 with M16/K1025, each with two independent A/B pairs. These use the same complete three-component source, frozen lambda, eta_n, vstar and unchanged engine. The original three rounds are not rerun. See `highdim_preregistration.json`, `highdim_summary.md` and `highdim_completion_receipt.json`. Δt is approximately .1733 and .3466, respectively; Γ1 and hardband aliasing remain possible, and coarse/fine differences are not quadrature error bounds. Two-pair t1 and Hoeffding intervals address finite-node MC targets only. Diagnostic exp(log winner mass) may underflow at n512; scoring, bands and microbox ratios stay logarithmic.

A documentation correction verified `validation.json` reports maximum winner error 2.220446049250313e-16, not zero. The original validation and production datasets were retained unchanged.
