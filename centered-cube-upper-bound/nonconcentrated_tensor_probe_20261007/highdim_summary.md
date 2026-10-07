# Independent high-dimensional controls

The original three rounds were not rerun. Source formula, all three overlapping tensor components, fixed lambda, eta_n and vstar are unchanged; controls use independent registered seeds. No dimension order is fitted.

| n | M / side | nodes | Δt | Γ1 screened cert | Θhard screened cert | Γvstar screened cert | Γvstar fullhard upper | cert fraction all queries |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 128 | 32 | 513 | 0.173287 | 0.0185302 | 0.0148095 | 0.0886741 | 0.112873 | 0.917717 |
| 512 | 16 | 1025 | 0.346574 | 0.00812282 | 0.00779491 | 0.089351 | 0.103566 | 0.975427 |

Screened cert is a subset passing the sufficient microbox count bound for the exact nonconcentrated residual, conditional on computed winner/counts. Failed certificates remain unknown. Fullhard upper omits the microbox gate and is not cooperative residual energy.

Each energy mean averages two independent A/B cross-product pair estimates. SE is the sample SD divided by sqrt2. Stored approximate t1 intervals and marginal bounded Hoeffding intervals address finite-node MC statistics only; two pairs give broad, fragile intervals. Γ1 and hardband discontinuities may alias the finite time grids, especially Δt=.34657 at n512. Stored coarse/fine differences diagnose sensitivity without giving a deterministic quadrature error bound. General floating comparisons are screened, without outward interval arithmetic.

At n512, diagnostic exponentiation of winner log mass can underflow to zero; winner scores, lambda comparisons and microbox share computations remain in log space. That mean-mass diagnostic is not a numerical zero-mass claim.

Full control JSON contains per-node screens, sources/seeds/hashes and all quantity/axis energies. Four NPZ files preserve full A/B profiles.
