# C048 pre-native acceptance report

The copied-node second derivative formula and exact physical chain are derived in theorem.md. New code retains all zeros and repeated nodes, extends the real analytic band jets to order7, and saves21 independent entries of the symmetric6×6 Hessian. It uses DD7 for repeated base nodes and an explicitly derived equivalent barycentric formula for distinct base nodes.

Registration input SHA field:
`316b359ae0f8ddb5cafe09a46a5d6edb25af18edd16813155fe864a978b8cd4d`.
The frozen C016 manifest/source/full-result, fan, C010 header archive and runtime dylibs were checked before writing the new source. ABI is64-bit ulong/fmpz,48-byte arb_struct, FLINT3.6.0, as checked in the native program. New source/build/binary hashes and full compile command are in build_record.json. No frozen file was changed, no environment installed, and no system path changed.

Compilation succeeded once with exit0; no compile fix. The observed source-generation-plus-build tool wall was0.738015209 seconds (pure compiler time was not separately instrumented). The binary is private scratch and will not be archived.

Exactly16 scalar Hessian quantities passed an independent algorithm: a degree7 Hermite interpolation matrix is formed at the augmented eight labelled nodes, and exact Fraction Gaussian elimination obtains the coefficient-of-degree7 functional. Normalized atan derivatives are calculated from exact complex rational powers, while only atan values use Arb1024. This differs from both native triangular DD and its fast barycentric formula. No finite differences are used. Configurations cover all-zero, mixed zero, four repeated zeros, distinct mixed signs, six identical nonzero nodes, and two nonzero repeated clusters. Both diagonal factor2 and mixed-label factor1 are exercised. Complete reference dyadics and exact interpolation weights are retained in scalar_validation.json.

**Full native Hessian executions so far:0.** The only permitted next run is the exact C016 oldc point, denominator1000000007, a=1/250,T=6/125 and the same fan. Its existing band, six-gradient and endpoint densities must match saved C016 complete dyadics. All967680 logical terms must be retained. A full Hessian independent reference is explicitly out of scope; passing this next run cannot be reported as full independent Hessian certification.

No new spatial point, global grid, grad(p), second-order Euler, MATLAB, subagent, stage/commit/push or Cn bound is included.
