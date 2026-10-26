# A057 report: positive-mixture log-curvature

Status: new analytic lemmas and an executable conditional one-box certificate interface completed; seven fixed symbolic/matrix/interface configurations passed. Parent mathematical review remains the integration gate. There is no new kernel observation, global g certificate, or C6 upper bound.

## Main result

For the original unit-weight body32 model, F=H_a or B_(a,T), a>0, and any legal six-column basis N,

    -2 ||N^-1 u||^2/a^2 <= D_u^2 log F <= (25/4)||N^-1 u||^2/a^2.

The proof independently pays positivity, finite mass, normalization, first/second derivative domination, and the full infinite tail. It is valid at all space points, including S3, zero/repeated projections, ties and zeros. The explicit determinant-128 basis has face quadratic Q6=I5/4-J5/32.

The same proof with a fixed full-frame right inverse gives the stronger result

    D_u^2 log F >= -||u||^2/(16a^2),
    log F(c+u) >= log F(c)+grad log F(c).u-||u||^2/(32a^2).

The fixed choice R=A^T/32 is Loewner-optimal among fixed right inverses of this body32 A. It is a single matrix for every direction and uses the same center jet as the six-basis result. The derivative argument uses an ordinary translation of the remaining 26 variables to a fixed integration space; no conditional density on a zero-probability fiber is used.

The exact 32-corner whole-box rule preserves gradient directions. With interval jets it also preserves the common F(c) denominator: at each corner d, form b=sum min(g_i^-d_i,g_i^+d_i), f=clamp(b,[L_F,U_F]), then lower-bound f exp(b/f-d^TQd/a^2). One center value and five matching face derivatives suffice. No 21-field Hessian evaluation or separate positive J2 bound is required.

## Actual validation

The registered run used the existing local Python 3.15.0, SymPy 1.14.0 and python-flint 0.9.0 at 384 Arb bits. It took approximately 0.067 seconds for the seven configurations. `results.json` records all exact algebraic results, outward exponential endpoints, versions, and source hashes. Arithmetic verifies the symbolic identities and finite implementation behavior; the general theorem is proved in `theorem.md`, not inferred from samples.

1. Three scalar rational gap identities certify the Cauchy score and both log-curvature constants.
2. One legal six-body basis: det N=128, its inverse Gram and face eigenvalues 3/32, 1/4.
3. Its exact general face quadratic and unit five-box corner range [15/32,39/32].
4. One anisotropic box has exact joint tangent exponent -9891209/65894400; separate gradient/quadratic payment loses 11/8960.
5. The shared-denominator minimizer's three branches and a synthetic passing one-box interface run.
6. A synthetic unresolved run; nonpositive a and wrong-model rejection.
7. AA^T=32I, R=A^T/32, Q32=I5/32, and the exact equal-halfwidth curvature-cost ratio Q6/Q32=39/5.

All production G/H/p/gradient/Hessian calls: **0**. Full-fan calls: **0**. MATLAB calls: **0**. The two interface jets are labelled synthetic and cannot be cited as body32 values. The batch stopped at its registered seven configurations; no coverage pilot or basis search followed.

## Files and replay

The necessary files are `registration.md`, `input_manifest.json`, `theorem.md`, `box_certificate.py`, `check_fixed.py`, `results.json`, and this report. No files outside this directory were edited. Nothing was staged, committed or pushed.

From this directory, the original run command was:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_fixed.py
```

For replay, use a clean copy without `results.json`; the script intentionally refuses to overwrite an accepted batch. Its read-only absolute source dependencies and hashes are recorded in `input_manifest.json`.

For a future genuine jet:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B box_certificate.py INPUT.json
```

The input fields are:

- `model`: `body32_unit_cauchy_v1`, and `model_definition`: the exact constant in `box_certificate.py`.
- `window`: `{kind:"H",a:"rational"}` or `{kind:"B",a:"rational",T:"rational"}`.
- `center`: six rational strings with first coordinate 1; `face_box`: five pairs of rational endpoint strings, containing the center.
- `quadratic`: `frame32` (default, stronger) or `basis6`; `target`: a positive rational string.
- `jet`: matching model, model_definition, window and center; `F_interval`; five `face_gradient_intervals`; and a nonempty `source_proof` identifying the external verified jet source.

The result `PASS_CONDITIONAL_ON_MATCHING_CENTER_JET` means the tangent bound passes for that whole box if the supplied jet proof is valid. It does not authenticate that proof, establish any global coverage, or imply a point estimate is a certified interval. Otherwise the checker returns `UNRESOLVED`; invalid metadata are rejected.

## Remaining work

The parent's proof review should focus on the fixed-fiber representation (theorem §2), its local domination on the infinite tail, and the shared-denominator whole-box minimization (§4). Those are the new interfaces worth independent review.

If accepted, a later finite pilot needs fresh matching cutoff/band values and five gradients at fixed rational centers in low-value regions, S3 away from the paid S2 tube, and boundary/seam boxes. Use Q32 first; it has the same jet cost as Q6 and uniformly smaller curvature payment. A cutoff ladder must be fixed before kernel execution, and changing a or T requires matching jets. No such numerical pilot is authorized or executed in this subtask.

Successful local boxes still must be assembled into an exact whole-root cover with the A040 tube payment, with all unresolved leaves and seams accounted for. Full G>=F then permits the genuine full-G symmetry/homogeneity transfer. The target g>=1/65000 remains unresolved.
