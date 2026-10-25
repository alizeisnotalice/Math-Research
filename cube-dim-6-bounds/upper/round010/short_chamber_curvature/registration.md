# A010 registration: analytic obstruction to short-potential chamber convexity

Registered before the exact structural verification. This work owns only
`upper/round010/short_chamber_curvature/` in the curvature-astra clone. It
does not alter another agent's files, stage, commit, push, invoke MATLAB,
or launch subagents.

The target is the body32 Cauchy short potential K_T at T=4/125 on the
strict ordered positive face. The proposed proof is stronger: for every
T>0, a negative second directional derivative occurs in this chamber.
Its mechanism is a rank-two singularity, regularity away from S3, and
permutation symmetry. It does not rely on subtracting infinite G and H.

The sole registered computational batch is exact finite algebra:

- the 32 sign-vector columns, their Gram identity, and one six-column basis;
- c=(1,1/2,1/2,1/2,1/2,1/2), z=(0,51,22,2,-24,-51)/200;
- the exact path inequalities, four positive partial sums, and the
  complementary-column bounds in a positive density lower estimate;
- the symbolic pairing of the eight sign patterns of a three-column span
  into four opposite pairs, and the scalar Cauchy derivative identities.

Before computation, `verify_structure.py` writes an input manifest with
exact rational inputs and SHA-256 hashes of itself, the registration,
the theorem, and read-only prior sources. Its wall-clock limit is 1200
seconds. It stops on a failed assertion or timeout; otherwise it stops
after the one structural verification. There is no spatial search.

New H point passes: **0**. New G point passes: **0**. New density point
passes: **0**. No Monte Carlo, grid, numerical minimization, or derivative
sampling is authorized or needed. The output is an analytic existence
proof of negative curvature, not a numerically located Jensen triple.

The theorem gives no global lower bound for K_T and no improvement of C6.
The established upper bound remains 64 unless another independently
proved global argument changes it.
