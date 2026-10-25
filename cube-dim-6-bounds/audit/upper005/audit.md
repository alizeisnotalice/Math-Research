# upper005 Cauchy-product face heat-kernel audit

## Verdict

Accepted as a proved improvement of the body32 fixed-face heat-kernel budget:

    p_t(1,y) <= [1/(243*pi^6)] * t^-4/(1+1024*t^2)
              <= [1/(243*pi^6)] * t^-4.

Thus the short-time positive-face L(6/5) coefficient is

    Knew = (6/7)*(1/243)^(1/6)*pi^(-11/6)
         = 2*3^(1/6)/(7*pi^(11/6)).

This is a local analytic estimate only. It supplies neither a completed
face integral certificate nor a C6 bound.

## Product inequality

For N=32 and real z_i with sum z_i=s, dividing the asserted Cauchy-product
inequality by the positive product reduces it to

    N*(N*t^2 + sum z_i^2) >= s^2 + N^2*t^2.

The difference is exactly

    N*sum z_i^2 - (sum z_i)^2
    = sum_{i<j}(z_i-z_j)^2.

So the inequality is correct for every real fiber point; no exceptional
zero/repeated configuration is involved.

On the x_1=1 fiber, all body directions have first component one, so
sum z_i=1. The factor is therefore

    32*t / (pi*(1+1024*t^2)).

## Deleting one factor and the quotient density

Using the delta/coarea form of the convolution makes the Jacobian explicit.
After deleting factor i, integrate z_i against the first-coordinate delta.
Its coefficient is v_i,1=1, hence there is no Jacobian multiplier. The
remaining five equations become

    y-v_i' = sum_{j != i} (v_j'-v_i') z_j.

Consequently the resulting function is precisely the five-dimensional
quotient density q_i,t at y-v_i'. It is finite: its 31 quotient directions
span R^5. The quotient symbol is bounded below by 12 times the l1 norm.
Hence Fourier inversion gives

    ||q_i,t||_infty
    <= (2*pi)^-5 * product_{k=1}^5 integral exp(-12*t*|u_k|) du_k
    = (12*pi*t)^-5.

The exact check found quotient rank five and the equal-coordinate Rademacher
absolute sum 60, consistent with the sharp l1 constant 12. The standard
convex symmetric argument gives the all-vector inequality from that
equal-coordinate value.

Summing the 32 deleted factors gives

    32*t/[pi*(1+1024*t^2)] * 32*(12*pi*t)^-5
    = [32^2/12^5] * t^-4/[pi^6*(1+1024*t^2)]
    = [1/243] * t^-4/[pi^6*(1+1024*t^2)].

Thus the stated new Cf is exactly 1/(243*pi^6), a factor 32 below the
previous 32/(243*pi^6).

The delta derivation initially holds almost everywhere. Both p_t and every
q_i,t are continuous because their Fourier multipliers are integrable under
the displayed l1 coercivity. The inequality therefore upgrades to every
fixed face point y, not merely an a.e. slice statement.

## Short-time norm budget

The x_1 marginal remains Cauchy of scale 32*t. At x_1=1 its total y-mass is

    32*t/[pi*(1+1024*t^2)].

Reflection symmetry in the five free coordinates gives one thirty-second of
this mass on the positive orthant, so the positive unit-face L1 bound is at
most t/pi. Interpolation with the new Linfinity bound, at q=6/5, gives

    ||p_t(1,.)||_Lq([0,1]^5)
    <= 243^(-1/6)*pi^(-11/6)*t^(1/6).

Integrating from zero to epsilon proves the displayed Knew times
epsilon^(7/6). The retained denominator (1+1024*t^2)^-1 could be useful in
a later full-time budget, but has not been used here to claim a new cost.

## Actual algebra check and scope

The independent SymPy check ran successfully:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python
      /Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/upper005/check_cauchy_face_budget.py

It verified the sum-of-pair-squares identity, quotient rank, the value 60,
and 32^2/12^5=1/243.

This audit does not provide vertex intervals, a continuous face box cover,
an upper bound for the full face integral I, a bound on g, a cost comparison,
or a C6 conclusion.
