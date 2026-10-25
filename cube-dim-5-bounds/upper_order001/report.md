# Maximum-order conditioning does not improve 11110

The proposed maximum-order partition is valid and has no entropy penalty.
Its 16 conditional contributions must be summed, and the local code does
so. Nevertheless, the complete resulting lower-bound formula at 11110 is
strictly between 4.20630244543289e−5 and 4.31012868806624e−5. This lies below
52e−6 and below the target 1e−4. It is an obstruction for this particular
formula at this point, not an upper bound on actual G.

`proof.md` supplies the conditional Jensen/least-squares derivation,
normalization, exact projection multiplicities, and rigorous stable
truncated-moment series. The first-moment series is proved by subtracting
the leading Laplace term before Hankel inversion, so no integration
constant is omitted. Inverse moments use the exact full inverse moment
minus an absolutely integrated tail series.

The registered finite computation uses 1024 geometric bins beginning at
z=1. Exact Stieltjes weights are differences of F^16 divided by16 per
winner; all 16 winner events are then summed. The bin kernels have interval
upper and lower bounds, based on monotonicity of the truncated mean and
harmonic mean. The omitted low interval is at most 1.17851e−7 and the
entire infinite tail at most 3.59695e−15. Therefore tail recovery and
quadrature refinement cannot raise this formula to the comparison value.

The remaining loss is in replacing the 15 residual variables by their
separate first and inverse moments. In the far maximum tail the first
moment grows like z^1/4 while the harmonic mean stays bounded; since all
projection residuals q are positive at11110, the resulting kernel decays
like z⁻¹. This is an exact description of this layer, not a claim about
uncomputed top-two or top-three constructions.

The actual local Arb256 run exited0 in0.35 seconds, using the existing
Python/Flint venv. No numeric stable CDF, spatial grid, MATLAB invocation,
or additional lower-bound computation was used. There are no running
handles and nothing was committed or pushed.

Reproduce from the independent clone:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_order001/order_certificate.py
```
