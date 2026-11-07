# r060 independent weighted-mesh audit

## Verdict: asymptotic construction is sound after one required hypothesis repair

The advertised family `R=ceil((log N)^2)/N`, `M=floor(N^(1/4))`, for sufficiently large `N`, has the stated crossing scale, strict pointwise gap, and normalized actual modular budget `1+o(1)`. It is correctly only a weighted pointwise/composed benchmark construction, not an endpoint lower bound.

### Required repair to the stated finite parameter range

The initial condition `0<R<1/(2(M-1))` proves only `p1>0`; it does **not** imply the subsequently used `p1>=alpha=2R`. The latter is needed for:

* the first recovery height to be at most `1+1/N`;
* the displayed uniform downshift bound; and
* the nonnegative first low cell needed by the budget argument.

For example, `M=3, R=.24` obeys the stated `R<1/4`, but has `p1=.04<R`, so its first low value is negative. Replace the finite hypothesis by `0<R<=1/(2M)` (or explicitly impose `p1>=2R`). This holds with room to spare for the advertised asymptotic family, so no asymptotic claim changes.

### Weighted rank and global rectangle reduction

With that repair, `D_i=R r_i/p_i` gives at prefix `r_j` the exact table `a_ij=1-R/p_i` for `i=j` and one otherwise. Weighting by `p_i` therefore loses exactly `R`; the within-cell weighted loss rate is also exactly `R`.

For a proper outer suffix starting at slab `l>=2`, an inner prefix before cell `l` makes every retained fibre equal one. Thereafter its aggregate loss is nonnegative. A partial outer slab produces a fractional-linear average for each fixed inner prefix, hence has an endpoint maximizer. The only competing branch that retains a deficient fibre is maximized by the full outer support. Thus the exhaustive common-action value is indeed `max{(1-R)/(1+eta), T/(T+eta)}`, with `T=1-p1`.

Solving equality gives `eta*=RT/(p1-R)=2(M-1)R^2/(1-(2M-1)R)`, so `eta*=O(MR^2)`. At `b=(1+eta*)(1+N^-2)`, the strict inequalities follow from `R>N^-2/(1+N^-2)`. The adaptive quantity `1/(1+eta)` is valid only for the stated adaptive/composed benchmark, where each outer fibre may select its own inner prefix of value one; it must not be renamed as the strong maximal value.

### Actual budget

The recovery mass calculation is correct:

`H=p1(Q-1)/Q^(M-1)+alpha(1-Q/Q^(M-1)) <= alpha+Q^-1` for `M>=3`.

After discarding the nonpositive low-cell contribution, this yields the displayed upper bound (even a slightly sharper version subtracts `H Phi_N(b)`). Here `N eta*=o(1)`, so `Phi_N(b)=1+o(1)` and `Phi_N(b(1+1/N))=O(1)`; hence the upper budget is `1+o(1)`. Jensen (or `Phi_N(t)>=t`) gives the matching lower `b(1-R)=1-o(1)`.

### Unit-cube halo check

Assuming the stated lift/support `Q=[0,1]^N` and `f_N<=cap_N 1_Q`, direct domination gives `M_S f_N <= cap_N M_S 1_Q`. Put `L=log(cap_N)`. In the region where exactly `k` coordinates lie outside the cube, product factorization of the uncentered strong maximal function and logarithmic coordinates give volume at most `cap_N L^k/k!`; the coordinate choice and two sides contribute at most `binom(N,k)2^k`. Therefore

`|{M_S f_N>1}| <= cap_N sum_{k=0}^N (2N L)^k/(k!)^2`.

Since `N log(cap_N)->1`, the right side is bounded by `cap_N I_0(2 sqrt(2N log(cap_N))) -> I_0(2 sqrt(2))`. This establishes a bounded halo-volume upper bound only. It does not convert the adaptive pointwise gap into any endpoint-constant lower bound.

The exact finite Wolfram rows and the high-precision `weighted_check60` rows are consistent with these algebraic conclusions; they are useful corroboration, not a replacement for the endpoint-reduction proof.

Integration: the required finite hypothesis is corrected in the final proof. The general rank range M q_N²/N→0 follows from the displayed exact crossing and budget bounds.
