# The Hessian of the same positive band

Fix0<a<T and f(b)=atan(b/a)−atan(b/T). This is the C016 band, not the full Green function. For a six-node labelled tuple z=(b0,...,b5), F(z)=DD5(f;z). With repeated entries interpreted as normalized Hermite derivatives,

    ∂_i F=DD6(f;z+bi),
    ∂_ij F=DD7(f;z+bi+bj)       (i≠j),
    ∂_ii F=2 DD7(f;z+bi+bi).                         (1)

For example take a rectangular contour enclosing all real nodes inside the analytic strip |Im z|<a. Its Cauchy divided-difference representation is

    F=(2πi)^-1∮ f(ζ) /∏_(k=0)^5(ζ−bk) dζ.

Differentiating a labelled node inserts one additional corresponding pole. Differentiating that same label twice differentiates both copies of its factor and gives the coefficient2. Distinct labels carry coefficient1, even if their numerical values coincide. This proves(1) for distinct real nodes and, by real analytic continuation in the labelled nodes, at all zero/repeated configurations. No sign(0) convention is involved.

The normalized jets needed through order7 are obtained by integrating the coefficient series of u/(u²+(b+z)²). If that series is Σg_k z^k, then

    g0=u/(u²+b²),
    g_k=−(2b g_(k−1)+g_(k−2))/(u²+b²),  g_-1=0,
    f^(k)(b)/k!=(g_(k−1),a−g_(k−1),T)/k.

At zero the value and even jets vanish, while odd jet k equals
(-1)^((k−1)/2)(a^-k−T^-k)/k. The program retains them. An eight-entry Hermite table evaluates DD7, including eight identical zero nodes.

For distinct base nodes the implementation uses an algebraically equivalent fast formula. Put w_i=∏_(j≠i)(bi−bj)^-1, c_i=f(bi)w_i, S_i=Σ_(j≠i)(bi−bj)^-1 and d_i=∂_iF. The divided-difference recursion gives

    H_ij=(d_i−d_j)/(bi−bj), i≠j.

Differentiating d_i=f'(bi)w_i−Σ_(j≠i)(c_i+c_j)/(bi−bj) gives

    H_ii=f''(bi)w_i−2f'(bi)w_i S_i+c_i S_i²
                  +Σ_(j≠i)(c_i+2c_j)/(bi−bj)².          (2)

Thus the fast distinct branch computes exactly(1). Every repeated base tuple instead uses(1) directly in triangular Hermite tables. There is no zero-excluding branch. If finite precision cannot separate two distinct accurate rational nodes, the complete invocation fails; no term is skipped or silently assigned zero.

## Physical signs, permutations and symmetry

In a cone, bj=Σ_k r_(j,k) s_k x_(σ(k))/q_j. Its Jacobian is J_(j,σ(k))=r_(j,k)s_k/q_j, independent of x. Therefore the physical Hessian is exactly Jᵀ Hnode J, with no second derivative of the node map. The Fourier cone weight and final factor(2π)^-6 are the same as C016.

For node indices i≤j and physical y indices k≤l, the chain coefficient before dividing by q_i q_j is

    r_(i,k)r_(i,l)                              if i=j,
    r_(i,k)r_(j,l)+r_(j,k)r_(i,l)                if i<j.

This includes the two off-node contributions; it works also for k=l. Then multiply by s_k s_l and send the pair(k,l) to the sorted physical pair(σ(k),σ(l)). The21 saved fields are ordered(0,0),(0,1),...,(0,5),(1,1),...,(5,5). Mirroring gives the complete symmetric6×6 Hessian; no diagonal or mixed field is omitted.

Since f is odd and DD order5 is odd, F(−b)=F(b). The node Hessian is even under global negation; both Jacobian factors change sign, leaving the physical Hessian unchanged. The existing exact global-negative pairing therefore still has weight2. The new code does not canonicalize by an additional negative reflection.

This mathematical derivation applies to the same matching band. The finite acceptance batch does not supply an independent complete full-group Hessian reference, a box oracle or a Taylor remainder. It cannot replace a full-G derivative by a band derivative. No new grad(p) is computed and no second-order Euler identity is claimed.
