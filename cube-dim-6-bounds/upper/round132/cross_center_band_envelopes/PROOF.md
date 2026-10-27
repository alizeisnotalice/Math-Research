# C132: cross-center complete-box positive-band envelopes

This is a local certification mechanism for the original body32 unit-weight model. Its inputs are only the frozen C127 recovery cache and static C129 recovery geometry seed. It does not prove a global Green floor or a C6 bound.

For each fixed positive band W=[a,T], 0<a<T, let F_W(x)=integral_a^T p_t(x)dt. A057's fixed-fiber argument proves F_W is finite, positive and C2 everywhere on R6 and Hess log F_W >= -I6/(16a²). The argument uses AAᵀ=32I6 and R=Aᵀ/32. Its differentiation/normalization is legitimate for a>0. Every real center c and every segment c+s(x-c), 0<=s<=1, lies in this global domain. Integrating the Hessian lower bound yields

    F_W(c+u) >= F_W(c) exp(grad F_W(c)·u/F_W(c) - ||u||²/(32a²)).

No part requires c to lie in the target box. That check in the old executable is an interface restriction. C132 does not defeat it with altered metadata: it authenticates each original jet against its original expected center/window, then applies the displayed global theorem in an independent function. On the face x1=1, u=(0,d), so the quadratic is sum d_i²/32. The target box and segment may cross ordered-face boundaries, repeated/zero Fourier projections, or the singular strata of the uncut Green function; the positive band remains globally smooth.

For any chosen source center c, its true source jet certifies 0<L<=F_W(c)<=U and g_i^-<=partial_(i+1)F_W(c)<=g_i^+. For each of the target box's32 displacement corners d, compute exactly

    b_d=sum_i min(g_i^- d_i,g_i^+ d_i),
    f_d=min(U,max(L,b_d)),
    C_d=f_d exp(b_d/f_d - ||d||²/(32a²)).

For a fixed positive f and gradient vector, the log tangent is concave in d; hence its minimum over the whole rectangle is attained at a corner. Minimizations over the jet rectangle and spatial rectangle are both minima and commute without a minimax assumption. Minimizing the gradient dot product gives b_d. The derivative of log f+b_d/f is (f-b_d)/f², giving precisely the displayed clamp. Arb512 encloses each exponential outward, and L_(W,c,B)=min_corner lower(C_d) is a reliable lower bound on the entire box B. The minimum of the outward upper endpoints is an upper enclosure of this robust tangent minimum, not an upper bound for F_W.

For a single band W, each center proves independently that F_W(x)>=L_(W,c,B) for every x in B. Therefore F_W(x)>=max_c L_(W,c,B) throughout B. This maximum concerns complete-box reliable bounds. We do not assert an equality between max-min and min-max, nor apply a corner rule to a pointwise maximum of tangents. The inequality max_c inf_x T_c(x)<=inf_x max_c T_c(x)<=inf_x F_W(x) suffices, with no saddle-point or center-interchange hypothesis. In particular, values from different centers of the same band are never added.

The five fixed positive time bands have disjoint interiors and shared endpoints of zero dt measure. Nonnegativity implies G(x)>=sum_W F_W(x)>=sum_W max_c L_(W,c,B). An absent band contributes0. Source window and center value/gradient are never transformed. Euclidean nearest16 selection is only a fixed proposal choice; correctness rests on each authenticated source jet and all32 corners. Selection uses exact squared Euclidean distances on an integer lattice, with semantic-key tie-breaks.

The experiment uses only the unresolved boxes of the first lexicographic unclosed original root from the frozen geometry seed, node IDs sorted, at most128. Existing passed boxes are not recomputed. New target1/63500 and reference1/65000 are compared against the same computed local sum. Even if all selected boxes pass1/63500, existing passed leaves have not been recertified at that stronger floor; neither the root nor the full domain has thereby been proved at that target. If more than128 unresolved boxes existed, the experiment would additionally leave unresolved boxes untested.

The saved full exact M6 fraction is separately checked with rational arithmetic:

    U^5=(9/8)^6 M6 / (384^5 (1/63500)^6) <64^5.

This check is conditional on a future full-domain Green floor1/63500 and all prior ingredients. It uses the complete source Fraction, never a displayed decimal. It supplies no missing global coverage. Its exact arithmetic evidence is conditional_M6_check.json.

No exponential/native call occurs in preparation or source-only authentication. Scientific execution requires both main and independent review approval in parent_release.json. Runtime is single CPU at nice>=10, Arb512, no time cap, with at most128*5*16*32 corner exponentials. The actual selected budget, missing bands and original source hashes are precommitted in the input packet.
