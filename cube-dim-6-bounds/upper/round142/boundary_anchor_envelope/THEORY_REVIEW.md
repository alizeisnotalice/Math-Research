#142: independent mathematical review of the proposed141 rational envelope

**PASS for the stated mathematical mechanism, conditional on authentic matching137 jets, exact implementation and finite registration.** No141 batch, corner table, exponential, native kernel or historical baseline is recomputed here. This is a theory review; the author's final frozen driver has not yet been audited. The reviewed degree set is the finite odd set{1,3,5,7,9,11,13,15} supplied by the141 author.

For one positive band w=[a,T], A057 proves on all R6

    log B_w(c+d)>=log B_w(c)+grad B_w(c).d/B_w(c)-||d||^2/(32a^2).

Thus c may be a boundary anchor and may lie outside the target box. The complete joining segment is paid by the global positive-cutoff log-Hessian bound, including ties and trueS2/S3. For the face x1=c1=1 use its five-coordinate displacement; no derivative term in coordinate1 is omitted because that displacement iszero. If those first coordinates differed, a five-gradient rule would not suffice. Keep the actual source center in metadata: do not replace it by a box midpoint to evade the old executable checker's center-in-box guard. The new driver needs this global theorem, not an unauthorized mutation of A057 or an n4 continuity argument.

Let0<L<=B_w(c)<=U and g_i^-<=partial_i B_w(c)<=g_i^+. At a displacement d define

    b(d)=sum_i min(g_i^- d_i,g_i^+ d_i),
    f(d)=clamp(b(d),[L,U]),
    r(d)=b(d)/f(d)-q(d)/a^2,   q(d)=sum_i d_i^2/32.

For fixed f>0 the exponent increases with g.d, so minimizing the interval gradient gives b. The derivative of log f+b/f is(f-b)/f^2; hence the exact minimizer over[L,U] is that clamp, including b<=0. The robust tangent value C(d)=f(d)exp(r(d)) is the minimum over the independent interval uncertainty set, retaining the same f in the value and gradient denominator. Actual within-jet correlations may improve the bound but cannot invalidate this conservative minimization. Use one complete authentic same-window137 accepted jet, without exchanging its fields with fullG or a different window/attempt.

The32-corner rule comes from the **original exponential tangent**. For each fixed admissible(f,g), its log is a concave quadratic in d. Consequently its minimum on any rectangular displacement box occurs at a vertex. Minimization over the uncertainty set commutes with the finite vertex minimum; equivalently log C is the infimum of concave quadratics and is concave. Therefore

    inf_targetbox B_w >= min_(32 vertices v) C(v).

The anchor need not be the box center, and the vertices must be those of the entire actual target rectangle shifted by that anchor. An ordered-domain intersection does not justify dropping raw rectangle vertices without a separate polytope argument.

For every real r and odd n, Taylor's Lagrange remainder gives

    exp(r)-P_n(r)=exp(xi) r^(n+1)/(n+1)!>=0,
    P_n(r)=sum_(k=0)^n r^k/k!,  xi between0 and r.

Thus E(v)=f(v)max(0,P_1(r(v)),P_3(r(v)),...,P_15(r(v))) satisfies0<=E(v)<=C(v). Clipping atzero is legitimate because exp(r)>0, including at negative r; every registered odd order is valid on the entire real line. Taking the maximum across finitely many such lower bounds is valid. It follows directly that min_v E(v)<=min_v C(v), so min_v E(v) is a valid whole-box lower. **No concavity of P_n or max(0,P_n) is used or implied.** Applying a corner rule directly to the polynomial/clipped spatial substitute would be an unproved different argument. A nonpass may reflect this rational relaxation, not failure of the actual band or exponential envelope.

For disjoint positive time bands the safe aggregate is sum_w[min_v E_w(v)], with each band's own source anchor, gradient interval, f-clamp and its actual a_w^-2 curvature. It is not generally min_v[sum_w E_w(v)]: a sum of log-concave exponentials need not be log-concave, and separate minimizing vertices may differ. Combining the band jets into one true prefix is possible only with its own valid positive-mixture curvature (at least the minimum time endpoint), not the first band's larger a. Multiple anchors for the **same** band can give alternative lower certificates whose maximum is valid; summing them double-counts that band. The proposed141 author formula explicitly uses per-band vertex minima before summing and is therefore safe.

Finally a saved old leaf certified only at g0=1/65000 cannot be relabeled as passing g*=1/63500. Its saved absolute lower must itself be shown>=g*, or a separate valid stronger certificate must be provided. Testing frozen C129 unresolved boxes does not upgrade old g0 leaves or certify a complete g* cover. All pass counts remain local until exact domain coverage and the target of every retained leaf are paid. No global g, actual model-cost lower bound or newC6 statement follows from this review.

Two bounded independent local SymPy proof sanity checks actually passed: d/df(log f+b/f)=(f-b)/f^2, and the five-dimensional original tangent log-Hessian is exactly−I/(16a^2). They read no141 data and evaluate no corners, numeric exponentials or native kernels; they do not replace the analytic remainder or global A057 proof.

Read-only mathematical sources: A057 positive_mixture_curvature/theorem.md §§2,4 and its box_certificate.py (whose center-in-box restriction is an API restriction); the already executed137 same-window jet/results identities; the141 author's stated formula supplied directly to this reviewer. No source or saved scientific artifact is modified by this file.
