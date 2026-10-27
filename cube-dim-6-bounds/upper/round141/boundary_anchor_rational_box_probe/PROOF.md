# C141: finite rational lower probes from five existing boundary-anchor jets

This round proposes128 fixed original unresolved terminal C129 leaf boxes, ranked only by their exact worst frame32 quadratic distance to the R137 boundary center c=(1,1,1,3/8,1/4,1/8). There is no spatial refinement or new point-window query. Only the five already accepted R137 jets are allowed. Targets g0=1/65000 and gstar=1/63500 are accounted separately. The source center may be outside a selected box: the A057 positive-band log-Hessian lower bound holds on all R6 and pays the entire center-to-point segment, independently of ordered-wedge boundaries or singular Fourier strata.

Fix one source band w=[a,T]. With its original certified center value interval0<L<=F_w(c)<=U and original face-gradient interval rectangle, define for a displacement d at a target point

    b(d)=sum_i min(g_i^- d_i,g_i^+ d_i),
    f(d)=min(U,max(L,b(d))), q(d)=sum_i d_i²/32,
    r(d)=b(d)/f(d)-q(d)/a²,
    T_w(d)=f(d) exp(r(d)).

This is the original exact robust lower tangent value, not the actual band value and not an upper bound on it. The shared f in the prefactor and gradient denominator is retained. For each fixed admissible positive center value f and fixed gradient g, log f+g·d/f-q(d)/a² is concave in d. Its infimum over the independent source interval rectangle is also concave, so log T_w is concave; equivalently, minimum over jet parameters and minimum over the spatial rectangle commute and each fixed-parameter minimum is paid at a corner. Consequently inf_B F_w>=inf_B T_w=min_(v in32corners) T_w(v). The corner rule belongs to this true robust tangent, never to the polynomial surrogate below.

For every real r and every odd n in the fixed set D={1,3,5,7,9,11,13,15}, let P_n(r)=sum_(k=0)^n r^k/k!. Taylor's integral remainder (or its Lagrange form) is nonnegative because the omitted degree n+1 is even and exp is positive. Thus P_n(r)<=exp(r) for every real r, including negative r. Define p(r)=max(0,P_1(r),P_3(r),...,P_15(r)). Since exp(r)>0, p(r)<=exp(r). This maximum need not be concave; none is assumed. At each original box corner v, C_w(v)=f(v)p(r(v))<=T_w(v). It follows purely from the true T_w corner rule that

    F_w(x)>= L_w(B):=min_(v in32corners) C_w(v),  for every x in B.

A zero p or L_w is a valid nonnegative bound and cannot be labeled a strictly positive certified floor. Exact Fraction recurrence computes only the listed finite polynomials. No exponential, native kernel, or box checker is called. The program records the selected odd degree or zero and all exact f-clamp/exponent/lower data needed to inspect each unique rational source-point evaluation.

The five fixed bands have disjoint interiors and matching adjacent endpoints of zero time measure. Therefore G(x)>=sum_w F_w(x)>=sum_w L_w(B). The minimum32corners is taken separately within each band before adding bands; no min(sum bands) corner rule is used, since a positive sum of log-concave functions need not be log-concave. Every band uses its own a in q/a². Source windows, centers, values and gradients are never changed. These are lower-only certificates: failure to reach either target is unresolved, not a genuine g counterexample, actual-G upper bound or proof of impossibility for the original robust exponential formula.

The fixed batch contains128 boxes, five bands and32 corners, at most20480 box-band-corner cases. One shared Taylor recurrence through degree15 per unique source-point identity needs at most307200 exact updates. No new leaf or partition is created. Native/exp/checker calls, new semantic queries, and subdivisions are all forbidden. Each selected box is accepted at g0 and/or gstar only by comparing its exact sum of bandwise lower bounds with the corresponding exact target. Any hash, schema or budget mismatch stops without retry. Execution remains UNRUN until main and independent exact-hash gate approval.

C129 terminal metadata records20 fully closed roots among the55 continued roots, but6 were already closed at start; C129 added14 root closures, moving the whole325-root forest from276 to290 closed roots. The remaining35 roots have10266 unresolved terminal leaves. Even if some C141 boxes pass gstar, neither those remaining leaves nor previously g0-only paid leaves are thereby paid at gstar. No global g0/gstar coverage or new C6 bound is claimed here. Separate mathematical rechecking of saved native/certificate artifacts remains the assigned independent audit's responsibility.
