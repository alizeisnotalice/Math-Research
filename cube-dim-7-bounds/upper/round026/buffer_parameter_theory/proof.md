# R026: unique minimum of the free-buffer limiting parameter family

2026-10-27. Pure Sol theorem/specification, no mathematical scalar/root calculation, MC, MATLAB, commit or push. Numerical execution needs a later explicit main release. Parent retains proof and integration responsibility. This is optimization of a proved upper-bound family, not of the true unknown operator constant.

## 1. Accepted family and operator scope

The accepted R025 free-buffer theorem proves, for every t>1, the unchanged Lebesgue centered cube maximal operator on generalL1and all positive common radii satisfies

    C7<=F0(t),
    F0(t)=[(1+1/t)^(7/2)+sqrt(1+7logt)]².                   (1)

Its limiting step first proves finite-parameter inequalities and then takes an infimum of their scalar coefficients as an auxiliary dilation tends to infinity. No selector, cover or input-dependent operator limit is needed. The frozen R025 theorem is an input to the present calculus proof, and its exact source bytes are locked. We do not reinterpret(1)as a pointwise kernel inequality or as an equality for C7.

Put H0(t)=(1+1/t)^(7/2)+sqrt(L(t)), L(t)=1+7logt>0. The F0andH0derivatives have the same sign because F0'=2H0H0'andH0>0.

## 2. Independent derivative computation

Direct differentiation gives

    H0'(t)=7/(2t)*[1/sqrtL-(t+1)^(5/2)/t^(7/2)].

Let W=t^(7/2)/(t+1)^(5/2)>0. Rationalizing the difference yields the exact positive-prefactor identity

    H0'(t)=7/[2t W sqrtL(W+sqrtL)] *G(t),
    G(t)=t^7/(t+1)^5-1-7logt.                              (2)

The squaring comparison here is legitimate only because W,sqrtLare positive. Differentiating the rational term once more gives

    G'(t)=E(t)/t,
    E(t)=t^7(2t+7)/(t+1)^6-7.                              (3)

There is no sign-changing omitted t factor in(2)-(3).

## 3. Global sign analysis for every t>1

Write J(t)=E(t)+7=t^7(2t+7)/(t+1)^6>0. Its logarithmic derivative is

    J'/J=7/t+2/(2t+7)-6/(t+1)>0,

since already7/t-6/(t+1)=(t+7)/[t(t+1)]>0. Hence E is strictly increasing for all t>0. At the endpoint

    E(1)=9/64-7=-439/64<0.

For t>=1, t+1<=2t and2t+7>=2t, so J(t)>=t²/32. Thus E tends to infinity and has exactly one zero e0>1. Equation(3)shows G strictly decreases on(1,e0)and strictly increases on(e0,infinity).

Also G(1)=1/32-1=-31/32<0. Its minimum is therefore negative. On the other hand t^7/(t+1)^5>=t²/32andlogt<=t-1for t>=1, so G tends to infinity. Consequently G has exactly one zero t_*>e0in the entire domain t>1, with negative sign before it and positive sign after it. Equations(1)-(2)now show

    F0 strictly decreases on(1,t_*),
    F0 strictly increases on(t_*,infinity).                 (4)

This proves a unique global minimum of this complete scalar family. It is not a global optimum or lower bound for the actual C7.

## 4. Exact analytic bracket, independent of any scalar run

The elementary exponential series gives5/2<e<3. The lower inequality follows by retaining1+1+1/2and positive remaining terms; the upper one follows by dominating1/k!for k>=2by the geometric sequence1/2^(k-1), strictly for later terms.

Because e<3and3^(4/3)<5(the latter is81<125), one has log5>4/3. Therefore

    G(5)<78125/7776-31/3=-2227/7776<0.

Because e>5/2, e²>25/4>6andlog6<2. Consequently

    G(6)>279936/16807-15=27831/16807>0.                    (5)

These are purely analytic rational endpoint bounds, not software results. In particular t_*is already proved to lie strictly in(5,6). Later outward scalar evaluation of the endpoint signs is an implementation check and a source for the bisection trace, not the missing existence proof.

## 5. Fixed finite bisection and legal candidate expression

After main release, use the exact initial bracket[5,6], certify its two Gsigns, and perform exactly80rational midpoint steps. The global sign pattern in(4)justifies replacing the negative endpoint if G(mid)<0, or the positive endpoint if G(mid)>0. A ball containing0does not establish equality: stop without guessing, choosing a new point or automatically raising precision.

After80successful steps the bracket width is2^-80. Its exact dyadic midpoint t_c is within2^-81of t_*. The theorem(1)gives C7<=F0(t_c), so one fresh outward evaluation and strict rational upper rounding gives a legal bound independently of the closeness estimate. The value interval for F0(t_c)does not by itself enclose the exact family-minimum value F0(t_*); no such numerical minimum-value claim is made here.

The previously evaluated t=11/2is a simple rational candidate and should be reused exactly from R025's completed result/terminal row, with manifest/source binding. Its F0value must not be recomputed merely to fill a new row. The first bisection midpoint is also11/2, but its newly evaluated Gsign is a different authorized functional; reusing the old F0row does not supply that sign.

Thus the proposed new scientific evaluations number83: two endpoint Gvalues,80midpoint Gvalues, and one F0candidate value. A reused F0(11/2)row is separately labelled, not a new computation. There are no Evalues, extra fixed t, optimization-library calls, interval-Fover-bracket evaluation or adaptive rows in this scope. Intermediate library operations and exact consistency checks must be reported honestly.

## 6. Limits of the conclusion

The unique minimum in(4)concerns only this explicit family supplied by the accepted free-buffer selection proof. Choosing any t>1preserves the original allL1/allcommon-radii operator theorem; no Green floor is involved. Taking an infimum of these upper bounds is legal, but even an exact lower bound for that family minimum would not imply a lower bound for C7. Other selection arguments or unrelated operator estimates may be better.

No mathematical scalar, endpoint sign, bisection step, candidate value or new numerical C7bound was evaluated by this theory owner. The existing R025 t11/2certificate retains its executing owner's provenance. The root-localization and strict candidate comparison are future finite certificate tasks only, pending main release and independent review.
