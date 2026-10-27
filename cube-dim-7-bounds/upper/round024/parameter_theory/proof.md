# R024: the unique minimum of the proved seven-dimensional parameter family

2026-10-27. Pure Sol derivation/specification, with no mathematical scalar computation, optimizer, MC, MATLAB, commit or push. Numeric work requires a later explicit main release after source locking. This studies the already valid upper-bound expressions, not the unknown best weak-type constant.

## 1. Eliminate tau exactly, including the fixed a36 expression

The accepted allL1/allpositive-common-radii centered-cube theorem gives, for every a>1,tau>0,

    C7<=K(a,tau)=(1+A(a)/tau)(tau+B(a)),
    A(a)=(1+2/a)^7,
    B(a)=(1+1/a)^7(1+7 loga).

Both A,B are positive. Expanding K=A+B+tau+AB/tau and applying AM-GM proves the unique minimizing threshold at each fixed a is

    tau_*(a)=sqrt(A(a)B(a)),
    F(a)=min_tau K(a,tau)=(sqrt(A(a))+sqrt(B(a)))².           (1)

No scalar approximation of tau is required to use this theorem: it is one fixed positive real parameter determined only by a and dimension, not by f,lambda or the selected covering family.

At a36 the exact new expression is

    F(36)=[(19/18)^(7/2)
               +(37/36)^(7/2)sqrt(1+7log36)]²,
    tau_*(36)=(703/648)^(7/2)sqrt(1+7log36).                 (2)

It is no greater than the accepted R023 a36,tau7 expression by(1). A strict numerical comparison of complete intervals, if achieved later, may quantify an improvement; no display or equality assumption establishes that comparison here.

## 2. Derivative factorization independently computed

Put H(a)=sqrt(A(a))+sqrt(B(a)), L(a)=1+7loga and N(a)=a-7loga. Direct differentiation gives

    (sqrtA)'=-7sqrtA/[a(a+2)],
    (sqrtB)'=7/(2a)*((a+1)/a)^(7/2)
                                  *N/[(a+1)sqrtL].

Thus

    H'(a)=P(a)(R(a)-1),
    P(a)=7((a+2)/a)^(7/2)/[a(a+2)]>0,
    R(a)=1/2*((a+1)/(a+2))^(5/2)*N(a)/sqrt(L(a)).            (3)

In particular F'=2H P(R-1), so F and H have the same derivative sign. The5/2exponent in R comes from cancelling the extra(a+2)/(a+1)factor against the7/2power; it must not remain7/2.

## 3. There is exactly one critical point on all a>1

On1<a<=7, N'=1-7/a<=0 and N(1)=1, so N<=1. L>1 and0<(a+1)/(a+2)<1. If N<=0 then R<=0; if0<N<=1 then R<1/2. Therefore R<1 throughout this first range.

For a>=7 define S(a)=N/sqrtL. Its derivative has the sign of

    T(a)=(1-7/a)L-7N/(2a)
        =7loga-5/2-7/a-(49/(2a))loga.

For a>=7,

    T(a)=7loga(1-7/(2a))-5/2-7/a
         >=(7/2)loga-7/2>0.                                (4)

The final strict inequality uses log7>1. For an elementary justification, e=sum1/k!<3<7, bounding the terms k>=2by the corresponding geometric series, so log7>1. Thus S is strictly increasing on[7,infinity).

Also N(7)=7(1-log7)<0, N'(a)>0 for a>7 and N tends to infinity. Hence it has a unique zero a0>7. For7<a<=a0, R<=0. For a>a0, S>0 and S'>0, while the positive factor((a+1)/(a+2))^(5/2)also strictly increases. It follows that R strictly increases on this latter interval. The product argument is used only where S is positive; multiplying two increasing functions would not by itself prove increasing R on the negative-S range.

Finally R(a) tends to infinity, since R is asymptotic to a/[2sqrt(7loga)]. This elementary limit follows, for example, from loga<=2(sqrta-1), obtained by integrating1/t<=1/sqrtt on[1,a]. Thus N grows at least linearly after a sufficiently large a while sqrtL grows at most as a constant times a^(1/4).

By continuity and strict increase on(a0,infinity), there is exactly one a_*>a0 with R(a_*)=1. Combining all ranges,

    F strictly decreases on(1,a_*),
    F strictly increases on(a_*,infinity).                  (5)

So a_* is the unique global minimum of this entire a-family. The jointly minimizing pair of this two-parameter family is(a_*,sqrt(A(a_*)B(a_*))). This statement is exact analytic theory; no numerical location or value was computed here.

## 4. A finite, auditable root localization and legal candidate bound

A later authorized numerical batch may first certify R(35)-1<0 and R(37)-1>0. These signs are proposed acceptance tests, not results assumed by this proof. If they hold, uniqueness places a_* in(35,37). At every step take the exact rational midpoint of the current bracket, evaluate the sign of R(mid)-1 by outward arithmetic, and replace only the endpoint of matching sign. The global sign pattern(5)justifies each replacement. If a sign enclosure contains0, stop without guessing or increasing precision automatically.

After exactly80successful midpoint steps, the bracket width is2/2^80=2^-79. Its exact dyadic midpoint a_c has distance at most2^-80from a_*. Evaluating F(a_c)produces a legal direct C7upper expression regardless of its proximity to the minimum, because a_c>1 and the original bound holds at every a. The numerical interval for F(a_c)is not by itself an interval enclosing the exact family minimum F(a_*). No numerical minimum-value enclosure is claimed without a separate certified interval/Taylor argument.

The finite contract evaluates only F(36), the two endpoint signs,80midpoint signs and one F(a_c). It does not scan a, run a numerical root finder with an undocumented stopping tolerance, change the bracket, or test other parameters. Positive square roots may be implemented as r³sqrtr for r^(7/2)and r²sqrtr for r^(5/2), avoiding an unnecessary general real-power operator. The exact rational a values must never be rounded before forming the balls.

## 5. Quantifiers and limits of optimization

The accepted theorem already holds for each parameter pair on every generalL1input and all positive common cube radii. Taking the infimum of these upper bounds, or choosing a_c and its exact positive tau_*, preserves that same operator and scope. It introduces no restriction on f or scale and no Green-floor assumption.

Conversely a lower bound on the parameter-family minimum is not a lower bound on C7. The family may be intrinsically slack relative to the true optimal operator inequality. The uniqueness in(5)means only one best a within these specific proved expressions after tau is eliminated. It is not global optimality among covering methods, not the best possible cube maximal constant, and not a C7counterexample. R023's finite grid numeral remains a separate already proved upper bound, with its actual numerical source retained.
