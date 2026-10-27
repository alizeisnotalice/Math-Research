# R026 independent pure audit: limiting free-buffer family

2026-10-27. **Accepted; no source mathematical erratum found.** The derivative chain, all-t>1unique minimum, analytic(5,6)bracket and fixed83new-evaluation contract are correct for the limiting free-buffer upper-expression family. No mathematical numerical or symbolic software, endpoint scalar, midpoint or candidate was evaluated by this auditor. Parent retains proof and numerical-integration responsibility.

## 1. Sources and inherited operator theorem

This audit binds to R026proof SHA256322fb89887e6612b83bf04dd29cf1735407fe667de9d06235c24ab1135774a6e, registration7809c7b57c23152cd9f2b2ca7ef43433a9819cced2ec8313dd31fc9939f47380 and parameterse0ca97979452e512faacd9f3073f357fdb5574c19f48eb4396b5623d71517c50. The whole frozen proof/spec/parameters, sources and receipt were read. The accepted R025independent free-buffer proof was also read in full, at SHA256614f18791d22a23018ea0767ce33c615ee273bd41ab3a4fdba0bddf24c8d5c87. Its prior independent audit is an accepted dependency, not reexecuted by this review.

R025already proves, for each fixed t>1and positive tau, an allL1/allpositive-common-radii centered-body weak bound. It first proves finite-ainequalities and takes the infimum of the scalar coefficients as the auxiliary dilation a tends to infinity. The input operator and left side are unchanged; selectors, increments and witness families need not converge. Optimizing tau after that valid scalar limit yields

    C7<=F0(t)=[(1+1/t)^(7/2)+sqrt(1+7log t)]^2.

The present calculus is on this **limiting scalar family**. It does not assume a finite selector at a=infinity, change the original centered cube, restrict f/radii, or interpret F0as an actual operator value. The independent clone HEAD isb135609efca22428ac8f6a9c9485d114a4a03846; only audit/round026/buffer_parameter was written. Author and lower/previous sealed files remain unchanged.

## 2. Independent derivative chain and safe squaring

Put H=(1+1/t)^(7/2)+sqrtL, L=1+7log t>0 and W=t^(7/2)/(t+1)^(5/2)>0. Direct differentiation gives

    H'=7/(2t)[1/sqrtL-1/W]
       =7/[2t W sqrtL(W+sqrtL)] *G,
    G=t^7/(t+1)^5-1-7log t.

The prefactor is strictly positive. Rationalizing is valid because W and sqrtLare positive; no sign information is lost by replacing their difference with the squared difference. F0'=2HH'has the same sign. Differentiating G independently gives

    G'=E/t,
    E=t^7(2t+7)/(t+1)^6-7.

The positive tdivisor cannot change the sign. No numerical differentiation or symbolic software was used.

## 3. Complete-domain uniqueness proof

Let J=E+7>0. Its logarithmic derivative is

    J'/J=7/t+2/(2t+7)-6/(t+1)
         =(t+7)/[t(t+1)]+2/(2t+7)>0.

Thus E is strictly increasing for every positive t. E(1)=-439/64<0, and for t>=1the rational inequalities t+1<=2t and2t+7>=2t give J>=t^2/32, so E tends to infinity. It has precisely one zero e0>1. Hence G strictly decreases on(1,e0)and strictly increases on(e0,infinity).

G(1)=-31/32<0, so the first decreasing portion cannot contain a zero and the minimum is negative. Meanwhile

    G(t)>=t^2/32-1-7(t-1)

by log t<=t-1. This tends to infinity. Therefore there is exactly one t_*>e0with G(t_*)=0; G is negative before it and positive after it on the entire domain t>1. H and F0are accordingly strictly decreasing then strictly increasing. This proves a unique **global minimum of F0over all t>1**, not just local uniqueness on a proposed bracket. No endpoint computation, numerical root finder or bounded scan is part of the argument.

The endpoint t→1cannot hide a better value because the function is strictly decreasing from that endpoint up to t_*. No assertion that E is positive on the whole domain is needed. Bisection can use the globally established sign pattern even if an initial interval included part of the earlier decreasing-Gbranch.

## 4. Analytic bracket is genuinely proved

The source's elementary exponential-series arguments5/2<e<3are correct, yielding log5>4/3and log6<2. For an independent lower-log route, phi(t)=log t-2(t-1)/(t+1)has derivative(t-1)^2/[t(t+1)^2]>0and phi1=0, so log5>4/3directly. An independent upper-log route is exp2>1+2+2+4/3>6, so log6<2. These are analytic proof steps, not software endpoint values.

Consequently the exact rational estimates in the frozen source are valid:

    G5<78125/7776-31/3=-2227/7776<0,
    G6>279936/16807-15=27831/16807>0.

Together with uniqueness, they already prove t_*in(5,6). This differs from a merely proposed unverified numerical bracket. The later outward evaluations of the same signs are implementation checks and trace inputs; they do not supply the missing theoretical existence argument, since none is missing here. This audit has proved the signs analytically but executed zero endpoint scalar evaluations.

## 5. Eighty-step localization and candidate value

Once separately released, the numerical owner must use the exact bracket[5,6], require outward strict endpoint signs, then perform exactly80exact rational midpoint Gtests. The global sign pattern justifies replacing the negative or positive endpoint according to a strictly separated sign. A zero-containing ball is unresolved, not proof of equality and not a reason to pick a sign, increase precision or alter the bracket.

The initial width1halves eighty times to2^-80. The final exact midpoint t_cis within2^-81of t_*. Its value F0(t_c)is a legal operator upper expression because every t_c>1is legal in the inherited theorem, independently of its approximation error. Its own scalar interval does **not** automatically enclose F0(t_*). The source correctly avoids a numerical minimum-value claim or a new interval-Fover-bracket evaluation.

The strict rational grid upper floor(U*10^12)+1over10^12is larger than the outward upper endpoint even if Uis already exactly on the grid. For a claimed strict improvement, compare the new upper endpoint to the actual saved old lower endpoint; a display or rounded grid numeral is insufficient. A nonseparation should be reported honestly without adding rows, while a sound candidate remains a legal bound in its own right.

## 6. Exact old-row reuse and functional accounting

The accepted R025result has winner t='11/2' and its terminal id1is the completed optimized-tau row at the same t. Their full interval and grid payload match. This is a metadata/source binding, not a re-evaluation of F0. The future source gate correctly requires those bytes and their manifest, and stops if the old row cannot be bound rather than silently recomputing it.

The first bisection midpoint is also11/2. Computing **G**there is a new functional, explicitly authorized as the first of the eighty sign tests; it neither repeats the old F0calculation nor obtains its sign from the old F0value. The new scientific count is therefore two endpoint Gvalues plus eighty midpoint Gvalues plus one candidate F0value, totaling83. The old F0row is reused with separate provenance, not counted as a new evaluation. No Evalues, alternate fixed t, additional F0at11/2, optimizer or adaptive parameter is within this contract.

Positive square-root evaluation r^3sqrt(r)for r^(7/2), r=1+1/t, is correct. The internal A+L+2sqrt(AL)cross-check is an equivalent algebraic form, not a new parameter row. Actual log/sqrt/arithmetic call counts and raw widths still require honest reporting by the executing owner. The author did not release or execute this batch through writing the theory files, and this audit does not trigger it.

The audit's read-only source checks include the R025result/terminal/manifest hashes and exact JSON payload binding of the winner to id1. No interval endpoints are converted to floating values and no logarithm, root or F0/Gexpression is evaluated for that metadata check.

## 7. Strict accepted range and stopping

Accepted: inherited limiting-family use; exact derivative factorization; Emonotonicity and global Gsign pattern; unique family minimum; pure analytic(5,6)bracket; exact80step localization and83new-functional contract; correct old-row reuse and candidate/strict-rounding quantifiers. No mathematical source erratum was found.

Not claimed: any software endpoint signs, successful bisection steps, exact t_*value, candidate interval/value, numerical minimum-value enclosure, new numerical C7bound, or optimality/lower bound for actualC7. An exact lower bound on this family minimum would still be a lower bound on upper-certificate expressions only. Other selection mechanisms or unrelated operator arguments may improve them. The lower branch stays paused; no new feature or analysis is introduced.

Only full source reading, independent mathematical reasoning and metadata/hash operations were used. No mathematical numerical or symbolic software, new endpoint/midpoint/candidate calculation, sample/cache functional, phase/DD, MATLAB, source edit, commit or push occurred. Source hashes, receipt and manifest accompany this review. The audit is complete and stops.
