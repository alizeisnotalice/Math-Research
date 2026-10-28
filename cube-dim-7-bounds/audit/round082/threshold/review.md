# R082 explicit threshold — independent ACCEPT with a tau-tail addendum

For every integer n with L=log n>=L0=2^24=16777216, the original
centered Lebesgue cube maximal constant satisfies

    C_n <= Phi81=F(A_n(t),B_n(t)-(13/10)n)
        < n log n-n log log n-(3/8)n,
    t=n/[(2/5)L], F(A,B)=(sqrt A+sqrt B)^2.

This is a NEW explicit-threshold acceptance of the already reviewed
R081 geometry and full operator framework. R081's original proof
identity remains
`1081f7052cd0fdea55186a71d062dd09137d4c9f6eaf31fe82c68871db64516d`.
It was originally an eventual theorem with uncomputed L_*. This audit
does not rewrite its proof, execution history or earlier acceptance.

All sources were read only. This is independent temporary staging,
not a Git clone. Only `audit/round082/threshold` was written. No
shared MATLAB, commit, push, producer replay, nearby-threshold test or
parameter search was used. The producer's single registered 54-gate
receipt is preserved and was not reexecuted.

## 1. Read set, identities and new actual arithmetic

Read in full: R082 registration.json, README.md, verify.py and
result.json. The registered L0 is 2^24, its source proof hash is the
frozen R081 final proof above, and its receipt links to the current
unchanged registration and script hashes. The author program computes
base envelopes, explicitly leaves mathematical acceptance separate,
and its actual PID 28988 records all 54 gates PASS. This review
supplies the all-tail argument and invokes the separately completed
R081 independent full-chain review, rather than inferring a theorem
from this finite PASS.

A NEW bounded independent script `check.py` actually ran, PID 29629.
It uses the weaker lower envelope N12(L)=L^12/12!<=exp L=n, exact
Fraction arithmetic, and the single positive Taylor term
T8(y)=y^8/8!<=exp y. It neither imports nor executes the producer
program. Eight grouped checks cover the new sensitive domains,
probabilities, complete payments, source identities and the tau-tail
polynomial. Exact values are preserved in this directory's result.json.

Display-only independent values at the sole registered base are

    N12>1.0382*10^78, x0=6710886.4,
    h0=1/281474976710656,
    rho_lower>3.8677*10^70,
    baseline mean upper<2.499*10^-62,
    baseline failure upper<3.998*10^-53,
    each count failure upper<3.853*10^-70,
    first-profile failure upper<6.165*10^-73,
    source-face Chebyshev failure<0.000354370<0.001,
    weighted cap/T8 payment<3.051*10^-27<e,
    tau/T8 payment<1.146*10^-21<e,
    complete root/T8 remainder<8.021*10^-46<1/200.

This independently confirms margins using deliberately weaker
envelopes than the original N16/S16 receipt. It does not replace the
analytic tail and original-operator arguments below.

## 2. Tau derivative checked first: explanation error and exact repair

The README's claim that the cubic alone dominates negative degrees
0--2 for all y>=4 is incorrect at y=4. The correct derivative
conclusion does hold. Keeping the quartic and setting u=y-4 gives

    yS16'-(5/2)S16
      >=53/6+(33/2)u+(27/4)u^2+(13/12)u^3+u^4/16>0.

All omitted coefficients are nonnegative. Hence
L^(5/2)/S16(L/50) decreases once L/50>=4, a condition satisfied by
L0 directly. The precise error, the -43/6 cubic countervalue and
the proof repair are separately documented in ADDENDUM.md. This is
a repair of a local explanatory sentence; no original receipt or
hash is edited. R081's exponential-limit proof did not use that
sentence. The full independent proof here can alternatively avoid
it altogether by using T8, as in Section 6.

## 3. Explicit active domain, local lemma and complete matching

Retain R081's exact s=4 tuple: qmin=1/4, qmax=.99, x=.4L,
h=L^-2, epsilon=1/(x n^2), delta=10^-5, rlo=1+delta,
rhi=3.96-delta, and every error e=epsA=10^-4.

At L0 the weaker N12 check proves

    N12>=300, x>=2, x/N12<1, x*N12^2>1,
    2h<delta, h<=.001, rlo-2h>1,
    1-h>0, 3/4-h>0, epsilon<1,
    qmax/(x*N12^2)<=h/4,
    N12/x>rhi, N12/(4x)>4rhi and >=12,
    rhi<4, 4rhi>5.

Thus the independent local lemma's explicit h<=.001 and rho>=12
requirements are covered, as are beta<2 and kappa>1. The matching
interval q_±=(r±h)/4 stays within [1/4,.99]. Every active source
has R in [r-2h,r], R>1, |w_i|<=r<=R+2h and

    g-3R/4=epsilon*q-h/4<=0.

The independent s=4 coordinate formulas and AC derivative bounds
apply throughout this interval, including |w_i|>R and its central
deletion. Sources unmatched at the band edge have zero columns and
receive no R>=1 local assertion. Full-band R,a0,G remain positive.
Both annuli lie in norm<t and are disjoint.

For ALL larger L the same guards persist: N12/x grows as L^11,
x/N12 decreases as L^-11, x*N12^2 grows as L^25, and the ratio
[qmax/(x*N12^2)]/h decreases as L^-23. Each trim/positivity guard
improves as h decreases. There is no uniformity left implicit in
this part of the previous eventual proof.

## 4. All five target-event and whole-matching payments

The s=4 full-target Fubini and bounded-log arguments give baseline
log-sum mean and variance envelopes

    Mlog=576x^2/n+48x/n,
    Vlog=2304x^2/n.

The new N12 base check gives Mlog<=e/2 and
4Vlog/e^2<=e/8. For larger L their adverse envelopes are positive
sums of L^-10 and L^-11, so both remain valid. Each endpoint count
failure is <=4/(e^2 n)<=e/8; its N12 envelope decreases as L^-12.

The target-cap positive-face fraction is at most
1+x/(2N12)-1/(2rhi)<=.9 and decreases, while negative pinned
faces already lie inside D because t>rhi. Since n>=N12>=300,
the cap failure is <=.9^(n-1)<=.9^299<=e/8. This is a per-face
bound and loses no hidden half-face factor.

The fixed strict inequality mu4=3-5log(8/5)<.65 belongs to the
frozen and independently reviewed R081 constant proof. The new pin
payment 4/N12<=1/200 gives mean <=.655x. The remaining .005x gap
and free-coordinate variance <=16x^2/n yield the actual first-profile
failure <=640000/n<=e/8. These N12 bounds decrease with L.

The five events need not be independent; their total failure is
<=5e/8<e. R081's fixed strict mass inequality then gives

    lambda(A)>=(1-e)n log(rhi/rlo)>1.37n.

For the whole matching interval, the independently derived local
constants imply product-log drift <=80/L. The new base has 80/L<=e,
and this decreases. Individual |J-1|<=6/rho<=e also improves as
rho_lower grows. Consequently the same P0=1-2e and Ulo/Uhi are
valid on every matching q, not only at the baseline.

Endpoint count calibration is paid by 2h+1/n<=3e/4, whose envelope
decreases as a sum of L^-2,L^-12. First-profile transport is paid by
80h<=.01 and also decreases. Thus sum J1<=.67x throughout matching,
and the fixed rounding .67/(1-e)<.671 stays valid. The true
inverse-source variance remains <=4x/[3(1-e)]<2x. Pins are
deterministic, and no mixed-face variance or post-cap independence
is substituted.

## 5. Entire source-retention tail and weighted-cap tail

Let a=.029 and p=.5. The source-face Chebyshev threshold is
u(x)=a*x-p, and its bad probability is bounded by

    f(x)=2x/(a*x-p)^2.

The new base check has u(x0)>0 and f(x0)<.001. Directly,

    f'(x)=-2(a*x+p)/(a*x-p)^3<0 for a*x>p.

Since x=.4L grows, the bound holds for EVERY larger L, every target,
matching q, omitted coordinate and sign. The retained fraction .999
is uniform before faces are mixed and before source-cap/soft
conditioning. This explicitly closes the previous eventual retention.

For the weighted source cap, the whole-matching count gives
nu=1/4-(5/4)e=.249875. The arbitrary positive source pin costs at
most one counted coordinate; the full captured lower-strip product
estimate gives

    caploss/k<=M exp[-K_actual],
    K_actual=(nu-1/n)*eta*x/(1+e), eta=1/198.

Replacing n by the weaker N12 gives

    K12(L)=C[nu*L-12!*L^-11],
    C=(2/5)*eta/(1+e)>0.

The base has nu-1/N12>0 and K12>4. Its derivative is
C[nu+11*12!*L^-12]>0. The cap payment is therefore bounded by
M/T8(K12(L)), which decreases and is already <e at the base.
The new weighted argument sums all actual pre-cap inverse face
coefficients, uses bsoft<=theta/G pointwise, and integrates only the
complete matched dq/q interval. Theta cancels against k exactly.
No unweighted cap mass or obsolete eta-chi condition is used.

The radial-width bound 1/[1-(h/rlo)^2]<=1+e is paid at L0 and
improves. The complete raw row and column estimates from R081
therefore remain valid. In particular all target pin signs are included
in the column coefficient and the complete column is <=2k<Mk.

## 6. Positive roots, complete tau payment and full scalar cost

The lower scalar domain is established for all L>=L0 directly:

    exp(L/2)>=(L/2)^2/2=L^2/8>=x=(2/5)L.

The last inequality holds already for L>=16, so it is valid here
without inferring a tail from a base S16 check. Hence log t>=L/2,
B>=nL/4>1.3n, and the SAME original-selection threshold
tau=sqrt[A(B-1.3n)] is positive for the entire tail. Also x>=2 and
n>=300 imply n log x>=n/2>1, so B<=1+n log t<=nL; A<=exp x.

Consequently

    tau/k<=rhi L^(5/2)exp(-L/50)
       <=rhi L^(5/2)/T8(L/50).

Its base value uses sqrt(L0)=4096 exactly, and was independently
checked <e. The right side is a fixed constant times L^(-11/2),
so it decreases for every larger L. This independent proof avoids
the flawed cubic explanation entirely; the repaired S16 derivative
argument in ADDENDUM.md also validates the producer's own S16 tail.

After weighted cap and THIS tau are paid, row-tau>=b*k with the
unchanged R081 fixed positive b. Both ambient degrees stay <=M*k.
The previously independently verified fixed constant inequality
(1-e)Gamma>1.3 remains strict and parameter independent.

The true complete root remainder is

    R(L)=1/n+exp(-.6L)+2sqrt L exp(-.3L).

The new independent base upper bound

    1/N12(L0)+1/T8(.6L0)+8192/T8(.3L0)<1/200

includes 8192=2sqrt(L0). Its corresponding tail envelopes have
powers L^-12, L^-8 and L^(-15/2), respectively. Thus R(L)<1/200
on the whole tail. Equivalently, the producer's S16 scalar envelopes
decrease because yS16'-(1/2)S16>0 for y>=2 and reciprocal prefixes
decrease. The domain .3L0>=2 is explicitly covered; 8192 is never
used as a uniform bound for 2sqrt L.

With log(5/2)<23/25 and 1.3-23/25-3/8=1/200, this proves the
strict full F bound, including positive A, cross and terminal self
terms, rather than just the deficit coefficient.

## 7. Every finite history, genuine deficit, general L1 and all scales

All numerical side conditions in the independently reviewed R081
finite-history bridge are now valid uniformly in n on this explicit
domain. In particular beta<2, t<n, rhi<4 and rho>4rhi>5. For
EVERY a_alg>=8n^4, the original cap factor ccap>1; any genuine
target item meeting A has radius<5 and any included source increment
has radius>5. The latter precedes that SAME target in the original
selection history.

The displacement inequality is verified algebraically for the ACTUAL
n, rather than deduced from a base substitution:

    n^-3-(beta+2)/(a_alg+1)
      > n^-3-4/(8n^4)
      =(n-1/2)/n^4>0,

because n>=N12>=300. This applies at every actual n on the tail.
Included pairs capture the true target in its original eligibility
test; bounded weights leave its actual row<=tau. No original item,
weight or selected history is replaced by a probabilistic model.

The finite envelope equals c_a*lambda on A times P, with c_a>=1.
The disjoint support sets make the two missing-mass terms cost at
most c_a*M*k times the TOTAL full-cap deficit. Energy, degrees and
target mass scale with c_a^2,c_a,c_a. The same finite-a inequality
therefore gives

    Delta_a >=[c_a^2 E-tau*c_a*lambda(A)]/(c_a*M*k)
       >=lambda(A)*b/M>=Gamma*n>1.3n

for every original finite history and all such finite dilations. No
gamma<1 restriction is required because the positive root was proved
directly. The complete finite cover factor and terminal self-weight
are retained; only their uniformly valid numeric a_alg limit is taken,
yielding the displayed F function.

The original general-L1 extension is unchanged: positive true
witnesses cover every strict level set with their open shrunk cubes,
compact subsets have finite subcovers, selected dilates remain true
witnesses, and every positive terminal radius is admissible under
normalization. Inner regularity and level homogeneity recover all
L1 inputs, levels and positive scales. There is no compact-support,
scale-cutoff, scale-lattice, actual-increment iid or selector-limit
assumption. These mathematical bridges were explicitly accepted in
the R081 full independent review; every eventual numerical condition
they need is covered above.

## 8. Final acceptance and scope

ACCEPT R082's explicit threshold log n>=2^24, with the precise
tau-tail explanation correction in ADDENDUM.md. The old 54-gate
PASS is finite arithmetic evidence, not the reason the dimension
tail or operator theorem holds. This review supplies the monotonicity,
uniform-domain and finite-history reasoning, with additional bounded
independent exact arithmetic. No claim of minimal threshold is made.

The complete new function is also strictly below accepted Phi80 on
this domain: the R081 full-function lower bound for Phi80 and the
strict 3/8 separator were already accepted, and 2^24>=4096 puts
the entire new domain inside the actual common domain. Other branches
retain their own validity domains and may still enter a minimum.

The original R081 eventual identity and the producer R082 execution
identity remain separate and unchanged. Source hashes are recorded
and rechecked before handoff; no original document is silently fixed.
