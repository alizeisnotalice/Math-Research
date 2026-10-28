# R071 full-chain independent review

Disposition: full-chain ACCEPT. The new J2 facts needed here were also
independently reconstructed directly from the frozen author text. This
conclusion is not conditional on a variance lower bound or on a parameter
model. The final external local-source binding is recorded below.

## Source and scope

Author source read in full:

`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round071/second_profile_cantelli/PROOF.md`

SHA256 `0e36488654bfa1413c19335a3eb9cf7b900e16b25cf9b73b2b8f36efc7381a10`.

This review accepts, for every integer n with L=log n>=8192, the genuine
original-operator statement

    C_n <= F(A_n(t), B_n(t)-(13/100)n)
        < n log n-n log log n-(7/200)n,
    t=n/((91/100)L).

The quantifiers include general L1 inputs and every positive terminal
witness radius. The leading coefficient remains one. No global
optimality, new lower bound, normal approximation, or statement about a
different maximal operator follows. Neither a successful tuple check nor
an independently checked profile alone establishes this conclusion.

Only this private audit directory was written. No author/main source,
Git state, accepted R068 source, MATLAB session, production search or n=7
calculation was changed or run.

## 1. Wide matching interval and the correct budgets

The wide h=1/(256L) is paid explicitly. At L0=8192 the trim is 4h;
thereafter h decreases. Thus q_-/q_+ and q_b=q_+ stay in the source
band, and active R=sq-h lies in [r-2h,r] with R>1. Inactive edge
sources may have R<1, but have zero final column and need no local R>=1
law. The source-radius and epsilon gates give rho>4rhi, kappa>1,
epsilon q<=2h/3, and consequently g-R/3<=0. This uses the relative
g<=R/3 domain; it does not import the old absolute d<=1/3 domain.

The source-count calibration keeps the sign ellpin'<0. Endpoint means
and the monotone ratio Nsum/(2n ellpin) give the whole-interval bounds
1+-eN once 2h+1/n<=3eN/4. With eN=2^-19, 2h=eN/2 at the base;
the old tolerance would have failed. The exact captured source-face
coefficient is product(J_i) sum_i N_i/J_i /(2ellpin).

The relative lemma's rigorous 3/rho log derivative, over matching width
4h/3, costs at most 6Lh=3/128. It would be invalid to keep the earlier
nearly-unit product factor, or to use the old coarse 256/rho derivative
with this h. Baseline log error eB and exp(-u)>=1-u justify exactly
P0=1-eB-3/128 <= product J <=1/P0. The pointwise J_i bounds separately
justify its denominators and produce the displayed Ulo/Uhi. No union
bound over individual J_i is needed.

## 2. New J2 identity and target law

Direct integration of phi^2/ell on the two tails gives

    J2=R[Qg^3+Qv^3]/(3rho).

The central removed segment has phi=0, so this remains valid when
|w|>R. The positive-part joins are locally Lipschitz; the derivative
bound integrates by absolute continuity despite onset kinks. On the
active domain Q<=1/5, |Q_R|<=1/(5R), Q_u<=3/(5R). With R'=3/2,
g'<1 and |v'|<=3/2, these imply

    |(Qg^3+Qv^3)'| <= 63/(250R).

Together with (R/(3rho))'=h/(3q rho), this proves |J2'|<.1/rho;
the author's weaker 1/rho is safe. The imported first-profile derivative
<.5/rho gives transport <=xh; the weaker second-profile derivative
gives <=2xh. At the base these are below 10^-5 x and 10^-6 x,
respectively, and improve thereafter.

At q_b, R=r. Fubini with the uniform TARGET coordinate integrates
capture probability ell_R, which cancels the inverse ell_R and yields
the full-SOURCE first and second phi moments in equation (8). It does
not identify those two probability laws. The coefficient c2 follows by
integrating log^2, with its stated sign. The baseline displacement
a-5/6 in [0,h/(3R)] and derivative -log^2 a >=-1/25 give the lower
coefficient c2-(c2+1/25)h. The free-coordinate count n-1 and a separately
nonnegative deterministic target pin give equation (9); the upper pin
is safely bounded by 1/rho. Independent alternating log8/9 enclosures
verify the resulting whole-face mean corridor [.005275,.005295]x.

The six target conditions are events under each signed ambient target
face, not hypotheses on actual increments. Their dependence on one
another is immaterial. The two-sided J2 failure is contained in ONE
event |S-ES|>=x/200000, so one Chebyshev estimate pays both sides.
The six bad rates are each <=epsA/8; the complement fraction on a
positive pinned face uses (3/4)^(n-1), without an erroneous half factor.
Negative pinned faces are already in the cap. Thus their union costs
less than epsA and the annular mass exceeds mL n.

## 3. Source laws, Cantelli and clipping

For every fixed captured signed SOURCE face, its pin contributes a
deterministic Qpin<.183. Its free coordinates are independent in the
full inverse-capture law, with exact means J1_i/J_i and second moments
J2_i/J_i. Hence the same target's whole-matching profiles give

    EV <= .183 + .04422 x/(1-eJ),
    Var V <= sum_(i != pin) J2_i/J_i
          <= .005301 x/(1-eJ) < .00531 x.

Nonnegativity permits adding the omitted source pin back only in this
upper bound. No target-law pin value is substituted for an arbitrary
source pin, and no unproved variance lower bound is used. The lower half
of the J2 profile is true but dispensable here.

The threshold gap is >x/1330. Consequently the per-face Cantelli tail
is <=v/(v+u^2)<=14/25, using .00531*1330^2/x<=14/11.
The proof applies this within each face and then mixes probabilities;
it never mixes variances or assumes equal face means. It retains at
least 11/25 of the inverse row. Since bsoft/(theta/G) equals one on
G>=theta and lies in [0,1] elsewhere, these probability estimates imply
the stated soft-row lower bound.

## 4. Radial row/column and cap removal

Both weights and the kernel are defined on the FULL source band before
its cap restriction. This preserves the product probability laws.
Radial calibration cancels dr/r against arad=r/rhi exactly. The only
radial-width correction is atanh(h/r)/(h/r), controlled by
1/(1-(h/r)^2). Thus the row normalizer really is k=nh theta/rhi.

For each source, all target signed faces are paid using
T_i/(2ell_i(r))<=1. Since r>=R on matching, product ell_i(r)<=G(z).
Multiplication by bsoft leaves at most theta; the matched radial width
is <=2h. Therefore the full ambient column is <=2k. M>2 controls both
degrees even after restrictions. No target face or source-pin factor
is discarded.

Only after these full-band calculations is the original source-cap
complement removed. Its UNWEIGHTED mass is bounded by
Csrc n exp(-(1-1/n)eta x); bsoft and arad are both <=1. Dividing by k
retains Csrc*rhi/h. The positive finite-n correction is bounded by two,
not silently dropped. Both cap and tau exponent rates equal 81/20000.
The single positive tau is sqrt(A_n(B_n-gamma n)); B_n>=nL/4 verifies
that its second factor is positive. The proof pays cap loss eS k and
the SAME tau eT k separately in b.

## 5. Entire dimension tail, not a sampled threshold

New local computation used the weaker independent n envelope L^32/32!
and exponential prefix P40, not the author's N64/S64 program. It
verified all delicate probability margins, cap and tau base losses and
the complete scalar buffer. For P40 as well as S64,

    yP'-P >= -1+y^2/2,
    yP'-(3/2)P >= -3/2-y/2+y^2/4.

The omitted terms are nonnegative. At the present y>4 both are
positive, so L/P(cL) and L^(3/2)/P(cL) decrease throughout the tail.
The rational 91 pays sqrt(8192) only at the base; the proof then uses
the actual L^(3/2), not the incorrect global bound sqrt L<=91.

The other errors are negative powers after substituting an n lower
envelope. The J2 lower mean bound improves as h and 1/N decrease, while
its upper bound decreases. The Cantelli gap rate increases and its
tail ratio decreases. Source lower radii increase; epsilon q/h
decreases; 6Lh is constant. Thus every finite condition extends from
8192 to all larger L. There is no residual unquantified eventual guard.

## 6. Same finite history and operator conclusion

The target and source point supports are disjoint. Their physical
item radii satisfy r_j<3 and r_i>=rho>4, so a source precedes the SAME
target item in the actual largest-eligible-radius history. For
algorithmic a>=8n^4, x>1 gives t<n, beta<2, ccap>1 and
zeta-(beta+2)/(a+1)>0. Hence every kernel-included increment piece
really enters that target's eligibility test. A cross-band item is
charged only for its included point-measure; the proof does not copy
whole item weights. With both weights <=1 the actual weighted row is
<=arad(w)tau<=tau.

On D the finite envelope is exactly lambda_a=ca lambda, ca=alpha^-n>=1,
and D is contained in the actual finite cap. The energy scales ca^2,
the degrees ca, and the target mass ca. The disjoint-support deficit
identity therefore gives

    Delta_a >= (ca^2 E-tau ca lambda(A))/(ca M k)
             >= (E-tau lambda(A))/(M k)
             >= Gamma n > gamma n.

This is uniform in the actual history and in algorithmic a. General
L1 input, strict-level positive witnesses and compact finite subcovers
give the original cover/load inequality; terminal-radius normalization
works for every positive radius. Inner regularity and only a scalar
a->infinity limit then give C_n<=F(A_n,B_n-gamma n). No ambient
probability law or compactness limit of selectors is imposed on mu_a.

Finally B_n<=T-n log(.91)+1 and A_n<=n exp(-.09L). The complete root
cost includes the positive A and cross terms. A new positive-series
log-price bound and P40 remainder verify

    -log(.91)<189/2000,
    1/n+exp(-.09L)+2sqrt L exp(-.045L)<1/2000.

Since .13-.0945-.035=.0005, strict slack proves the stated final
root-cost inequality. This does not compare only a displayed gamma.

## Actual new computation and reproducibility

* preflight.py: twelve new Fraction checks, actual PID 99742, PASS.
* critical_tail.py: twelve additional compressed groups, actual PID 219,
  PASS; N32, alternating log8/9, positive log-price8 plus a geometric
  tail and full positive exponential P40. No author helper was imported.

The author 43-gate verifier and older passed programs were not rerun.
The analytic arguments above, rather than the finite checks alone,
establish the infinite dimension tail and the geometry/probability law
transfers. Receipts retain each execution's actual source identity.

Local imported accepted identities: R064 relative-coordinate lemma
`c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`
and R066 first-profile lemma
`67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030`.
Their active-domain hypotheses were checked anew, including the much
wider h; their previous tuple thresholds were not imported.

## Final external local-source binding

After the preceding full reconstruction, the R070 local source was
actually read in full at

`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/upper/round070/second_profile/PROOF.md`

SHA256 `4d2f4fc2d0f3353ba82e32487dba68804d3021d4c0d477137c11f6c1a004e2b8`.

Its closed form, finite-h baseline mean, derivative/AC bridge and
per-source-face second-moment upper bound agree with the independently
reconstructed R071 chain. At the present base h=2^-21<2^-20, n>=2000
and x/n<1/1000, all its relevant guards hold. R071 uses weaker profile
probability and derivative bounds, but explicitly pays them at 8192.
This audit does not import or need R070's source variance lower bound,
its proposed R072 application, or any mixture variance conclusion.
The author's R071 source remained frozen at 0e364886 throughout this
binding; none of the passed programs was rerun for this metadata step.
