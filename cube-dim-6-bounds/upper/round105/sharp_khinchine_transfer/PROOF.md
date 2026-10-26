# 105 prepared source-qualified sharp Khintchine and angular transfer

**Status:** theory submission and unrun exact script, awaiting primary static release
and later independent review. This is105's new upper-theory work, not an audit of104
or108. No kernel, native, old grid,104/100 program or new Cn result is claimed.

## 1. Public mathematical input and its qualification

Use Xinyuan Xie, *Stability of Szarek's Inequality with Best Constant*,
[arXiv:2610.07526v1](https://arxiv.org/abs/2610.07526v1),2026-10-05 preprint,
not a journal article. The supplied original8-page PDF has SHA-256
5731a47fcbd91d116778df9fa0033353a82f9e24db3672e6cdb48f6552bc31d4;
the original TeX is794b7efaead1ea2288baac88ef0d6c7fd87812e9225abd001439f155bd41f074.
The complete TeX mathematical proof was read. The PDF's bytes/page extraction are
bound by the main acquisition receipt; no new independent visual8-page read is
claimed here. No conversation or LLM operational text is a mathematical input.

The external input used is Theorem1: for a real unit vector a sorted in descending
nonnegative coordinates,

    E|sum a_j eps_j| >=1/sqrt2+c_*||a-(e1+e2)/sqrt2||,
    c_*=(2-sqrt2)^(3/2)/8.

The project's use does not require the optimality claim about c_*. The source's
proof was checked beyond its abstract: its large a1+a2 branch invokes FW26
Theorem1; the smaller branch uses its proved Theorem2, convexity of F(1/s), the
cosine-tail lemma and a special gamma-ratio bound attributed to Kershaw. The
case split, positive/negative Fourier contributions, Gaussian-tail absorption,
concavity endpoint reduction, and zero-coefficient cases have been analytically
checked with no defect found in the used chain. This remains a source-qualified
public theorem input, not a claim that the project independently established the
optimal stability theorem or read every original source in its bibliography.

[Fang--Wang arXiv:2608.27908v1, Theorem1 and proof](https://arxiv.org/html/2608.27908v1)
were read from the official HTML. Its proof conditions on the first two signs,
uses the classical sharp Khintchine inequality on the remaining signed sum, and
then a quadratic identity. This matches the branch quoted by Xie. The whole FW26
paper was not independently reviewed. The saved HTML is3317d13947b2b3e69f63ef7f5f5edc01476450eca2cc44cbc445123d7ca7845b.

The classical Szarek/Haagerup original PDFs were not read in full; main acquisition
records403/truncation. They are explicit external inputs in that proof chain,
not silently upgraded to independently read original proofs. The official AMS
Kershaw page also returned403. AppendixA supplies an independent proof of the
specific gamma-ratio estimate needed here, without claiming to have read the
1983 paper. External theorem acceptance remains a final-review qualification.

## 2. Same full body32 symbol and the equality directions

The model is unchanged:32 unoriented columns (1,±1,...,±1), unit weights,
alpha1, common time, Fourier convention(2pi)^-6, Lebesgue measure. Its full symbol
is psi(xi)=32 E|sum_(j=1)^6 xi_j eps_j|. It is invariant under physical-coordinate
signs and permutations; this statement concerns the full symbol, not an omitted
column law.

Let

    E={ (s e_i+t e_j)/sqrt2 :1<=i<j<=6, s,t in{-1,1} }.

These are60 distinct unit vectors, including antipodes. Sorting |theta| for a
unit theta minimizes distance to E by choosing the two largest absolute entries
and matching their signs. Hence d(theta)=min_(u in E)||theta-u|| is exactly the
distance in Xie's proper-chamber statement. Its theorem transfers to

    psi(theta)>=c0+32c_*d(theta), c0=16sqrt2.

Dropping the nonnegative deficit also gives the classical sharp coercivity
psi(xi)>=16sqrt2|xi|. The former098/104 moment coefficient is c_old=32/sqrt3.
The difference is in a genuine full-frame input, not a pointwise radial
heat-kernel assertion.

## 3. Uniform chordal-cap area on every three-dimensional subspace

Fix ANY three-dimensional linear subspace S in R6 and its unit sphere S_S^2.
For u in E put w=P_Su, r=|w|<=1 and a=1-delta^2/2, with delta=1/8. The condition
||theta-u||<delta, theta in S_S^2, is <theta,w>>a. If r<=a, the cap is empty up
to a possible point. If r>a, its axis is w/r and the condition becomes

    <theta,w/r>>a/r.

The spherical-cap area fraction is(1-a/r)/2<=(1-a)/2=delta^2/4=1/256.
Thus the estimate is uniform even if u is outside S; projection cannot enlarge
the cap. It does not assume that S contains any equality direction.

The union of60 such caps has fraction at most60/256=15/64, with no disjointness
assumption. Its complement has fraction at least49/64 and d(theta)>=1/8.
Boundary circles/points have surface measure zero; using >= on the complement
also makes the estimate valid without discarding them.

Two finite exact square comparisons certify the following coarse constants:

    c_*>1/20,
    32c_*(1/8)/c0 =sqrt2*c_*/8 >sqrt2/160 >1/114.

For the first, c_*^2=(20-14sqrt2)/64 and sqrt2<248/175 follows from
248^2>2*175^2;20-14*(248/175)=4/25. For the second,
2*114^2>160^2. All relevant quantities are positive, so squaring preserves
the intended directions. No numerical estimate of a distance or cap is used.

Therefore psi>=c0 everywhere on S_S^2 and psi>c0*(115/114) on its complement
of fraction at least49/64. Positivity and the worst allowed union area give

    average_(S_S^2) psi(theta)^-3 <= eta*c0^-3,
    eta=15/64+(49/64)*(114/115)^3<1.

The excess is an angular-average improvement; it is NOT a new pointwise
coercivity constant larger than16sqrt2, which equality directions would forbid.

## 4. Actual three-dimensional quotient peaks and exact coefficients

Let Q3(S) be the time-one peak of the actual full body32 law projected onto S.
Coercivity and ordinary positive radial Fourier integration give

    Q3(S)=(1/(4pi^3))int_(S_S^2) psi(theta)^-3 dS
          <=eta/[pi^2(16sqrt2)^3].

This holds for all3D S, in particular every actual projection span(i,j,e)^perp
that occurs after restoration in104. Neither an omitted-frame law nor a
selected body-triple peak is substituted for this actual quotient.

The sharp-only alternative has eta replaced by1. Relative to the moment-bound
Q3 coefficient, the common linear factors are

    r_K=(2/3)^(3/2)=sqrt(8/27), r_105=eta*r_K.

For the098 principal branch these give K_K=sqrt2/(192pi^4), K_105=eta*K_K.
For104's Sigma=180/sqrt20+450/sqrt32+300/6 and
D_old=3sqrt6 Sigma/(64pi^4), they give

    D_K=Sigma/(16pi^4), D_105=eta*Sigma/(16pi^4).

No Sigma, peak, source convolution, geometry, time integral or spatial grid is
recomputed. These are exact coefficient identities. Every fixed104 P and B norm
is linear in that same uniform Q3 upper coefficient, so a common factor r scales
all positive norm budgets by r and the associated I allocations by r^(6/5).
For r_K this I factor is(2/3)^(9/5); for105 it is eta^(6/5)(2/3)^(9/5).
This is a relation between upper-budget expressions. Their numerical lower
endpoints would not be lower bounds on actual G or actual I.

The104 data are completed saved dyadics, but at preparation no final104 manifest
was supplied and its proposal text still says prepared.105 does not certify
104's independent proof acceptance; any use as an actual paid budget remains
conditional on that separate acceptance. No display decimal is an input.

## 5. Pre-run steering: conditional108 dominance and no duplicate104 recalculation

The primary has proposed a separate108 mechanism. Its proof is being handled
elsewhere and is NOT assumed proved or audited in105. For ONE conditional
algebra comparison only, suppose its claimed uniform peak bound is valid:

    Q3(S)<=C_BL:=1/(4096pi^3).

Relative to the moment coefficient this would have factor
r_BL=8/(3sqrt3*pi). Compared with105's cap coefficient,

    C_BL / C_105 =2sqrt2/(pi*eta).

The fixed rational eta satisfies eta>31/32. With pi>3 and
93^2>8*32^2, one has pi*eta>93/32>2sqrt2, so C_BL<C_105<C_K.
The registered finite exact program will certify this inequality and a
separated coefficient interval. It will not prove108's premise.

All original104 budget expressions are strictly increasing in their positive
common Q3 coefficient. Consequently, IF the108 premise is later accepted,
its norm and I allocations are strictly smaller than ALL planned105 sharp-only
and cap allocations on exactly the same regions and same positive inputs.
No extra integration or old-program run is needed to establish that conditional
dominance. The105 plan is therefore amended BEFORE execution: certify the
source-qualified backup coefficient and this one conditional comparison, and
omit full104 budget recomputation if the comparison passes. There is no scan,
parameter adjustment, or optimization of caps.

## 6. Operator scope and remaining gaps

The genuine098 time-first identity and104 restoration use the same full symbol;
this transfer changes only their uniform Q3 upper coefficient. Its effect stays
inside the original body32 Cauchy model. It does not alter the common radius,
Lebesgue measure, general-L1 input class, or all-positive-scale quantifiers.

This is neither a global g_L certificate nor complete same-model face I_U.
Far and other unpaid spatial allocations remain unpaid. Conditional plan S or U
arithmetic could only be interpreted after those pieces and global g are paid;
no remaining positive core allowance is a complete positive S. There is0 new
Cn result, and no C6 lower-bound work. The105 backup theorem/source acceptance
must undergo later independent review.

## AppendixA. Independent proof of the needed gamma-ratio special case

This proves the exact bound attributed to Kershaw in Xie's cosine-tail proof;
it is not a claim to have read Kershaw's original paper. For x>0 set

    R(x)=Gamma(x+1/2)/Gamma(x+1), B(x)=(x+1/4)R(x)^2.

The convergent Euler-product/digamma difference, or its positive Laplace
integral, gives

    digamma(x+1)-digamma(x+1/2)
      =(1/2)int_0^infinity exp(-(x+1/4)t)sech(t/4)dt.

Thus (log B)'=1/(x+1/4)-int_0^infinity exp(-(x+1/4)t)sech(t/4)dt>0.
The strict inequality uses0<sech(t/4)<1 for t>0. Xie's elementary cosine
integral and dominated convergence give F(s)->sqrt(2/pi), and
F(2x)=sqrt(2/pi)*sqrt(x)R(x). Therefore B(x)->1 as x->infinity.
Monotonicity yields B(x)<1 and hence

    F(s)<=sqrt(2/pi)*sqrt(s/(s+1/2)), s>0.

This is exactly the special estimate used on F(s-1/2) for s>1/2. It uses only
ordinary convergent gamma-product differences and the stated elementary limit,
so no uniform endpoint, sign, or unread general gamma inequality is hidden.
