# A027: original body32 positive-part mechanism, restored exactly

## Definitions and normalization

The model is V={(1,s2,...,s6):sj=+/-1}, with32 independent Cauchy variables
of the SAME scale t and all weights1. The symbol is psi=sum_v|v.xi|,
L=-psi(D), and p_t is the positive probability density with Fourier
transform exp(-t psi). No axis or pair directions occur in this audit.

Let G(x)=integral_0^infinity p_t(x)dt, an extended positive function,
and g=essinf_{[-1,1]^6}G. Let q=6/5 and

    I=integral_[0,1]^5 G(1,y2,...,y6)^q dy2...dy6,
    B=I^(5/6)=||G(1,.)||_{L^(6/5)([0,1]^5)}.

I is a five-dimensional Lebesgue integral over the entire positive face
of volume1. It is NOT an integral over the whole boundary, an ordered
chamber average, or a time integral. With ordered free coordinates it is
120 times the ordered-chamber integral. The whole cube boundary integral
J of G^q is384 I, accounting for6 coordinates,2 signs and32 face orthants.
The sign/permutation multiplicities in the kernel have already been paid.

Use R_c=(G-c)_+ for the SPATIAL positive part, and H_epsilon=integral_epsilon^infinity
p_t dt for the TIME-truncated potential. These are different functions.
G=H_epsilon+U_epsilon by positive Tonelli, where U_epsilon is the short-time
integral. In particular short-potential failures at a few fixed times do
not imply failures of the R_c family.

## Source theorem and all-input/all-scale quantifiers

The original complete finite-source text is
`Math-Research-finite-dim/cube-finite-dim-bounds/audit/second_review/theory/paper_finite_source.md`,
Sections A–D, pinned in scalar_plan.json. It states its applications for
n=3,4; the n=6 adaptation actually rechecks the required steps in
upper/theory/001_n6_body_frame_and_singular_strata.md, including its
expanded proof, not merely its summary. The relevant requirements are:

1. A true positive, conservative, symmetric semigroup; a spanning basis
   gives p_t<=Ct^-6 and psi comparable to Euclidean norm. This yields
   Dom_2 L=H1, and n>1 gives the integrable large-time bound.
2. For bounded compactly supported f>=0, the discrete obstacle construction
   yields f=sigma-Lu, u>=0 in L1 intersect L-infinity intersect H1,
   0<=sigma<=kappa, sigma=kappa on D={u>0}, integral sigma=integral f.
   The large-cutoff bound ||L chi_R||infinity=O(R^-1) recovers mass.
3. Weak-L^(6/5) control from short/long time splitting makes G locally L1
   and R_c globally L1, including all singular strata. Semigroup Jensen
   gives L R_c=mu_c-delta_0 on the WHOLE space, with mu_c a positive
   probability measure. The small-cutoff L1 scaling epsilon^5 removes
   a possible positive atom at zero. This never deletes S2.
4. For a certified g0<=g and 0<c<g0 set

       K_c=R_c/[64(g0-c)].

   Then K_c>=1_Q/64 a.e. The Fourier finite-source identity gives, for
   K_{c,r}(x)=r^-6 K_c(x/r),

       K_{c,r}*Lu=[mu_{c,r}*u-u]/[64(g0-c)r].

   On the same D-complement it is nonnegative for every fixed r. Rational
   r give one common null set; L1 dilation continuity and bounded f extend
   to all positive real r. Taking kappa=lambda/integral K_c controls the
   level set by D. Finally min(|f|,j)1_[-j,j]^6 increases to |f|; positive
   convolution and commuting suprema extend the inequality to general L1.

Thus the original mechanism does have the required all-L1/all-scale
interface. This audit checks the proof chain and its dimensional inputs;
it does not claim a new obstacle implementation or a fresh full formal
verification of every earlier theorem.

## Exact cost and optimal truncation

G has degree -5. Cube-radial coordinates x=ry have Jacobian r^5 dr dS(y).
For a fixed face value a=G(y), positive radial integration gives

    integral_0^infinity (a r^-5-c)_+ r^5 dr
        =(5/6)c^-1/5 a^6/5.

Tonelli and J=384I yield integral R_c=320 c^-1/5 I. Therefore

    integral K_c=5 I c^-1/5/(g0-c).

Its logarithmic derivative is -1/(5c)+1/(g0-c), so the minimum is exactly
c=g0/6. The resulting sufficient bound is

    C6 <= 6^(6/5) I0/g0^(6/5) = (6 B0/g0)^(6/5),

provided both global certificates I<=I0 and g>=g0 hold for this same G.
For the best ideal positive-part family using the true g and I, write
C_model=(6B/g)^(6/5). This is a family cost, not a lower bound on C6.
Competing with64 is exactly equivalent to

    B < (16/3)g,  or  I < ((16/3)g)^(6/5).             (C)

## What the new exact point excludes

B023 and independent F023 agree on an exact rational Gamma at
x0=(50,50,50,45,27,0)/50, with G(x0)=Gamma/(128 pi^5).
Our Arb512 postprocessing encloses it near1.75407245845747e-5 and proves

    1/60000 < G(x0) < 1/55000 < 1/50000.

G is continuous near this point (it is off S2). The cube contains
positive-volume neighborhoods approaching it from inside. Hence

    g <= G(x0),

not g>=G(x0). In particular the old g>=1/50000 goal is strictly false;
the goal g>=1/60000 is merely not contradicted by this point.
Any competitive ideal model must satisfy the exact necessary condition

    B < Gamma/(24 pi^5),
    I < [Gamma/(24 pi^5)]^(6/5).                       (N)

These are competition ceilings, not upper bounds on actual B or I.
To exclude the family rigorously, an actual lower bound B_L>=Gamma/(24pi^5)
(or the equivalent I_L) would suffice. No audited existing result supplies
it. An upper-envelope value larger than this number does not supply it.

## Explicit surviving certificate targets

Choose g0=1/60000. Then (C) asks for

    B<1/11250,
    I<(1/11250)^(6/5),

the latter threshold enclosed near1.37599537155658e-5.
The already proved whole-face short-time estimate at epsilon=1/1000 is

    ||U_epsilon||q <= t0,
    t0=(6/7)243^-1/6 pi^-11/6 (1/1000)^7/6,

near1.33048729415938e-5. It includes all S2 neighborhoods and intersections.
Thus a less demanding alternative to the old H-norm goal1/16000 is

    g>=1/60000,      ||H_(1/1000)(1,.)||q<=1/14000.   (T)

Minkowski would then give B<=1/14000+t0, I<=that quantity^q, and a cost
enclosed near60.4266939170373. The corresponding I goal is near
1.29916954919234e-5. For this allocated B budget alone the required g is
strictly greater than(3/16)(1/14000+t0), near1.58875208194060e-5.
All quantities were evaluated with outward Arb balls in results.json.
Neither global premise in (T) has been proved.
