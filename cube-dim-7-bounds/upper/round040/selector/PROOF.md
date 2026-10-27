# R040 — a fresh-witness-mass selector and its causal-spill interface

## Outcome and scope

This is a different selection mechanism, rather than a new bound on the old
pointwise load. It freezes disjoint portions of the original input mass, not
small-cube volume increments. The resulting finite-to-operator theorem is
proved below. It removes the R037 common-spike saturation mechanism on that
diagnostic, but no improved general-dimensional upper bound is proved.

The precise remaining problem is a uniform estimate for the one causal-spill
operator in Section 2. Section 5 gives an all-dimensional, actual-history
counterexample to closing that problem by assigning each rejection to one
source having a fixed share of its witness mass. Arbitrarily many earlier
sources can be needed simultaneously, even with strictly decreasing source
radii and all original witness tests valid.

Scope throughout: Lebesgue measure on R^n; centered axis-parallel cubes with one
common coordinate radius; arbitrary L1 input; all positive radii. This is not
lower-bound work or a special-dimensional optimization. All new files are in
`upper/round040/selector/`. No shared MATLAB session or Git synchronization was
used. The Python checks use exact rational arithmetic, not floating evidence.

## 1. The selector

Write V_i=Q(y_i,R_i) for a true witness at level lambda>0:

    integral_{V_i} f > lambda |V_i|,       f>=0, f in L1.

Fix a>1, put Q_i=Q(y_i,R_i/a), and process the finite input family in
nonincreasing R_i order, with a fixed tie order. Fix 1/2<theta<1. Maintain H,
the union of previously selected witness cubes, initially empty. Select V_i
if and only if

    integral_{V_i minus H} f >= theta lambda |V_i|.              (1)

For its selected index j, freeze

    F_j=V_j minus union_{ell<j} V_ell,
    g_j=f 1_{F_j},      m_j=integral g_j,

then add V_j to H. Rejected witnesses are never added to H. Thus the F_j are
actual, pairwise disjoint measurable mass-ownership regions, and

    m_j>=theta lambda |V_j|,
    sum_j m_j<=||f||_1,
    sum_j |Q_j|<=||f||_1/(theta lambda a^n).                    (2)

At a rejection, with H_pre the then-current selected union, (1) fails while
the original witness is strict. Consequently

    integral_{V_i intersect H_pre} f
      >(1-theta)lambda |V_i|.                                 (3)

Every source in H_pre has radius at least R_i. This is a simultaneous
prefix-mass test. There are no coefficients w_j, no old eligibility test,
and no universal pointwise weighted-load estimate in this argument.

## 2. One explicit operator interface

An admissible source datum at level lambda is a finite list

    V_j=Q(y_j,R_j), R_1>=...>=R_N>0,
    H=union_j V_j,       g>=0 in L1, supported in H,
    F_j=V_j minus union_{ell<j}V_ell,       g_j=g1_{F_j},

satisfying, for every j,

    m_j=integral g_j>=theta lambda |V_j|,
    integral_{V_j}g>lambda |V_j|.                             (4)

The original selector produces exactly such data by taking g=f1_H. In
particular restriction to H does not change any selected witness integral,
since every selected V_j is contained in H. Conversely, a source-only family
satisfying (4), with input g, is selected in the given order by (1).

Put c_a=1+1/a. Define the causal-spill maximal operator

    T_{D,a}g(x)
      =sup_{R>0} [1/|Q(x,c_a R)|]
          integral_{Q(x,c_a R)} sum_{j:R_j>=R} g_j.             (5)

Every source piece has its own upper scale R_j. Smaller-scale sources cannot
pay larger-scale tests. All sources allowed at a test scale are summed
before testing; the operator does not choose one source independently.

Define the single constant needing an estimate:

    Psi_{n,a}(theta)
      =sup_D lambda |{T_{D,a}g>(1-theta)lambda/c_a^n}|/||g||_1, (6)

over all finite admissible data (4), with positive total mass. Translation,
homothety, normalization of lambda, and finite source cardinality are
unrestricted. Supremizing only a bounded range of radii would not suffice.

**Interface theorem.** For each n>=1 and 1/2<theta<1,

    C_n <= Psi_{n,a}(theta)     for every a>1,
    hence C_n<=liminf_{a to infinity} Psi_{n,a}(theta).         (7)

Thus an explicit uniform estimate Psi_{n,a}(theta)<=K_{n,a} with
liminf_a K_{n,a} below the current accepted bound would give a genuine
general-L1/all-positive-scales improvement. No such estimate is supplied
here. A finite family with a small value in (6) does not bound its supremum.

**Finite-family proof.** Let U=union_i Q_i for any original finite witness
family and apply the selector. If V_i was rejected, (3) is supplied by sources
having R_j>=R_i. For x in Q_i,

    V_i subset Q(x,c_a R_i).

The test R=R_i in (5) therefore gives

    T_{D,a}g(x)>(1-theta)lambda/c_a^n.                         (8)

If V_j was selected, its own g_j lies in V_j and is allowed at R=R_j. For
x in Q_j the same containment and (4) yield

    T_{D,a}g(x)>=theta lambda/c_a^n
                  >(1-theta)lambda/c_a^n.                    (9)

The condition theta>1/2 is used precisely here. Consequently the whole U,
including selected and rejected small cubes, lies in the one halo (8).
Equations (6), ||g||_1<=||f||_1 give

    lambda |U|<=Psi_{n,a}(theta)||f||_1.                       (10)

Measurability of (5) needs no selector limit: for a finite datum, the active
source sum is constant between the finitely many R_j. Integrals over moving
cubes are continuous in center and positive radius for L1 functions. A
countable dense collection of radii, together with one-sided approximation
at each R_j, has the same supremum. Thus the halo is measurable.

**Passage to arbitrary L1 and all scales.** For an arbitrary f in L1 replace
it by |f|. At every x in {M_n f>lambda} choose any actual positive witness
radius R_x. For a compact subset K of this level set, the open cubes
Q(x,R_x/a) cover K and have a finite subcover. Apply (10) to that finite
family, obtaining lambda|K|<=Psi||f||_1. Inner regularity gives the same
inequality for the full level set. The level set is open, since translation
of each fixed-radius cube average is continuous for L1 input. No support
truncation, lower radius cutoff, upper radius cutoff, fixed scale lattice,
or extension from numerical samples is used. This proves (7).

## 3. Baseline checks, and the precise unclosed obligation

The interface is finite. If T_{D,a}g(x)>0, some active source g_j contributes
at some R<=R_j. Since its support lies in V_j,

    x in Q(y_j,(1+c_a)R_j)=Q(y_j,(2+1/a)R_j).

Using (4), this gives the explicit crude bound

    Psi_{n,a}(theta)<=(2+1/a)^n/theta.                         (11)

For every datum, T_{D,a}g<=M_n g, so the old operator constant also gives

    Psi_{n,a}(theta)<=c_a^n C_n/(1-theta).                     (12)

Equations (7), (11) recover the familiar exponential estimate C_n<=2^n by
letting a tend to infinity and then theta increase to one. They do not
improve the accepted n log n bound. Equation (12) is a check on direction,
not a useful recurrence: its factor 1/(1-theta) is greater than one.

In particular, merely replacing the rejection in (3) by the unrestricted
maximal function of g loses a factor 1/(1-theta), while ||g||_1 can equal
||f||_1. There is no uniform mass contraction to compensate. The following
source constraints must be used together for a stronger estimate: exact
fresh regions F_j; strict selected-witness inequalities in (4); the
individual upper scales R_j; and the sum over every allowed previous
source. The current proof does not bound that coupled halo by a constant
with a smaller n log n leading coefficient.

## 4. What happens to the R037 common-spike diagnostic

Suppose a finite original witness family all contain the entire support of
f, which lies in Q(0,epsilon). Let R_max be the largest witness radius.
The first witness is selected and captures all input mass. Every later
fresh mass is zero, so there is exactly one source.

For this source, any positive value of (5) must have R<=R_max and a cube
Q(x,c_a R) meeting the support. Therefore its halo lies in
Q(0,c_a R_max+epsilon). The first witness satisfies
||f||_1>lambda(2R_max)^n. Hence its actual halo ratio is bounded by

    lambda |halo|/||f||_1
      <(c_a+epsilon/R_max)^n.                                (13)

At fixed n,R_max this approaches c_a^n as epsilon decreases to zero. This
is a fixed datum/level diagnostic. It is not a claim about the whole weak
norm of shrinking homothetic inputs. Unlike the old R037 weighted-load
history, this selector does not carry n log n worth of selected payments
at the common spike. That establishes a real change of mechanism, but
does not control other common witnesses for that geometry or the supremum
in (6).

## 5. An exact many-source rejection obstruction

Fix n>=2, N>=1 and 1/2<theta<1, and normalize lambda=1. Let the final
witness be V_*=(0,2)^n, of radius 1 and volume v_*=2^n. Put

    A=(1-theta/2)v_*,       B=(3theta/4)v_*.

For j=1,...,N set

    R_j=3-j/(2N),    delta_j=j/[2(N+1)],    delta_0=0,
    V_j=(delta_j-2R_j,delta_j) x (-R_j,R_j)^(n-1).

The R_j are strictly decreasing and greater than 1. Define pairwise
disjoint rectangles

    P_j=(delta_{j-1},delta_j) x (0,1)^(n-1),
    B_j=(delta_{j-1},delta_j) x (-1,0) x (0,1)^(n-2),
    P_*=(1,2)^n.

Put a constant density on P_j with total mass A/N, a constant density on
B_j with total mass 2|V_j|-A/N, and a constant density on P_* with total
mass B; set f=0 elsewhere. Every prescribed mass is positive and f is a
single nonnegative L1 function.

Each P_j and B_j lies inside V_j and outside every earlier V_ell, because
its first coordinate exceeds delta_{j-1}. Every earlier support rectangle
is inside V_j, and every later support rectangle is outside V_j; the
P_* rectangle is outside all V_j. Thus, in the forced nonincreasing-radius
order,

    integral_{V_j}f=sum_{ell<=j} 2|V_ell|>|V_j|,
    integral_{V_j minus H_pre}f=2|V_j|>=theta|V_j|.             (14)

All sources are genuine witnesses and are selected. Their actual frozen
pieces contribute precisely A/N each to V_*. The target has

    integral_{V_*} f=A+B=(1+theta/4)v_*>v_*,
    integral_{V_* minus H} f=B=(3theta/4)v_*<theta v_*,
    integral_{V_* intersect H} f=A>(1-theta)v_*.              (15)

It is therefore a genuine rejected witness, of strictly smaller radius.
Yet each individual earlier source supplies only 1/N of its captured
mass, or the fraction

    (1-theta/2)/[(1-theta)N]                                  (16)

of the rejection threshold. For any proposed positive fraction independent
of N, choose N sufficiently large and no source meets it. All small cubes
Q_i=V_i/a are defined by the usual concentric homothety, so this is an
actual finite input at every a>1, not an abstract allocation relaxation.

This disproves only a single-source mass-share claim. It does not disprove
the selector, the simultaneous halo estimate, or a fractional allocation
that combines sufficiently many sources. The bulk mass outside the target
is deliberately retained and fully paid; no independent worst histories
are spliced together. The example has a large original input mass and is
not presented as a difficult operator ratio.

## 6. Actual checks, provenance, and status

Run from the dedicated checkout:

    python3 cube-dim-7-bounds/upper/round040/selector/verify.py

The script uses standard-library `fractions.Fraction` on the local Python
3.15.0 runtime. It computes exact rectangular intersections and differences,
integrates the one consistent piecewise-constant f, checks every witness
and every selector decision, and checks each actual residual-source
contribution and total mass conservation. The three bounded diagnostics
are (n,N)=(2,8),(2,32),(3,8), theta=3/4. They give target mass 19/16,
fresh mass 9/16, and captured mass 5/8 times target volume, and individual
source shares 1/N. These are minimal interface checks, not dimensional
optimization or evidence for a uniform halo theorem.

`verification.json` records the actual results and SHA256 hashes of the
R032 joint-cover proof, R035 payment-compression proof, and frozen R037
proof plus clarifications. R032 changes weights inside the old selector;
R035 reassigns its scalar payments afterwards. The present selector
instead chooses on exact fresh f-mass before creating payments, and its
remaining bound is the causal halo (5), not an old W bound or an
unrestricted power-scale supremum.

Strict results: (7), the exact finite/general mapping, (11)-(13), and the
all-n>=2 obstruction (14)-(16). Open: an explicit uniform estimate of (6)
that improves the currently accepted upper bound. No new C_n improvement
is asserted. No fetch, rebase, commit, push, production scan, n=7 tuning,
or old low-dimensional project change was made.
