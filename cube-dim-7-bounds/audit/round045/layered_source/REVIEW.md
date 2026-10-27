# Independent Sol audit: R045 layered_source

**Verdict: PASS for the proposed sufficiently-large-n theorem, relative to
the accepted original finite cover/envelope framework. No blocking gap was
found. No explicit finite domain is certified by this audit.** The
coordinator retains final proof integration and acceptance responsibility.

Audited source:
`/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round045/layered_source/PROOF.md`.
Its actual SHA256 was checked as
`22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac`.
The source was not edited. This is an audit of that inverse-source-weight
proof, not of the separate R045 geometry/softcap proof or its 1D audit.

The full proposed chain was checked, including both final operator
statements. The ten numerical integration basepoints and the author's
arithmetic receipt were not rerun. A new symbolic identity check and exact
rational constant check were actually run locally; they are supporting
checks, not substitutes for the inequalities below.

## 1. Parameter domains and source normalization

For fixed source radius rho=tq, q=1-v, the definitions give exactly

    q0=(kappa-1)rho=q-zeta rho,
    R=q+d-h,
    g=R-q0=d-h+zeta rho.

Matching implies R<r<R+2h, including at either truncated endpoint of I.
Eventually h<d/16 and zeta t<h/2, so R>1, r<2,
d/2<g<d and q0>1/2. In particular kappa rho>rho>4>r.

At a fixed source coordinate s, the fraction of target coordinates in
[-R,R] captured by the kappa rho interval is exactly the proposed ell(s).
Only one endpoint can be lost; the lost tail has height g/(2R). Therefore
1-d/2<=ell<=1, and omega is exactly the reciprocal product of these
reference target capture fractions.

For an actual matched target radius r>=R, its coordinate fractions
ell_i(r)=min(1,(1+(kappa rho-|z_i|)/r)/2) are nonincreasing in r and at
least 1/2. The complete target-face sum is bounded by twice n times their
product: each deleted-coordinate product is at most twice the full
product, while the two signs and the face factor 1/2 are both included.
The matched logarithmic target interval is at most 2h/R. Consequently

    C0(z)<=omega(z) (4nh/R) product_i ell_i(R)<=4M.

The omega cancellation here is exact, pointwise in the source, and occurs
before any high-Phi truncation. No density or concentration of an actual
history measure is used in this column estimate.

## 2. The coordinate row law and its two-radius moment

For |w|<=q0 the captured source interval is the whole source interval,
so J(w)=F=E[1/ell(Z)]. Beyond q0, increasing |w| removes a one-sided
interval of a positive symmetric density. Thus J is even and
nonincreasing in |w|, as required by the proof.

The exact free-coordinate loss is

    E(1-ell(Z))=g^2/(4R rho).

The pointwise inequality 1/ell-1<=2(1-ell) therefore gives the stated
bound on F-1. The one-sided excluded source interval has length at most
g+2h<=3d and inverse weight at most two. These establish the two-sided
pointwise bounds (9); neither bound ignores a source face.

The reference mean E_{[-R,R]}J=1 follows directly by Fubini and inverse
normalization. It was also independently verified through the exact
formulas

    F-1=[2R log(2R/(2R-g))-g]/rho,
    E(b)=(R/rho)log[(2R-g+b)/(2R-g)],  0<=b<=g,

where E(b) is the excluded one-sided mass when |w|-q0=b.
The identity R^-1 integral_0^g E(b)db=F-1 gives reference mean one.
The symbolic check verifies these primitives by derivatives and boundary
values, avoiding a numerical quadrature or an unchecked logarithm branch.

For a free actual target coordinate, the expansion from R to r is exactly

    E_r J=1+(1/r)integral_R^r [J(w)-1]dw.

Since J(R) is at most its reference mean, the extra integral is
nonpositive; its magnitude is at most 6hd/rho. The probability of the
nonconstant boundary region is at most 3d, and the pointwise deviation
bounds give

    E_r(J-1)^2<=4d^4/rho^2+27d^3/rho^2<=31d^3/rho^2.

For two different matched source radii, the same W is indeed shared.
Writing J=1+delta, both delta means are nonpositive. Cauchy--Schwarz for
their product, using the last displayed bounds at the same actual target
radius, gives

    E_r[J_rho J_rho']<=1+31d^3/(rho rho').

No independence between source radii is assumed. Both of their reference
radii lie below that same actual r by matching, so the two nonpositive-
mean bounds apply simultaneously.

## 3. All pinned source faces and the row second moment

At pinned source coordinate +rho capture is equivalent to w_j>-q0;
at -rho it is equivalent to w_j<q0, since the opposite distant endpoint
lies beyond the whole target cube. Thus, up to null boundaries,

    N_j(w)=1+1_{|w_j|<q0}.

The inverse pinned-coordinate factor is 1/ell_* with
ell_*=1-g/(2R). The complete source-face sum is consequently exactly

    S_rho=(1/2)sum_j [N_j/ell_*] product_{i!=j}J_rho(w_i).

Using 1<=N_j<=2, 1/2<=ell_*<=1 and eventually 1/2<=J<=2 gives the
two bounds n/4 and 4n times the full product. These estimates sum all
2n source faces, including every possible omitted-coordinate term.

On any fixed target face, exactly one target coordinate is pinned. Its
J factor is at least 1/2 for the first moment and its product of two J
factors is at most four for the second. The other n-1 coordinates are
independent uniforms on that target face. Since n/rho<=2x and
xhd=epsilon, their product mean is at least exp(-24 epsilon), with the
factor 1/2 retained. This is greater than 1/4 after the pinned factor.

The matching source interval has logarithmic length between h and 4h,
even at the endpoints of I. This proves conditional row mean at least
M/16. For the second moment one integrates both source radii before
target expectation. The shared-W crossmoment from Section 2 yields

    E[product_i J_rho(w_i)J_rho'(w_i)]
      <=4 exp(31 n d^3/rho_min^2)<=12

eventually, because rho_min>=t/2 and x^2/n tends to zero. Multiplying
the complete source-face bounds and both logarithmic lengths gives
3072M^2<=4096M^2, exactly in the required upper direction.

The full target mass lies between nd/4 and nd. Therefore E0>=ndM/64
and the integrated second moment is at most 4096ndM^2. These are bounds
for the complete row, not for a chosen source face or a chosen radius.

## 4. Source high-Phi cutoff and target row cutoff

On each fixed maximal source face, the pinned contribution to Phi is at
most d. For the free coordinates, -log ell<=2(1-ell), so their mean
contribution is at most d^2/(2rho) each. Uniformly in the source radius,

    E Phi<=d+xd^2<=2xd^2

when xd>=1. Markov under the ORIGINAL full source shell law therefore
removes mass at most 2nd/K. The correct removed weighted energy is
integral_bad C0(z)d lambda0(z), and Section 1 already bounded C0,
including its omega, by 4M. Thus the loss is 8ndM/K=2^-17 ndM.
There is correctly no Q multiplier in this mass deletion estimate.
Using a source weight ceiling in place of the normalized column bound
would give a different, inadequate estimate; the frozen proof does not
make that replacement.

Deleting targets with R0>H_R M costs at most integral R0^2/(H_R M),
which is 2^-12 ndM. This bounds removed energy, not just removed target
mass, and uses the complete row second moment established above.

## 5. Both original corner-cap complements

For target radius r<t the original cap-complement condition is that all
coordinates are at least 1-beta r. The negative pinned face is excluded;
each free coordinate on a positive pinned face has retained fraction

    ell_cap(r)=1+1/(2t)-1/(2r)<=3/4

eventually. Its mass is bounded by nd(3/4)^(n-1). Cauchy--Schwarz with
the full row second moment gives the stated energy loss
64ndM(3/4)^((n-1)/2), which is eventually at most 2^-12 ndM.

For source rho=t(1-v), the cap threshold is -rho+v. Again the negative
pinned face is excluded; free-coordinate retained fraction is exactly
1-v/(2rho). With v>=d/4 and rho<=t the complement mass is at most
nd exp(-xd/16). Its removed energy is at most 4M times that ORIGINAL
mass, using the same inverse-normalized column estimate, and eventually
at most 2^-12 ndM.

All deletions are from the same full weighted kernel. Union-bounding
their energies is valid even where the deleted sets overlap. The exact
remaining lower coefficient is

    1/64-3/4096-1/131072=1951/131072>=1/128.

The final kernel has every ambient row at most H_R M and every ambient
column at most 4M<=H_R M. Source and target distributions of an unknown
actual measure have not been substituted for these ambient bounds.

## 6. Weighted energy lemma and the actual original history

For a nonnegative bounded weighted kernel, the stated expansion

    E_lambda,lambda-E_mu,mu
      =E_(lambda-mu),lambda+E_mu,(lambda-mu)

is exact. Since mu<=lambda, each term is bounded using the appropriate
ambient degree. The disjoint target/source annuli ensure that their two
deficit masses sum to at most the TOTAL envelope deficit. The energy
lemma therefore applies to omega>1; no K<=1 assumption is needed.

The original finite-a bridge was checked against the definitions in the
accepted R033 actual-increment/envelope proof and the original paper.
Its critical points remain valid for this new weighted kernel:

* Actual scaled increment points have r_i>=m(z); the increments are
  disjoint and their measure integrates to the genuine original weights.
* D is contained in H_a because the displayed threshold difference is
  positive for m>=1. Thus lambda_a=alpha^-n lambda EXACTLY on A and P.
* Every actual target item meeting A has c r_j<=||u-w||<3 and hence
  r_j<3, since c>1. Every source item meeting P has r_i>=rho>4.
  The source is therefore earlier in the same radius-ordered history.
* The center displacement is less than
  (kappa+2/(a+1))r_i<alpha beta r_i. For a>=8n^4 the positive margin
  is zeta-(beta+2)/(a+1)>0. Thus every contributing source item really
  captures that one target center in the original eligibility test.
* Integrating only its actual increment pieces gives an UNWEIGHTED row
  at most tau. On the retained P, omega<=Q, so the actual WEIGHTED row
  is at most Q tau. The proof explicitly pays this factor; neither
  inverse normalization nor mean source weight removes it.

Writing c_a=alpha^-n, energy scales by c_a^2, degrees by c_a, and target
mass by c_a. Hence the finite-a deficit bound is

    [c_a^2 E-Q tau c_a lambda(A)]/[c_a H_R M]
      =[c_a E-Q tau lambda(A)]/[H_R M]
      >=[E-Q tau lambda(A)]/[H_R M].

This direction is correct even though the kernel was designed using the
limiting corner measure. With Q tau<=M/256, the numerator is at least
ndM/256, giving Delta>=2^-32 nd uniformly over all original finite
histories and all a>=8n^4. No limit of actual histories is invoked.

## 7. Paying Q tau and optimizing the COMPLETE scalar cost

The old cap integral bound gives B_n>=nL/2 eventually, while
B_n<=nL and A_n<=exp x=n exp(-Bd^2 L). Hence B_n-delta0>0 for
delta0=2^-32 nd. The chosen actual threshold tau is well-defined.

Since Q=exp(Kxd^2)<=exp(KLd^2) and B=4K,

    Q tau/M
      <=(sqrt L/h) exp[-(B/2-K)d^2 L]
      =(sqrt L/h) exp(-Kd^2 L) ->0.

The prefactor is polynomial in L; Kd^2 is a FIXED positive constant.
Consequently the actual Q tau<=M/256 guard holds eventually. This
condition does not use the desired new deficit to prove its own
positivity or existence.

The final finite-a COVER FACTOR is retained. Letting a tend to infinity
in the already uniform numeric bound gives exactly

    (1+A_n/tau)(tau+B_n-delta0).

At tau=sqrt[A_n(B_n-delta0)] this equals
F(A_n,B_n-delta0). The new deficit is therefore inserted into the
complete cover/load cost, not just into one row term. The optimization
identity was independently checked symbolically.

Finally x=(1-Bd^2)L and 0<Bd^2<1/2 give

    B_n<=T+2Bd^2 n+1,
    2Bd^2 n=delta0/512.

The full positive remainder A_n+2sqrt(A_n B_n)=o(n)=o(nd), because d
is fixed here. More explicitly, eventual sufficiency for the claimed
half-deficit is

    1+A_n+2sqrt(A_n B_n)<=(255/512)delta0,

which follows from that limit. Thus the two audited conclusions are

    C_n<=F(A_n(t),B_n(t)-2^-32 nd),
    C_n<=T-delta0/2=T-2^-97 n

for one common sufficiently-large-n domain. The exact Fraction check
confirmed the ratio 1/512 and the final exponent 97.

## 8. General L1, all positive scales, and eventual-domain status

At each point of an original strict level set, choose any true positive
witness radius. Its a-shrunk open cube covers that point, and a compact
subset has a finite subcover. Every selected a-dilate remains one of
these true witnesses. Apply the uniform finite-a bound and then
Lebesgue inner regularity. This supplies arbitrary L1 input and every
positive common radius without a support truncation, bounded-scale
reduction, convergent selector sequence, or sampled-family assumption.

Every eventual guard used above has a common domain, with the one fixed
d=2^-64 and the fixed rational constants. In particular x is a fixed
fraction of L, h/d tends to zero,
zeta t/h=d/(epsilon n^2) tends to zero,
n d^3/rho_min^2=O(x^2d^3/n) tends to zero,
the source-cap exponential tends to zero, and Q tau/M tends to zero.
Only finitely many such guards are used. This verifies the qualitative
existence of one threshold; it does NOT compute or certify that threshold,
and no previously certified R042 threshold transfers to these parameters.

## 9. Local evidence and files

Actually run:

    python3 cube-dim-7-bounds/audit/round045/layered_source/symbolic_check.py

Local Python 3.15.0 and existing SymPy 1.14.0 were used. The check verifies
the loss integral, inverse-coordinate mean, excluded-tail primitive,
reference mean J=1, complete scalar optimization identity, and exact
energy/deficit/buffer constants. It contains no integration basepoints,
numerical quadrature, dimension scan, or finite-threshold certification.
Results are in `symbolic_result.json`.

All new audit files are confined to this owned audit directory. No
author/main checkout was changed; no shared MATLAB, Git synchronization,
old low-dimensional project, or n=7 work was used.
