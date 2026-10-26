# Full-time high-endpoint budget from joint radial and angular distances

Fix the original body32 frame, q=6/5 and rho=1/64. Read the complete
C037 intervals M1,M2,M3, without evaluating their defining integrals.
For A containing1 and proper, put k=|A|, m=6-k,
U_A=span(1_A,1_B), T_A={dist(x,U_A)<rho} on F={(1,y):y in[0,1]^5}.
The target H is any subset of

    union_A {x in T_A: mean_B(x)>3/4} intersect {dist(x,span ones)>=rho}.

In particular the proposed disjoint high cell, after removing cap, low
and middle cells, is such a subset. No claim that the tubes are disjoint
is made. All estimates are for the full positive-time Green G.

## 1. Geometry and the endpoint-ones incident terms

On a high tube write x_A=alpha1_A+e and x_B=u1_B+z, with zero-sum
residuals and r^2=|e|^2+|z|^2<rho^2. Since x1=1 and x_i<=1,
e1 is maximal and alpha=1-e1 lies in[1-rho,1]. Thus all B coordinates
are >3/4-rho=a=47/64, and A coordinates are >=1-sqrt(2)rho>a.
Also alpha,u belong to[3/4,1]. Orthogonal decomposition gives

    s^2:=dist(x,span ones)^2=r^2+(km/6)(alpha-u)^2
      <=rho^2+3/32=R^2=385/4096.                       (1)

C037 gives G<=Phi=sum_(ordered v!=w) M_j/(|x|^2 d_v^2 rho_vw).
There are exactly31 pair planes containing ones, each represented by
the pair (ones,w) and its reverse. To prove completeness, the span of
two body vectors has agreement and disagreement supports; a body vector
in that plane must take coefficients with |a+b|=|a-b|=1, hence is one
of the original two body lines. Thus a plane containing the body line
ones must have ones as one column. Its other column w!=ones determines
a distinct unsigned partition plane; no two different w give the same
plane. The three orbit counts for w are

    n=(binom(5,1)+binom(5,5),
       binom(5,2)+binom(5,4),binom(5,3))=(6,15,10).

These give62 ordered incident terms. The remaining counts are
(192,480,320)-2n=(180,450,300)=30n, so no term is discarded and no
full direction enumeration is required. Set W=6M1+15M2+10M3.

On H, |x|^2>=r0^2=1+5a^2. Every non-ones body line has

    d_w^2>=D^2=(10/3)a^2.                              (2)

Indeed write x=a ones+delta with delta>=0. For a body vector with h
negative entries, 1<=h<=5, the quadratic distance at a ones equals
(2/3)h(6-h)a^2; its linear increment has coefficients
2a[1-w_i(6-2h)/6]>=0, and its quadratic increment is nonnegative.
The minimum of h(6-h) over this integer range is5.

Every pair plane not containing ones has mixed signs on at least one
of its agreement/disagreement support blocks. Two opposite signs in
that block give residual squares at least (x_i+x_j)^2/2>=2a^2.
Therefore its actual distance is >=g=sqrt(2)a. This covers all930
nonincident ordered terms, including planes mixed on both blocks.

## 2. Exact five-dimensional projected coordinates and symmetry

Let z=P_(ones perp)x and s=|z|. Because x1=1, the inverse is
x=ones+z-z1 ones. The map from the five face coordinates to z has
Gram matrix I5-11^T/6, determinant1/6. Hence dy=sqrt(6) dz.
The condition x_i<=1 is exactly the maximal-coordinate cone
C1={z1>=z_i for every i}. The six maximal cones partition ones perp
up to null ties. This is a fraction1/6 only for permutation-invariant
integrands and sets; no individual-plane cone symmetry is assumed.

For the31 planes above, let ell_w=P_(ones perp)span(ones,w), a line.
The true distances are d_ones=s and rho_(ones,w)=dist(z,ell_w).
Define S(z)=sum_(w!=ones) M_class(w)/dist(z,ell_w). The full collection
of lines, with class weights, is invariant under all coordinate
permutations. Thus S and every union of equal-radius tubes about these
lines are invariant, and their full integrals on C1 equal1/6 of their
full-space integrals. Combining all992 terms yields on H

    G(x)<=r0^(-2)[s^(-2)S(z)+D^(-2)S(z)+30W/(D^2 g)].  (3)

This retains the actual joint line-plane singularity s^(-2)/dist(z,ell)
in the first31 ordered terms. It does not replace s by rho or assign
all angular peaks their maxima simultaneously.
The uniform positive-coordinate bounds are applied only on the true H.
After (3), its nonnegative right side is integrated on an enlarged cone
and angular cover; (2) and the mixed-plane gap need not hold at the added
points, and are not asserted there.

## 3. Angular covering and rearrangement, including all overlaps

Write z=s omega on S4, whose area is A4=8pi^2/3. Since H lies in at
least one unsigned tube, its angular section lies in

    E_s=union_(31 lines ell) {omega:dist(omega,ell)<rho/s}.

For a fixed axis, put t=sin(theta), c=sqrt(1-t^2). Uniform sphere measure
has the exact cumulative distribution

    P(sin(theta)<t)=1-(3/2)c+(1/2)c^3
      =(1-c)^2(c+2)/2 <=t^4, 0<=t<=1.                 (4)

The difference factor is 2(1+c)^2-(c+2)=3c+2c^2>=0.
A positive union estimate, with no disjointness assumption, therefore
gives |E_s|/A4<=mu(s):=min(1,31(rho/s)^4).

For f=sin(theta)^(-q), (4) implies P(f>lambda)<=lambda^(-4/q).
The decreasing rearrangement consequently satisfies f*(v)<=v^(-q/4).
For every measurable angular set of measure fraction <=mu,

    (1/A4) integral_E sin(theta)^(-q) <=
          mu^(1-q/4)/(1-q/4).                         (5)

One can obtain the same inequality by layer-cake integration and the
bound min(mu,lambda^(-4/q)). It includes zero-angle singularities as
null sets. Minkowski on E_s gives

    [integral_Es S(omega)^q]^(1/q)
      <=W [A4 mu(s)^(1-q/4)/(1-q/4)]^(1/q).            (6)

There is no assumed independence between lines or between E_s and the
fixed angular singularity. Equation(5) is valid for every measurable E.

## 4. Fixed radial integrals and the complete high budget

Set s0=31^(1/4)rho. The exact gates rho<s0<R are checked by powers.
Put b=7/5=5-3q, K=A4/[sqrt(6)(1-q/4)]. The factor1/sqrt(6) is
sqrt(6) times the legitimate permutation-cone fraction1/6.
Integrating the two branches of mu at s0 gives

    J3=K/b [2s0^b-rho^b-s0^(2b)/R^b],                (7)
    J1=K [s0^(2b)R-(14/19)s0^(2b+1)
                         -(5/19)rho^(2b+1)],          (8)
    V=A4/sqrt(6) [31rho^4(R-(4/5)s0)-rho^5/5].        (9)

These are upper-integral formulas, not exact integrals on clipped H.
Before simplification, their positive radial integrands are

    J3=K integral_rho^R s^(4-3q) mu(s)^(1-q/4) ds,
    J1=K integral_rho^R s^(4-q) mu(s)^(1-q/4) ds,
    V=A4/sqrt(6) integral_rho^R s^4 mu(s) ds.

They are finite and positive. Minkowski applied to (3), using (6), proves

    ||G||_Lq(H) <= Uhigh
      :=W/r0^2 [J3^(5/6)+J1^(5/6)/D^2
                            +30 V^(5/6)/(D^2 g)].     (10)

All positive times and all992 ordered terms are paid; all31 possible
incident planes and their overlaps are retained. Cap subtraction is
used only through s>=rho, and removal of low/mid regions is safely
relaxed. The radiusR and transition s0 are forced by (1) and covering,
not tuned numerical parameters.

If this fixed formula exceeds a desired whole-face budget, that is a
failure of this particular certified upper expression to meet that
target. It is not a lower bound on G, not an impossibility theorem for
sharper angular integration, and not an exclusion of the original Rc
route. No global g, middle long-time tail, far region or whole-face I
is established here; there is no new C6 bound.
