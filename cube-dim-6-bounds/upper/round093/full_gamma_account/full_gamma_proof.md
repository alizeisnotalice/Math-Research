# A093: complete local Green upper-Hessian account and its estimate-only obstruction

The object is the original unit-weight body32 Cauchy Green function G, with
the exact normalizations accepted in076/085. Fix

    x*=(1,1,1,1,1/2,0), rho=1/512, R=1/41472, T0=1/64,
    Q={(1,x2,...,x6): |(x2,...,x6)-(1,1,1,1/2,0)|_infinity<=h},
    h=1/2048,                  g0=1/65000 (candidate only).

This note pays every term in the076 decomposition, including the virtual
normal-deficit upper tail. It gives a finite **complete chosen upper bound**
Gamma_env for D²G on the rho ball. It does not compute the smallest valid
semiconcavity constant. The chosen bound is much too large; the resulting
vertex certificate is provably unusable on Q, even with exact vertices.
That failure excludes this budget expression, not the actual Green function,
the semiconcavity theorem, the original model, or every possible curvature
estimate. Actual G/F vertices, matching positive-window jets, global g_L,
complete I_U, and a positive final margin S remain UNPAID.

086's aggregate theorem was recorded as author-sealed/pending independent
review when this arithmetic batch was frozen. The agent subsequently received
and read its independent PASS audit; the parent retains final acceptance. Original source pins and
the registration are not retroactively edited to claim earlier acceptance.
The new dependence status is recorded separately.

## 1. Geometry and notation inherited without rerunning085

Use085's orthogonal coordinates x=(y,z), with U={z=0}, y*=(2,1/2,0),
and its complete four-column in-plane frame (2,+/-1,+/-1). The independent
primary appendix identifies U as the unique incident rank3 span at x*.
All nonincident rank3 spans are at distance at least1/32, and the true S2
distance at x* is sqrt(1/8). Thus the076 decomposition holds throughout
the closed rho ball, with the fixed scalar cutoff R above. In particular,
dist(x,S2)>1/4 there. The entire symmetric five-dimensional box Q lies in
that ball since5h²<=rho².

The full box is in the hyperplane x1=1, but extends beyond the positive unit
face at x2,x3,x4=1 and x6=0. Its volume is(2h)^5; its intersection with
[0,1]^5 has volume2h^5. The full-box certificate below is deliberately
stronger than the clipped-face requirement. No out-of-face volume is charged
to the actual positive-face integral I.

Let A_k, k=0,...,4, be rational upper bounds for ||D^k a_U|| on085's natural
tangential box. Let E_k be the global C4 extension bounds, and Z_k,
k=0,1,2, the global C2 bounds for the in-plane correction r_U in

    p_t^U=t³ a_U+t⁵ r_U,                0<t<=T0.

These are exactly the085 stored scaled rational bounds divided by27, using
pi³>=27. They are extracted in `upstream_constants.json`, not recalculated.
The extension equals the original amplitude for |y-y*|<=3rho, so the natural
A_k apply to coefficients evaluated at y in our ball; E_k are needed only
where the full outside fiber may move y-Cw arbitrarily far.

The fixed quotient differentiation basis is B_J=2I3, with tangent transport
||C_J B_J^-1||<=1. Its scalar shift alpha=B_J^-1 q satisfies, for |q|<=1,

    |D P_t|<=P_t/t,        |D² P_t|<=(5/4)P_t/t².          (1)

Indeed ||alpha||1²<=3/4 and2||alpha||2²<=1/2, using the Cauchy log-score
bounds. These inequalities concern the pure selected Cauchy product, without
cutoffs. The28-label coarea moments from085 become the rational bounds

    M0<=m0 t^-3, M1<=m1 t^-2, M2<=m2 t^-1,
    (m0,m1,m2)=(1/216,11/108,11/54).                       (2)

Here M1 integrates |w| and M2 integrates |w|² over Bw=z against the original
28-variable law. They have not been replaced by a seven-variable weighted
law. Also ||C||²=28 and ||C||<=6.

## 2. Exact decomposition and completeness of the ledger

The fixed parameter partition of076, the unique-U qualification, and its
amplitude expansion give the local identity

    G = R4 + R5 - Rcut + Rquad - Rsigned
        -a_U(y)N(z) + Rnormal + H_T,                      (3)

where

    R4 = integral_0^T p_t^(rank>=4)dt,
    R5 = integral_0^T t⁵ integral_(Bw=z) r_ext(t,y-Cw) product(c_t chi) d_Bw dt,
    Rcut = integral_0^T t³ integral_(Bw=z) a_ext(y-Cw)
                                      [1-product chi(w_j)] P_t(w) d_Bw dt,
    Rquad = integral_0^T t³ integral_(Bw=z) E(y,w)P_t(w)d_Bw dt,
    E(y,w)=a_ext(y-Cw)-a_ext(y)+Da_ext(y).Cw,
    Rsigned=sum_j [Da_U(y).c_j] L_j(z),
    L_j(z)=integral_0^T t³ m_j(t,z)dt,
    Rnormal=a_U(y)[T q_1(0)+D_T(z)],
    D_T(z)=integral_T^infinity t³(q_t(0)-q_t(z))dt,
    H_T=integral_T^infinity p_t^(original6D)dt.              (4)

P_t is the original28-variable product law; q_t its3D quotient density;
m_j its labelled signed first moment. In particular Rnormal and H_T are
different functions and both must be paid. Rnormal appears because the
negative deficit norm uses all positive times, whereas the in-plane expansion
is truncated at T. It cannot be discarded as an irrelevant constant: its
factor a_U varies tangentially.

All signs in(3) are retained. For the signed remainder pieces we pay an
absolute Hessian bound, so their signs do not require positivity. For the
negative cusp we use only the valid upper bound from085. Each displayed
integral is interpreted in the fixed-fiber coordinates of076. No moving
event indicator or selector is differentiated.

## 3. Rank>=4: a finite body32 worst-case bound without enumerating J

Every invertible6x6 matrix of actual body columns has determinant divisible
by2^5=32, hence absolute determinant>=32. To see divisibility for any +/-1
matrix, first sign-normalize its first row, then subtract one column from
the other five and factor2 from each. Every5x5 cofactor is divisible by16;
Hadamard gives its squared absolute value<=5^5<64², so its absolute value
is at most48. Therefore every inverse entry is at most3/2. The sum of the
squared row norms is at most54 over four tail rows and27 over two completion
rows. These facts hold for any086 fixed completion, not just a chosen sample.

There are at most binomial(32,4)=35960 independent four-subsets. From086,
write

    a1=4/R, a2=20/R², b1=1+4/R, b2=2+8/R+12/R²,
    d_tail=(a2+3a1²)T³/3+a1b1 T²,
    d_comp=(b2+b1²)T+2a1b1 T².

Its uniform physical PSD trace bound is dominated by

    Gamma4 =35960/[32*729*R^8] [54d_tail+27d_comp].          (5)

The729 uses pi^6>=3^6. The coefficients of T,T²,T³ in the bracket are,
respectively,27(b2+b1²),108a1b1,18(a2+3a1²). This pays every J, even though
no J list is generated. It uses a theorem covering the entire sum, not the
partial-record API whose completeness would remain an input obligation.
No2^32 patterns or actual body inverse matrices are enumerated.

## 4. Unexpanded cutoff-removal payment

Use

    1-product_j chi(w_j)=sum_j eta(w_j)product_(ell<j)chi(w_ell).

Each j term is differentiated in its registered quotient basis containing j.
The amplitude remains unexpanded. The large scalar factor and its two
derivatives have bounds t/pi times(R^-2,4R^-3,20R^-4). For another selected
factor, either c_t or c_t chi, write

    s1=1+2/R,       s2=2+4/R+6/R²,
    H1=max(4/R,s1), H2=max(20/R²,s2,H1²).                  (6)

The basis norm is at most Lambda=9/4. The tangent columns each have norm
at mostsqrt3, and the registered inverse Frobenius norm is at mostsqrt2.
Hence ||C_J B_J^-1||<=3sqrt2<9/2. For a unit physical direction h=(p,q),
the derivative of y-Cw has norm at mostsqrt(1+(9/2)²)<V=5. This estimates
the whole physical direction, not only a normal derivative.

The sum over all28 labels of D_j=|det B_(J(j))|^-1 is D=11/2. After solving
the basis coordinates and integrating the residual subprobability factors,
the three selected scalar factors contribute at most D_j/(pi³R²)t^-1.
Their first/second derivatives contribute additionally Lambda H1/t and
Lambda²H2/t². Retaining all amplitude product-rule terms gives

    Gamma_cut = D/(27R²) [V² E2 T³/3
                         +V Lambda H1 E1 T²+Lambda² H2 E0 T].     (7)

The first term comes from two amplitude derivatives, the second from the
two mixed hits, and the last from two density derivatives. Since t<=1,
replacing a large-factor derivative's absent t^-1 by t^-1 is a safe upper
bound. Its large support is retained. No weighted index must share a basis
with j, because the amplitude has not yet been Taylor expanded.

## 5. Quadratic amplitude error and the t5 correction

For a unit physical direction let p=Yh, q=Zh, alpha=q/2 and e=C_J alpha.
Then |p|,|e|<=1. Put d=Cw. Differentiating the Taylor error at a fixed fiber
gives the bounds

    |E|<=E2 |d|²/2,
    |E'|<=E3 |d|²/2+E2 |d|,
    |E''|<=E4 |d|²/2+2E3 |d|+E2.                          (8)

The first term of E'' is the second y derivative of the Taylor remainder,
the middle is its mixed y,d derivative, and the last is its d,d derivative.
Thus C4 is genuinely paid. Using |d|²<=28|w|², |d|<=6|w|, (1), and(2),
the integrand of the Hessian of Rquad is bounded by

    14 E4 m2 t² +(12 E3 m1+28 E3 m2)t
                     +E2 m0+12 E2 m1+(35/2)E2 m2.

Consequently

    Gamma_quad=14E4m2 T³/3+(12E3m1+28E3m2)T²/2
                         +(E2m0+12E2m1+(35/2)E2m2)T.       (9)

The fixed basis used for differentiation is distinct from the per-label
bases used to bound the positive weighted integrals; their coarea values
are equal, so this does not differentiate a basis choice. Fixed-t continuity
and the integrable bounds justify the time differentiation exactly as in076.

For R5 retain the cutoffs and use absolute selected-factor peak bounds,
rather than relative scores of a cutoff which may vanish. With

    h2=max(s2,s1²),

the selected product has peak m0 t^-3, first derivative bound m0 s1 t^-4,
and second derivative bound m0 h2 t^-5. The derivative of y-Cw has norm<=2,
so the complete physical Hessian payment is

    Gamma5=m0[4 Z2 T³/3+2 Z1 s1 T²+Z0 h2 T].              (10)

All three terms of the product rule are present. The frozen085 Z_k are
uniform for every0<t<=T0, so(10) also holds with any smaller T without
recomputing the old table.

## 6. Signed moments, their linear term, and their upper tail

The marginal quotient coercivity is16sqrt2>=c_q=22. Keep all28 original
weighted labels. For their containing-j bases use Lambda_j<=1 for12 axis
labels and Lambda_j<=9/4 for16 tetrahedral labels. Their determinant sums are

    D=sum D_j=11/2,
    D1=sum D_j Lambda_j=21/2,
    D3=sum D_j Lambda_j³=753/16.                           (11)

076's exact renormalization is

    L_j(z)=F_j(z)+T Dm_j(1,0).z-Q_j(T,z),
    F_j=-3/(4pi²) integral_S² sign(b_j.theta)
                         (theta.z)|theta.z|/psi(theta)^5 dtheta. (12)

The two terms in the all-time definition of F_j are not separately
integrable; the linear subtraction must stay inside that definition. The
following bounds pay both the pure z Hessian and the coefficient derivatives:

    ||D²F_j||<=F2=6/(3*c_q^5),
    sum_j ||D²Q_j(T,.)|| <=rho/T *48 D3/27,
    sum_j |L_j|<=D T²/(4*27),
    sum_j |DL_j|<=3D1 T/(2*27).                            (13)

The first two follow from076's angular formula and its third-derivative
bound K_(j,3)<=48D_jLambda_j³/pi³. The last two follow directly from the
absolutely finite coarea first moment: |m_j|<=D_j t^-2/(2pi³) and
|Dm_j|<=3D_jLambda_j t^-3/(2pi³). They therefore include the finite-time
linear term in(12); that term has zero pure z Hessian but nonzero mixed
terms after multiplication by a varying coefficient.

Let b_j(y)=Da_U(y).c_j. Since |c_j|<=2, its value, gradient and Hessian are
bounded by2A1,2A2,2A3. The natural A_k apply because y stays in the ball.
The full product-rule payment for Rsigned is

    Gamma_signed =2A1[28F2+rho/T *48D3/27]
                        +4A2[3D1T/(2*27)]+2A3[DT²/(4*27)].      (14)

The computation records these four pieces separately. This is a bound for
the entire finite-time signed layer with its coefficients, not only for
the angular Hessian of F_j.

## 7. Negative norm profile and normal truncation remainder

085's selected negative-profile semiconcavity bound was a rational
coefficient C_cusp times pi^-4. Here use

    Gamma_cusp=C_cusp/81.                                 (15)

The coefficient in the extraction is
3465893065395297320960/65143564997786663569923. It pays only -a_U N.

For q_t in dimension3, radial Fourier integration and angular averaging give

    q_1(0)<=Q0=1/(9*c_q³),
    ||D²q_t||<=Q2 t^-5,        Q2=4/(9*c_q^5).             (16)

The factor4 uses integral_S²(theta.u)²dtheta=4pi/3 for a unit u, followed by
integral r^4 e^(-ct r)dr=24/(ct)^5. Evenness gives Dq_t(0)=0. Therefore on
|z|<=rho,

    |D_T|<=Q2 rho²/(2T), |DD_T|<=Q2 rho/T, |D²D_T|<=Q2/T.

Rnormal in(4) has the complete payment

    Gamma_normal=A2 TQ0
                 +Q2/T [A2 rho²/2+2A1 rho+A0].           (17)

The first term is the variable-amplitude constant Tq_1(0); the second
contains its normal tail's value, gradient and Hessian. Both are recorded.

## 8. Original six-dimensional large-time Hessian tail

The accepted original frame identity AA^T=32I gives

    psi_6(xi)>=sum_v(v.xi)² / max_v|v.xi|
               >=32|xi|/sqrt6>=13|xi|.                    (18)

The zero frequency follows by continuity. No new frame enumeration is needed.
For a unit u, use area(S5)=pi³ and its angular mean (theta.u)²=1/6. Fourier
inversion with the original factor(2pi)^-6 gives

    |D_uu p_t(x)|<=105/[8pi³*13^8] t^-8,
    Gamma_infinity=15/[8*27*13^8] T^-7.                  (19)

This pays H_T for every x, with the correct6D normalization. It is distinct
from the artificial3D normal tail in(17). It does not extend a divergent
positive J2 moment down to t=0.

## 9. Complete upper bound and the observed main loss

Sum(5), (7), (9), (10), (14), (15), (17), (19). All are nonnegative rational
upper bounds. Subject to the accepted source hypotheses and086 final
acceptance, the complete distributional statement is

    D²G <= Gamma_env(T) I,   |x-x*|<=rho, 0<T<=T0,
    Gamma_env(T)=C0+C1T+C2T²+C3T³+C_(-1)/T+C_(-7)/T^7.   (20)

Every coefficient is strictly positive. The exact12-component/20-coefficient
ledger and total coefficients are in results G007. No analytic component
of(3) is omitted. At the only evaluated time T=1/64, the exact rational is

    624890693102357658534530657446023861791036558637538365390129657767989292076124686052975125519829712340930023398539001856
    /2200689204522269050109470405752843608840214497771433470971763107552359375,

displayed as2.83952269052008e47. The rank4 worst-case term dominates; the
next largest term, cutoff removal, is about1.73703926815951e19. This is a
diagnosis of estimate loss caused primarily by R^-8 and the derivative
cutoff powers in the aggregate majorant. It is not evidence that the true
curvature has this size.

## 10. Exact obstruction for this fixed vertex certificate and its time family

The face map x=(1,y_face) is an isometric coordinate embedding. With32 true
vertex lower bounds l_v, the semiconcavity certificate would require

    min_v l_v >=g0+P_env(T),
    P_env(T)=(5h²/2) Gamma_env(T).                         (21)

At T=1/64 the exact penalty is displayed as1.69248741300111e41. To interpret
this without importing unknown vertex values, we separately derive a coarse
but valid upper bound for the true G everywhere in the rho ball.

There dist(x,S2)>1/4. Set r_value=1/640. Since sum|a_i|<80, any parameter
fiber reaching this ball must have at least three independent columns whose
scalar parameters have modulus>r_value. Each chosen scalar density is at
most t/(pi r_value²). The selected3-column coarea Jacobian is at least1
because its positive Gram determinant is an integer. Its3D quotient has
coercivity at least13 by restriction of(18), so its peak is at most
1/(pi²13³t³). A positive union bound over at most4960 triples gives

    p_t(x)<=M=4960/[3^5*13³*r_value^6].                    (22)

No event is differentiated. Also original6D Fourier inversion gives
p_t<=P t^-6, P=15/(8*27*13^6). Splitting this value estimate at the fixed
T0=1/64 yields the true uniform bound

    G(x)<=Gbar=M T0+P/(5T0^5)
           =11700693505150503878656/1172914587
           <1.0e13.                                      (23)

The positive-density union estimate first holds almost everywhere; the
strict distance margin holds on a neighborhood, and continuity of the
fixed-time full density extends the constant bound to every point there.

The exact arithmetic verifies P_env(T0)>Gbar. Since any genuine lower
vertex endpoint is at most its actual G value, it is impossible for(21)
to pass using this Gamma_env, even if all32 vertices were known exactly.
There is no reason to spend new kernel calls on those vertices for this
specific certificate.

Merely tuning T in this same frozen account cannot rescue it. For all T>0,

    Gamma_env(T)>=C1T+C_(-7)T^-7
       >=(8/7) C1 [7C_(-7)/C1]^(1/8).                    (24)

The eighth power of this lower bound is exactly
8^8 C1^7 C_(-7)/7^7. The program compares its induced box-penalty eighth
power with Gbar^8 using integers/Fraction, and certifies that the penalty
lower bound is at least10^35, hence greater than Gbar. This applies
to every0<T<=T0 for which the frozen r bounds are valid. No second T was
evaluated and no time scan occurred. Strict convexity of the positive
expression, Gamma_env''=2C2+6C3T+2C_(-1)/T³+56C_(-7)/T^9>0, also gives a
unique balance, but computing its location is unnecessary after(24).

Crucially, (24) is a lower bound on the **chosen upper-bound expression**
Gamma_env, not on the actual smallest semiconcavity constant of G. Neither
(21)-(24) nor the huge value upper bound(23) excludes the original model,
the target g0, a sharper local representation, or a new joint estimate.

## 11. Same-box comparisons with a positive window and the original mixture

Fix the explicit comparison object F=integral_(a)^b p_t dt with a=1/250,
b=1/2, at the same x*, same body32 model and same full box Q. This window is
only a registered comparison object: no F value, gradient, or vertex was
computed here, and no old061 point is substituted.

The057 all32-column fixed-fiber right inverse R_A=A^T/32 has
||R_Au||1<=|u| and ||R_Au||2²=|u|²/32. The scalar Cauchy log upper second
derivative is1/(4t²). It follows for unit u that

    D_uu p_t <=129/(128t²) p_t,
    D_uu F <=129/(128a²) F,
    |grad log F|<=1/a.                                    (25)

Let F_U be a genuine matching center upper endpoint. On Q,
F<=F_U exp(sqrt5 h/a). Use sqrt5<=9/4 and exp(s)<=1/(1-s) for0<=s<1.
Here s=1125/4096, so a certified window upper-Hessian constant gives the
vertex penalty

    P_F <=[10078125/194707456] F_U
         approximately0.051760344503705084 F_U.            (26)

Thus a positive-window32-vertex certificate can use
min_v F_L(v)>=g0+(10078125/194707456)F_U. Every required scalar datum remains
UNPAID. Equation(26) is a finite comparison coefficient, not a passed box.

The original positive-mixture lower-log-curvature interface instead uses
one matching center value and its five face derivatives:

    log F(c+u)>=log F(c)+grad log F(c).u-|u|²/(32a²).

Its whole-box quadratic loss is exactly

    q_box=5h²/(32a²)=78125/33554432
                       approximately0.0023283064365386963.      (27)

For example a conservative sufficient interval condition is
log F_L-h D1_upper/F_L-q_box>=log g0, where D1_upper bounds the l1 norm of
the matching face gradient. The sharper shared-denominator corner interface
of057 may also be used. The gradients cannot be assumed zero. If they were
zero, purely for scale comparison, the required multiplicative factor
exp(q_box) lies between1+q_box and1/(1-q_box). This hypothetical case is
not claimed for the actual x*.

The data cost therefore favors testing a newly qualified matching positive-
window center jet over buying32 complete-G vertices for the excluded account.
That is a conditional future pilot decision, requiring separate parent
release and valid production runtime. No body/native call is authorized by
this result. F's fixed time window has no full-G spatial homogeneity; only
the eventual genuine G floor may enter the global symmetry/homogeneity
transfer.

## 12. Final cost S and precisely what remains unpaid

The adopted guidance uses q=6/5 and

    S=64(g_L/6)^(6/5)-I_U,
    I=integral_[0,1]^5 G(1,y)^(6/5)dy, B=I^(5/6).

For the **candidate** g0=1/65000, the conditional competitive norm gate is
B<4/48750; equivalently I_U^5<64^5(g0/6)^6. These are exact planning gates,
not an already certified g_L or I_U. This local curvature account improves
neither global endpoint and does not establish S>0. The current certified
upper bound for C6 remains64; the accepted lower bound remains
2073015624830/10^12.

The analytic account(3) is complete; its only source-acceptance condition
at this note's initial seal is the parent's final acceptance of independently
passed086. Actual vertices, a positive-window center jet, full-domain g_L,
complete I_U, and final S are UNPAID. The most important new information for
S is negative but actionable: the present worst-case aggregate/cutoff
expression cannot certify the prescribed box for any allowed T, so spending
more native calls or scanning T on this unchanged expression cannot reduce
the global g deficit. A sharper joint estimate must first remove its main
loss, or the positive-window comparison must be supplied real matching data.
