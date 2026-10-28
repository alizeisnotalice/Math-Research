# R078 theory: bounded-coordinate Bernstein retention in the R073 kernel

Status: proposed replacement lemma and full-cost consequence for the main
agent's independent review. The R073 theorem is accepted in the source
report; this note does not itself update the accepted ledger or paper.
Only this owned directory is written. No commit or push was made.

The fixed source is `upper/round073/cantelli_frontier/PROOF.md`, SHA256
`50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b`.
All its geometry, target set, kernel and parameters remain fixed. The
change is solely its Section 6 retention proof and the resulting gamma.
The proof source and R072 report were read, together with the available
R075 empirical-transfer proof. R074/REPORT.md and R075/REPORT.md were
absent from the specified original checkout. Neither is needed below.

## 1. Precisely applicable law and uniform coordinate bound

Fix an integer n, L=log n>=8192, x=(231/256)L>=7392, any target
w in the R073 good set A, any source q in its entire matching interval,
and any captured signed source face j. Work before source-cap deletion
and before bounded soft clipping. This is the R073 auxiliary inverse-
captured fullsource law. It is not the actual increment measure.

The active guards give R>1, rho>4r, a0>0 and 0<g=R-a0<=R/3.
For every source coordinate z_i in [-rho,rho],

    kappa rho=rho+a0,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R))
               >=(R+a0)/(2R)=1-g/(2R)>=5/6.

Thus all coordinate losses satisfy

    0<=phi_i=-log ell_R(z_i)<=Q0=log(6/5)<b=183/1000.       (1)

This holds on the full source cube, hence on every captured portion and
on every active matching q, not just at the baseline. The source-face
pin has deterministic loss Qpin=-log(1-g/(2R))<=Q0<b.
The free coordinates are independent in the inverse-capture law:
both the capture indicator and product inverse weight factor over
coordinates, and their normalizations are the individual J_i. A face
may omit any coordinate; it need not omit the target pin.

R073 Section 5 consequently supplies for V_j=Qpin+sum_(i!=j)phi_i

    E V_j<=p+m x, Var(V_j)<=v x,
    p=b=183/1000,
    m=2211/[50000(1-2^-20)], v=531/100000.                (2)

The variance is the true sum of centered free-coordinate variances,
with no random-pin term. No lower variance or empirical-profile gate
is needed. Write Y_i=phi_i-E phi_i. Then E Y_i=0 and
|Y_i|<=b, since both phi_i and its mean belong to [0,b]. The single-
coordinate centering does not subtract the larger total face mean.

## 2. Self-contained Bernstein lemma and finite choice of lambda

For independent centered Y_i with |Y_i|<=b and variance sum<=V,
the inequalities E|Y_i|^k<=b^(k-2)E Y_i^2 and
k!>=2*3^(k-2) (k>=2, by induction) give for 0<=lambda<3/b

    E exp(lambda Y_i)
       <=1+lambda^2 E Y_i^2/[2(1-b lambda/3)],
    log E exp(lambda sum_i Y_i)
       <=lambda^2 V/[2(1-b lambda/3)].                   (3)

The series can be integrated termwise because the variables are bounded;
absolute moments control all remainders. Independence gives the product
of the individual moment generating functions. Markov's inequality and
the admissible choice lambda=u/(V+bu/3), for u,V>0, imply

    P(sum_i Y_i>u)<=exp[-u^2/(2(V+bu/3))].                (4)

The lambda lies strictly below 3/b because V>0. Here V=vx>0 is a
deterministic upper bound, so (3) remains valid even if the actual
variance is zero. Apply (4) with

    chi=63/1280, a=chi-m=127067153/25420000000>0,
    u(x)=a x-p,
    H(x)=u(x)^2/[2(vx+(b/3)u(x))].                      (5)

Since chi*x-E V_j>=u(x), the event V_j>chi*x is contained in
the event sum Y_i>u(x). Thus, for every captured signed face,

    P_inv(V_j>chi*x | j)<=exp[-H(x)].                     (6)

All quantifiers precede this inequality; it is uniform in the retained
target, matching q, face and sign. It is therefore legitimate to mix the
facewise probabilities after establishing (6). No variance formula for
the face mixture is claimed.

## 3. One exact base check and proof over the entire dimension tail

At x0=7392, exact Fraction arithmetic in verify.py gives

    u0=14603570859/397187500>0,
    H0=71088093944611332627/4364039646326901875>16,
    S32(16)=sum_(k=0)^32 16^k/k!>1000000.                (7)

The positive exponential series proves exp(16)>S32(16)>10^6;
no floating exponential or numerical logarithm is used as a certificate.
For x>=x0, u>0 and direct differentiation gives, with d=b/3,

    H'(x)=u[a v x+p v+d a u]/[2(vx+d u)^2]>0.            (8)

Every term in the numerator is positive. Hence H(x)>=H0>16 on
the whole tail, and the NEW uniform retained fraction is

    P_inv(V_j<=chi*x | j)>Ret78=999999/1000000.           (9)

This is a finite-tail inequality, not a central-limit approximation,
compound-Poisson limit, dimension grid, or assertion about actual mu.
Equation (9) replaces only R073 Section 6 Cantelli's Ret73=9717/10000.

## 4. Unchanged bounded kernel and complete energy transfer

Use the exact R073 kernel with bsoft=min(1,theta/G), theta=exp(-chi*x).
The ratio bsoft/(theta/G)=min(1,G/theta) equals one on
V=-log G<=chi*x and is nonnegative everywhere. Therefore (9) raises
its full inverse-row lower coefficient from Ret73 to Ret78. The upper
coefficient is unchanged. Source cap deletion is then paid by the SAME
unweighted cap budget eS; independence is not asserted after deletion
or clipping. Both weights remain at most one. The complete column
bound remains 2k; target set and its mass are unchanged.

In the exact notation of R073, fix

    eB=eN=eT=2^-19, eJ=eRad=eS=eD=2^-20,
    P0=1-eB-6/768,
    Ulo=P0(1-eN)/(1+eJ), Uhi=(1+eN)/[P0(1-eJ)],
    M=2(1+eRad)Uhi, mL=3043/10000,
    brow78=2 Ret78 (1-eRad)Ulo-eS-eT,
    Gamma78=mL*brow78/M.

The exact new deficit gate is

    Gamma78=14891706016275974635743948569258721 /
            49711977879893038028841549824000000,
    gamma78=599/2000 < (1-eD)Gamma78 <1.                (10)

All original finite geometry, target-probability, source-cap and radial
gates remain unchanged. B_n>=nL/4>gamma78*n. The NEW threshold

    tau78=sqrt[A_n(B_n-gamma78*n)]

is positive and <=tau73, so the old upper bound on tau/k is still
valid. The largest-eligible-radius algorithm is used with this selected
numeric threshold. R073's same-history proof uses no lower bound on
tau: true target items still have radius<3, source increments radius>4,
the eligibility test still bounds the actual weighted row by tau78,
and weights are unchanged and <=1. The same finite-a cap equality,
ambient degrees, energy difference, target mass and finite-family bridge
give Delta_a>=Gamma78*n>gamma78*n. These statements must be reviewed
with the unchanged R073 Sections 7--10; (9) is not by itself an operator
theorem.

## 5. Complete root cost, display coefficient and comparison scope

Retain the entire original cover factor and terminal self-weight. The
R073 finite-a-to-numeric-limit argument, with tau78, gives the candidate
full bound

    C_n<=Phi78=F(A_n(t),B_n(t)-(599/2000)n),
    F(A,B)=(sqrt A+sqrt B)^2,
    t=n/[(231/256)log n], log n>=8192.                  (11)

The old scalar price -log(231/256)<257/2500 and its exact entire-
tail remainder budget 1/5000 are unchanged. Indeed

    599/2000-257/2500-393/2000=1/5000.

Thus the same scalar proof gives

    Phi78<T-(393/2000)n, T=n log n-n log log n.          (12)

There is also a direct comparison of full functions, not only displays:
on this domain A_n>0, B_n-gamma78*n>0 and gamma78-gamma73=17/2000>0.
F(A,B) is strictly increasing in B>0. Consequently

    Phi73-Phi78=(17/2000)n
      +2 sqrt(A_n)[sqrt(B_n-(291/1000)n)
                    -sqrt(B_n-(599/2000)n)]
      >(17/2000)n>0,
    Phi73=F(A_n(t),B_n(t)-(291/1000)n).                 (13)

This strictly improves the R073 branch on its entire valid domain.
It does not claim to improve every other accepted branch or the lower-
dimension domain. The R069 branch remains stronger in its stated larger
domain, and the accepted envelope must still take the minimum using each
branch's own valid domain. No paper edits or global-optimality claim.

## 6. Search, actual evidence and finite acceptance object

A read-only search of all original upper/*.md for Bernstein/Bennett
found only the R035 typical-distribution discussion: changing its old
Markov estimate to Bernstein may alter constants but does not break its
old O(n/sqrt L) bottleneck. That statement concerns a different old
one-shot typical geometry and does not exclude the present soft inverse-
source retention proof. No previous rejection was found in R060--R077.

The finite acceptance object is the self-contained lemma (1)--(9), the
new exact scalar gates (10), the unchanged-bridge audit in Sections 4--5,
and comparison (13). verify.py checks only new finite arithmetic; it
does not rerun the old 53 gates or certify their entire original proof.
verification.json records exact fractions and source identity. Python
was actually run once for the final script. No MATLAB session,
parameter grid, CP-tail calculation, n=7 or lower-bound work was used.

The independent local clone initially failed checkout because the source
partial clone lacked unrelated objects. A fresh remote sparse clone was
created at `cube-r078-theory-bernstein-sparse-20261028`, HEAD
`a2a41644d82c97ecd973ad8c96104a565aad67f3`; only the frozen R073 source
directory was checked out, followed by these owned R078 files. Original
checkout HEAD at read time was `b5fcd9480a8ddd08968c1f4f60cd7cd50c47ca10`.
