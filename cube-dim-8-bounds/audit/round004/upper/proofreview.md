# R004: independent audit of shared-time 4+4 Green anchors

2026-10-27. **Verdict: PASS for the three exact point values and the exact global peak/tail normalization.** These point certificates do not prove a continuous cube minimum, a numerical \(g_0\), a face integral \(I_0\), or a new bound for \(C_8\). The accepted general Green/source/obstacle theory is unchanged. The parent retains responsibility for the final proof and project integration.

For \(A=(1,1,1,0)\), \(B=(1,1,0,0)\), the verified results are

\[
G_8(A,A)=\frac{10495703317}{2149908480000}\pi^{-7},\quad
G_8(A,B)=\frac{15547486331}{1920360960000}\pi^{-7},\quad
G_8(B,B)=\frac{143}{10584}\pi^{-7}.
\]

The independently verified auxiliary improvement is

\[
\|p_t^{(8)}\|_\infty=p_t^{(8)}(0)
=\frac{121}{9437184\pi^8}t^{-8},\qquad
\int_T^\infty p_t^{(8)}(x)dt
\le\frac{121}{66060288\pi^8}T^{-7}
\]

uniformly in \(x\). This gives a tighter tail budget for the same model, not an operator-constant improvement by itself.

All authored deliverables are confined to `audit/round004/upper/`; older audit packages and source-author files are unchanged. No author script was executed. The independent implementation reconstructs the full density functions, rather than checking only finitely many times or running the author's checker. It shares the exact flint polynomial arithmetic library with the successful producer, a disclosed lower-level dependency, but uses a different real-variable assembly and a different partial-fraction algorithm. Fraction arithmetic handles the local expansions and integrations. SymPy is used only to parse the saved author polynomial strings for exact comparison, not to produce a cone sum or perform partial fractions.

## 1. Four-dimensional geometry and Fourier normalization

The model uses eight **unnormalized** directions
\(W=\{(1,\epsilon_2,\epsilon_3,\epsilon_4)\}\), each of length 2 and coefficient 1. Thus \(\psi_4(\xi)=\sum_{v\in W}|v\cdot\xi|\), with Fourier inverse factor \((2\pi)^{-4}=1/(16\pi^4)\).

The locked n4 geometry uses rays

\[
r_A=(1,0,0,0),\ r_B=(1,1,0,0),\ r_C=(1,1,1,0),\
r_D=(1,1,1,1),\ r_E=(2,1,1,0),\ r_F=(3,1,1,1),
\]

and bases \((r_A,r_B,r_E,r_F)\), \((r_E,r_B,r_C,r_D)\), \((r_E,r_B,r_D,r_F)\), of determinant magnitudes 1,1,2. The audit independently verifies their inverse rows. On the sorted chamber \(a\ge b\ge c\ge d\ge0\), their nonnegative coordinates are respectively

\[
(a-b-c-d,b-c,c-d,d),
\]
\[
(a-b,b-c,b+c-a-d,d),
\]
\[
(c-d,b-c,(-a+b+c+d)/2,(a-b-c+d)/2).
\]

The three regions are \(a-b-c\ge d\), \(a-b-c\le-d\), and \(|a-b-c|\le d\), including the chamber ordering. They cover the chamber and have disjoint interiors. All 384 signed permutations consequently give 1152 cones, whose overlaps occur only on null boundaries. The script independently generates all of them using a group-action convention different from the producer's coordinate lookup, checks uniqueness, and finds 240 rays.

For every body direction and cone, the direction's four ray projections have one weak sign, verified exactly. Hence \(\psi_4\) is linear on the entire cone; it is not inferred from an interior sample. The six ray rates are exactly \(8,8,12,12,16,24\).

For any \(t>0\), Fourier inversion is absolutely convergent. On a cone write \(\xi=\sum_j\lambda_jr_j\), \(\lambda_j>0\). Integrating four elementary exponentials gives the actual probability density

\[
p_t^{(4)}(x)=\frac1{16\pi^4}
\sum_{R}|\det R|\prod_{j=1}^4(tq_{r_j}-ix\cdot r_j)^{-1}.
\]

This establishes the prefactor and the meaning of the full sum, including negative real parts of individual terms. No cone or individual pole may be discarded.

## 2. Independent full rational identity, including zero and repeated projections

The producer groups signed projection multisets and accumulates a rational function in a complex variable \(z\). This audit instead sets \(u=t^2\), keeps every one of the 1152 cones, and uses one common real denominator from the start.

For a cone let \(b_j=x\cdot r_j/q_{r_j}\), and let \(e_k\) be the elementary symmetric functions of the four signed \(b_j\). Then, for \(t>0\),

\[
\Re\prod_j(t-ib_j)^{-1}
=\frac{u^2-e_2u+e_4}{\prod_j(u+b_j^2)},
\]
\[
\Im\prod_j(t-ib_j)^{-1}
=t\frac{e_1u-e_3}{\prod_j(u+b_j^2)}.
\]

The cone weight, including Fourier's rational factor, is \(|\det R|/[16\prod_jq_{r_j}]\). For each squared projection, the common denominator uses the largest multiplicity appearing in any raw cone denominator. Every raw zero factor and repeated factor is retained. Each raw numerator is multiplied by the exact polynomial quotient from that fixed common denominator, then **all 1152** real numerators and imaginary numerators are summed. The imaginary polynomial is identically zero. Only after this complete sum is formed does the audit take a numerator/denominator gcd.

The actual exact computation records:

| Point | Raw common denominator degree in \(u\) | Raw zero multiplicity | Cones containing zero projection | Cones with repeated signed projections | Reduced degrees \(N/D\) |
|---|---:|---:|---:|---:|---:|
| A | 16 | 2 | 528 | 816 | 6/8 |
| B | 12 | 3 | 672 | 864 | 2/4 |

The reduced denominators have no zero root. This certifies the cancellation globally; it does not assign an artificial value to zero projected nodes or delete a cone's divergent time contribution. Repeated projections are preserved as multiplicities throughout. The full function identities, not merely sample values, are verified against the saved author functions by exact cross multiplication.

Writing \(p_t^{(4)}(x)=\pi^{-4}P_x(t)\), \(u=t^2\), the resulting formulas can be expressed compactly as

\[
P_A(t)=\frac{N_A(u)}{48(144u+1)(64u+1)^3(16u+1)^4},
\]
\[
\begin{aligned}
N_A(u)={}&9+4080u+878336u^2+80203776u^3\
&+3149660160u^4+57965281280u^5+425201762304u^6,
\end{aligned}
\]
\[
P_B(t)=\frac{1+250u+25344u^2}{3(64u+1)^2(36u+1)(16u+1)}.
\]

All coefficients in these displayed numerators and denominators are positive, consistent with the true positive density. The reduced squared-pole multiplicities are

\[
A:\ (1/144,1),(1/64,3),(1/16,4),
\qquad B:\ (1/64,2),(1/36,1),(1/16,1).
\]

Their positivity and nonzero constant terms are exact, not inferred from high-precision samples.

## 3. Shared-time 4+4 normalization and regularity of the points

The eight-dimensional frame is the direct sum of two independent four-dimensional frames. At each fixed time,

\[
p_t^{(8)}(x,y)=p_t^{(4)}(x)p_t^{(4)}(y),\qquad
G_8(x,y)=\int_0^\infty p_t^{(4)}(x)p_t^{(4)}(y)dt.
\]

There is one common time integral. Multiplying \(G_4(x)G_4(y)\) would introduce two independent times and be a different quantity. The product has prefactor \(\pi^{-8}\), and its exact time integral supplies one factor of \(\pi\); this explains \(\pi^{-7}\) in each final answer.

Both A and B require at least two body directions in their own block: neither is a multiple of a body direction, since each has a zero coordinate. Two directions suffice (take the average of suitable body directions). The blocks are disjoint coordinate subspaces, so each paired point requires at least four independent active directions in total. It is consequently outside the eight-dimensional \(S_3\) singular set. The following direct rational convergence check independently establishes finiteness without relying solely on that classification.

After the full cancellation,

\[
P_A(0)=3/16,\qquad P_B(0)=1/3.
\]

The three complete product limits at zero are therefore \(9/256,1/16,1/9\), multiplied by \(\pi^{-8}\). The integrands are bounded and continuous near zero. These statements are made only for the fully combined density functions; they do not justify integrating individual raw zero-pole terms from zero.

For each raw cone and \(t>0\), \(|t-ib_j|\ge t\). Summing its absolute weight gives

\[
|P_x(t)|\le\left(\frac1{16}\sum_R\frac{|\det R|}{\prod_jq_{r_j}}\right)t^{-4}
=\frac{11}{3072}t^{-4}.
\]

The script exactly verifies the parenthesized coefficient. Hence every paired integrand is bounded by \(121/(9437184\pi^8)t^{-8}\), integrable at infinity. Together with the small-time result, the three integrals are proper and all subsequent **finite** partial-fraction terms are integrable.

## 4. Independent partial fractions and exact time integration

For each full rational product \(R(u)=N(u)/D(u)\), all denominator factors are \((u+a^2)^m\) with positive rational \(a\), and the function is proper. The audit obtains partial fractions without using the author's list or SymPy's `apart`/`residue`.

At the root \(u=-a^2\), divide out its full multiplicity and write

\[
H(v)=\frac{N(-a^2+v)}{[D(u)/(u+a^2)^m]_{u=-a^2+v}}
=\sum_{j\ge0}h_jv^j.
\]

The numerator and denominator Taylor coefficients are computed by the exact binomial theorem. If the denominator coefficients are \(d_j\), the quotient recurrence is

\[
h_j=\frac{n_j-\sum_{\ell=1}^jd_\ell h_{j-\ell}}{d_0}.
\]

Then \(h_j\) for \(0\le j<m\) is the partial-fraction coefficient of \((u+a^2)^{-(m-j)}\). This preserves all confluent multiplicities. After doing every root, the program reconstructs the entire common-denominator numerator and checks equality identically. The independent local expansions give 16 terms for AA, 12 for AB and 8 for BB. The saved author lists are separately cross-multiplied against the same reconstructed products and also agree identically.

For \(a>0\), the elementary integral is

\[
\int_0^\infty(t^2+a^2)^{-k}dt
=\pi\frac{\binom{2k-2}{k-1}}{2\,4^{k-1}a^{2k-1}}.
\]

One may derive it from \(I_1=\pi/(2a)\) and differentiation under the integral, giving \(I_{k+1}=(2k-1)I_k/(2ka^2)\). All \(a\) are positive rationals with exactly square rational \(a^2\); every coefficient is a Fraction, so summing these integrals produces exact rational coefficients of \(\pi\). Negative partial-fraction coefficients present no convergence issue because each term is individually integrable and the sum is finite.

The resulting coefficients of \(\pi^{-7}\) are exactly

| Point | Exact coefficient |
|---|---:|
| AA | \(10495703317/2149908480000\) |
| AB | \(15547486331/1920360960000\) |
| BB | \(143/10584\) |

As an independent dimensional-normalization consistency check, the same local-expansion integration applied to the individual density functions gives \(G_4(A)=17/(640\pi^3)\) and \(G_4(B)=1/(24\pi^3)\). These known four-dimensional identities are checks, not the input from which the eight-dimensional values were inferred.

## 5. Exact global peak and improved uniform tail budget

The new density identities both have leading coefficient \(11/3072\) in their \(t^{-4}\) asymptotics. This equals the raw weight sum computed independently above. There are two compatible strict explanations of the peak normalization.

First, at \(x=0\), every raw cone projection is zero and the complete positive-time density sum is directly

\[
p_t^{(4)}(0)=\frac{11}{3072\pi^4}t^{-4}.
\]

Second, the actual scaling \(p_t^{(4)}(x)=t^{-4}p_1^{(4)}(x/t)\) and continuity of \(p_1^{(4)}\) imply that either computed leading coefficient is \(\pi^4p_1^{(4)}(0)\). Both checks agree, without a new numerical experiment.

The Fourier integrand \(e^{-t\psi_4(\xi)}\) is nonnegative and integrable, so its inverse Fourier transform obeys

\[
|p_t^{(4)}(x)|\le(2\pi)^{-4}\int e^{-t\psi_4(\xi)}d\xi=p_t^{(4)}(0).
\]

The actual probability density is nonnegative, and equality is attained at the origin. The same-time product therefore gives the exact eight-dimensional global peak coefficient

\[
C_{\rm exact}=\frac{121}{9437184\pi^8},\quad
\|p_t^{(8)}\|_\infty=C_{\rm exact}t^{-8}.
\]

Relative to the earlier basis bound \(1/(256\pi^8)\), this is smaller by the exact factor

\[
\frac{1/256}{121/9437184}=\frac{36864}{121}
\]

(approximately 304.66). Integrating from \(T\) gives the uniform budget

\[
R_T(x)\le\frac{121}{66060288\pi^8}T^{-7}.
\]

For any finite-measure seven-dimensional face domain \(E\), this additionally implies \(\|R_T(1,\cdot)\|_{8/7,E}\le|E|^{7/8}121/(66060288\pi^8T^7)\). It does not alone give a finite full-plane norm by integrating a uniform constant over an infinite plane. The existing Green theory already explains the separate whole-plane interpolation budget; this audit does not alter that theory.

These are auxiliary bounds for the fixed density model. Without compatible continuous \(g_0\) and \(I_0\) certificates they do not produce a numerical improvement to \(C_8\).

## 6. Execution history and provenance

The author registration specifies only AA, AB and BB, the full 1152-cone treatment, exact rational integration after cancellation, a 200000-term stopping cap, and 80/160-digit numerical quadrature checks. The successful source implements these operations without a \(1152^2\) cone-product expansion, and its saved functions and integrals match the independent full identities above.

The initial direct SymPy calculation was **interrupted**, with `KeyboardInterrupt`, exit code 1 and wall time **66.53380966698751 seconds**. Both `run_interrupted.py` and `receipt_interrupted.json` are retained and hashed by the author manifest. At that time the receipt used the path `run.py` for those older bytes; this audit correctly matches its hash to the preserved `run_interrupted.py`, not to the later successful script. This recorded interruption is not discarded, described as success or conflated with the later run.

The subsequent exact flint accumulation completed with exit 0 in **2.0598439580062404 seconds**; its receipt hashes match the successful `run.py`, registration and locked n4 sources. The author's saved check reports about **2.9496 seconds**, and the manifest records check exit 0. All manifest file hashes and external source hashes were independently matched to actual files. Saved receipts are self-recorded execution evidence, not a separate immutable preregistration timestamp.

The 80- and 160-digit mpmath values agree with the exact forms in the saved result, but they are numerical checks without outward-error enclosure. This audit does not use them as proof or repeat quadrature. Its strict acceptance rests on full rational density identities, cancellation, properness, exact local partial fractions and exact integration. The actual independent run exits 0 and records its timing/interpreter/library versions/source hash in `actualreceipt.json`; its exact algebra and point results are in `result.json`.

The current scope remains three point anchors and a tighter exact global peak/tail budget. No global minimum, minimum location, continuous face integral, or numerical \(C_8\) claim is accepted from these points. Earlier Green theory and older audit materials remain unchanged.

## Reproduction

From the dedicated clone root, use the existing mathematical Python environment:

```text
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-8-bounds/audit/round004/upper/independent_check.py
```

This reconstructs geometry, every density term, all cancellations and the exact local partial fractions, and verifies saved-package hashes. It does not import or execute either author run/check script. The hash manifest records author files before and after the final freeze and verifies all earlier audit deliverables remain unchanged. No unrelated historical repository blobs were repaired or declared valid.
