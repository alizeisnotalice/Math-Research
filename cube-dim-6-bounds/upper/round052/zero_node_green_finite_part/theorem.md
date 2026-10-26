# Scalar Green finite part at zero Fourier nodes in dimension six

This is a new mathematical interface, with twelve formal scalar checks. It is
not a new complete point evaluation, numerical engine, Green-floor certificate,
or bound on the centered cube maximal constant. All sources in
`registration.json` are read-only. The historical F023 full symbolic result is
not recomputed here.

## 1. Frame, density, and the true exceptional set

Let $a_1,\ldots,a_N\in\mathbb R^6\setminus\{0\}$ span \(\mathbb R^6\).
Positive spectral weights may be absorbed in the columns: a term
\(w_j|v_j\cdot\xi|\), \(w_j>0\), is represented by \(a_j=w_jv_j\).
Let

\[
 \psi(\xi)=\sum_j|a_j\cdot\xi|,\qquad
 c_t(s)=\frac{t}{\pi(t^2+s^2)}.
\]

The law of \(\sum_j a_j Z_j\), with independent density-\(c_t\) variables,
is a positive conservative convolution semigroup. Its Fourier transform is
\(e^{-t\psi}\), its smooth density is $p_t$, and
\(G(x)=\int_0^\infty p_t(x)\,dt\) is allowed to be infinite. Define

\[
 S_2=\bigcup_{i<j}\operatorname{span}(a_i,a_j).
\]

This is also the union of spans of at most two columns, since the frame spans
six dimensions. Parallel pairs and repetitions are permitted. The set is
closed, being a finite union of subspaces. In this theorem the frame may be
arbitrary: no coordinate, sign, permutation, or generic-position symmetry is
assumed.

Choose a nonsingular six-column matrix $B$, set
\(C=\|B^{-T}\|_F^2>0\), and take

\[
 c=(1+C)^{-1}>0.
\]

Since \(\|\xi\|\le\sqrt C\|B^T\xi\|_2\le\sqrt C\psi(\xi)\)
and \(1+C\ge\sqrt C\), we have \(\psi(\xi)\ge c\|\xi\|\).
Fourier inversion consequently gives, for every $t>0$,

\[
 0\le p_t(x)\le\frac{15}{8\pi^3c^6}t^{-6}.                 \tag{1}
\]

Indeed the area of \(S^5\) is \(\pi^3\) and
\(\int_0^\infty r^5e^{-ctr}dr=120/(ct)^6\). Polynomial multiples of
\(e^{-t\psi}\) are integrable, proving the asserted fixed-time smoothness.
In particular the large-time part of $G$ is finite uniformly in $x$.

## 2. Independent-large-triple coarea estimate off $S_2$

Fix \(\delta>0\) and $x$ with \(\operatorname{dist}(x,S_2)\ge\delta\).
Set

\[
 S=\sum_j\|a_j\|_1,\qquad R=\frac{\delta}{2S}.
\]

For every representation $x=\sum_j a_j z_j$, call $j$ large if
\(|z_j|>R\). If the large columns had rank at most two, their whole sum
would lie in some component of $S_2$: choose a maximal independent subset
of the large columns and, if necessary, extend it to two frame columns.
The remaining sum has norm at most
\(R\sum_j\|a_j\|_2\le RS=\delta/2\), a contradiction. Thus the large
columns have rank at least three. Selecting a basis greedily from them gives
an independent triple $J$ all of whose three parameters exceed $R$.
This does **not** say that every three large columns are independent.

For each independent triple $J$, put $U_J=\operatorname{span}(a_j:j\in J)$
and \(D_J=\det(A_J^TA_J)>0\). The orthogonal projection onto $U_J^\perp$
of the sum from the remaining columns has a three-dimensional density
\(q_{t,J}\). Its symbol on $U_J^\perp$ is exactly the restriction of
\(\psi\), because the selected columns vanish on this subspace. It is
therefore at least $c\|\eta\|$, and

\[
 \|q_{t,J}\|_\infty\le
 (2\pi)^{-3}4\pi\int_0^\infty r^2e^{-ctr}dr
 =\frac1{\pi^2c^3t^3}.                                  \tag{2}
\]

The map from the selected three scalar parameters to $U_J$ has Jacobian
\(\sqrt{D_J}\). Disintegrate the remaining independent parameters by their
orthogonal projection, and bound each selected density on its large event by
\(t/(\pi R^2)\). Positive coarea then gives the following bound for the
pushforward density restricted to that event:

\[
 \frac1{\sqrt{D_J}}\frac{t^3}{\pi^3R^6}
 q_{t,J}(\operatorname{proj}_{U_J^\perp}x).
\]

Every parameter fiber over such $x$ is covered by the union of these events.
Summing over independent triples only and using (2) yields

\[
 p_t(x)\le \frac1{\pi^5c^3R^6}
          \sum_{J:\,|J|=3,\,D_J>0}D_J^{-1/2}
 =:M_\delta<\infty,\qquad t>0.                           \tag{3}
\]

The pushforward/coarea argument initially gives the density inequality almost
everywhere. For a fixed $x$ as above, apply it with any \(0<\delta'<\delta\)
in a neighborhood of $x$, use the fixed-time continuity of $p_t$, and then
let \(\delta'\uparrow\delta\). This proves the pointwise version (3), including
the boundary of the distance condition.

For rational columns this constant has a completely rational overestimate,
without enumerating triples. Let $d_A\ge1$ be a common denominator of all
column coordinates. For an independent triple, the Gram determinant of the
integer matrix $d_AA_J$ is a positive integer. Hence
\(D_J\ge d_A^{-6}\), and

\[
 M_\delta\le
 \frac{\binom N3d_A^3}{\pi^5c^3R^6}
 \le\frac{\binom N3d_A^3}{32c^3R^6}.                     \tag{4}
\]

Here only the elementary inequality \(\pi>2\) is used; for example
\(\pi/4=\int_0^1(1+s^2)^{-1}ds>1/2\). Combining (3) for $0<t\le1$
with (1) for $t\ge1$ proves $G(x)<\infty$ off $S_2$. These bounds are
uniform on any set a positive distance from $S_2$; dominated convergence
also proves that $G$ is continuous on \(\mathbb R^6\setminus S_2\).

The dimensional balance in (3) is three selected Cauchy tails $t^3$ against
a three-dimensional quotient peak $t^{-3}$. It is a six-dimensional bound;
it cannot be copied without changing exponents in another dimension.

## 3. Positive lower bound on $S_2$: mandatory rejection

Suppose $x\in S_2$. Choose two independent columns whose span contains $x$,
extending a one-dimensional span if necessary, and complete them to a basis

\[
 B=(a_{j_1},\ldots,a_{j_6}).
\]

Then $b=B^{-1}x=(b_1,b_2,0,0,0,0)$. For each row put
\(C_i=\sum_{\ell\notin B}|(B^{-1}a_\ell)_i|\).
Restrict all remaining scalar parameters to \(|z_\ell|\le t\); this event
has product mass $2^{-(N-6)}$. For $0<t\le1$, the first two solved basis
parameters have magnitude at most \(|b_i|+C_i\), so their Cauchy densities
are at least \(t/[\pi(1+(|b_i|+C_i)^2)]\). The other four magnitudes are
at most $C_it$, and their densities are at least
\(1/[\pi t(1+C_i^2)]\). The exact full-basis change of variables gives

\[
 p_t(x)\ge c_x t^{-2},\qquad
 c_x=\frac{2^{-(N-6)}}{|\det B|\pi^6}
 \prod_{i=1}^2\frac1{1+(|b_i|+C_i)^2}
 \prod_{i=3}^6\frac1{1+C_i^2}>0.                         \tag{5}
\]

Thus $G(x)=+\infty$ on all of $S_2$, including the origin. This proves the
exact scalar finiteness dichotomy for the present frame class:

\[
                    G(x)<\infty\quad\Longleftrightarrow\quad x\notin S_2.
                                                                    \tag{6}
\]

## 4. Exact all-pairs certificate for rational inputs

For rational $x,a_i,a_j$, write $X=x\cdot x$, $g_{uv}=a_u\cdot a_v$,
\(u_i=a_i\cdot x\), and \(d_{ij}=g_{ii}g_{jj}-g_{ij}^2\).
For every unordered pair $i<j$, compute exactly

\[
 \rho_{ij}=\begin{cases}
 X-\dfrac{g_{jj}u_i^2-2g_{ij}u_iu_j+g_{ii}u_j^2}{d_{ij}},&d_{ij}>0,\\[6pt]
 X-u_i^2/g_{ii},&d_{ij}=0.
 \end{cases}                                                   \tag{7}
\]

These are squared distances to the true pair spans. Nonzero columns ensure
\(g_{ii}>0\); a negative Gram determinant or residual is an input/arithmetic
failure. A zero residual rejects finite Green evaluation by (5). All
\(\binom N2\) strictly positive rational residuals certify off-$S_2$, with
\(\rho=\min_{i<j}\rho_{ij}>0\). For example the rational choice
\(\delta=\min(\rho,1)/2\) satisfies \(\delta^2\le\rho\) and can be inserted
in (3) or (4). This includes dependent pairs, not just pairs from a preferred
orbit or generic rank class. A certificate must bind the point and every
column, record the pair count and exact data, and witness the frame rank six.
No such full point certificate is computed in this round.

## 5. Complete finite-fan formula at positive cutoff

Assume a finite complete simplicial fan of \(\mathbb R^6\), with disjoint
interiors up to null boundaries, on each cone of which \(\psi\) is linear.
For one full-dimensional cone choose its rays $r_1,\ldots,r_6$, and set

\[
 q_i=\psi(r_i)>0,\quad
 b_i=\frac{r_i\cdot x}{q_i},\quad
 \alpha=\frac{|\det(r_1,\ldots,r_6)|}{\prod_iq_i}.
\]

For rational columns, rays, and point, these data are rational. Cone
integration in the absolutely convergent Fourier integral gives

\[
 p_t(x)=(2\pi)^{-6}\sum_{\mathcal C}\alpha_{\mathcal C}
                  \operatorname{Re}\prod_{i=1}^6(t-i b_i)^{-1}.       \tag{8}
\]

A symmetry-orbit representation may be used only after verifying the actual
frame symmetry and its precise multiplicities. In a grouped representation
the stored weight includes those multiplicities once; no extra factor is
introduced here. An ordered-chamber fan alone is not a complete fan for an
arbitrary frame.

Let $[b_1,\ldots,b_6]f$ denote the fifth confluent divided difference,
using normalized repeated jets $f^{(k)}(b)/k!$. The resolvent identity is

\[
 [b_1,\ldots,b_6](t-i b)^{-1}
       =i^5\prod_{j=1}^6(t-i b_j)^{-1}.
\]

Since $i^{-5}=-i$ and
\(\operatorname{Re}[-i/(t-i b)]=b/(t^2+b^2)\), it follows that

\[
 \operatorname{Re}\prod_j(t-i b_j)^{-1}
      =[b_1,\ldots,b_6]\frac{b}{t^2+b^2}.                        \tag{9}
\]

For a fixed \(\epsilon>0\), all required jets of the right-hand integrand
are integrable on \([\epsilon,\infty)\). Also
\(\int_\epsilon^\infty b/(t^2+b^2)dt=\arctan(b/\epsilon)\), including
negative and zero $b$. There are finitely many cones and finitely many jets,
so (8) and (9) imply

\[
 H_\epsilon(x):=\int_\epsilon^\infty p_t(x)dt
   =(2\pi)^{-6}\sum_{\mathcal C}\alpha_{\mathcal C}
      [b_1,\ldots,b_6]\arctan(b/\epsilon).                     \tag{10}
\]

This argument only integrates at a positive cutoff. It does not assert that
an individual cone has a convergent Green integral as \(\epsilon\downarrow0\).

## 6. Laurent coefficients and the scalar finite-part theorem

At a fixed nonzero $b$,

\[
 \arctan(b/\epsilon)=\frac\pi2\operatorname{sign}b
                         -\arctan(\epsilon/b).
\]

Its value has constant term \((\pi/2)\operatorname{sign}b\), and every fixed
positive-order normalized $b$-jet is $O(\epsilon)$. At $b=0$, the
normalized jets are exactly

\[
 \frac1{k!}\partial_b^k\arctan(b/\epsilon)\big|_{b=0}
 =\begin{cases}
 (-1)^{(k-1)/2}/(k\epsilon^k),& k\text{ odd},\\
 0,&k\text{ even},
 \end{cases}                                                  \tag{11}
\]

where the even case includes $k=0$.

For a fixed six-node multiset, DD5 is a fixed linear functional of the
normalized jets at its distinct nodes, through order at most five. Define
\(\Gamma_{\mathcal C}\) to be this functional applied to the algebraic data

* value \(\operatorname{sign}b\) at a nonzero node;
* value zero at the zero node;
* every positive-order jet zero, at every node.

Equivalently it is the leading coefficient of the unique degree-at-most-five
Hermite interpolating polynomial for those data. This is a constant-coefficient
functional, not a claim that sign is differentiable at zero. Equations (11)
and the nonzero-node expansion give

\[
 [b_1,\ldots,b_6]\arctan(b/\epsilon)
 =\frac\pi2\Gamma_{\mathcal C}
  +L_{1,\mathcal C}\epsilon^{-1}
  +L_{3,\mathcal C}\epsilon^{-3}
  +L_{5,\mathcal C}\epsilon^{-5}+O(\epsilon).                   \tag{12}
\]

For rational nodes all four coefficients apart from the displayed \(\pi/2\)
are rational. There is no logarithmic term and no other negative power.

If $x\notin S_2$, Sections 1–2 prove $H_\epsilon(x)\to G(x)<\infty$.
Insert the finite sum of (12) into (10). Successively eliminating its largest
negative power proves the exact identities

\[
 \sum_{\mathcal C}\alpha_{\mathcal C}L_{5,\mathcal C}
 =\sum_{\mathcal C}\alpha_{\mathcal C}L_{3,\mathcal C}
 =\sum_{\mathcal C}\alpha_{\mathcal C}L_{1,\mathcal C}=0.        \tag{13}
\]

Therefore the scalar Green value is

\[
 \boxed{\quad G(x)=\frac{\Gamma(x)}{128\pi^5},\qquad
 \Gamma(x)=\sum_{\mathcal C}\alpha_{\mathcal C}\Gamma_{\mathcal C},
 \qquad x\notin S_2.\quad}                                    \tag{14}
\]

The normalization is \((2\pi)^{-6}(\pi/2)=1/(128\pi^5)\). No numerical
small-cutoff cancellation is needed. Computing the three complete principal
channels can be an additional diagnostic, but cannot replace the complete-fan
and off-$S_2$ premises. Conversely, a constant term by itself proves nothing
on $S_2$. For the six-zero formal multiset, DD5 is exactly
\(1/(5\epsilon^5)\), while its constant coefficient is zero. It would be
incorrect to certify $G(0)=0$; in fact (5) gives $G(0)=+\infty$.

## 7. Residue form and historical F023 normalization

For node multiplicities $m_b$, define

\[
 a_{b,\ell}=[u^\ell]\prod_{c\ne b}(b-c+u)^{-m_c}.
\]

Residues of $f(z)/\prod_b(z-b)^{m_b}$, or Hermite interpolation, give

\[
 [b_1,\ldots,b_6]f
   =\sum_b\sum_{k=0}^{m_b-1}a_{b,m_b-1-k}\frac{f^{(k)}(b)}{k!},
 \qquad
 \Gamma_{\mathcal C}=\sum_{b\ne0}\operatorname{sign}(b)a_{b,m_b-1}.
                                                                    \tag{15}
\]

This identifies the constant functional with the order-one pole coefficient
used in B023/F023. Under negation of all six nodes the DD5 of the odd atan
function is unchanged, so a premerged global-negation weight must not be
multiplied by two again. The F023 report pinned in `registration.json`
independently checked its historical complete 40320-tuple Gamma and the three
principal cancellations. That old complete calculation is cited only as
historical consistency evidence; it is not rerun or claimed as a new result.

Read-only C010 source inspection confirms the normalized confluent jets and
the \((2\pi)^6\) denominator for $H$. C013 confirms the generic scalar
\(128\pi^5\) normalization, but its zero-node rejection and gradient domain
are not changed here.

## 8. Oracle contract, rejection, and derivative limitation

A future scalar implementation may evaluate (14), with exact rational
arithmetic or rigorous Arb enclosure of the finite functional, only after:

1. binding the exact nonzero rational columns, a rank-six witness, point, and
   complete rational fan with actual weights/multiplicities;
2. checking all pair distances (7) strictly positive;
3. keeping every zero and repeated node, using exact equality and exact signs;
4. applying normalized Hermite jets to (14) and enclosing its normalization.

Any failed frame/fan validation, missing pair, rank failure, or unresolved
interval division must return unsupported/unresolved rather than a value.
A point on true $S_2$ must return infinite/rejected by (5). An interval that
cannot distinguish two distinct exact nodes is a precision problem, not
permission to merge them or discard a term. Zero columns, if present in an
external input, must first be removed with an exact proof of unchanged symbol;
they are outside the stated input contract.

The interface is **scalar only**. The artificial zero higher jets in the
constant functional are coefficient extraction at a fixed point; they are
not derivatives of $G$ with respect to that point. When a point moves,
nodes collide or cross zero and the discarded Laurent terms may enter
differentiation. Neither copied DD6/DD7 rules nor the existing generic
gradient routine are certified at zero nodes by this argument. Separate
uniform regularity and cancellation proofs would be required. The continuity
proved off $S_2$ does not assert a global gradient or Hessian theorem.

For the new sparse-four frame this theorem supplies an eligible mathematical
scalar interface once the full-fan certificate and the particular point's
all-pairs certificate are bound. No such full point evaluation, new engine,
or derivative implementation is included in B052.
