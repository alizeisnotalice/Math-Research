# Dimension-free first angular response, and a second-response obstruction

Let c range in a fixed compact interval I⊂(0,∞), and let μ solve μ=1+e^(−μ). Write p_c=(1+e^c)^−1 and

\[
 \psi_c(t)=\frac{1+e^{-c+ict}}{(1+e^{-c})(1-it)},\qquad
 P_c(t)=\prod_{j=1}^b\psi_c(t_j).
\]

## First derivative

There is a constant C_I, independent of b, such that

\[
 \sup_{c\in I,\,\sum_jt_j=0}|\partial_cP_c(t)|\le C_I. \tag{1}
\]
At c=μ the sum restriction may be removed.

Indeed |ψ_c(t)|≤(1+t²)^(−1/2), and the nonvanishing numerator makes the log derivative well defined:

\[
 a_c(t)=p_c+\frac{e^{-c+ict}(-1+it)}{1+e^{-c+ict}}.
\]
The mean m(c)=1+cp_c has m'(μ)=0. Taylor's formula on |t|≤1 and the denominator bound |1+e^(−c+ict)|≥1−e^(−inf I) give

\[
 |a_c(t)-im'(c)t|\le C_Ih(t),\quad
 h(t)=\begin{cases}t^2,&|t|\le1,\\|t|,&|t|>1.\end{cases} \tag{2}
\]
For large t this follows directly from |a_c(t)|≤C_I(1+|t|). The linear sum cancels on Σt_j=0, or vanishes at μ.

Put S=Σ_{|t_j|≤1}t_j² and let m be the number of remaining frequencies. The small-frequency part is bounded by C_I S exp(−S/4). The large-frequency part is bounded by

\[
 C_I\sum_{|t_i|>1}\frac{|t_i|}{\sqrt{1+t_i^2}}
 \prod_{j\ne i,\,|t_j|>1}(1+t_j^2)^{-1/2}
 \le C_I m2^{-(m-1)/2}.
\]
Both expressions are bounded uniformly. Multiplying (2) by |P_c| proves (1), including b=1 in the unrestricted assertion at μ.

## Operator consequence

Let H_b={x∈R^b:Σx_j=0}, let Π be orthogonal projection onto H_b, and define

\[
 T_cf(x)=E f(x+\Pi(Z_{c,1},\ldots,Z_{c,b})).
\]
Plancherel and (1) imply

\[
 \|T_c-T_d\|_{2\to2}\le C_I|c-d|. \tag{3}
\]
The bounded multiplier derivative gives strong L² differentiation. The fundamental theorem of calculus in L², with a measurable absolutely continuous representative in the parameter, and Minkowski give

\[
 \left\|\sup_{c\in J}|T_cf-T_{c_0}f|\right\|_2
 \le C_I|J|\|f\|_2,\qquad c_0\in J\subset I. \tag{4}
\]
The supremum can first be taken over a countable dense set. This is for the projected **unconditioned** law, not the law conditioned on the radial sum.

## Second derivative is not dimension free

For an integer M take b=M²+1 and frequencies M followed by M² copies of −1/M. Their sum is zero. Put p=p_μ and a=p exp(−μ/2)>0. Then

\[
 P'_\mu\longrightarrow-a,\qquad
 \frac{P''_\mu}{M}\longrightarrow-ia. \tag{5}
\]

To see this, set B_c=ψ_c(−1/M)^(M²). Uniform Taylor expansions of log ψ and its first two c derivatives on a compact neighborhood of μ give

\[
 B_\mu=e^{-i\mu M-\mu/2}(1+O(M^{-1})),\quad
 B'_\mu/B_\mu=O(1),\quad B''_\mu/B_\mu=O(M).
\]
Here Var(Z_μ)=μ, the potentially O(M) term in B'/B vanishes because m'(μ)=0, and B''/B can retain −im''(μ)M. For the large factor,

\[
 \psi_\mu(M)=O(M^{-1}),\quad
 \psi'_\mu(M)=-p e^{i\mu M}+O(M^{-1}),\quad
 \psi''_\mu(M)=-ipM e^{i\mu M}+O(1).
\]
In the product rule, ψB'' and 2ψ'B' are O(1); the displayed leading large-factor derivatives therefore prove (5). Thus sup_b,t |P''_μ(t)|=∞. A first-derivative operator bound cannot justify a dimension-free second-order angular Taylor theorem.

## Conditional radial law: the remaining step

For Σξ_j=0, the angular characteristic function conditional on S=s is the Fourier-inversion ratio

\[
 \frac{\int e^{-i\lambda s}\prod_j\psi_c(\xi_j+\lambda)d\lambda}
 {\int e^{-i\lambda s}\psi_c(\lambda)^b d\lambda}.
\]
Inside the integrals the shifted frequencies sum to bλ. The zero-sum cancellation in (1) therefore cannot be applied unchanged; differentiating also changes the denominator. Numerical conditional checks do not prove a uniform multiplier bound for this ratio.

Even a conditional version of (4) would not alone settle the desired physical upper bound. A rare deficit of probability O(ρ), ρ=√(logN/N), may have L² size O(√ρ). A parameter interval of length O(ρ) then only supplies an O(ρ^(3/2)) response bound, which is larger than the required o(logN/N)=o(ρ²). The source's actual modular, the conditional law and the reuse of its angular mass across receiver cones must still be controlled together.

No statement here is an upper bound for the general strong maximal constant C_N.
