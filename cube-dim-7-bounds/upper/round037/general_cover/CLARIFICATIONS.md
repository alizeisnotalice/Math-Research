# Clarifications to the frozen R037 proof

The frozen `PROOF.md` SHA256 is `6261341afbdafe868616e3e365f3484e0dab7e9e8e6ca09f83c56c93751df90e`. This file adds notation and scope clarifications only; it does not alter that proof or its conclusions.

## 1. The cap integral in Section 4

The integral denoted I_n(t) there is specifically

\[
I_{\rm cap}(t)=\int_1^t\ell_t(\rho)^{n-1}\frac{d\rho}{\rho},
\qquad\ell_t(\rho)=1+\frac1{2t}-\frac1{2\rho}.
\]

Its accepted bound is

\[
I_{\rm cap}(t)\le\log_+(t/n)+8+\log2.
\]

Source: `cube-dim-7-bounds/revision_v3/asymptotic/proof.md`, Section 1, equation (1), and the corresponding integral-uniform lemma in `cube-dim-7-bounds/paper/main.tex`. These are the cap-shell integral and its bound used in the original scalar budget Bn(t). They are not the two-point capture-capacity integral also denoted I_n(t) in other R035 documents. The two objects must not be interchanged.

## 2. The spike diagnostic is at a fixed level

In Section 5, fix n, total mass Mf>0 and a level λ>0. Let f≥0 have that mass and support in (−εf,εf)^n, and put R0=(Mf/λ)^(1/n)/2. Then

\[
(-[R_0-\varepsilon_f]_+,[R_0-\varepsilon_f]_+)^n
\subset\{M_nf>\lambda\}
\subset(-R_0-\varepsilon_f,R_0+\varepsilon_f)^n,
\]

up to irrelevant boundary choices. For the lower inclusion, at ||x||∞<R0−εf choose a radius strictly between ||x||∞+εf and R0, which contains the entire support and has average >λ. The upper inclusion is the one proved in Section 5. Therefore

\[
\left(\frac{[R_0-\varepsilon_f]_+}{R_0}\right)^n
\le\frac{\lambda|\{M_nf>\lambda\}|}{M_f}
\le\left(\frac{R_0+\varepsilon_f}{R_0}\right)^n
\longrightarrow1.
\]

This is a fixed-level limit with n,Mf,λ held fixed. It is not a claim that the whole weak-(1,1) norm of a family tends to 1 under homothetic support shrinkage. That norm is dilation invariant. The fixed-level observation is sufficient for the intended diagnostic: a common-star high-W history can have genuine common witnesses without its large pointwise W implying a comparable Cn lower bound.

