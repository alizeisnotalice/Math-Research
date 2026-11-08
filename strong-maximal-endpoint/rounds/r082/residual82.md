# r082: residual column covers, complete duality, and the degeneracy test

## 1. Complete finite column LP

Fix a cap-two source \(F\), receivers \(x_1,\ldots,x_m\), receiver masses
\(\mu_i>0\), and demands \(0\leq s_i\leq1\).  A real source column is
\[
 q=(a,G),\qquad 0\leq G\leq F,\qquad
 A_{iq}=\mathbf1_{\{H_aG(x_i)>1\}},
\]
with modular charge \(J(G)\).  Given a fixed-corner conversion factor \(K\),
the exact soft-cover-plus-residual LP is
\[
\begin{split}
 {\cal C}_{F,K}(s;\mu)
 =\inf_{w_q,r_i\geq0}\ &K\sum_qJ(G_q)w_q+\sum_i\mu_ir_i\\
 {\rm subject\ to}\quad&
 \sum_qA_{iq}w_q+r_i\geq s_i\qquad(1\leq i\leq m).
\end{split}                                             \tag{1}
\]
Equivalently,
\[
 {\cal C}_{F,K}(s;\mu)
 =\inf_{w\geq0}\left\{
 K\sum_qw_qJ(G_q)+
 \sum_i\mu_i\left(s_i-\sum_qA_{iq}w_q\right)_+\right\}. \tag{2}
\]
Thus (2) is the literal finite version of \(KB+D\).

For a finite complete menu, LP duality gives
\[
\begin{split}
 {\cal C}_{F,K}(s;\mu)
 =\max_{\phi}\ &\sum_i s_i\phi_i\\
 {\rm subject\ to}\quad&
 0\leq\phi_i\leq\mu_i,\\
 &\sum_iA_{iq}\phi_i\leq KJ(G_q)
 \qquad\hbox{for every column }q.
\end{split}                                             \tag{3}
\]
The upper cap \(\phi_i\leq\mu_i\) is exactly the polar of direct residual
payment.

For a fixed anchor \(a\), define the complete set price
\[
 c_a(S)=\inf\{J(G):0\leq G\leq F,\ H_aG(x_i)>1\ (i\in S)\}.
\]
Then the full constraints in (3) are equivalently
\[
 \sum_{i\in S}\phi_i\leq Kc_a(S)
\quad\hbox{for every feasible }(a,S),                   \tag{4}
\]
or, in the only form suitable for full certification,
\[
 \boxed{\
 \sup_{a,\ 0\leq G\leq F}
 \left[\sum_{i:H_aG(x_i)>1}\phi_i-KJ(G)\right]\leq0.
 \ }                                                     \tag{5}
\]
Checking (5) only on columns returned by an optimizer is not a dual
certificate: a restricted primal is an upper bound, but a restricted dual
need not be feasible for the complete LP.

Same-anchor fractional consolidation remains valid for full hard coverage,
through the set price \(c_a(S)\).  It does not permit one to replace a
fractional load \(wG\), \(0<w<1\), by a genuine source of cost \(wJ(G)\):
the latter must still have actual average above one.  Thus (1) makes source
reuse explicit without silently identifying fractional accounting with a
physical shared subsource.

## 2. The decisive degeneracy for physical-volume residuals

Suppose now that \(\mu\) is physical Lebesgue measure and that the
fixed-corner weak charge is valid for every real column:
\[
 |\{x:H_aG(x)>1\}|\leq KJ(G).                            \tag{6}
\]
Then the unconstrained optimization (2) is completely degenerate:
\[
 \boxed{\qquad {\cal C}_{F,K}(s;dx)=\int s\,dx. \qquad}  \tag{7}
\]
Indeed, for every weighted family, with
\(\Gamma=\sum_qw_q\mathbf1_{\{H_{a_q}G_q>1\}}\),
\[
 \int s-\int(s-\Gamma)_+
 =\int\min\{s,\Gamma\}
 \leq\int\Gamma
 \leq K\sum_qw_qJ(G_q).                                 \tag{8}
\]
Hence the objective in (2) is at least \(\int s\), while the zero-column
choice attains equality.  Dually, every \(0\leq y\leq1\) automatically
satisfies the column inequalities in (3) by (6), and \(y=1\) attains
\(\int s\).

Therefore the raw expression \(KB+D\), with residual paid at unit physical
volume and the same \(K\) used to convert columns to volume, cannot itself
be a new source-budget upper-bound mechanism.  It merely repackages the
soft volume \(\int s\).  This is an exact obstruction, not a numerical
heuristic.

## 3. Two nondegenerate formulations

### A. Capacity-constrained source allocation

Keep residual in physical volume but impose an independently justified raw
source capacity \(L\), rather than charging \(KJ\) in the optimization:
\[
 {\cal R}_{F}(L;s)=
 \inf_{\substack{w,r\geq0\\ \sum_qw_qJ(G_q)\leq L}}
 \left\{\sum_i\mu_ir_i:
 \sum_qA_{iq}w_q+r_i\geq s_i\right\}.                  \tag{9}
\]
For a finite complete menu, its exact dual is
\[
 {\cal R}_{F}(L;s)=
 \max_{\substack{0\leq\phi_i\leq\mu_i\\\tau\geq0}}
 \left\{\sum_i s_i\phi_i-\tau L:
 \sum_iA_{iq}\phi_i\leq\tau J(G_q)\quad\forall q\right\}. \tag{10}
\]
This is nontrivial because \(L\) is fixed externally.  A useful theorem
would be an independently proved capacity allocation
\[
 L\leq\Lambda_NJ(F),\qquad
 {\cal R}_{F}(L;s_F)\leq R_NJ(F),                        \tag{11}
\]
which would yield the physical estimate \(KL+{\cal R}_F(L;s_F)\).
If instead one minimizes over \(L\) with coefficient \(K\), (7) returns.

The independently justified source-capacity estimate in (11) remains missing.
The unconstrained physical-volume LP does not supply it. No general endpoint
upper bound is claimed by this note.
