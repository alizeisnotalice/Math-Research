# r034 — matched reference and the affine drift compensator

Date: 2026-11-06.  This note audits the proposed matched-reference
comparison for the **max** Bellman terminal from r032.  It is an exact
one-step identity plus a local expansion; it is not a Bellman upper bound.

Put

\[
m(c)=1+\frac c{1+e^c},\qquad \mu=1+e^{-\mu},\qquad
\delta_a=\mu-m(a).
\]

Here \(m\) has its unique global maximum at \(c=\mu\).  For the fixed
action operator \(K_c\), let \(v_t^a=K_a^t\phi_N\), where

\[
 \phi_N(r)=\begin{cases}e^{-r},&r\geq0,\\1-r/N,&r<0.\end{cases}
\]

The correctly drift-compensated reference is

\[
 \widehat v_t^a(r)=v_t^a(r)+\frac {t\delta_a}{N}.
\]

## Exact one-step formula

Writing \(T v=\sup_c K_c v\), constants pass through every \(K_c\), so
for every `c` and every `t>=1`,

\[
 K_c\widehat v_{t-1}^a-\widehat v_t^a
   =(K_c-K_a)v_{t-1}^a-\frac{\delta_a}{N}. \tag{1}
\]

This is the requested exact ``gain minus running drift'' formula.  It is
not legitimate to omit the last term when comparing with the true MAX
Bellman operator.

There is a sharper cancellation.  Let

\[
 L_{t-1}^a(r)=1-\frac rN+\frac{(t-1)m(a)}N,
 \qquad w_{t-1}^a=v_{t-1}^a-L_{t-1}^a.
\]

Since \(K_cL_{t-1}^a-L_{t-1}^a=m(c)/N\), (1) becomes the exact identity

\[
 K_c\widehat v_{t-1}^a-\widehat v_t^a
   =(K_c-K_a)w_{t-1}^a+\frac{m(c)-\mu}{N}. \tag{2}
\]

Thus the entire *affine/mean* advantage is removed.  In particular, at the
mean-maximising action `c=mu`,

\[
 K_\mu\widehat v_{t-1}^a-\widehat v_t^a
      =(K_\mu-K_a)w_{t-1}^a, \tag{3}
\]

with no residual mean term; at every other `c` the displayed mean term is
non-positive.  This is independent of any Taylor expansion or stochastic
ordering claim.

On the whole left ray the identity is even more concrete.  Since all
increments are nonnegative,

\[
v_t^a(r)=1-\frac rN+\frac{t m(a)}N,\quad r\leq0,
\qquad
\widehat v_t^a(r)=1-\frac rN+\frac{t\mu}N. \tag{4}
\]

The latter is exactly the true MAX iterate on that ray.  Hence any putative
deterministic improvement obtained purely from its affine left continuation
is an artefact of using the un-compensated reference.  A remaining gain must
come from the kink/non-affine residual `w`, i.e. from actual localization.

## Size under the matched choice

Let `q=N mu-d` and choose `a=a_N` from

\[
2a-m(a)=q/N.\tag{5}
\]

Set \(x=\mu-a\).  Wolfram 15.0 was actually used to differentiate the
explicit `m`; its numerical values are

\[
 \mu=1.2784645427610737951\ldots,\qquad m'(\mu)=0,
 \qquad m''(\mu)=-\frac{\mu-1}{\mu}
 =-0.2178117057198001\ldots . \tag{6}
\]

The stationarity \(m'(\mu)=0\), rather than a merely bounded derivative,
is decisive.  From (5),

\[
 \frac dN=2x-\delta_a,
 \qquad
 \delta_a=\frac{\mu-1}{2\mu}x^2+O(x^3),
\]

and hence

\[
 x=\frac d{2N}+O\!\left(\frac{d^2}{N^2}\right),\qquad
 \delta_a=\frac{\mu-1}{8\mu}\frac{d^2}{N^2}
       +O\!\left(\frac{d^3}{N^3}\right). \tag{7}
\]

At the critical scale \(d^2=\mu N\log N\,(1+o(1))\), this gives

\[
 \delta_a=g\frac{\log N}{N}+o\!\left(\frac{\log N}{N}\right),
 \qquad g=\frac{\mu-1}{8}. \tag{8}
\]

So the final-time compensator is indeed only `O(log N/N)` in **payoff
units**.  At a terminal density of order `log N/N`, it can move a root by
`O(1)`, and cannot by itself change a `log N` root coefficient.

## Consequence for a proposed `g log N` source

Equations (1)–(4) settle the bookkeeping question but not the sign or the
size of the non-affine source.  They imply:

* the compensation does remove a linear mean term exactly, including at
  `c=mu`;
* because of the stationary identity `m'(mu)=0`, that removed final amount
  is only `g log(N)/N+o(log(N)/N)`, hence corresponds to at most an `O(1)`
  terminal-root shift;
* therefore it cannot cancel a genuine horizontal
  \(b_Nd=\Theta(\log N)\) contribution.  If a calculation replaces such a
  term by \((1-1/N)b_Nd\), the change is only `O(log(N)/N)` at this scale,
  again an `O(1)` root effect.

Accordingly, a `g log N` coefficient may still arise from the localized
quadratic/kink residual `(K_mu-K_a)w` (or an equivalent Poisson source), but
it cannot be justified from the affine mean transport.  Any sharp argument
must retain (2), prove a global localized bound for that residual, and carry
the `O(log N/N)` compensator when claiming the constant term.  This note
does **not** prove that the residual has coefficient `g`, `2g`, or either
sign after all Bellman iterations.

## Conditional local square completion (audit of the r034 algebra)

This paragraph records what the proposed local calculation does prove once
its stated density/curvature estimates are supplied.  Set
\(b=m(a)-a>0\), and write the local quadratic numerator/denominator in the
form

\[
 A=p_a\alpha+E_A,quad B=-p_a d+E_B,qquad
 \alpha=bd+\tfrac12\mu\eta f. \tag{9}
\]

The running compensator in (1) is \(-\delta_a/N\).  Wolfram series
expansion (actually run, with `a=mu-x`) gives

\[
 \delta_a=\frac{p_a}{2}b^2+O(b^3),\qquad a\to\mu, \tag{10}
\]

where replacing `p_a` by its limit changes only the cubic remainder.  Thus,
if `E_A=E_B=0`, the compensated quadratic gain is exactly, to the displayed
order,

\[
 \frac{p_a}{2}\left\{\frac{(bd+\frac12\mu\eta f)^2}{d}
                  -\frac{b^2}{N}\right\}
 =\frac{p_a}{2}\left\{b^2(d-N^{-1})+b\mu\eta f+
                    \frac{\mu^2\eta^2f^2}{4d}\right\}
   +O\!\left(\frac{b^3}{N}\right). \tag{11}
\]

In particular the cross term has its **signed** value; it must not be
replaced by an absolute value.  Under the local comparisons asserted in the
proposal,

\[
 d-N^{-1}\le f(1+C\epsilon),\qquad
 \frac{f^2}{d}\le f(1+C\epsilon),qquad
 \epsilon=|\eta|+t^{-1}, \tag{12}
\]

(11) yields

\[
 \mathrm{gain}\le
 \frac{p_a}{8}f(2b+\mu\eta)^2+
 C f(b^2+\eta^2)\epsilon+O(b^3/N). \tag{13}
\]

So the advertised square and its centering
\(2b+\mu\eta\simeq r/t-(2a-m(a))\) are algebraically correct, conditional
on (9), (12), and a denominator of the stated sign.

There is one indispensable accuracy audit.  For
\(|E_B|\le p_a d/2\), elementary expansion of
\(-A^2/(2B)\) first gives, without replacing `d` by `f`, an additional
absolute error bounded by

\[
 C\left(\frac{|\alpha||E_A|+E_A^2}{d}
       +\frac{\alpha^2|E_B|}{d^2}\right),
 \qquad \alpha=bd+\tfrac12\mu\eta f. \tag{14}
\]

Thus one must keep the `d/f` ratios through the entire moderate band; it is
not valid to assume globally that `d` and `f` are comparable.  A bare
`O(f b)` error in `A` is indeed not harmless, but the stronger statement
that *every* `b`-dependent error needs an extra `1/t` was incorrect.

In the r033 centre bookkeeping, the weaker estimates proposed by the main
calculation are already sufficient after its bridge/localization relations
are inserted:

\[
 E_A=O\!\left(f\{\rho^2+t^{-1}\}\right),\qquad
 E_B=O(\rho d+f/t),\qquad \rho^2\le C\log(t)/t. \tag{15}
\]

For the `A` terms this gives
\(\rho E_A+E_A^2/f=O(f\{(\log t/t)^{3/2}+(\log t/t)^2\})\), which is
summable in the normalized bridge calculation.  Terms of type `f b eta` or
`f b^2` therefore belong to the acceptable cubic error, unlike a bare
`f b` term.  The `rho d` component of `E_B` must be settled separately by
the r033 relation involving `Q/N`; it cannot be reduced by asserting
`d~f` on the full moderate band.  This paragraph records the required error
interface only, not a proof of those estimates.

This remains a local module.  It neither supplies (12), proves (15) for the
actual iterates, nor controls exits from the bridge band, so it is not a
global MAX-Bellman estimate.  The only conclusion used outside this audit is
the exact compensated-drift identity (1)--(4).
