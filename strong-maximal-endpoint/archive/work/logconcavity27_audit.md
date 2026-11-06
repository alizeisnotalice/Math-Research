# Block Gamma--binomial log-concavity certificate: audit

Date: 2026-11-06. Verdict: accepted for the stated common-parameter density
family; no blocking issue found.

## Exact polynomial certificate

On `j <= x <= j+1`, the script constructs
\[
 P_{b,j}(x)=\sum_{k=0}^{j}\binom bk(x-k)^{b-1},\qquad
 Q_{b,j}=P_{b,j}'{}^2-P_{b,j}P_{b,j}''.
\]
The power-basis construction expands each `(x-k)^(b-1)` correctly.  The
`shift(q,a,width)` routine is the exact coefficient expansion of
`q(a+width*t)`.  For a degree-d polynomial `sum a_i t^i`, the conversion
used is
\[
 \beta_k=\sum_{i=0}^k a_i\binom{k}{i}/\binom{d}{i},
\]
which is the correct power-to-Bernstein formula on `[0,1]`.  Nonnegative
Bernstein coefficients therefore certify `Q>=0` on each listed interval.

There are 4+5+6+7=22 finite pieces in the JSON, covering `[0,b]` for
each b=4,5,6,7.  The reported tail polynomial is `Q(x)` shifted by x=b+u;
its all-nonnegative power coefficients certify the whole `[b,infinity)`
tail, not merely one extra unit interval.  Shared integer endpoints are
covered from both sides.  The added term at an integer knot is a multiple
of `(x-j)^(b-1)`, so for b>=4 both P and P' glue continuously.  Also
P(x)>0 for x>0 by the k=0 term.  Hence Q>=0 implies `(log P)''<=0` on
each positive piece and the continuous first derivative joins the pieces
into a globally concave log P; multiplying by exp(-x) only adds a linear
term to the log density.  The zero left support is compatible with the
usual extended log-concavity convention.

## Independent checks

Wolfram 15 `Resolve` independently returned `True` for every finite
interval and every tail at b=4,5,6,7.  It returned `-4` for the stated
b=3 right-piece defect at x=1.  This agrees with the JSON certificate.

## Scope of the induction

Provided the common-c parameter family is explicitly normalized so that
the order-b density is closed under convolution in b, the base orders
4--7 imply every integer b>=4: b=4,...,7 are base cases, b=8,...,11 are
4 plus a base case, and induction adds an order-4 factor.  Convolution
preserves log-concavity.  This proof does not extend to unequal c_i; that
restriction should remain part of the theorem statement.
