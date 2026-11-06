# r023 independent certificate audit — 2026-11-06

Read-only audit of r023 report and certificate code; no repository/Git change.

Accepted. The certificate proves the stated bounded-amplitude comparison only.

The two-level input has values 4 and 16, mass 244, and cost `244+968 log(2)`.
Each receiver average is at least rational `tau=10001/10000`.

The strip candidate formulas are the valid `[l,r]`, `[l,y]`, and `[y,r]`
averages, and the corner formula is exactly the stated mass/rectangle bound.
On every 1/64 cell the chosen Fraction Möbius function is pole-free and
monotone, so its endpoint minimum supplies a continuous lower receiver area.
The interior, strip, and corner regions are disjoint modulo boundaries.

`math.ceil` on `Fraction` uses an exact integral protocol; all geometry and
sums are Fraction operations, and static compilation passed. The
constant-rectangle ratio is increasing for A>1, so the cap maximum is A=16.
Wolfram independently verified its derivative identity and the atanh log2
identity. The positive log2-series bound has the right directed inequality
orientation, and the recorded rational gap is positive.

No blocker found in continuous receiver validity, strip partial intervals,
orientation/disjointness, or logarithm directions.

## Addendum: near-unit section

Accepted. For `f_eps`, the stated right-strip mass divided by its rectangle
area is `(4+eps G(y))/x`, with
`G(y)=4-(1-y)/(4-y)` below one and `G=4` above one. Taking
`x-4=(1-eps^2)eps G` gives excess average
`eps^3 G/(4+(1-eps^2)eps G)`, at least
`3 eps^3/(4+4 eps)`; the interior average is larger. Reflections/exchange
give four disjoint strips.

Their total integral is `15+3 log(4/3)`, so the announced volume and the
actual cost `1+15(1+eps)(1+log(1+eps))` are correct. Symbolic Wolfram checks
gave zero residual both for this integral and for the first-variation gap
`(6 log(4/3)-1)/8`. The Fraction certificate uses a lower log for witness
volume, upper log for cost, and upper/lower logs in the monotone constant
rectangle ratio in the correct directions; its positive epsilon=1/100 gap is
therefore rigorous. No scope issue beyond fixed dimension two.
