# r076: arbitrary anchor count from separated witness copies

## Result

There can be no finite-corner count bound depending only on N and the
near-corner enlargement threshold. Starting from one compact strict cap-two
witness, arbitrarily many far-separated copies preserve its score but force
at least one genuinely distinct anchor per copy. This does not make
source-budget allocation equally expensive: the natural localized allocation
has a reuse factor independent of the number of copies. Thus raw anchor
count is not the appropriate global objective; weighted source allocation is.

## Translated copies preserve the score

Fix a compact cap-two strict witness \(F_0\) and a compact strict receiver
set \(E_0\), with a finite exact-corner cover of cardinality \(M_0\). Let L
be the length of an interval containing the first-coordinate projection of
\(\operatorname{supp}F_0\), and let W be the diameter of the
first-coordinate projection of \(E_0\).

For an arbitrary positive integer K choose a separation D and put

\[
 F_j(y)=F_0(y-jD e_1),\qquad E_j=E_0+jD e_1,\qquad
 F_K=\sum_{j=1}^K F_j,\quad E_K=\bigcup_{j=1}^K E_j.       \tag{1}
\]

Take D large enough that the supports and receiver sets are disjoint.
Then \(0\le F_K\le2\), it is compact, and modular additivity gives

\[
 J(F_K)=KJ(F_0),\qquad |E_K|=K|E_0|.                       \tag{2}
\]

Strict successful rectangles translate with the copy, so every point of
\(E_K\) is still strict-successful. In particular the score
\(|E_K|/J(F_K)\) is exactly the score of the original witness. Thus any
near-optimal or prescribed packet efficiency of the initial witness is
preserved.

## A projection bound for every orientation

Let \(S\) be the first-coordinate projection of \(\operatorname{supp}F_K\).
It has measure at most KL. For any nondegenerate oriented fixed-corner
rectangle having a as one vertex and x as the opposite vertex, write
\(\ell_1=|x_1-a_1|\). Independently of its orientation,

\[
 \begin{aligned}
 H_aF_K(x)
 &\le {2\,|S\cap [\min(a_1,x_1),\max(a_1,x_1)]|\,
             \prod_{i=2}^N\ell_i
          \over \ell_1\prod_{i=2}^N\ell_i}\\
 &\le {2KL\over |x_1-a_1|}.                               \tag{3}
 \end{aligned}
\]

The first inequality only uses \(F_K\le2\) and the fact that its support is
confined to the first-coordinate set S. It is therefore valid for every
signed orientation; no positivity assumption on x or a is being used.

If this fixed-corner average exceeds \(e^{-\tau}\), then

\[
 |x_1-a_1|<2e^\tau KL.                                    \tag{4}
\]

Choose D further so that

\[
 D-W>4e^\tau KL.                                          \tag{5}
\]

An anchor satisfying (4) for one point of \(E_j\) and one point of \(E_k\),
\(j\ne k\), would put those two first coordinates within
\(4e^\tau KL\) of each other. This contradicts (5). Hence one fixed anchor,
even when all its orientations are free, can serve points from at most one
copy \(E_j\).

It follows that every near-corner cover of all of \(E_K\), at threshold
\(e^{-\tau}\), requires at least K genuinely distinct anchor locations.
Since K is arbitrary, there is no count bound depending on N and tau alone.

## Localized allocation does not grow with K

Translate the \(M_0\) exact corners of the original cover along with every
copy. For each such translated corner assign its own copy source \(F_j\).
The total assigned modular, allowing the M_0 repetitions within a copy, is

\[
 \sum_{j=1}^K\sum_{\ell=1}^{M_0}J(F_j)
 =K M_0J(F_0)=M_0J(F_K).                                  \tag{6}
\]

Thus the source reuse factor is M_0, independent of K. The copies prove
that unweighted global anchor count is not the right quantity to control.
They do not remove the intrinsic, possibly exponential, lower bound on
M_0 or on an allocation factor for one full high-efficient witness. Nor do
they prove a general upper bound.

## Exact elementary packet model

Take F_0=(3/2)1_[0,1]^N and E_0=[1/4,3/4]^N. The corner at zero has average exactly 3/2 on E_0. For a rational enlargement A>=1 and K copies, choose D=8AK+2 along the first coordinate. The projection separation margin is D−1/2−4AK>0, so the lower bound above gives at least K corners even at threshold 1/A. The K translated zero corners give an exact cover with no enlargement. Thus the minimum count is exactly K in this model, whereas using each packet as its assigned source gives reuse factor 1. This packet has low score; it is a transparent certificate model, not a near-extremal C_N construction. The general translated-copy theorem above separately applies to high-score packets.
