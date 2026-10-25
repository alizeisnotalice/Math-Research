# C2/012: the augmented frame survives two exact point tests

The sole registered new candidate

    psi_1(xi)=sum_body32 |v.xi|+sum_(i=1)^6 |xi_i|,
    T=2/45

passes the necessary pointwise tests at both already known B010 points.
Both certified K_(1,T) values are strictly above 1/90000, the proposed
sufficient global target for 62.5. This is **not** a global lower
certificate, a global-minimum computation, or a new C6 bound. The
established upper bound remains **64**.

The batch is complete and stopped: two existing spatial points, two new
G passes plus two new H passes, no other beta or time, and no further
search. All four passes retain the full 967680 group/cone contributions.

## Exact candidate results

Let x^(1),x^(2) be the precise rational B010 points stored verbatim in
`input_manifest.json` and `results.json`. Point 1 has denominator 10^12;
point 2 has denominator 50000000350000000000000. No input was rounded.

|Quantity|Point 1|Point 2|
|---|---:|---:|
|G_beta, beta=1|1.7545531224790340e-5|1.7495289527864908e-5|
|H_beta,2/45|6.3519154102431369e-6|6.3392272306223105e-6|
|K_beta,2/45|1.1193615814547203e-5|1.1156062297242598e-5|
|K minus 1/90000|8.2504703436092e-8|4.4951186131486e-8|
|K minus 1/92160|3.4292137010276e-7|3.0536785279815e-7|
|Formal point cost T/(64 K)|62.0393317003024|62.2481683896735|

The displayed decimals summarize outward certificates; they are not
the certificate inputs. The independent saved-data check proves the
short rational enclosures

    1119361581454720/10^20 < K(x^(1)) < 1119361581454721/10^20,
    1115606229724259/10^20 < K(x^(2)) < 1115606229724260/10^20.

The strict below-64 necessary threshold is T/4096=1/92160. The separate
global target m0=T/4000=1/90000 would imply 125/2 only if established
almost everywhere on the whole cube. Both points surviving these tests
does not supply that quantifier. In fact m<=K(x) at a continuous face
point, so the displayed formal point costs are lower bounds for the
cost based on the true essential minimum, not operator upper bounds.

## What was proved for the new model

`theorem.md` checks the positive 38-variable Cauchy construction, the
symmetric conservative Markov semigroup, mass ||K_beta,T||_1=T,
distributional identity L_beta K_beta,T=p_beta,T-delta_0, and the exact
dilation/time-average identity. Its direct maximal-ergodic argument
gives weak (1,1) factor one for general L1 and all positive common cube
radii. Consequently a global K>=m gives C6<=T/(64m).

The proof uses a general greedy-rank argument to establish G_beta
finiteness and continuity off the new union S2_beta of pair spans.
It explicitly allows dependent distinct triples such as v-w=2e_i.
The six distinct positive coordinates of both candidate points exclude
body/body, body/axis and axis/axis pair spans. Hence their G-H subtraction
is legitimate and their point values control the essential minimum.
No body32-only g, I, singular-set classification, or old Green point
certificate is inherited as a new-frame numerical bound.

F012 separately reported PASS for the general positive spanning-frame
theory, including the full S2 singularity statement. Its report hash is
df66b3907fecf8f0f31a7b9a41636a449d0703627acf4573cb435484694bd87a.
F also read the new fan structure and hashes, but did not take ownership
of a second kernel replay. The parent retains final integration and
proof responsibility.

## Exact fan transfer and finite audit

The original rational chamber-cell certificate and its independent
verifier are the geometric source. `audit_frame.py` actually replayed
that verifier into this owned directory: all 18 arrangement vertices,
4096 sign words, 21 cells and section volume 1/120 agreed exactly.
It then checked all 126 ray occurrences (18 unique primitive rays):

- each ray is nonnegative and ordered and matches its original rational ray;
- each old rate equals the explicit body32 symbol;
- each new rate equals old rate plus the coordinate L1 norm and, independently,
  the explicit sum over all 38 columns;
- all 38 linear forms have constant weak sign on each base cone;
- the determinants remain exact and unchanged; the new frame Gram is 33I.

The same 46080 signed permutations therefore give the same complete
geometric cover with new rates and weights. The new fan SHA-256 is

    450caf04d124334474b91f3d9fdabce4eadf3259ea47def9094ae5eb056b6c4e.

The original explanatory scaling parenthesis was arithmetically wrong:
.048/(1+1/12)=72/1625. The explicit registered T=2/45 was kept throughout.
The corrected exact explanation is .048/(1+2/25)=.048/1.08=2/45.
This correction changed no parameter, point, result, or batch history.
The time choice is a scale heuristic, not an optimization claim.

## Actual arithmetic, range handling and limits

G was computed with the complete exact fmpq Green/Laurent implementation
from B010. Only its five function definitions were loaded through the
AST; its executable top level was never run. Generic, zero-node,
confluent, and three principal-part channels were retained. Scalar
confluent anchors were checked independently by Hermite interpolation.
The two actual points are generic for this new fan: zero and repeated
node counts are zero, and the full summed principal coefficients are
exactly zero. Thus the new point computations do not themselves test a
nontrivial principal-part cancellation.

Point 1's H used the existing C010 native Arb512 engine. Source, binary,
runtime-library hashes and all input integer bounds were checked before
the call. Its maximum node denominator is 170000000000000 and the
maximum dot numerator bound is 12000000000000, safely within int64;
node cross-products fit int128. The complete exact dyadic output is
retained in `point1_native_H.json`.

Point 2 does not fit the native integer limits. Its H used the existing
Python arbitrary-integer node and confluent Arb evaluator at 512 bits,
with the new fan. Its exact input was neither truncated nor rounded.
The evaluator's legacy returned field is named H_arb384; the actual
context was reset to 512 after import and is explicitly recorded. Its
100-digit enclosing serialization is the source of the wider, still
decisive final point-2 interval.

The registered 1200-second kernel stage completed normally in
79.440336875 seconds, Python PID 28399, terminal session 27644, exit 0.
Point-1 G/H took 22.6482/1.27366 seconds; point-2 G/H took
45.8118/9.58629 seconds. The finite fan audit took 0.47350 seconds.
There were no failed kernel calls, skipped terms, replacement points,
extra times, Monte Carlo evaluations, grids, or MATLAB calls.

`check_results.py` then ran at 768 bits with **zero kernel calls**:
it reconstructed both exact rational/pi^5 Green values, independently
decoded the native dyadic H endpoints, used the full serialized
point-2 H enclosure, checked 30 input hashes, checked the complete group
counts and all target signs, and proved the short rational intervals
above. This is an independent saved-data check, not a second full
kernel certification.

## Exact old/new comparison and remaining gap

At the same exact two points, the old body32 model at T=4/125 had K
approximately 7.8151629735836953e-6 and 7.7858856970127892e-6, with
formal point costs 63.9781923537702 and 64.2187696374524 respectively.
The second point excluded that old fixed-time candidate. The new
beta=1,T=2/45 candidate instead survives both points. This comparison
changes both the frame and the time and does not isolate the effect of
the added axes. All old values were reconstructed from saved exact
certificates; no original-frame kernel was recomputed.

The remaining requirement is a global essential lower certificate for
the new frame at its fixed time. Neither survival at the two points
nor the new cone formula supplies that certificate. No conclusion is
made about other beta values, other times, or a global minimizer.

## Reproduction and retained material

With the recorded Python environment, run in this directory:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python audit_frame.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python evaluate.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_results.py

The first command is finite algebra; the second is the bounded kernel
batch and should only be rerun when a replay is explicitly needed; the
third rechecks saved results without kernels. Manifests lock precise
inputs and original sources before each stage. The existing native
engine build record and dependency hashes are referenced rather than
duplicating its binary or headers. Only the proof, registration, final
code, exact new fan, manifests, complete results, saved-data check and
delivery receipt are retained. No scratch, compiled binary, or Python
cache was created here. Nothing was staged, committed, or pushed.
