# A008: G is not convex even inside one strict ordered chamber

The proposed chamber-convexity shortcut is **false** for the existing body32
Green model. A strict midpoint counterexample occurs entirely inside

    x1=1 > x2 > x3 > x4 > x5 > x6 > 0.

The batch stopped at the first negative diagnostic: one center, two directions,
and **five** new generic rational G points. The negative triple was then
evaluated independently with exact fmpq barycentric divided differences and
Arb384. No further directions or centers were tested, no MATLAB was called,
and round007 was not changed. This is a rigorous obstruction to the new
convexity hypothesis, not a C6 bound.

## Explicit counterexample

Let D=100000000700 and

    c=(100000000700,85327182900,69128743300,
       51738912700,30728176300,11736192700)/D,
    v=(0,0,0,1,-1,0),  h=1/100,
    x_minus=c-hv,  x_plus=c+hv.

Thus only the fourth and fifth coordinates vary. Their numerators are

|Point|x4 numerator|x5 numerator|
|---|---|---|
|x_minus|50738912693|31728176307|
|c|51738912700|30728176300|
|x_plus|52738912707|29728176293|

Every point has x1=1 and strictly decreasing positive remaining coordinates.
The chamber is convex, so the whole segment stays strictly inside it. The
midpoint relation c=(x_minus+x_plus)/2 and all ordering inequalities were
checked with exact rationals. This counterexample does not cross a permutation
boundary or S2 singular stratum.

The actual Arb384 point values are near

    G(x_minus) = 4.42986822336336743436088480564e-5,
    G(c)       = 4.41190672415134693304461916148e-5,
    G(x_plus)  = 4.39329352480233763027576886798e-5.

The convexity residual has the wrong strict sign:

    J := [G(x_minus)+G(x_plus)]/2-G(c)
       = -3.25850068494400726292324666610...e-9,

with the exact-rational enclosure

    -163/50000000000 < J < -13/4000000000 < 0.

Consequently G(c) is strictly greater than the endpoint average, contradicting
convexity on this one ordered chamber.

## Independent exact arithmetic

The C++ diagnostic uses the established full generic cone formula, a triangular
DD5 recurrence, and this platform's binary64 `long double`. Its values were
used only to select the candidate. The first direction e6 had a positive
diagnostic Jensen gap approximately 5.46805e-9. The next direction e4-e5 gave
approximately -3.25850069e-9, causing an immediate stop at five points.

For independent certification, `certify_counterexample.py` evaluates

    G(x) = [sum_(complete group and 21 cones) alpha_R DD5(sign;b_R)]
             /(128*pi^5)

with exact fmpq arithmetic. Unlike the C++ recurrence, it uses the barycentric
identity

    DD5(sign;b) = sum_j sign(b_j)/product_(k!=j)(b_j-b_k).

The identity sum_j 1/product_(k!=j)(b_j-b_k)=0 allows summing twice the positive
class or minus twice the negative class. All node signs, nonzero tests, and
distinctness tests are exact. Each of the three certified points uses the
complete **967680** group/cone terms, including zero contributions; no
structurally failed term is skipped.

Each point therefore has G(x)=R_x/pi^5 with an exact rational R_x. The
certificate stores all three rational numerators/denominators and the exact
negative numerator/positive denominator of

    (R_minus+R_plus)/2-R_center.

Thus the strict negative sign is already an exact rational result; Arb is
needed only for pi and the displayed magnitude. `check_certificate.py`
independently reconstructs this rational residual, verifies geometry, and
checks three scalar barycentric/triangular DD anchors. It does not claim a
second full-kernel replay.

Exact certification completed normally in **55.929112958 seconds**, execution
session **86533**, Python PID **19132** (also recorded in the certificate).
Individual point times were 19.09, 16.62, and 20.20 seconds. The certified
Jensen ball is

    [-3.258500684944007262923246666104441680077647091112842275249184242320961707124745597479711881935095043e-9 +/- 1.74e-109].

Input fan, C++ source, reference formula, point table, and scripts have locked
SHA-256 values. The adaptive registration records that additional direction
definitions were appended after the first nonnegative diagnostic and before
the second direction was evaluated; no center or previously evaluated point
was changed.

## Scope and stopping status

This excludes convexity of the full body32 G on the entire ordered positive
chamber. It does not exclude smaller convex subregions, a different potential,
or a different global lower-bound argument. It does not improve global g, I,
or either C6 endpoint. The finite test is complete and all calculation has
stopped; no grid or additional search was started after certification.

## Reproduction and retained artifacts

Use the existing Python environment. From this owned directory:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python prepare_probe.py
    c++ -O3 -std=c++20 /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round003/direct_kernel/direct_kernel.cpp -o direct_probe

Run the diagnostic binary from the original read-only direct-kernel directory,
passing the absolute path to this directory's points.txt, redirecting output
to this directory's probe_results.json. Then run from this directory:

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python probe_more.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python certify_counterexample.py
    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_certificate.py

`probe_registry.json` is the final diagnostic ledger; `negative_triple.json`
is the exact certification input. `exact_certificate.json` contains the
three actual fmpq results, and `checked_certificate.json` contains geometry,
exact residual reconstruction, scalar anchors, and hashes. The compiled
binary, adaptive scratch point file, initial duplicate diagnostic output,
terminal log and Python caches are unnecessary after verification and removed.
No commit or push was made.
