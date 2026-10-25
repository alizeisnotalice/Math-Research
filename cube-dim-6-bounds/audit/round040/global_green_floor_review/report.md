# A040 independent review: PASS within the stated scope

The positive coarea S3 obstruction, the paid closed S2 tube, the exact
compact-domain reduction, and the strict-margin finite-band existence
argument pass this independent review. The global floor 1/65000 remains
unproved. There is no new C6 bound, dimension reduction, paid global cover,
or assertion that the actual Hessian of G diverges.

The source under review is the theory clone's
`upper/round040/global_green_floor_strategy/`. Its locked report, theorem,
and manifest SHA256 values are respectively

    4fc6f58973cedc177dc82db375772ae1509925f7c3171cfff7df032a3772aa87
    a23bcd3b84b1789fc63b30a6afe5fbec8b9f30053802a542755b21f5852ba47e
    ed4606e24e34c98c0e37289561974d7b7cf1a05d583af90ab126950adc0a36ef

The 28 read-only source paths and hashes, including A009, A012, A015 and
all three A038 exact inputs/raw outputs, are in `input_manifest.json`.
No author writer was executed or imported. The review recomputes the
eight registered algebra/endpoint configurations using Python Fraction,
independently of the author's SymPy/Arb calculation. It inherits the
previously reviewed full-group engine validity of the A038 outputs; it
does not rerun or independently reconstruct those full kernels.

## 1. Positive coarea, the explicit family, and all S3

For a six-column body basis A and the other columns R, the ordinary
change of variables from the six selected scalar parameters gives

    p_t(x)=|det A|^-1 integral product_i f_t((A^-1(x-Rz))_i)
                                      product_j f_t(z_j) dz.

Every factor is nonnegative. For a fixed t>0 this identity is pointwise:
the selected density is bounded and continuous, so its expectation over
the other parameters is continuous. There is no conditional probability
on a zero-measure fiber.

Our rational Gauss-Jordan calculation, without importing author code,
verifies |det A|=32 for the columns allones, flip5, flip6, flip2, flip3,
flip4. It verifies the affine identity for every b,c (by separately
checking constant, b and c coefficients), and obtains row absolute sums
(49,15,15,15,15,15) from the remaining 26 actual body columns.

On x=(1,1,1,1,b,c), 0<=c<=b<=1, restrict |z_j|<=t. Its product probability
is exactly 2^-26. For 0<t<=1/100 the three nonvanishing solved coordinates
are bounded by 149/100,13/20,13/20, and the last three by 15t. The resulting
six density factors, determinant and probability give exactly

    c*=5000^3 / (2^31 pi^6 226^3 11101 2113^2).

We checked all rational denominators. Our strict comparison does not
reuse the author's Arb pi enclosure: pi=4 integral_0^1 (1+s^2)^-1 ds<4,
and the rational number obtained by replacing pi^6 with 4^6 is already
strictly greater than 1/10^20. Thus the explicit p_t>=c*>1/10^20 is
independently certified on the whole stated triangle and time interval.

For any x in S3, choose three independent body directions spanning x
(extending a smaller independent set when necessary), then complete to
a body basis. The last three coordinates of A^-1 x vanish. On the same
remaining-parameter box, their solved coordinates are bounded by C_i t,
so their densities are at least [pi(1+C_i^2)t]^-1. The first three solved
coordinates remain bounded for 0<t<=t0, hence their densities are at
least c_i t. Multiplication cancels t^3 against t^-3 and leaves a positive
constant depending on x and the basis. This argument needs no uniform
constant over the unbounded union S3.

The independence statement is valid for body directions: for a distinct
pair v,w both agreement and disagreement coordinate blocks are nonempty.
If a third body line were av+bw, its absolute coordinate values force
|a+b|=|a-b|=1, so ab=0 and it would repeat a pair line. The frame spans
R6, as the displayed six-column basis also proves. Lower-rank points can
therefore be extended to a triple.

Consequently both integral_0^T t^-2 p_t(x)dt and integral_0^T t^-1 p_t(x)dt
are infinite on all S3 for every T>0. On the explicit triangle the band
moment is at least c*(1/a-1/min(T,1/100)). The A015 negative-Hessian
inequality is independently consistent with AA^T=32I and right inverse
A^T/32: log f_t has second derivative >=-2/t^2, giving
D_h^2 p_t >= -|h|^2 p_t/(16t^2). Integrating this particular lower bound
over all time has an infinite positive moment on S3. It yields no finite
curvature constant there. Signed derivatives may still cancel; neither
the argument nor the computations imply divergence of an actual Green
derivative or Hessian.

## 2. The tube payment and the exact remaining domain

The A012 distance proof is valid including ties and zeros. Signed block
sums are maximized by unsigned sums. For a fixed block size k, the
quantity s^2/k+(S-s)^2/(6-k) is convex in the block sum s, so one of its
extreme subset sums maximizes it. This gives the five contiguous splits
on the sorted positive face, without empty-cluster exceptions.

Our exact substitution at x=(1,1,1,1,1/2,0) verifies min q_k=1/8. Thus the
explicit S3 obstruction is not hidden within the paid small S2 tube.

The inherited A009 barrier is a lower bound on the entire cube,

    G >= K_(4/125) >= (1/(5*10^9)) (1/d2-40000)_+,

with infinity on S2. Its full positive coarea/truncation proof and saved
independent rational-corollary check were read and hash-locked; the old
eight parameter budgets were not rerun. Substitution of the new closed
radius delta0=1/125000 gives exactly 17/10^6. The direction is essential:
smaller d2 gives a larger lower bound. The boundary d2=delta0 is included.
No pair multiplicity or disjoint-tube premise is introduced.

We reconstruct A038 point3's raw dyadic G enclosure and verify

    1/65000 < G(point3) < 100033915510893/6250000000000000000
                           < 17/10^6.

The exact point also satisfies min q_k>delta0^2. Hence the paid tube is
strictly above a known retained witness and contains no global minimum.

The precise compact set is

    D={1=x1>=x2>=...>=x6>=0; q_k(x)>=delta0^2 for k=1,...,5}.

It is nonempty, closed and bounded. The rational S3 example above belongs
to D; small perturbations into the strict chamber keep all five strict
distance inequalities, so D is still five-dimensional. It includes all
unpaid boundary pieces, ties, zeros, S3 strata and intersections. The
source's phrase about retaining the complete face boundary must be read
as retaining those boundary pieces outside the paid tube; it does not
mean that paid boundary points belong to D. This is a clarification of
the exact displayed set, not an exception to its proof.

## 3. Off-S2 continuity and essential infimum

Here is an independent justification of the regularity needed, without
asserting higher Green regularity. If dist(x,S2)>=delta and x=sum a_v v,
at least three |a_v| exceed R=delta/(60 sqrt6). Otherwise subtract the
at most two large directions: the other at most 30 terms have norm sum
<=30 sqrt6 R=delta/2, contradicting the distance. Cover the parameter
fiber by the 4960 possible selected-triple events. All triples are
independent. Coarea along a triple contributes 1/sqrt(Gram)<=1. On its
event, its three density factors are <=t^3/(pi^3 R^6).

The remaining projected process in the orthogonal three-space has symbol
equal to the full symbol restricted there, since the three removed
directions vanish there. The A009 moment coercivity psi>=8 sqrt6 |xi|
therefore remains valid. Fourier integration bounds its density by
1/[pi^2(8 sqrt6)^3 t^3]. This proves the finite bound

    p_t(x) <= 4960/[pi^5(8 sqrt6)^3 R^6],

uniform in t>0 on this off-tube set. Coarea first gives the inequality
almost everywhere; to include an exact distance boundary, apply it with
delta'<delta on an open neighborhood, use fixed-time continuity, and
let delta' increase to delta. This supplies a detail implicit in A040.

The same Fourier coercivity in dimension six gives O(t^-6) globally.
Combining the constant short-time bound with the integrable large-time
bound gives an integrable dominator on every compact set avoiding S2.
Since p_t is continuous for fixed t>0, G is finite and continuous there.
This includes D's S3 pieces.

G attains m=min_D G. Signed permutations and G(ry)=r^-5 G(y) for r>0
reduce the cube to the ordered positive face; G(0)=infinity. On the
removed tube G>=17/10^6>G(point3)>=m, and elsewhere the normalized face
point lies in D. Therefore every cube point has G>=m. Conversely, a
minimizer in D is off S2 and G is continuous in its full-space
neighborhood. Every small neighborhood of any cube boundary point has
positive six-dimensional intersection with the cube. Thus for every
eta>0 a positive-measure cube subset has G<m+eta. This proves exactly

    ess inf_[-1,1]^6 G = min_D G.

No differentiation, neighborhood deletion, or lower-dimensional search
is used in this argument.

## 4. What finite-band termination does and does not prove

For n>=2, B_n=integral_(1/n)^n p_t is smooth in all space. Positivity
gives B_n increasing strictly to G on D, and both B_n and G are
continuous on compact D. Dini's argument therefore gives uniform
convergence. It can alternatively be proved by the finite subcover of
the open increasing sets {B_n>G-eta}.

If m>g*=1/65000, one fixed n has B_n>g* on D with a uniform margin.
For this fixed finite window, matching band values/gradients and its
positive moment U2 are finite on compact neighborhoods. The A015 Taylor
remainder, with loss |h|^2 U2/32, tends to zero as cell diameter tends
to zero; its translation comparison tends to one. Sufficiently accurate
outward enclosures then preserve the margin. Compactness supplies a
finite passing cover. This is an existence/conditional termination
statement, provided the value/gradient/moment evaluations can be made
rigorously accurate, the time windows are eventually increased, and a
fair rational refinement yields diameters tending to zero. It is not
an estimate of cost, a completed run, or permission to launch a grid.

For whole-root coverage, points with d2<delta0 admit a sufficiently small
cell having a single fixed closest-plane witness for all its vertices.
Distance to that one plane is convex, so the whole cell is paid. One
cannot select different planes at different vertices. The seam d2=delta0
is in D and B_n has strict margin there under the hypothesis above; band
neighborhoods cover cells crossing it. These neighborhoods and the
interior tube neighborhoods cover the compact root. All leaves must be
paid by one valid rule, with exact union and no discarded volume.

If m=g*, a minimizer is outside the higher-paid tube. At that point
every finite window B_n is strictly below g*, because the density is
positive and the omitted time intervals have positive mass. Thus a
pure finite positive-band lower certificate cannot reach g*. This is
not merely a possible exception at a measure-zero point: any finite
collection of such windows has its maximum strictly below g* there,
and by continuity it fails on a relative cube neighborhood of positive
measure. Infinite limiting arguments or new positive omitted-tail
payments are separate mechanisms. If m<g*, the target is false.
Current evidence distinguishes none of these alternatives.

Full G lower values cannot be converted to band lower values from
G>=B_n, and full grad G cannot replace the matching band gradient. No
such new band data were computed in this audit.

## 5. Saved directions: exactly three point exclusions

For each saved A038 point c we reconstruct all six derivative dyadic
intervals and multiply by w-c with w=(1,1,1,1,1,0); the first displacement
is exactly zero. The resulting independent rational outer intervals are

| Saved point | Lower bound | Upper bound |
|---|---|---|
| 1 | -1377853212569/10^18 | -172231651571/(125000000000000000) |
| 2 | -1095169392327/10^18 | -547584696163/(500000000000000000) |
| 3 | -189151799551/(200000000000000000) | -472879498877/(500000000000000000) |

All three upper bounds are strictly negative. The points and w are in
the convex ordered face, so the right-hand direction is feasible. The
saved generic nonzero-node engine hypotheses give differentiability at
each point. Each point is therefore not a local minimum on that face.
Full term counts, precision, zero-node status, all seven interval widths
and the independently reconstructed Euler interval were also checked
as saved-output integrity checks, not substitutes for the inherited
engine proof.

At a supported generic point, analytic continuity would imply some
existential neighborhood on which a fixed directional sign persists.
These three saved intervals supply no quantitative neighborhood radius,
uniform certified cell sign, or paid flow to a face. They do not justify
deleting any proposed neighborhood/cell, following a sign through S3,
or locating the global minimum on first4=1. The source's negative
neighborhood claim is understood in this certificate sense.

## 6. Actual run, limitations and reproduction

Eight configurations were registered before arithmetic. The sole
recording run used PID47091 and returned exit0, tool chunk795c52. Exact
checks and result recording took 0.013170 seconds; the printed full
process work including runtime metadata hashing was 0.031308 seconds.
Python3.15.0 and standard-library Fraction were used; executable and
Fraction-source hashes are in `run_receipt.json`. No third-party math
library, child process, G/H/p/native evaluation, optimizer, grid, MC,
MATLAB call, author writer, commit or push was used. All handles closed.

The audit's minimal read-only checker is:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/audit/round040/global_green_floor_review/check_saved.py
```

It verifies all input hashes and recomputes the same eight configurations,
then compares against `results.json`; it launches no kernel. `--record`
refuses to run again once results exist. No further configuration is
authorized or performed. The parent retains final proof integration.
The read-only saved-result replay passed as PID47296, exit0, in
0.009901 seconds (tool chunk298082), with no new configurations.

No fatal proof gap was found in the reviewed scope. The unresolved work
is the actual global floor and a finite affordable certificate with
matching band inputs or additional tail information. Public and internal
C6 upper bounds both remain 64.
