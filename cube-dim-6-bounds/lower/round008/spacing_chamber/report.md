# Round 008 D: exact second-spacing chamber and event candidate

The first examined future endpoint event produced a strictly improved
primary candidate. With positive u^6=2,v^6=3, set

```
a1=1+u^4/2,
a2=u^2*v^2/3+u^2*v^4/6,
a3=a4=a5=a6=(1+u)/2.
```

The complete static witness ratio is

```
R_event=2.067397333267381527719117937221768242786600566397580681... .
```

The full rational lower endpoint is strictly greater than E007's complete
Rideal upper endpoint. The certified lower margin begins
`0.0004731950858685122932502951912473057661876653160581554728...`.
These decimals are displays of the actual exact computation; the full
36 field coefficients, rational interval and Fraction comparison are in
`chamber_results.json`.

Search stopped immediately after this one promoted primary candidate.
E independently passed polynomial-union certification and the internal
rational-parameter static cross-check for the supplied family. General L1
acceptance and final C6 registration remain with E/main/F. This report does
not announce a completed new C6 bound on its own.

## Current input and precise family

The comparator is E's `audit/lower007/candidate_number_field/
candidate_results.json`, SHA-256
`b5661c0eb23c264286085053b3c59ab9f992918dcb1dc4a79fdb01e1ddcb6b81`.
Its complete coefficients and exact interval are retained in this round's
result. The rounded strict lower bound `2.06692413818` and the old main
ledger were not used as the candidate comparator. The already accepted
1e-6 step was not re-evaluated with a union calculation.

Write s=a2 and start from

```
s_current=(1000001/1000000)*(u^3+u^4*v^3)/4
         =1.040919250010784167967237171657210153855196010194997393... .
```

All other spacings and all atom masses stay fixed. The family consists
of every standard central consecutive block with lengths in {1,2,3,4}^6.
For each length vector, M=product(lengths) and the one common radius is
`r=M^(1/6)/2`. The positive-half-cell center intervals come directly from
intersecting the cube-containment slabs with `[0,a_j/2]` in every axis.
Six reflections give the usual ratio `64*half_cell_union_volume` with
period source mass one. No independent rectangle radii are introduced.

All 4096 length vectors were examined. Of these, 3866 have a nonpositive
width in a fixed coordinate and are empty for the whole parameter line;
103 have a nonpositive moving-coordinate width at s_current; 127 are
nonempty raw witnesses. All 127 raw boxes participate in event analysis
and the interface. No reduction to the 96 contained-maximal boxes is used.

There are no missed births when s increases. In coordinate 2 an odd
length l has interval

```
[0,min(s/2,r-(l-1)*s/2)],
```

and an even length l has interval

```
[max(0,l*s/2-r),s/2].
```

For l>=2 and s>0, either interval has positive width exactly when
`(l-1)*s<2r`. For l=1 the width is always positive. Thus a currently
empty moving-coordinate interval remains empty for larger s, while fixed
coordinate emptiness cannot change. This analytically accounts for all
necessary births in the registered family; no length-5/6 extension was
repeated.

## Exact endpoint chamber

The raw family supplies 37 distinct unclipped affine endpoint functions,
including the zero and s/2 cell faces. Their pairwise crossings have
rational slopes in the absolute parameter s, so each crossing lies in
Q(u,v). A finite algebraic inventory contains 217 distinct locations to
the right of s_current. This inventory is not a 217-candidate volume scan.

Only the nearest event was examined. Lexicographic left/right evaluation
of the actual min/max clipped endpoints shows that it changes the
effective raw endpoint order; it is not a crossing of two unused clipped
branches. The event is

```
r_324-s = 2s-r_2304,
r_324=u^2*v^4/2,
r_2304=u^2*v^2,
s_star=(r_324+r_2304)/3=u^2*v^2/3+u^2*v^4/6.
```

For example, the odd upper endpoint belongs to length vector
`(2,3,2,3,3,3)` and the even lower endpoint to `(3,4,3,4,4,4)`;
the other three tail permutations of each are also present. Both sides
of the event have 127 positive raw boxes. There is no box birth/death
or clipping-branch switch at this event, only the endpoint order exchange.

This is the first actual **raw-arrangement** event to the right of the
current parameter. It has not been proved to be the first kink in the
final union-volume function, since overlap can hide a raw arrangement
change. It is not asserted to be a local or global optimum. The stop rule
requires stopping at the first accurately improving candidate, rather
than continuing across the event to seek an optimum.

The exact chamber used for certification is

```
s(t)=s_current+(s_star-s_current)*t,  0<=t<=1,
s_star=1.042497096645529320820680798849419545849308095334345178...,
s_star-s_current
 =0.001577846634745152853443627192209391994112085139347784934... .
```

`chamber_interface.json` exports the 127 boxes, their lengths, complete
36-rational-coefficient affine endpoints, and a per-axis ordered partition
of all used endpoint IDs. Identical affine functions are globally
deduplicated; no different affine functions are called equal. For each
of the 120 adjacent ordered pairs, the gap is nonnegative at t=0 and t=1
and is positive at at least one endpoint. Hence it is strictly positive
throughout 0<t<1. Exactly one adjacent pair touches at the right boundary.
All orders hold on the closed interval with boundary equality allowed.

The static candidate calculation reconstructs the source-containment
boxes at t=1 and checks 3048 center/source endpoint inequalities. It finds
127 positive raw boxes and 117 distinct geometric boxes. A recursive exact
slab union of all distinct boxes gives the displayed result. No overlapping
raw volumes are added independently.

Only coordinate 2 varies, so on a fixed-order chamber each slab width in
that coordinate is affine while all cross-sectional volumes are constant.
Consequently the full ratio is affine on the chamber closure. The static
event answer was checked coefficient by coefficient against

```
R_event=R_current+g*(s_star-s_current),
g=R007_directional_derivative / ((u^3+u^4*v^3)/4)>0.
```

The division is performed by an explicitly verified field inverse. Writing
`w=u*v^3`, w^6=54 and

```
((u^3+u^4*v^3)/4)^(-1)
 =-(2*u^3/53)*sum_{k=0}^5 (-w)^k.
```

Thus the ratio is strictly increasing throughout this certified chamber;
the event endpoint is the largest value on this chamber, without any
claim about parameters beyond it. The complete static certificate does
not depend on trusting the inherited derivative answer: their equality
is an additional consistency check.

## Whole-interval common-radius check

`check_interface_geometry.py` independently reads the exported affine
interface. At both t=0 and t=1 it checks all 127 boxes against their
source-block endpoint indices, their half-cell faces, and the same
radius r in all six coordinates: **6096 exact source endpoint checks**
passed. Every center/source coordinate difference and every cell face is
affine in t, so validity at the two boundary parameters implies validity
for every t in [0,1]. This check makes zero union calls.

The independent E interface is supplied with normalization 64, candidate
t=1 and an internal static validation parameter t=1/100. That validation
is a cross-check of the same parameter family, not another search point
or a second promoted input.

E subsequently reported that the independently reconstructed polynomial
has degree one, its constant coefficient exactly equals E007, all 120
order gaps pass with the one expected right-boundary tie, and its t=1/100
value agrees with a separate static union in all 36 field coefficients.
The event value and comparison with E007's full upper endpoint also pass.
The separate evidence is in E's `audit/lower008/parametric_union/`; these
independent calculations were not repeated by D.

## Finite resources, provenance and remaining acceptance

The candidate run used PID **21034**, exec session **21654**, exited **0**,
and took **9.841197415982606 seconds**, with **3579** cached union states.
It examined **1** event candidate and computed **1** complete new static
union. Registered limits were 12 examined events, one full candidate union,
1,000,000 states and 1200 seconds. The promotion stop was reached before
any resource limit. No event after s_star was evaluated.

The separate whole-interval geometry checker used PID **21356**, exited
**0**, and took **0.2883514999994077 seconds**. All processes have ended.
There were no MATLAB calls, child agents, commits, pushes, or ledger edits.
Earlier round007 and round006 files were read-only inputs.

Key SHA-256 values:

* Current E007 certificate: `b5661c0eb23c264286085053b3c59ab9f992918dcb1dc4a79fdb01e1ddcb6b81`.
* Exact field source: `9072275d47bd66c3ba9ec7de0c435efdba0785c39a02b122ea1ebb960403ec82`.
* Executed event/certificate source: `bbff141d3c2deed2a165786c885752942aa3ef5a1065c7a2b19a87de268befb2`.
* Exported chamber interface: `99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933`.
* Geometry checker: `601465d11f88ca8e56f14849c6932cd62d3d8e6a9d8880795962027dd837be31`.

Reproduce from this directory:

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B chamber.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_interface_geometry.py
```

Still pending at D handoff: explicit finite-L1 truncation/smoothing,
F's independent acceptance, and main's final threshold registration.
No analysis past the first event is required for this handoff or claimed
as completed.

Under the latest user steering, upper-bound work has priority. D has
finished only this current round008 candidate and report; all lower-bound
processes remain stopped. No new scan, direction, round or child agent
will start without a new low-frequency batch from main.
