# C055 stable precheck: one sparse4 joint scalar probe, not launched

**PRECHECK ONLY; FULL NATIVE CALLS0.** The exact new-frame point is
x0=(50,50,50,45,27,0)/50, beta=1/80 and the only proposed T is32/825.
No G/H/K value or positive-time kernel has been computed. The current compiled
driver has run only ABI and complete-input validation modes. A separate main
release after B055 acceptance is mandatory before its sole full invocation.

## New exact true-S2 certificate

This is a new actual152-direction certificate, not the unexecuted053 plan or
the original body32 point's geometry. All152 nonzero distinct unoriented
directions and11476 unordered pairs were checked exactly. Their positive
spectral scaling does not change the spans. The exact results are

    min line distance²=2991/10000,
    min pair distance²=363/4375>0.

The unique pair minimum is between sparse4 directions(1,1,1,1,0,0) and
(1,1,1,0,1,0). A second fmpz3×3/2×2 Gram determinant calculation independently
checks this witness. There are no zero distances. A full-rank six-body-column
witness has determinant−32. Thus the point is strictly off the **new true S2**.
The complete pair-distance stream SHA and all witnesses are retained in
S2_certificate.json; primary checking took0.02196287500555627seconds.

B052 proves for arbitrary spanning finite positive Cauchy frames that scalar
G is finite exactly off true S2, and that on a complete finite fan its value
is the complete epsilon^0 coefficient of H_epsilon. F050's complete fan audit
and F052's final theory PASS are now pinned in dependency_pins.json. In
particular formal constants on S2 must NOT be interpreted as finite Green
values; the negative fake finite part in the F052 axis example is not a
permitted fallback. The native input is hard restricted to the x0 signed
orbit, excluding such different points; launch_once.py additionally checks
the pinned152/11476 certificate and strict positive minimum before creating
the native release token.

## Complete input and orbit accounting

native_input.txt binds478 integer-scaled rays, all13229 cones and1920 signed
image representatives. It is448661bytes, SHA
`610b5e73a07e354ab644ef27c2653bdf39ec9341a04dc3a2577cb502c8034174`.
Its canonical data digest is exactly the frozen053 converted-data digest
`bbe9274e5b235d5c7f94b8c226ca124bd2c1f778ba012a124440cdf711cd47b3`.

The new writer generates distinct absolute multiset permutations and32 signs
of the five nonzero coordinates, rather than replaying053's46080group loop.
It obtains3840 images and1920global-negative representatives, each with weight24
already including stabilizer12. Negation uses the lexicographic first NONZERO
coordinate, not the first coordinate: zero-first images are retained.

Each ray with q=p/d is stored as R=d r,Q=p. Nodes remain exactly(r·x)/q,
and scaled determinant D'=D∏d gives exactly the old per-cone D/∏q weight.
The paired-orbit factor24 is applied exactly once in the native alpha.
There is no additional factor12 or2.

    paired image/cone tuples=13229·1920=25399680,
    logical contributions for EACH of G,H=609592320,
    actual G DD5 evaluations=25399680,
    actual H DD5 evaluations=25399680,
    total actual DD5 field evaluations=50799360.

Every tuple evaluates both fields, including zeros/repeats. Zero and repeated
tuple counts are recorded separately and multiplied by24 only in the output
ledger. Node cache hits do not suppress or miscount either DD field.

## Actual compiled outer driver

One new build succeeded, zero compile fixes. ABI checks require FLINT3.6.0,
8-byte ulong/fmpz and48-byte arb_struct. Mathematical operations are Arb512
and exact rational/integer parsing, with int128 intermediate checks; no
binary64/longdouble mathematical substitution. Source and binary SHA are in
build_record.json. The entire normalized input was parsed by native validate
mode, which checks all sizes, IDs, six distinct cone rays, int64 ranges,
orbit uniqueness, absolute multiset and lexicographic negative parity, but
does not construct a time ball or call the DD routine.

The persistent joint node cache holds real atan(b/T) normalized jets through
order5, with genuine zero-node odd jets. Simultaneously it evaluates G's
constant functional: value(π/2)sign(b) at nonzero nodes, value0 at zero,
all higher constant jets0. Those latter zeros are finite-part coefficient
extraction under the accepted theorem, not derivatives of physical G.
Repeated nodes use Hermite insertion; distinct exact nodes are sorted and
must be interval-separated. Neither field has a zero/repeat rejection or a
term-dropping branch. The G functional already includesπ/2, so **both fields
are divided by(2π)^6=64π6**, not G by another128π5. K is their interval
difference.

No gradient/Hessian is present. No old author full or05316scalar batch was
executed. There are six new input/structure tests: full S2 certificate;
independent minimum witness; full driver input and multiset orbit; ABI;
native complete stream validation; exact rank-six witness. Limit8 is not
exceeded.

## Release, checkpoints and stopping

launch_once.py requires a new full_release.json explicitly marked main release
after accepted B055, with pinned B055 and this precheck manifest hashes.
It checks every frozen local/external input, exact eligibility, input/source/
binary identity, and F050/F052 pins. Only then does it write the native token
and an irreversible attempt1 marker. None of those full files exists yet.

The native full loop has no artificial memory/CPU shutdown. It checkpoints
every32 completed images, preserving completed tuple counts and both DD
field counts; it writes a complete result only after all1920images and13229
cones match the exact required totals and G/H/K are finite. Checkpoint and
launch/output files are append-only/new phase artifacts, not replacements
for this frozen precheck. No restart/continuation of a partial sum is allowed.

The launcher limits the one actual native process to1800seconds. On timeout
it terminates only that child PID, waits, escalates to kill only if needed,
then waits again. It records PID, command, source/input/output hashes, wall
time, exit state and checkpoint. Failure/ABI/hash/count/unresolved distinct
nodes/nonfinite output means stop without retry; no partial or fake finite
G is returned. The caller will additionally retain its actual tool session
and chunk IDs when the future run is launched.

For a complete result, the exact comparison threshold is

    T/4096=(32/825)/4096=1/105600.

Kupper<threshold excludes only this fixed frame/time necessary condition;
Klower>threshold passes only this one point's necessary condition. An
overlapping interval remains unresolved, with no automatic re-run. No
classification gives a whole-cube mT, global g/I or Cn bound.

Current state: stable precheck/input/core supplied to B055 for independent
audit; waiting for B055 PASS and the main's separate one-native release.
No old053/source/runtime files were modified, no environment was installed,
and there is no MATLAB, expansion, push or shared stage.
