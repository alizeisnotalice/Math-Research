# C022: one exact reduced point, structure only

This directory is the sole writable ownership. Earlier rounds, including
017, are frozen. No stage, commit, push, subagents, MATLAB, native kernel,
new G/H/K/p value, time evaluation, optimization, or Sturm calculation.

The sole input point is x0=(50,50,50,45,27,0)/50. It is new; no B015
numerical value is inherited. The original body32 integer fan must have
SHA-256 7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812.

Before enumeration, save a manifest with exact point, source, proof,
registration, and original-fan/geometry hashes. The one batch has a
1200-second wall limit, at most 967680 logical group/cone terms, and a
1-GiB process memory ceiling. A conservative RSS guard stops at 768 MiB;
check it regularly and record actual peak RSS. Any failed assertion,
time cap, or memory guard stops without further computation.

Enumerate all 46080 signed permutations to form exact point-orbit
multiplicities and the stabilizer. Enumerate each distinct orbit image
against the 21 cones with its full multiplicity. This may reduce actual
work but must account for precisely 967680 logical terms. Count exact
b and b^2 values, zero/repeated-node terms and multiplicities, unordered
signed node tuples, and tuples modulo global negation. Merge their
positive rational weights exactly and verify the total weight against
the old geometric certificate. Save the compressed exact tuple data.

Collect the candidate density denominator factors z+b^2, z=T^2, with
each exponent the maximum occurrence multiplicity in any term. Report
their raw degree only. Do not form or evaluate a density numerator,
expand a common denominator, perform rational-function cancellation,
or isolate roots. The prior 512-degree gate is only a comparison with
this raw upper bound; it cannot establish the reduced degree.

Also verify an exact off-S3 certificate using coordinates 1,4,5,6:
their absolute values are distinct and no signed sum is zero. The
proof explains why this excludes every three-body-direction span.
This is finite geometry, not a potential value.
