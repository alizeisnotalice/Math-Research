# C022: exact structural compression at one new rational point

The new point x0=(50,50,50,45,27,0)/50 admits substantial exact symmetry
compression, but its **raw** density denominator still exceeds the prior
512-degree gate. No G, H, K, p, or new-time value was computed. No
rational-function cancellation, common-denominator expansion, or Sturm
calculation was started. No value or inequality from B015 was transferred
to this different point.

## Complete group accounting

|Structural quantity|Exact result|
|---|---:|
|Signed-permutation group elements|46080|
|Point stabilizer|12|
|Distinct signed/permuted point images|3840|
|Multiplicity of every point image|12|
|Actual point-orbit/cone pairs processed|80640|
|Logical full-group/cone terms retained|967680|
|Distinct rational b values, including zero|1583|
|Distinct b^2 values, including zero|792|
|Unordered signed node tuples|80640|
|Tuples after global-negation pairing|40320|
|Self-paired tuples|0|
|Distinct squared-node multisets|39994|
|Minimum nonzero absolute b|1/5800|

The stabilizer is exactly permutation of the three equal positive
coordinates and the sign flip of the zero coordinate, of order 3!*2.
The enumeration independently checked all 46080 transformations and
their uniform orbit multiplicities. It then processed each of the 3840
images against all 21 original cones with multiplicity 12. Thus the
reduction does not omit any original group contribution.

All positive rational weights were merged exactly. Their total equals
520066141/48769479475200, the original fan's full exponential-integral
coefficient. This is a geometric coefficient checksum, not a density
value at x0.

## Zero and repeated nodes are material

There are 88320 logical terms containing a zero node, 201408 containing
a repeated signed node, and 12864 containing both. The zero-node
multiplicity histogram is

    0 zeros: 879360 terms
    1 zero:   79104 terms
    2 zeros:   9216 terms.

The signed-node multiplicity patterns are retained in `results.json`;
the maximum multiplicity is four. Consequently the C014 generic-only
engine, which rejects exact zero nodes, is not usable at this point.
Any later G calculation must retain the whole Laurent principal-part
sum; any H/band calculation must handle zero and confluent nodes.
The present calculation evaluates neither.

## Exact singular-set separation

The four coordinate absolute values 1,9/10,27/50,0 already exclude S2.
The stronger off-S3 claim is also proved, not assumed from the count of
four values. For any combination of three sign vectors, four selected
coordinates either repeat an absolute value or satisfy a signed sum
relation with coefficients +/-1. At coordinates 1,4,5,6 of x0, all
absolute values differ and the nonzero absolute signed sums are
11/25,16/25,34/25,61/25. Thus x0 is off S3.

The minimum pair-form absolute value is 1/10 and the minimum signed-sum
absolute value is 11/25. Distances to these containing hyperplanes give

    dist(x0,S3)^2 >= 1/200.

This is a certified lower bound, not the exact nearest-subspace distance.
The argument and its finite exact checks are in `theorem.md` and
`results.json`. Off-S3 local regularity is compatible with zero Fourier
nodes: individual cone singular terms can cancel in the complete sum.

## Candidate common denominator and the explicit gate

With z=T^2, the candidate denominator is

    D_raw(z)=product_beta (z+beta)^m_beta,

where beta ranges over the 792 rational b^2 values and m_beta is the
largest occurrence multiplicity in any single term. Its exponent counts
are

    exponent 1: 716 factors
    exponent 2:  56 factors
    exponent 3:  19 factors
    exponent 4:   1 factor.

Thus deg_z D_raw=716+112+57+4=889, equivalently degree 1778 in T.
The beta=0 factor has raw exponent two. All factors and multiplicities
are listed exactly in `denominator_factors.json`.

**889 is not the reduced degree.** No numerator was formed and no
factor cancellation was performed. The prior reduced-degree-at-most-512
gate therefore remains unresolved. Raw structure does not fit that
gate, so a global denominator expansion or root isolation is not
authorized automatically. The regularity theorem predicts cancellation
of the full zero-pole principal part, but does not certify any calculated
pole coefficients or the number of other surviving poles.

A later bounded exact cancellation stage can avoid a huge common
denominator initially. From a stored tuple form
P(z)=z^3-e2(b)z^2+e4(b)z-e6(b). For each raw factor (z+beta)^m, remove
that factor from the tuple denominator to obtain E(z), then compute

    a_(beta,k)=alpha/(m-k)! *
      [d^(m-k)/dz^(m-k) (P(z)/E(z))]_(z=-beta),  k=1,...,m.

Sum these coefficients over all stored weights in exact rational
arithmetic. There are at most 889 raw pole slots; derivatives have
order at most three. If r_beta is the highest index with a nonzero
summed coefficient, the actual reduced denominator is
product_beta(z+beta)^r_beta, because distinct linear factors are coprime.
This would compute the reduced degree without first expanding D_raw.
It still requires a separate bounded stage and exact identity checks;
it was **not executed here**. If its reduced degree exceeds 512, the
earlier root-isolation route must stop or receive a newly designed
certificate, rather than silently enlarging its budget.

## Frozen tuple interface and normalization

`merged_tuples.json.gz` is a lossless compressed JSON object with a
`tuples` array. Each row contains:

- `nodes`: six sorted exact rational strings, one representative modulo
  negating all six nodes;
- `weight`: the **full-group summed alpha**, already including the
  point stabilizer multiplicity and both global signs;
- `logical_terms`, `point_orbit_cone_terms`, and
  `global_negation_orbit_size`: exact accounting fields.

Do not multiply the stored weight by 12 or 2 again. For the density,
the whole compressed formula is precisely

    p_T(x0)=(2pi)^-6 sum_rows weight *
        Re product_(j=1)^6 (T-i nodes_j)^-1.

Equivalently it is (2pi)^-6 times the weighted even rational function
of z=T^2 above. No factor (2pi)^-6 is included in `weight`. Time
integration is with respect to T, not z. The global-sign merge is valid
because the two products are complex conjugates and have equal real
parts. All partner weights and counts were checked before merging.

The compressed archive SHA-256 is

    562e3703286fc9a9712efdad30bbebb71433c4f3f15f867c985ee41b80c6f719.

The uncompressed JSON content SHA-256, also recorded in `results.json`, is

    7deb058449e11fa74d7e1a4658311cb738a38fbbb57ab03118b494138cd39390.

## Execution, provenance, and retained files

The registered structural batch completed normally in 2.705768917
seconds, PID 35975, terminal session 1384, exit 0. Peak RSS was
189104128 bytes, below both the conservative 768-MiB stop guard and
the 1-GiB ceiling. The full logical count was reached exactly, and the
group enumeration was not repeated.

The original integer fan hash was checked against
7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812.
The input manifest fixes that fan, original geometry/verification sources,
the read-only 017 proof, and the current registration, proof and code.
The final 11.48-MB tuple JSON was losslessly compressed to 331552 bytes;
its restored byte hash was verified before removing the redundant local
plain copy. `archive_receipt.json` records both hashes and sizes. No
mathematical enumeration or kernel call was performed while archiving.

To reproduce the structural batch, run `python3 preflight.py` in this
directory. It creates the plain tuple JSON; it does not execute a native
binary or call any kernel evaluator. The retained compressed tuple file
can be read with `gzip.open(..., 'rt')`. The proof, registration, script,
manifests/results, exact factors, compressed weighted tuples, archive
receipt and delivery receipt are the necessary final artifacts.

There is no all-T conclusion and no C6 endpoint change. Nothing was
staged, committed, or pushed. All older directories remain frozen.
