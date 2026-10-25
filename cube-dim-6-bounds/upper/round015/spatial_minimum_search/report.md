# B015: a strict point obstruction for the fixed cutoff T=6/125

The sole released B015 batch stopped at its **second native evaluation**.
For the exact point z below, the complete native512 interval proves

\[
K_{6/125}(z)<\frac{111494916581}{10^{16}}
<\frac3{256000}=\frac{T}{4096}
<\frac3{250000}=\frac{T}{4000}.
\]

Thus the body32 short-potential model at this fixed cutoff cannot certify
C6<64, and also fails the sufficient 62.5 target. The point is in the strictly
ordered positive face and strictly outside S2. Continuity makes this an
essential-infimum obstruction, not just a measure-zero exception. The actual
model cost is at least 67.2676408037809295..., a **lower bound on this model's
cost**, not an upper or lower bound on C6. The existing C6<=64 is unchanged.

This conclusion uses the C014 generic full-group interval engine reviewed by
the main agent. Its independent F014 audit is a separate responsibility; B015
does not turn software-validation status into a new global certificate. No
other cutoff, time mixture or spectral frame is excluded by this report.

## Exact point and interval

Let D=100000000000000. The stopping point is

```
z = [100000000000000,
      99999998999999,
      99999997999988,
      90268521812170,
      53714088833376,
        157107163435] / D.
```

In decimals its last five coordinates are
`(.99999998999999,.99999997999988,.90268521812170,.53714088833376,.00157107163435)`.
The complete native value interval, retained as exact binary integers in
`native_outputs/02.stdout.json`, has midpoint approximately
`0.0000111494916580729044644601264373075060523722867918354`
and width below 9.81e-120. The purely rational postprocessor proves the simpler
outward enclosure

\[
\frac{5574745829}{500000000000000}
<K_T(z)<\frac{111494916581}{10000000000000000}.
\]

The exact upper endpoint lies below T/4096 by more than 5.69258e-7. No printed
floating-point midpoint is used for any threshold comparison.

The first evaluated point is the registered new denominator-1e14 rounding of
B010 point2, not B010 point2 itself. Its K value remains above both thresholds:

\[
\frac{124982932789}{10^{16}}<K_T(s)
<\frac{12498293279}{10^{15}},
\qquad K_T(s)\approx1.24982932789951082975\cdot10^{-5}.
\]

The second registered seed, the old A007 reference c, was never evaluated.

## Candidate construction and zero-node handling

SLSQP was used only to suggest candidates in the registered local box of
radius1/20 around the first seed, with x1=1 and all six gaps at least1e-8.
Function and derivative midpoints were both scaled by1e5 internally. There
were exactly two objective/gradient callback requests and two native passes.
No finite-difference gradients, grid, MC, new time, separate G/H pass or other
engine was used.

The second binary64 suggestion, after its prescribed exact denominator-1e14
quantization and feasibility rounding, was

```
[100000000000000,99999999000000,99999998000000,
 90268521812303,53714088834840,157107002384] / D.
```

It hit exactly one primitive zero form, x1-2x2+x3. **It was not sent to the
native engine.** Its six gap numerators were

```
[1000000,1000000,9731476187697,36554432977463,
 53556981832456,157107002384].
```

The registered gap-space repair chose the largest gap j=4 (zero-based).
The first allowed repair k=1 added the exact gap deltas

```
[1,11,121,1331,-162515,161051] / D.
```

Their sum is zero. Every repaired gap is at least1e-8, the point remains in
the original local box, and all9308 primitive signed/permuted Fourier zero
forms are nonzero. This produces exactly z. It is a new exact input, not a
tolerance replacement or skipped-zero evaluation. No k=2,...,6 repair was
needed. Both actual native passes retained all967680 group terms (483840
computed sign pairs), with zero_node_terms=0 and confluent_terms=0.

## Distance and continuity for the essential-infimum conclusion

For V={(1,s2,...,s6):s_i in{-1,1}}, the postprocessor checked every496 pair.
Writing a=z.v, b=z.w and g=v.w, the squared distance to span(v,w) is exactly

\[
|z|^2-\frac{6a^2-2gab+6b^2}{36-g^2}.
\]

The minimum over the pairs is

\[
\operatorname{dist}(z,S_2)^2
=\frac{6020805515960314567115330493}
       {40000000000000000000000000000}>0,
\]

so its distance is approximately0.38796924865124. The first seed's squared
distance is4381753622995444834749720723/40000000000000000000000000000.

There is also a direct open generic-domain certificate. For each primitive
integer form a in the9308-form list, all x with
`||x-z||_infinity < |a.z|/(2||a||_1)` retain its nonzero sign. The exact minimum
radius is1/80000000000000. C014 represents K_T there as a finite sum of smooth
divided differences of atan(T/b) at nonzero nodes (with the Hermite extension
for repeated nonzero nodes). Consequently K_T is continuous in an open
neighborhood of z. Since K_T(z)<T/4096 strictly, a smaller neighborhood has the
same inequality. Its intersection with the cube has positive six-dimensional
Lebesgue measure, even though z lies on the face x1=1. Therefore
essinf_cube K_T <=K_T(z)<T/4096. No convexity or coordinate monotonicity is
assumed, and no explicit continuity radius for the threshold itself is needed.

## Saved-output diagnostics for later topic selection

The following computations use only the two saved K, p_T and gradient
intervals, with exact Fraction endpoint arithmetic:

\[
q(T)=K_T/T,\qquad q'(T)=(T p_T-K_T)/T^2,
\qquad\text{point cost}=T/(64K_T).
\]

| Existing evaluated point | q(T), display midpoint | q'(T), display midpoint | Point-cost interval lower endpoint, display |
|---|---:|---:|---:|
| First seed s |0.0002603811099790647562|−0.0005254734001215743209|60.0081933795285156|
| Stopping point z |0.0002322810762098521763|−0.0004131987380349032385|67.2676408037809295|

Both q' intervals are strictly negative at this single T. This supplies only
a local time-derivative diagnostic: it does not identify an optimal cutoff,
give a time interval of exclusion, or prove all time mixtures infeasible.

At z the free gradient coordinates x2,...,x6 have display values

```
[-1.9016632670088636e-5, -1.9016631358056882e-5,
 -8.8440179587562348e-6, -1.6977746396651534e-6,
 +1.6921065402923250e-7].
```

Each sign is rigorously separated from zero. The coordinatewise derivative
descent signs are(+,+,+,+,−). However, the remaining single-axis allowances
inside the registered ordered local domain are only
`(1e-14,1.1e-13,1.33e-12,1.464e-11,1.61051e-9)`.
They explain that z is near search-box boundaries, not an asserted stationary
point or global minimum. No finite-step decrease is inferred from a derivative
alone, and no new point is evaluated for these diagnostics. The complete
intervals and analogous first-seed diagnostics are in `checked_results.json`.

## Execution, provenance, and stop receipt

Preparation and its runtime release refusal are preserved in
`preparation_report.md`, `checked_preparation.json`, and the four immutable
release-pinned files. The main agent then created release.json with SHA256
`01d229c70878f3a3ce4d9a4a4236c0a8d358f2a661250b25cec5a552c62e6d25`.
No native process was started before this explicit release.

The only `search.py --run` had PID32328 and tool session86155 (start chunk47dbad,
final chunkb2ba32), exit0. Native children were PID32338 and PID32371, both
exit0, approximately3.976s and4.081s outside the kernel. Total batch wall time
was8.07353575s, versus the registered1200s cap. The driver stopped immediately
after validating the second result. There was no continuation, second seed,
precision retry or additional kernel. All children were waited and the tool
session is closed.

Independent pure postprocessing ran twice, without kernel calls: initial
PID32898/chunk7e2102, then the parent-requested time/gradient diagnostics at
PID33017/chunka5738e (about0.057s), both exit0. It verified exact inputs,
repair reconstruction, range/genericity, source/raw-output hashes, all group
counts and endpoints, threshold residuals, every496-pair distance, and the
independently reconstructed interval identity x.grad K+5K+T p_T containing0.

Pinned mathematical-engine provenance:

- C014 final report:94352cfdd767f21a64678e19c0ff231b42ffada06bc19c5e596aea4b8b2699fc.
- C014 source:247a1f6c6967ad52f34320b999818214488700651b46cb8e5feb3eac7064953b.
- C014 binary:67119aa799041ed75b66f043444640a2a62556f58517cbdf3bd6f7cccd564f04.
- Integer fan:7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812.
- First native stdout:fbd3ac97abc87c6335ea8b9b5672e59702bab2b804c354ad0cefcf6198896a7c.
- Stopping native stdout:a053a55720bfa6545a11f833f725c03bbbf7cb39f9682030374ef500a5550397.

All exact paths, input hashes and runtime-library hashes are in
`preparation.json`; all output hashes and process records are in
`checkpoint.json` and `run_receipt.json`. The per-file allowlist is
`manifest.json`. The registered environment was Python3.15.0, SciPy1.18.1,
NumPy2.5.4, python-flint0.9.0, with native512 FLINT3.6.0. Prior rounds and C014
were read only. No MATLAB, subagent, stage, commit, or push was used.

The B015 batch is **completed and stopped**. Its sole mathematical exclusion
is the fixed body32 cutoff T=6/125 at the specified targets. Any further time,
frame, mixture or point computation requires a separately registered batch.
