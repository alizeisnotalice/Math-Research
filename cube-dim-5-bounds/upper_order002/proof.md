# Top-two order conditioning at 11110

Use β=3/4 and the 16 unnormalized body vectors in {1}×{−1,1}⁴. Their
frame sum is 16I. The stable variables S_j are iid with Laplace transform
exp(−λ^β). The Green normalization and exact time integral are those in
../upper_ridge001/proof.md, giving

    κ=3Γ(7/4)/(2^(7/2)π^(5/2)).

We retain the unordered pair {v,w} of directions whose stable variables
are largest. The stable law is continuous, so these 120 events partition
the probability space up to a null set. Conditional on retained values
z1,z2>0, the remaining 14 variables are iid conditioned on S≤t, where
t=min(z1,z2). Their conditional mean and harmonic mean are m(t),a(t) from
the exact series in ../upper_order001/proof.md. There is no change of
measure and no entropy penalty. The joint probability element for a
fixed unordered pair, integrated over the full positive square, is

    f(z1)f(z2)F(min(z1,z2))^14 dz1 dz2.

No extra factor 1/2 is present. The individual retained directions label
the two coordinates even though the event itself is an unordered pair.

## Conditional matrix and inverse identities

Let r=v·w∈{−3,−1,1,3} and

    B_c=c(16I−vvᵀ−wwᵀ)+z1vvᵀ+z2wwᵀ,
    D_c=(121−r²)c²+(55+r²)c(z1+z2)+(25−r²)z1z2.

Then det B_c=(16c)³D_c. The remaining 14 directions span R⁵, since the
largest eigenvalue of vvᵀ+wwᵀ is 5+|r|≤8<16. Exactly the same conditional
log-determinant Jensen and deterministic least-squares bounds as in
order001 apply, with finite first and inverse moments for every t>0.
After exact time integration,

    G(x) ≥ H_top2(x)
      = κ Σ_{unordered {v,w}} ∫∫ f(z1)f(z2)F(t)^14 K_vw(z1,z2) dz1 dz2,
    K_vw=(16m(t))^−3/2 D_m(t)^−1/2 [xᵀB_a(t)^−1x]^−7/4.

At x=(1,1,1,1,0), set B=(v·x)², E=(w·x)², C=r(v·x)(w·x), and R=r².
The inverse quadratic has a cancellation-free polynomial representation

    xᵀB_c⁻¹x = N_c/(16cD_c),
    N_c=h0 c²+h1 c z1+h2 c z2+h3 z1z2,
    h0=484−4R+11(B+E)+2C,
    h1=220+4R−11B+5E−2C,
    h2=220+4R+5B−11E−2C,
    h3=100−4R−5(B+E)+2C.

This follows directly from the two-column Woodbury identity. All four
coefficients are nonnegative for the classes used here. Hence direct
Arb evaluation avoids cancellation of large inverse terms. Local SymPy
checked the full 5×5 determinant and inverse against these formulas for
an integer representative of every class.

The kernel depends only on the four invariants (R,B,E,C), apart from
interchanging z1,z2 when B,E are interchanged. We take B≥E and integrate
both labels, which reduces the 120 pairs to these 11 exact classes:

| R | B | E | C | Multiplicity |
|---:|---:|---:|---:|---:|
|1|0|0|0|12|
|1|4|0|0|24|
|1|4|4|−4|12|
|1|4|4|4|12|
|1|16|0|0|12|
|1|16|4|8|8|
|9|0|0|0|3|
|9|4|0|0|24|
|9|4|4|12|4|
|9|16|4|24|8|
|9|16|16|48|1|

The code derives this table from all 120 integer pairs; multiplicities
sum to 120. Equality of the invariants directly proves equality of the
integrands, without needing an unverified geometric symmetry claim.

## Exact box weights

Let I_i=[b_i,b_(i+1)] with increasing endpoints. For i>j, the box with
z1∈I_i,z2∈I_j has the exact weight

    [F(b_(i+1))−F(b_i)] [F(b_(j+1))^15−F(b_j)^15]/15.

The other label orientation has the same weight but generally a different
kernel, so the code explicitly adds both. On a diagonal box I_j×I_j,
write L=F(b_j), H=F(b_(j+1)). Its exact weight is

    2∫_L^H (H-u)u^14 du
      =2[H(H^15−L^15)/15−(H^16−L^16)/16].

Every power difference means H^k−L^k, never (H−L)^k. The factor2 already
includes the two triangles. Kernel bounds over the whole diagonal box
are applied once to this weight.

We use 256 bins b_j=(1/4)(17/16)^j. The functions m(t),a(t) are increasing.
Also B_c increases in each of c,z1,z2, so its inverse quadratic decreases,
while D_c increases. Upper endpoint m,z1,z2 bound the determinant factor
below, and lower endpoint a,z1,z2 bound the inverse factor below after
raising to −7/4. Reversing these choices gives the upper kernel bound.
This covers each full box, not just corner samples of an unproved profile.

All probability and moment enclosures use the exact stable series from
order001 with 64 terms in each nonzero residue class modulo4. The first
moment series requires care: its k=1 and k=5 terms have the same sign.
The code sums all 64 terms of each class explicitly and applies an
alternating-remainder bound **only after the first omitted term**, whose
index is 257,258,or259. It does not apply the alternating test to the
finite prefix. For the first-moment multiplier the tail ratio is bounded
by the stable-tail ratio once k≥2. For the inverse-moment tail, its
additional multiplier ratio is below (k+4)/k. Thus the common sufficient
condition is z⁻³/k<1 at every omitted index. On z≥1/4 this is at most
64/257<1 and decreases thereafter. This explicitly handles both issues
without importing an unproved remainder bound.

As a normalization check, the sum of all box weights times120 agrees with

    F(U)^16−16F(U)F(ℓ)^15+15F(ℓ)^16,
    ℓ=1/4, U=b_256.

This is exactly the probability that all 16 variables are at most U and
at least two are at least ℓ. The Arb intervals overlap identically with
this closed probability formula.

## Uniform integrable complement bounds

Let s=max(z1,z2), t=min(z1,z2), and let v denote the direction of s.
Because a(t)≤t,

    B_a ≤ t(16I−vvᵀ)+s vvᵀ.

At 11110 the squared distance from any body line is at least4/5. Inversion
reverses the preceding matrix order, so

    xᵀB_a⁻¹x ≥ (4/5)/(16t)=1/(20t).

The determinant polynomial gives D_m≥16(m+s)(m+t)≥16st. Consequently,

    K_vw(s,t) ≤ C0 m(t)^−3/2 t^(5/4) s^−1/2,
    C0=20^(7/4)/(4·16^(3/2)).

This is independent of pair geometry and controls both omitted regions.
For t<ℓ, use s≥t, m≥a≥F(t)/n0, and n0=E S⁻¹=Γ(1/β)/β. After integrating
the larger variable by a probability at most1, the entire low omission
per unordered pair is bounded by

    2 C0 n0^(3/2) ℓ^β F(ℓ)^(27/2)/(27/2).

For the high complement, first consider exactly one retained variable
above U and t∈[ℓ,U]. Define

    A_U=∫_ℓ^U f(t)F(t)^14 m(t)^−3/2 t^(5/4) dt.

An upper bound for A_U is computed with the same bins: each exact weight
[F_hi^15−F_lo^15]/15 is multiplied by m(lo)^−3/2 hi^(5/4). This yields
the per-pair upper bound 2 C0 U^−1/2 T(U) A_U, where T=1−F.

If both retained values exceed U, let W=min(S1,S2) for two independent
original stable variables. The extra F(W)^14 factor can be dropped for
an upper bound, and m(W)≥m(U). The elementary Laplace estimate

    T(t)≤C_tail t^−β,      C_tail=(1−exp(−1))⁻¹

gives

    E[W^β 1_{W>U}]
      =U^βT(U)²+∫_U^∞ βt^(β−1)T(t)²dt
      ≤U^βT(U)²+C_tail²U^−β.

Thus the remaining high contribution is at most C0 m(U)^−3/2 times this
last expression per pair. Multiplying all complement bounds by120κ
bounds the complete omitted contribution to H_top2. Overlap with the
low omission is harmless. No infinite moment or extrapolated tail remains.

## Certified result and scope

The single registered batch produced the strict outward-rounded bounds

    0.0000608940491318362 < H_top2(11110)
      < 0.0000692276885126884 < 0.0001.

The low omission is less than7.226e−52 and the entire high omission less
than6.763e−7. In particular this layer strictly improves the previous
single-maximum formula, whose full upper bound was0.0000431012868806624.
That improvement is established from these computed intervals; monotonic
improvement between conditional Jensen layers was not assumed.

The upper endpoint still lies below the target1e−4, so recovering discarded
probability or refining the registered bins cannot make this specific
top-two formula attain the target at11110. This is not an upper bound on
actual G and does not exclude finer conditional constructions. No top-three
layer, spatial grid, full-cube cover, or new C5 bound is asserted.
