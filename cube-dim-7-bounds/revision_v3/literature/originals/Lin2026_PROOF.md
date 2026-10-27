# Informal proof

This document is the mathematical account of the Lean development. Section numbers are cited by the
planning documents in `.mathlib-quality/`. Every step below corresponds to a named Lean declaration.

## 1. Definitions

For `x ∈ ℝᵈ` and `r > 0` let `Q(x, r) = ∏ᵢ [xᵢ - r, xᵢ + r]`, the closed axis-parallel cube of side
`2r` centred at `x`; in Lean this is `Metric.closedBall x r` in `Fin d → ℝ`, which carries the sup
norm. For `f : ℝᵈ → ℝ`,

    M f (x) = sup_{r > 0} |Q(x, r)|⁻¹ ∫_{Q(x, r)} |f|        (values in [0, ∞]).

`C ∈ [0, ∞]` is a *weak type bound* in dimension `d` if `α |{M f > α}| ≤ C ‖f‖₁` for every integrable
`f` and every `α ∈ [0, ∞]`; `c_d` is the infimum of the weak type bounds. Closed cubes and the strict
level set are the conventions of Aldaz (2000); closed versus open cubes make no difference for
integrable `f`, and the strict versus non-strict level set gives the same constant.

`|Q(x, r)| = (2r)ᵈ` (`volume_closedBall_eq`).

## 2. The upper bound `c_d ≤ 2ᵈ`

**2.1 (dimension zero).** `ℝ⁰` is a point of measure `1`, every cube is the whole space, and
`M f = ‖f‖₁`. Hence `α |{M f > α}| ≤ ‖f‖₁` (`isWeakTypeBound_zero_one`).

**2.2 (bounded radii).** Let `d ≥ 1`, `0 < α < ∞`, `K = ‖f‖₁ < ∞`. If `α |Q(x, r)| < K` then
`(2r)ᵈ < K/α`. If `2r ≥ 1` then `(2r)ᵈ ≥ 2r`, so in all cases `r ≤ max(1, K/α)`
(`radius_le_of_mul_volume_lt`).

**2.3 (Vitali with an epsilon of room).** Let `E = {M f > α}` and `τ > 1`. For each `x ∈ E` choose
`r_x > 0` with `α < |Q(x, r_x)|⁻¹ ∫_{Q(x, r_x)} |f|`, i.e. `α |Q(x, r_x)| < ∫_{Q(x, r_x)} |f| ≤ K`; by 2.2
the radii are bounded. The Vitali lemma with enlargement `τ`
(`Vitali.exists_disjoint_subfamily_covering_enlargement`) gives a subfamily `u ⊆ E` whose cubes are
pairwise disjoint and such that every `Q(a, r_a)`, `a ∈ E`, meets some `Q(b, r_b)`, `b ∈ u`, with
`r_a ≤ τ r_b`. A common point `z` gives `|a - b| ≤ r_a + r_b ≤ (1 + τ) r_b`: the *centre* `a` lies
in `Q(b, (1 + τ) r_b)`. The disjoint family `u` is countable (each cube has nonempty interior and
`ℝᵈ` is separable). Therefore

    |E| ≤ Σ_{b ∈ u} |Q(b, (1 + τ) r_b)| = (1 + τ)ᵈ Σ_{b ∈ u} |Q(b, r_b)|,
    α Σ_{b ∈ u} |Q(b, r_b)| ≤ Σ_{b ∈ u} ∫_{Q(b, r_b)} |f| = ∫_{⋃ Q(b, r_b)} |f| ≤ K,

so `α |E| ≤ (1 + τ)ᵈ K` (`mul_volume_le_of_one_lt`). Letting `τ → 1` gives `α |E| ≤ 2ᵈ K`
(`isWeakTypeBound_two_pow`), hence `c_d ≤ 2ᵈ`. This is Tao's Exercise 42 (245A Notes 5); its hint
says the centres, not the cubes, are covered and that "one may need to first create an epsilon of
room", which is the role of `τ`.

## 3. The constants

Let `u = (2 + √22)/3` (`root`), the positive root of `3u² - 4u - 6 = 0` (`root_quadratic`). Put

    w = u² - 1 = (17 + 4√22)/9      (heavy),
    h = (1 + u)/2 = (5 + √22)/6     (hgap),
    V = h + 1 = (11 + √22)/6        (vgap),

and the witness sides

    1,   u = √(1 + w),   √w (sideH1),   √2·u = √(2(1 + w)) (sideLH2),
    √(2(2 + w)) (sideLHL2),   2h + 1 = √(2(1 + 2w)).

The identity `(2h + 1)² = 2(1 + 2w)` is equivalent to `3u² - 4u - 6 = 0`
(`two_mul_hgap_add_one_sq`). Numerically `u ≈ 2.2301`, `w ≈ 3.9735`, `h ≈ 1.6151`, `V ≈ 2.6151`,
`√w ≈ 1.9934`, `√2 u ≈ 3.1539`, `√(2(2 + w)) ≈ 3.4565`.

The slot has width `a = 2h - √(2(2 + w))/2 - u/2 ≈ 0.3868` (`slotW`) and height
`b = V - √2 u/2 - √w/2 ≈ 0.0414` (`slotH`). Rewriting `u`, `w`, `h`, `V` in terms of `√22`
(and `√2 · √22 = 2√11`) gives

    Φ = (2hV - 4ab)/(1 + w)       (phi_eq),

which is the closed form in `Challenge.lean`.

The following inequalities are used; each follows from `4.69 < √22 < 4.691` and
`1.414 < √2 < 1.4143`:

    V ≤ √2 u,   2V ≤ u + √2 u,   2h ≤ √(2(2+w)),   V ≤ √(2(2+w)),
    2V - √(2(2+w)) ≤ √w,   4h - √(2(2+w)) ≤ √2 u,   a > 0,   b > 0,   1 ≤ √w.

## 4. The weighted lattice and its witnesses

The atom with index `(c, r) ∈ ℤ²` sits at `(c h, r V)` and has mass `1` if `c` is even and `w` if
`c` is odd (`colWeight`). A *witness* at `(x, y)` is a side `L ≥ 1` and a finite set `A` of atoms,
each within sup-distance `L/2` of `(x, y)`, of total mass at least `L²` (`IsWitness`).

**4.1 (six witnesses).** For `(x, y)` with `0 ≤ x ≤ h`, `0 ≤ y ≤ V/2`:

| name | atoms `(c, r)` | side `L` | mass | valid for |
|---|---|---|---|---|
| light | `(0,0)` | `1` | `1` | `x ≤ 1/2`, `y ≤ 1/2` |
| hlh2 | `c ∈ {-1,0,1}`, `r ∈ {0,1}` | `2h + 1` | `2(1 + 2w)` | `x ≤ 1/2`, `1/2 ≤ y` |
| lh1 | `(0,0), (1,0)` | `u` | `1 + w` | `1/2 ≤ x ≤ u/2`, `y ≤ u/2` |
| lh2 | `c ∈ {0,1}`, `r ∈ {0,1}` | `√2 u` | `2(1 + w)` | `1/2 ≤ x ≤ √2u/2`, `V - √2u/2 ≤ y` |
| h1 | `(1,0)` | `√w` | `w` | `h - √w/2 ≤ x`, `y ≤ √w/2` |
| lhl2 | `c ∈ {0,1,2}`, `r ∈ {0,1}` | `√(2(2+w))` | `2(2 + w)` | `2h - √(2(2+w))/2 ≤ x`, `V - √(2(2+w))/2 ≤ y` |

In each case the mass equals `L²` exactly, and the distance conditions reduce to the stated ranges
together with the inequalities of Section 3. For example, for hlh2: the columns `-h, 0, h` are
within `h + 1/2 = L/2` of `x` because `0 ≤ x ≤ 1/2`; the rows `0` and `V` are within
`V - 1/2 = L/2` of `y` because `1/2 ≤ y ≤ V/2 ≤ V - 1/2`.

**4.2 (coverage of the quarter cell).** Every `(x, y)` with `0 ≤ x ≤ h`, `0 ≤ y ≤ V/2` outside the
open slot `(u/2, 2h - √(2(2+w))/2) × (√w/2, V - √2u/2)` has one of the six witnesses:

* `x ≤ 1/2`: light if `y ≤ 1/2`, otherwise hlh2;
* `1/2 < x ≤ u/2`: lh1 if `y ≤ u/2`; otherwise `y > u/2 ≥ V - √2u/2` (as `2V ≤ u + √2u`) and
  `x ≤ u/2 ≤ √2u/2`, so lh2;
* `u/2 < x ≤ h`: if `y ≤ √w/2`, h1 (as `h - √w/2 ≤ u/2` because `h = (1 + u)/2` and `√w ≥ 1`).
  If `y > √w/2` and `x ≥ 2h - √(2(2+w))/2`, lhl2 (as `√w/2 ≥ V - √(2(2+w))/2`). Otherwise the point is
  not in the slot only if `y ≥ V - √2u/2`, and `x < 2h - √(2(2+w))/2 ≤ √2u/2`, so lh2.

All witnesses use atoms with `c ∈ {-2, …, 2}` and `r ∈ {-1, 0, 1}` (`nearBox 0 0`).

**4.3 (symmetry).** The masses are invariant under `c ↦ -c` and `c ↦ c + 2k`. Reflecting the atom
indices in either axis turns a witness at `(x, y)` into one at `(-x, y)` or `(x, -y)`, and
translating them by `(2k, l)` turns it into one at `(x + 2kh, y + lV)`. Hence every point of
`cell = [-h, h) × [-V/2, V/2)` outside the four open slots `{u/2 < |x| < …, √w/2 < |y| < …}`, and
every translate of such a point by `(2kh, lV)`, has a witness using atoms of `nearBox k l`
(`exists_isWitness_of_abs`, `exists_isWitness_of_mem_goodCopy`).

## 5. From the lattice to integrable functions

Fix `N ∈ ℕ` and `ε > 0`. Keep the atoms with `|c| ≤ 2N + 2`, `|r| ≤ N + 1` (`atomBox N`) and replace
each kept atom of mass `m` at `p` by `m ε⁻²` times the indicator of the closed square of side `ε`
centred at `p` (`smeared N ε`). This is a finite sum of bounded functions with compact support, so
it is integrable, and `‖smeared N ε‖₁` is the total kept mass
`(2N + 3)((2N + 3) + (2N + 2)w) ≤ (2N + 3)²(1 + w)` (`lintegral_smeared`, `sum_colWeight_Icc`).

If `(L, A)` witnesses `z` and `A` is kept, the square of side `L + ε` centred at `z` contains the
whole smeared square of every atom of `A`, so its integral is at least `L²` and

    M(smeared N ε)(z) ≥ L²/(L + ε)² ≥ 1/(1 + ε)² > 1 - 2ε

(`lt_maximalFunction_smeared`; the middle inequality uses `L ≥ 1`, the last one
`(1 - 2ε)(1 + ε)² = 1 - 3ε² - 2ε³ < 1`).

## 6. The lower bound `Φ ≤ c₂`

`goodSet = cell \ slots` has area at least `|cell| - 4ab = 2hV - 4ab` (`ofReal_le_volume_goodSet`). Its
translates `goodCopy k l = goodSet + (2kh, lV)` are pairwise disjoint (the cells are half-open) and
have the same area. For `|k|, |l| ≤ N` the atoms `nearBox k l` are kept in `atomBox N`, so every
point of `goodCopy k l` lies in the level set `{M(smeared N ε) > 1 - 2ε}`. Hence

    |{M(smeared N ε) > 1 - 2ε}| ≥ (2N + 1)² (2hV - 4ab)      (ofReal_le_volume_levelSet).

If `C` is a weak type bound, apply it with `α = 1 - 2ε` and `ε = 1/(2N + 3)`, so that
`1 - 2ε = (2N + 1)/(2N + 3) =: q`:

    q (2N + 1)² (2hV - 4ab) ≤ C (2N + 3)² (1 + w),   i.e.   C ≥ q³ (2hV - 4ab)/(1 + w) = q³ Φ

(`ofReal_mul_phi_le`). As `N → ∞`, `q → 1`, so `C ≥ Φ` (`ofReal_phi_le`). Taking the infimum over
`C` gives `Φ ≤ c₂`.

## 7. Numerics

`1.685 < Φ < 1.686` follows from rational enclosures of `√2`, `√11`, `√22`, `√(70 + 8√22)` and
`√(17 + 4√22)` by interval arithmetic on the closed form (`phi_mem_Ioo`). High-precision evaluation gives `Φ = 1.68550999335552518466528…`.

## 8. Relation to the exploratory computation

The configuration was found by a numerical search over periodic weighted product measures and
certified independently in exact rational arithmetic (see `docs/HISTORY.md`). That computation used
twelve witness types on the full cell and an exact union-area computation; the Lean proof uses the
reflection symmetry to reduce to six witnesses on a quarter cell, and only needs the lower bound
`|goodSet| ≥ 2hV - 4ab`, not the exact area of the level set.
