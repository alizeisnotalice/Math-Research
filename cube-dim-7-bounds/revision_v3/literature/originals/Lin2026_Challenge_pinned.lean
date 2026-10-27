/-
Copyright (c) 2026 Yongxi Lin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yongxi Lin
-/
module

public import Mathlib.Analysis.SpecialFunctions.Sqrt
public import Mathlib.Analysis.SpecialFunctions.Pow.Real
public import Mathlib.MeasureTheory.Function.L1Space.Integrable
public import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
public import Mathlib.MeasureTheory.Measure.Lebesgue.VolumeOfBalls

/-!
# Challenge: the weak type constant of the centred maximal operator over cubes

For `f : ℝᵈ → ℝ` the centred Hardy–Littlewood maximal function over axis-parallel cubes is
`M f (x) = sup_{r > 0} |Q(x, r)|⁻¹ ∫_{Q(x, r)} |f|`, where `Q(x, r) = ∏ᵢ [xᵢ - r, xᵢ + r]` is the
closed cube of side `2r` centred at `x`. Its weak type `(1, 1)` constant `c_d` is the least `C` with
`α · |{M f > α}| ≤ C ‖f‖₁` for all integrable `f` and all levels `α`.

This file states:

* `weakTypeConstant_le_two_pow`: `c_d ≤ 2 ^ d` in every dimension (the Vitali covering argument
  with centred cubes);
* `ofReal_phi_le_weakTypeConstant_two`: `Φ ≤ c₂`, where `Φ = 1.68550999…` is the explicit
  algebraic number `phi`; the previously published lower bound is `c₂ ≥ 1.62119…` (Aldaz, 2000);
* `weakTypeConstant_two_le_upper`: `c₂ ≤ 3.879`, below the classical `4`, by comparison with an
  explicit kernel for the Cauchy generator;
* `lt_phi` and `phi_lt`: `1.685 < Φ < 1.686`, so that the size of the explicit constant `Φ` is
  part of the statement.
* `ballWeakTypeConstant_two_le_exp`: the centred weak type constant for Euclidean discs is at
  most `e`;
* `ballWeakTypeConstant_le_rpow`: in dimension `n ≥ 3`, the centred weak type constant for
  Euclidean balls is at most `(n / 2) ^ (n / (n - 2))`.

The type `Fin d → ℝ` carries the sup norm in Mathlib, so `Metric.closedBall x r` is exactly the cube
`Q(x, r)`, and `volume` is Lebesgue measure.
-/

@[expose] public section

noncomputable section

open MeasureTheory Metric
open scoped ENNReal

namespace CenteredMaximal

/-- The centred Hardy–Littlewood maximal function of `f : ℝᵈ → ℝ` over axis-parallel cubes:
`M f (x) = sup_{r > 0} |Q(x, r)|⁻¹ ∫_{Q(x, r)} |f|`, with values in `[0, ∞]`.

Since `Fin d → ℝ` carries the sup norm, `closedBall x r` is the closed cube `∏ᵢ [xᵢ - r, xᵢ + r]`
of side length `2r` centred at `x`; `volume` is Lebesgue measure. -/
def maximalFunction {d : ℕ} (f : (Fin d → ℝ) → ℝ) (x : Fin d → ℝ) : ℝ≥0∞ :=
  ⨆ (r : ℝ) (_ : 0 < r), (volume (closedBall x r))⁻¹ * ∫⁻ y in closedBall x r, ‖f y‖ₑ

/-- `C` is a weak type `(1, 1)` bound for the centred maximal operator in dimension `d`: for every
integrable `f : ℝᵈ → ℝ` and every level `α ∈ [0, ∞]`, `α · |{x : M f (x) > α}| ≤ C · ‖f‖₁`. -/
def IsWeakTypeBound (d : ℕ) (C : ℝ≥0∞) : Prop :=
  ∀ f : (Fin d → ℝ) → ℝ, Integrable f → ∀ α : ℝ≥0∞,
    α * volume {x | α < maximalFunction f x} ≤ C * ∫⁻ x, ‖f x‖ₑ

/-- The weak type `(1, 1)` constant `c_d` of the centred Hardy–Littlewood maximal operator over
axis-parallel cubes in `ℝᵈ`: the least weak type bound. -/
def weakTypeConstant (d : ℕ) : ℝ≥0∞ :=
  sInf {C | IsWeakTypeBound d C}

/-- The centred Hardy–Littlewood maximal function over **Euclidean balls**. The ambient space
`EuclideanSpace ℝ (Fin d)` has the Euclidean norm, unlike `Fin d → ℝ` above, whose sup norm
gives cubes. Open and closed balls yield the same averages for integrable functions. -/
def ballMaximalFunction {d : ℕ} (f : EuclideanSpace ℝ (Fin d) → ℝ)
    (x : EuclideanSpace ℝ (Fin d)) : ℝ≥0∞ :=
  ⨆ (r : ℝ) (_ : 0 < r),
    (volume (ball x r))⁻¹ * ∫⁻ y in ball x r, ‖f y‖ₑ

/-- `C` bounds the centred weak type `(1, 1)` maximal inequality over Euclidean balls. -/
def IsBallWeakTypeBound (d : ℕ) (C : ℝ≥0∞) : Prop :=
  ∀ f : EuclideanSpace ℝ (Fin d) → ℝ, Integrable f → ∀ α : ℝ≥0∞,
    α * volume {x | α < ballMaximalFunction f x} ≤ C * ∫⁻ x, ‖f x‖ₑ

/-- The optimal centred weak type `(1, 1)` constant for Euclidean balls in `ℝᵈ`. -/
def ballWeakTypeConstant (d : ℕ) : ℝ≥0∞ :=
  sInf {C | IsBallWeakTypeBound d C}

/-- The constant `Φ = 1.68550999335552518…`, an algebraic number of degree 16:
`Φ = ((77 + 16√22)/2 - (8 + √22 - √(70 + 8√22)) (11 + √22 - 2√2 - 2√11 - √(17 + 4√22)))
/ (26 + 4√22)`.

It is the covered area per unit mass of a periodic measure: unit masses and masses
`w = (17 + 4√22)/9` alternate along the columns `x = i h`, `h = (5 + √22)/6`, and the rows are
`y = j V`, `V = (11 + √22)/6`. -/
def phi : ℝ :=
  ((77 + 16 * √22) / 2 -
      (8 + √22 - √(70 + 8 * √22)) * (11 + √22 - 2 * √2 - 2 * √11 - √(17 + 4 * √22))) /
    (26 + 4 * √22)

/-- `1.685 < Φ`. -/
theorem lt_phi : (1685 / 1000 : ℝ) < phi := by
  sorry

/-- `Φ < 1.686`. -/
theorem phi_lt : phi < 1686 / 1000 := by
  sorry

/-- **Upper bound.** In every dimension `d`, `c_d ≤ 2ᵈ`: the Vitali covering argument, which for
centred cubes only needs to cover the centres, gives the factor `2ᵈ` instead of `3ᵈ`. -/
theorem weakTypeConstant_le_two_pow (d : ℕ) : weakTypeConstant d ≤ 2 ^ d := by
  sorry

/-- **Upper bound from kernel comparison.** `c₂ ≤ 3.879`, which is below the classical `2ᵈ = 4`.
The bound comes from an explicit comparison kernel for the Cauchy generator `-|D_u| - |D_v|`
in the diamond coordinates `r = |u| + |v|`. -/
theorem weakTypeConstant_two_le_upper : weakTypeConstant 2 ≤ ENNReal.ofReal (3879 / 1000) := by
  sorry

/-- **Lower bound.** `Φ ≤ c₂`. The previously published lower bound was
`c₂ ≥ 3/4 - √2/4 + √6/2 = 1.62119…` (Aldaz, 2000, Proposition 1.4 with `n = 2`). -/
theorem ofReal_phi_le_weakTypeConstant_two : ENNReal.ofReal phi ≤ weakTypeConstant 2 := by
  sorry

/-- **Planar disc bound.** The optimal weak type constant for centred Euclidean discs is at
most `e = exp 1`. -/
theorem ballWeakTypeConstant_two_le_exp :
    ballWeakTypeConstant 2 ≤ ENNReal.ofReal (Real.exp 1) := by
  sorry

/-- **Higher dimensional ball bound.** For `n ≥ 3`, the optimal weak type constant for centred
Euclidean balls is at most `(n / 2) ^ (n / (n - 2))`. -/
theorem ballWeakTypeConstant_le_rpow (n : ℕ) (hn : 3 ≤ n) :
    ballWeakTypeConstant n ≤
      ENNReal.ofReal (((n : ℝ) / 2) ^ ((n : ℝ) / ((n : ℝ) - 2))) := by
  sorry

end CenteredMaximal
