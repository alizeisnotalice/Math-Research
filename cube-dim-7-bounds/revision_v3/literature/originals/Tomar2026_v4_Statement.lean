/-
Copyright (c) 2026 Tejasvi Singh Tomar. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Tejasvi Singh Tomar
-/
/-
Statement of the bound `c₂ ≤ 3.59123` for the weak type (1,1) constant of the centred
Hardy–Littlewood maximal operator over axis-parallel squares (constant 47a of the
optimization-constants repository). Mathlib-only definitions.

The three definitions are those of Y. Lin's Lean development `centered-maximal-constant`
(`Challenge.lean`, commit 72c022b, Apache-2.0), so that the statement can be compared token by
token with Lin's `weakTypeConstant_two_le_upper`.
-/
import Mathlib.MeasureTheory.Function.L1Space.Integrable
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic

namespace Hl47a

open MeasureTheory Metric
open scoped ENNReal

/-- The centred Hardy–Littlewood maximal function of `f : ℝᵈ → ℝ` over axis-parallel cubes,
`M f (x) = sup_{r > 0} |Q(x, r)|⁻¹ ∫_{Q(x, r)} |f|`. Since `Fin d → ℝ` carries the sup norm,
`closedBall x r` is the closed cube `∏ᵢ [xᵢ - r, xᵢ + r]` of side `2r`; `volume` is Lebesgue measure. -/
noncomputable def maximalFunction {d : ℕ} (f : (Fin d → ℝ) → ℝ) (x : Fin d → ℝ) : ℝ≥0∞ :=
  ⨆ (r : ℝ) (_ : 0 < r), (volume (closedBall x r))⁻¹ * ∫⁻ y in closedBall x r, ‖f y‖ₑ

/-- `C` is a weak type `(1, 1)` bound: `α · |{x : M f (x) > α}| ≤ C · ‖f‖₁` for every integrable
`f` and every level `α ∈ [0, ∞]`. -/
def IsWeakTypeBound (d : ℕ) (C : ℝ≥0∞) : Prop :=
  ∀ f : (Fin d → ℝ) → ℝ, Integrable f → ∀ α : ℝ≥0∞,
    α * volume {x | α < maximalFunction f x} ≤ C * ∫⁻ x, ‖f x‖ₑ

/-- The weak type `(1, 1)` constant `c_d`: the least weak type bound. -/
noncomputable def weakTypeConstant (d : ℕ) : ℝ≥0∞ :=
  sInf {C | IsWeakTypeBound d C}

/-- The claim of the paper's Theorem 1.2, as a weak type bound. -/
def BoundClaim : Prop := IsWeakTypeBound 2 (ENNReal.ofReal (359123 / 100000))

/-- The claim of the paper's Theorem 1.2: `c₂ ≤ 3.59123`. -/
def ConstantClaim : Prop := weakTypeConstant 2 ≤ ENNReal.ofReal (359123 / 100000)

/-- The paper's form with non-strict level sets and real levels `λ > 0`. -/
def NonstrictClaim : Prop :=
  ∀ f : (Fin 2 → ℝ) → ℝ, Integrable f → ∀ l : ℝ, 0 < l →
    ENNReal.ofReal l * volume {x | ENNReal.ofReal l ≤ maximalFunction f x} ≤
      ENNReal.ofReal (359123 / 100000) * ∫⁻ x, ‖f x‖ₑ

/-- The weak type bound gives the bound on the constant. -/
theorem constantClaim_of_boundClaim (h : BoundClaim) : ConstantClaim :=
  sInf_le h

end Hl47a
