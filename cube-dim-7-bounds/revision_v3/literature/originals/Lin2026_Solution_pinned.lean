/-
Copyright (c) 2026 Yongxi Lin. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Yongxi Lin
-/
module

public import CenteredMaximal

/-!
# Solution: weak type constants for centred maximal operators

The statements of `Challenge.lean`, proved from the development in `CenteredMaximal`.
-/

@[expose] public section

open scoped ENNReal

namespace CenteredMaximal

/-- `1.685 < Φ`. -/
theorem lt_phi : (1685 / 1000 : ℝ) < phi :=
  phi_mem_Ioo.1

/-- `Φ < 1.686`. -/
theorem phi_lt : phi < 1686 / 1000 :=
  phi_mem_Ioo.2

/-- **Upper bound.** In every dimension `d`, `c_d ≤ 2ᵈ`: the Vitali covering argument, which for
centred cubes only needs to cover the centres, gives the factor `2ᵈ` instead of `3ᵈ`. -/
theorem weakTypeConstant_le_two_pow (d : ℕ) : weakTypeConstant d ≤ 2 ^ d :=
  weakTypeConstant_le (isWeakTypeBound_two_pow d)

/-- **Upper bound from kernel comparison.** `c₂ ≤ 3.879`, which is below the classical `2ᵈ = 4`.
The bound comes from an explicit comparison kernel for the Cauchy generator `-|D_u| - |D_v|`
in the diamond coordinates `r = |u| + |v|`. -/
theorem weakTypeConstant_two_le_upper : weakTypeConstant 2 ≤ ENNReal.ofReal (3879 / 1000) :=
  Cauchy.weakTypeConstant_two_le

/-- **Lower bound.** `Φ ≤ c₂`. The previously published lower bound was
`c₂ ≥ 3/4 - √2/4 + √6/2 = 1.62119…` (Aldaz, 2000, Proposition 1.4 with `n = 2`). -/
theorem ofReal_phi_le_weakTypeConstant_two : ENNReal.ofReal phi ≤ weakTypeConstant 2 :=
  le_weakTypeConstant fun _ hC => Lattice.ofReal_phi_le hC

/-- **Planar disc bound.** The optimal weak type constant for centred Euclidean discs is at
most `e = exp 1`. -/
theorem ballWeakTypeConstant_two_le_exp :
    ballWeakTypeConstant 2 ≤ ENNReal.ofReal (Real.exp 1) :=
  Ball.ballWeakTypeConstant_two_le_exp_of_ae_direct_certificates
    Ball.hasRealAEDirectObstacleCertificates_planarKernel.toENNReal

/-- **Higher dimensional ball bound.** For `n ≥ 3`, the optimal weak type constant for centred
Euclidean balls is at most `(n / 2) ^ (n / (n - 2))`. -/
theorem ballWeakTypeConstant_le_rpow (n : ℕ) (hn : 3 ≤ n) :
    ballWeakTypeConstant n ≤
      ENNReal.ofReal (((n : ℝ) / 2) ^ ((n : ℝ) / ((n : ℝ) - 2))) :=
  Ball.ballWeakTypeConstant_le_rpow_of_ae_direct_certificates n hn
    (Ball.hasRealAEDirectObstacleCertificates_newtonianKernel n hn).toENNReal

end CenteredMaximal
