# R016 origin-density and radial Jensen budgets

The complete original-process integral has the strict rational upper budget

    I < 4018421141119 / 10^18 = 0.000004018421141119.

This is the **Rao–Blackwellized** bound; its underlying expression is approximately 4.01842114111870165e-6. The raw origin-Jensen expression is approximately 4.12700249882942961e-6, with rational upper endpoint 4127002498830/10^18. Both strictly improve the R006 deterministic budget expression. They are upper budgets, not lower bounds on the unknown I.

The origin density is p₁(0)=5040 A_sum/π⁷≈5.38187647332190493e-14, where A_sum is recomputed exactly by summing the 212 previously certified simplex weights. Positive Tonelli for the original Gaussian mixture gives

    E det(Σ)^(-1/2) = (4π)^(7/2) p₁(0).

For Z=det(Σ)^(-1/14), this is E Z⁷; Jensen gives E Z≤2√π p₁(0)^(1/7). Substitution into the accepted complete-boundary interface simplifies to B_fan=p₁(0)^(1/6)/(840√π). Here p₁(0) is a time-one density, not G(0); the latter is infinite.

Conditioning away the independent Gamma(32) radial variable before Jensen yields the sharper multiplier

    η = [Γ(65/2)/Γ(32)] / [(65/2)(67/2)(69/2)]^(1/6)
      ≈ 0.973690018908026917 < 1,
    B_RB = η B_fan.

The raw seventh radial factor is 128 Γ(71/2)/Γ(32), not D⁷. The independent sharp-floor calibration ψ(ξ)≥32√2|ξ| gives p₁(0)≤6/[π⁴(32√2)⁷] and B_sharp≈4.94077714476538487e-6; the exact fan p₁(0) enclosure lies strictly below that bound. The floor's existing 21-orbit receipt was hash-verified; its enumeration was not repeated.

Using the rounded B_RB upper endpoint, the threshold for the existing comparison expression to be at most 128 is approximately g_required=2.59312073254e-6. **If a separate proof establishes the global face floor g≥1/320000**, the exact sixth-power comparison is strict and the resulting conditional expression is approximately 102.96228, below 128. The analogous raw-fan conditional expression is approximately 105.74441. Neither condition proves the required floor, and this batch does not establish a new C7 bound. R015's single-point result rules out its larger floor targets but does not prove or disprove this smaller global floor.

The registered run made ten distinct named 512-bit scalar enclosures, with no phase, new shape or Monte Carlo evaluation. It re-summed 212 exact α weights, verified both raw formulas overlap, checked the Gamma recurrence and η<1, separated the fan and sharp density intervals, and used exact sixth-power comparisons for the conditional costs. The prior primary scalar pilot overlaps the independent raw expression and retains separate provenance. Exact outward dyadic endpoints and rational rounded budgets are in result.json.

The run took 0.0117 s wall and 0.00533 s CPU, with individual process peak RSS 32,440,320 bytes; it was a single-process finite scalar task, not a CPU saturation test. The portable checker repeats the same ten scalar inputs as verification, recomputes the exact weight sum and rational comparisons, and evaluates no phase or sample. No MATLAB, commit or push was used.
