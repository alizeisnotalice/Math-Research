# R006 read-only implementation review

2026-10-27. Reviewer n7_theory (the theory author, not a fresh independent numerical auditor). Reviewed without execution:

    /Users/zhengzhihao/Desktop/T/Math-Research-cube7-radau4-r006/cube-dim-7-bounds/upper/round006/radau4_mellin/run.py
    /Users/zhengzhihao/Desktop/T/Math-Research-cube7-radau4-r006/cube-dim-7-bounds/upper/round006/radau4_mellin/assemble.py

No mathematical script, checker, panel, interval scalar assembly or MATLAB was run by this reviewer. Main/independent audit owns acceptance of the actual evidence.

Read-only findings:

- New exponents, group counts2048/1024/512/512, a=r/2, exact centres and halfwidths match the frozen proof. Assembly computes rho=3a/4 and each actual eta, its own Gamma normalisation, Cauchy error and recurrence tail. Series order129 and512-bit precision match the contract.
- R004/R005 complete certificate bytes are SHA256 checked against their copied manifests, and the correct old q endpoint fields are reused. No q0 Mellin call occurs.
- R1,...,R4 and all four explicit moment coefficients match the analytical expansion. D, alpha,beta, discriminant, stable root2 formula and root1=beta/root2 give enclosures of the true moment-defined nodes. Positivity and separated-root gates are explicit.
- Weight formulas, positivity, total mass and degree0...4 residual checks match the contract. These diagnostic overlaps are correctly supplemented by the theoretical exactness/global-majorant proof; no midpoint polynomial is used as an unproved global enclosure.
- `correction(d,4)` constructs1-T16(1+d)^(-1/4) by positive difference/doubling and a fourth root. Lower weights and lower nodes make the sum a lower correction. The complete-I denominator and q_a upper usage are correct.
- The R005 true-moment comparison reconstructs mu_old=R2-1,nu_old=R4/2-2R2+3/2. Its upper correction uses mu_upper,nu_lower, which is safe by the proved monotonicity, and its budget lower uses q_a_lower and1-delta_old_upper. Comparing the new interval upper to the lower endpoint of this old budget expression is the correct direction. It is not a lower bound on actual I. The assertednu_lower>mu_upper² is stronger than the theorem requires; if it passes, it is safe.
- The first4 panels of each group are included once and the remaining jobs start at index4. Full unique coverage is checked in assembly. No automatic expansion is in the reviewed code.

No mathematical or endpoint-direction obstruction was found in the reviewed implementation. This does not assert execution success, independent recomputation of coefficients or a numerical value. Any actual run/assembly failure belongs in Terra's receipt.
