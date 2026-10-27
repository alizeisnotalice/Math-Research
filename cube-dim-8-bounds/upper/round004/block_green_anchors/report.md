# 4+4 Green point anchors

For the accepted 16-direction 4+4 product model with shared t, let A=(1,1,1,0) and B=(1,1,0,0). With p4(t,x)=pi^-4 P_x(t), the exact rational P_A and P_B are recorded in result.json. They were obtained from all 1152 signed-permutation cones, grouped into 86 and 34 distinct node multisets respectively, and combined exactly before the final simplification. The intermediate calculation uses gcd cancellation of the accumulated full sum; it does not remove any cone or individual pole.

The three point integrals are exactly:

* G8(A,A) = (10495703317/2149908480000) pi^-7 ≈ 0.00000161637626420145318805610743267.
* G8(A,B) = (15547486331/1920360960000) pi^-7 ≈ 0.00000268057650309721195394310120012.
* G8(B,B) = (143/10584) pi^-7 ≈ 0.00000447339344911111557998972511141.

For each product, result.json records its exact partial-fraction identity in u=t². Applying ∫₀∞(t²+a²)^-k dt = π binom(2k−2,k−1)/(2·4^(k−1) a^(2k−1)) proves the stated coefficients. All a are positive rational. The simplified individual denominators are nonzero at t=0, so each product is locally bounded and integrable there. From the raw cone terms, |P_x(t)|≤(11/3072)t^-4 for t>0, hence |p4(t,x)p4(t,y)|≤pi^-8 (121/9437184)t^-8; the tail is integrable.

The check.py independently re-evaluates the complete unsimplified cone sum at four exact rational times per anchor, verifies each product's partial-fraction identity and coefficient, and independently integrates both 4D densities to recover the known coefficients 17/640 and 1/24 times pi^-3. mpmath quadrature at 80 and 160 decimal digits agrees with all three closed forms; it is numerical cross-checking, not an interval certificate. The exact point evaluations depend on the accepted model and cone geometry, and do not prove a cube minimum, an I certificate, g0, or a new C8 bound.

A first direct SymPy sum was interrupted after 66.53 seconds; its script and receipt are retained. The successful exact polynomial accumulation and validation are in run.py and check.py, with actual timings and hashes in receipts and manifest. Reproduce with the local Python environment and run both scripts; no locked source producer is executed.
