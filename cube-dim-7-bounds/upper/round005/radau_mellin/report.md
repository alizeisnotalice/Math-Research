# R005 two-moment Radau Mellin certificate

The complete original-process integral budget is now

    I < 7730422817512 / 10^18 = 0.000007730422817512.

The underlying deterministic upper expression is approximately 7.7304228175112773e-6, strictly below the R004 pair expression approximately 8.1799500703079561e-6. These are upper expressions for the whole positive-time/full-boundary integral; no lower bound on the actual I is implied. This calculation supplies no global g and no C7 bound.

Let U=T1+T2 for independent Hadamard blocks and q(r)=E U^(-r). R004's complete q(1/8) enclosure and its certificate hash are reused; it is not evaluated again. Three new positive Mellin integrals were enclosed, yielding q(1/16)≈0.795842739003644221, q(3/16)≈0.506748476644451785, and q(1/4)≈0.405408285241186356. These are pair moments q, not four-block moments.

Under independent U^(-1/8)-tilted laws, set W=cosh((log U1−log U2)/16)−1. The exact relations R1=q(1/16)q(3/16)/q(1/8)^2, R2=q(1/4)/q(1/8)^2, μ=R1−1 and ν=R2/2−2R1+3/2 give μ≈0.001745268318486792 and ν≈1.0320272152673e-5. Interval arithmetic proves μ>0 and ν>μ². A global Hermite majorant at nodes 0 and ν/μ, justified by complete monotonicity of g(w)=T8(1+w)^(-1/4), gives a correction δ=μ²/ν[1−g(ν/μ)]. Using a lower μ and upper ν is safe by monotonicity; the computed rigorous δ lower endpoint is approximately 0.02393287874463318418. Consequently M ≤ (q(1/8)^4/2)(1−δ_lower)^2 ≤ 0.077206558253009624, and the original whole-boundary interface gives the I endpoint. Full derivation, admissibility, and interval directions are in the frozen companion proof merge4_proof.md; exact dyadic intervals appear in certificate.json.

For the three new exponents 1/16, 3/16 and 1/4, the fixed interval computation used 1024, 512 and 256 panels on [0,32], respectively. The 12 pilot panels, four per exponent, are included in those counts. Every panel enclosed Taylor coefficients of orders 0–128 at 512 bits with acb_series; the pilot also checked direct/series overlap and Gamma recurrence. All 1792 unique terminal panel records are positive. Parameter-specific Cauchy envelopes bound the raw Taylor errors by approximately 3.509e-18, 1.520e-19 and 1.220e-19; the three normalised tail bounds are below 1.053e-25. Exact outward bounds and the complete q intervals are stored, not just these rounded displays.

Run wall time was 13.878 s and the sum of panel CPU time 75.640 s, approximately 5.45 core-equivalents averaged across wall time. Two sampled worker snapshots showed 15 processes summing to 625% and 528% CPU, and summed worker RSS 502080 and 360112 KiB; summed RSS can include shared pages. Individual peak RSS was 38,813,696 bytes, and sampled system swap use was 1147.56 MiB. Other machine work ran concurrently. These are observations, not a claim of sustained full CPU or a system-wide peak.

The portable checker verifies file hashes, all terminal panel IDs, parameter-specific tails and Cauchy envelopes, R004 certificate hash, independently reconstructs the scalar Radau assembly, and checks the rational endpoint and strict comparison. It does not rerun all acb_series coefficient calculations; run.py is the reproducible panel generator. The calculation has no unregistered extra panels, optimizer, MATLAB, commit or push.

The first checker attempt incorrectly required the conservative dyadic lower endpoint to overlap the interval it bounded. The corrected checker verifies the lower-bound inequality directly and passes; no panel was rerun.
