# Independent verification of three Phase A LP certificates

All three frozen cases pass. The checker reads the already audited case and asserts that its SHA-256 matches both the producer LP certificate and the independent Phase A result. Those prior results established completeness of the frozen maximal-clique list; no optimization, step replay, or historical atomic decomposition is performed here.

For each case, exact Fraction arithmetic verifies nonnegativity of m and z, each vertex covering inequality Σ_{C∋i} m_C≥h_i, each maximal-clique packing inequality Σ_{i∈C}z_i≤1, and both objective values. Weak duality gives Mlo=Σh_i z_i≤M≤Σm_C=Mhi, and hence Klo=U/Mhi≤Lopt≤U/Mlo=Khi. The checker verifies the exact recorded endpoints, a nonclique graph, the eligibility margin L0−Khi−10^−4 L0≥0 and the strict failure margin 2 Lfinal−L0−Khi>0. Exact fractions, source hashes and execution resources are in `audit_lp.json`.

| Case | Stop | Eligibility margin / L0 | Failure margin / L0 |
|---|---|---:|---:|
| C/2/05 full | nonedge coordinate stop | 0.6874 | 0.6875 |
| C/1/06 full | 128-step cap | 0.777677777726 | 0.111111111059 |
| A/3/15 partial | 128-step cap | 0.6999 | 0.5 |

The first stop rules out improving redistributions over individual nonedges; it does not rule out adjacent-vertex pair moves or general simultaneous moves. The capped cases are not certified as coordinate-stationary. The audited certificates establish statements for these finite frozen payment instances and algorithm outputs, not a new C_n bound or registered batch frequency.
