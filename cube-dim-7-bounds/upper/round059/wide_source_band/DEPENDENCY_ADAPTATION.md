# R059 dependency adaptation and definition-only revision

This revision changes source binding and makes the kernel interface explicit.
No estimate, parameter, scalar coefficient, numerical gate, or proof conclusion
was changed. The proposed result remains subject to independent review.

## Historical run preserved

* Original proof SHA256: `06a3ea9ac8e754ca0f3e63664aea99de20cd49a9c9080c49c8ceb3312e35f8ff`.
* Original and UNCHANGED `certificate.json` SHA256: `291fc41ee7fa40294992fe4e77ef011ebfa1f6657d8e0846278a8f3286a41421`.
* Original and UNCHANGED `certify.py` SHA256: `332cd3ae1cc7a2094a5ffa61d9cbeebd7cd0a9ea319c7410bad8532fb93fea28`.
* The original certificate honestly records local-lemma hash
  `6d1742b7be914df8e7a26fd0a3f8942ae5275ab0923ce6598d87ad7fd3634800`.
  That was the observed source identity at its actual run. No verified
  original-text copy of that identity was retained. In particular, this
  note does not claim that the old external path reproduces that source
  or that this author preserved a verified old byte snapshot.
* No certificate calculation, 64-term atanh check, old grid, or baseline
  integration was rerun for this revision. The historical receipt has
  not been falsely relabeled as a run against the final dependency.

## Current proof dependency, actually read and frozen

The external final lemma was read in full at
`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md`.
Its actual SHA256 is `97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec`.
This existing final local-lemma source is the current dependency. No
additional full duplicate of it is needed for this adaptation.

The local-lemma author independently reported that 6d174 was a genuine
intermediate identity and that reversing only its header/Section-8 source
identity wording in the final text reconstructs 6d174 exactly; the
mathematical core, guards, formulas and constants were unchanged. That
separate provenance check does not by itself prove which file this
R059 author's original run read. The original run receipt remains its
own historical record, and the current proof binding is based on the
actual final text read above.

The current R059 proof is bound to this final local interface as follows:

| Final lemma interface | R059 application checked directly |
|---|---|
| Section 2: d<=1/3, h<=d/8, q>=1/2, R>=1, rho>=4r, a0>0, 0<g<=min(d,R/2), positive zeta term | R059 Section 2 proves all these eventually on the whole matched interval; active R>=Rbase>1 and g=d-h+zeta*t*q<=d are retained. |
| Section 3: a.e. derivative and absolute continuity on [q_b-2h,q_b] | R059 Sections 3--4 supply the new K_w binomial event, positive x<=L, and integrate the bound over the whole interval. The numerical shift is recomputed as 1280 L d h=1280*2^-48; the final lemma's illustrative h=2^-24 value is not imported. |
| Section 4: baseline R=r under the full target free-coordinate law; 3/6/12 moments; full-matching bound 5d/(8rho) | R059 q_b=r-d+h gives R=r exactly. Baseline probability is estimated before cap restriction. Pointwise J bounds are uniform before the later good-target/cap selection. |
| Sections 5--6: full-source inverse moments and fixed captured-face law | R059 Section 5 rederives the two exact integrals, uses J>=1-eps_J, then applies Chebyshev per captured signed face and averages probabilities. No mixture variance or cap-conditioned independence is asserted. |

The final lemma's Section 8 identifies its separate R056 full-proof source.
Its R056-specific finite gates, d=3/10 tuple, and explicit threshold are
NOT imported into R059. This adaptation binds only the local interfaces
listed above; R059 retains its own whole-chain proof and eventual guards.

## Explicit kernel definitions

R059 now defines the full band `mathcal P_0` separately from the scalar
product bound `P0`, defines `G,b_soft` over that FULL band including D's
complement and inactive sources, and defines `P=D intersect mathcal P_0`.
After good targets A are defined, it gives exactly

    K(w,z)=1_A(w)1_P(z)1_matching 1_capture b_soft(z),
    R_lambda(w)=integral K(w,z)d lambda(z),
    C_lambda(z)=integral K(w,z)d lambda(w).

It explicitly states A intersect P is empty eventually, and inactive
source columns vanish. Thus pre-cap probability laws and the final
energy kernel require no implicit convention. No numerical bound changed.

## Final identities for independent audit

* Revised `PROOF.md` SHA256: `fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f`.
* Current external final local-lemma SHA256: `97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec`.
* The exact proof-only diff is supplied in `REVISION.patch`.

Only this owned R059 directory was written. No author/main/other-agent
file was changed, and no Git operation or shared MATLAB was used.
