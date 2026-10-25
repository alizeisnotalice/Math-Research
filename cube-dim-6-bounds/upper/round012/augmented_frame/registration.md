# C2/012: one augmented positive frame, registered before calculation

Owner: curvature-astra. Only this directory may be written. Round010 is
frozen. No agents, MATLAB, stage, commit, push, spatial search, or new
physical closest-S2 enumeration.

The new symbol is psi_beta=psi_body32+beta sum_i |xi_i|. The sole numeric
candidate is beta=1, T=2/45. This adds six independent coordinate-axis
Cauchy variables. The old geometric fan may be reused only after exact
ray/rate, determinant, sign-linearity and original coverage checks.

Stage 1 (1200 seconds maximum; zero kernel passes): replay the existing
exact fan verifier into this directory; compare every integer ray to its
original rational ray; check all 126 ray occurrences satisfy nonnegative
ordered coordinates, q_old=psi_body32(ray), and
q_new=q_old+sum(abs(ray))=psi_full38(ray). Check determinants and the
38-column Gram matrix. The input manifest is saved before this stage.

Stage 2 (1200 seconds maximum total; at most four full kernel passes):

1. At the existing B010 first point with denominator 10^12 and numerators
   (1000000000000,999999990000,966375016425,851798538049,
   486751693199,52461218363), compute G_beta once by the full exact
   fmpq Laurent route, then H_beta,T once by the C010 Arb512 engine.
2. If the upper enclosure for K=G-H is at most T/4096=1/92160, stop
   immediately: this candidate cannot give a bound below 64 through the
   short-potential domination interface.
3. Only if the first point passes that necessary condition, use the
   already exact B010 second point, denominator 50000000350000000000000,
   numerators (50000000350000000000000,49999999849999996500000,
   48392387637670880863725,42634261204591078083493,
   24357044587919354472043,2578553519241654955591).
   Compute G_beta once and H_beta,T once. Its denominator is outside
   int64: use arbitrary-integer fmpq/Python Arb, never rounding the point
   or passing it to the C++ engine. Stop after this second point regardless.

All zero/repeated-node terms and all principal-part coefficients are
retained. No original-frame G or H is recomputed. Complete-group counts
must be 967680 per pass. The native binary's exact source/binary/runtime
hashes and input range bounds are checked before use. Failures and timeout
stop the batch. No replacement point, beta or T is authorized.

The separate sufficient target m0=T/4000=1/90000 would give 62.5 only
with a global essential lower certificate. Pointwise survival never proves
global m, optimality, or any C6 improvement; the established bound is 64.
