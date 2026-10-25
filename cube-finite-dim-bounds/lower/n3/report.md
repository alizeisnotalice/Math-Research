# N3 two-orbit batch 001–002 (2026-10-25)

**No new bound.** The exact inherited R3 remains the comparison target; its original algebraic certificate is read directly. This work changes periodic positions, the mass split and all three periods together, rather than extending the already audited same-lattice counts.

## Construction and complete finite-witness certification

Take period vector P and two orbits P*(k+phi_b), k in Z^3, with masses w_b per atom; multiplication is coordinatewise. The certified rationalized survivor is

- P=(3.275875394988,1.072284832359,1.123649780712), exact terminating decimal rationals;
- phi_0=(0,0,0), phi_1=(0.498788116727,-0.000110812266,0.000382337232);
- w_0=0.983352091252, w_1=254161977187/250000000000, and w_0+w_1=2 exactly.

`certificate_input.json` specifies all 143 integer block-pair witnesses. Its SHA-256 is `a56bd687f588a95a10415739c8554b691ba6c024c91bdd5450a982f3a47b1a25`. Each key gives two triples of inclusive lower indices followed by two triples of upper indices. An orbit with any nonpositive side count contributes no atoms. For every witness, compute its selected mass M=sum_b w_b product_j max(0,h_bj-l_bj+1). All selected atoms fit into a cube centered at x exactly when each coordinate lies between

max_b P_j(h_bj+phi_bj)-r and min_b P_j(l_bj+phi_bj)+r,

where only nonempty blocks enter and r=M^(1/3)/2. Intersect each resulting receiver box with the **full** cell product[0,P_j]. The whole three-dimensional union is evaluated jointly; no reflection or product of marginal volumes is used.

`certify.py` uses Fraction throughout. For every rational mass, it brackets its positive cube root by q/10^24 and (q+1)/10^24 and checks both cubed integer inequalities. Shrinking/expanding every radius gives nested inner/outer receiver-box unions. A recursive exact coordinate sweep computes both volumes. Their mass-normalized enclosure is

[896511535967395056495143017663436374008572775667352485988199650559033313 / 500000000000000000000000000000000000000000000000000000000000000000000000,
 896511535967395056495144326227891673272183365588325794424157321414799229 / 500000000000000000000000000000000000000000000000000000000000000000000000].

This is approximately 1.793023071934790113. Its upper endpoint is strictly smaller, by a positive exact Fraction recorded in `certificate_result.json`, than the original R3 certificate's lower endpoint. Certification takes about 0.125 s. This is a complete enclosure of the stated finite witness union, **not** an upper bound on the complete maximal-function superlevel set at this point.

For strict superlevels and general L1, contract every retained receiver box by tau on each face, replace its threshold radius r by r-tau/2, and then smooth each selected atom by a same-mass nonnegative bump supported in infinity-radius less than tau/2. The original farthest distance is at most r-tau, so containment survives and (2r-tau)^3<M. There are finitely many witnesses, so tau can be common and their contracted union volume tends to U. Exactly two inner radius bounds are critical equalities `(2rlo)^3=M`; neither closed box is claimed to lie in a strict superlevel without this contraction. To make input mass finite, use receiver-cell translation indices q in [-N,N]^3, so each physical receiver cell is P*(q+[0,1]^3). The 143 keys have active-atom index_min=(-1,-1,-1), index_max=(1,2,2); thus b=max_j(max(abs(index_min_j),abs(index_max_j)))=2. A witness in translated cell q uses atom indices in q+[-1,-1,-1] through q+[1,2,2], all contained in [-N-2,N+2]^3. Retain both original orbit representatives P*(k+phi_0), P*(k+phi_1) for every integer k in that latter index block. Every index pair has mass w_0+w_1=2, so total finite input mass is exactly 2(2N+2b+1)^3=2(2N+5)^3; receiver volume is (2N+1)^3 U_tau. The negative transverse coordinate of phi_1 is kept: index pairs are grouped by their orbit label k, not by physical membership in the receiver cell. Thus an atom can lie slightly outside its named physical cell without affecting this exact mass count or containment. The ratio is at least [(2N+1)/(2N+5)]^3 U_tau/2. Hence N→infinity then tau→0 yields C3 >= U/2. `index_collar` in certify.py computes and asserts these bounds from all active key blocks, and certificate_result.json records index_min/index_max/b, pairing convention and both finite formulas. The actual endpoint check covered every active block and receiver extrema for N=0,1,7; no search or union-volume rerun was needed. The same argument applies to the inner rational union, whose certified lower endpoint is therefore a rigorous lower bound, albeit weaker than inherited R3.

## Numerical discovery and retained near candidates

Batch 001: 1124 seven-parameter evaluations, 4096 scrambled Sobol receivers (seed 202610253), then 131072 independent scrambled Sobol receivers (seed 202610254) for eight geometrically separated retained candidates plus the baseline. Periods were Pbase*exp(z_0,z_1,z_2), Pbase=(2a1,a2,a3) from the original formula, relative phase=(1/2,0,0)+(z_3,z_4,z_5), masses=(1+z_6,1-z_6). Every sample used actual common-radius cube counts in both periodic orbits, across all entry radii <=2.8. Bounds were z_0:2∈[-.035,.035], z_3:5∈[-.08,.08], z_6∈[-.18,.18]. The ±3 entry-index window covers those radii throughout these period bounds; radii beyond 2.8 were not searched. Wall time was 11.95 s.

The apparent sampled best 1.794035720467965 fell to 1.7887069308891985 on independent receivers and to 1.7886170081363033 in its extracted full joint witness union. Therefore it was not treated as a gain. Witnesses were extracted using a third seed 202610255 and 8192 receivers. Across the nine retained records, the eight distinct finite-union values were:

|Construction|Joint finite-union ratio (floating)|
|---|---:|
|Inherited baseline|1.7939332277947604|
|Perturbed 1|1.7886170081363033|
|Perturbed 2|1.7895924820873836|
|Perturbed 3|1.7889105716650595|
|Perturbed 4|1.7886446839449135|
|Perturbed 5|1.7886195996897893|
|Perturbed 6|1.7871302433611782|
|Perturbed 7|1.7880289795886637|

The duplicated baseline was retained in raw validation selection; its full union agrees with the original algebraic formula to 9e-15. All nonbaseline parameters and validation counts remain in the compact JSON results. The sampling gaps are material, and no sample-based exclusion is made.

Batch 002 was registered after this validation and before running. The union of all extracted keys gives a finite family of 143 witnesses. Six starts jointly optimized all seven parameters by Nelder–Mead, exactly 300 evaluations each (1800 total), using full floating coordinate-decomposition union volumes at every evaluation. Wall time was 7.02 s. Every run stopped at its evaluation cap; none establishes local optimality. Best nontrivial value was 1.7930230719346778, subsequently rationalized and certified above. The search includes receiver shifts across full-cell boundaries and unequal masses.

## Independent local software checks

Python actually ran both batches, extraction, the Fraction verifier and a repeat check after packaging. Python uses the provided local venv with NumPy/SciPy; certification uses only the standard library. MATLAB actually ran the original algebraic formula and returned 1.7939332277947688. `crosscheck_matlab.m` independently decoded the rational survivor, formed all 143 boxes and computed their full joint union by a three-dimensional difference grid; it returned 1.7930230719347879. These MATLAB values are floating checks, not the exact comparison proof. MATLAB variables use n3 prefixes and the shared session was neither cleared nor repathed.

## Proof-status separation and next gap

- **Numerical candidates:** both search batches and the seven initial perturbed unions. No surviving numerical value exceeds R3. Their full superlevel sets remain incompletely explored.
- **Local sufficient conditions:** the explicit selected blocks plus the max/min endpoint condition above certify containment for each witness; the uniform contraction and collar argument gives strict L1 transfer.
- **Complete certification:** exact enclosure and strict comparison for the one rational survivor's stated 143-witness union. Baseline comparison uses its original exact rational interval, never reconstructed decimals.
- **Excluded:** only the proposition that this particular certified finite union exceeds R3. Neither this configuration's complete bad set nor any whole two-orbit family is excluded. No claim of local or global optimality.

A next distinct finite experiment would use at least three shifted layers in a three-cell superperiod, including nonuniform layer masses and independent transverse phases; enumerate receiver arrangement cells or iteratively add uncovered cube witnesses before optimizing. This could change nearest-neighbor geometry beyond a single alternating pair. No next batch is running.

## Reproduction and file scope

Use `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python` with `search.py`, `witness.py`, `refine.py`, `certify.py` in that order. Search and refinement regenerate compact numerical summaries. By default `certify.py` reads the saved self-contained `certificate_input.json` directly and verifies it without importing or rerunning search/refinement. The optional `--rebuild-input` flag reconstructs that input from the saved floating survivor and witness family; it is not required for verification. Run `crosscheck_matlab.m` in MATLAB for the independent floating check.

Integration whitelist: `registration.md`, `report.md`, the five Python/MATLAB code files (`search.py`, `witness.py`, `refine.py`, `certify.py`, `crosscheck_matlab.m`), `results.json`, `witness_results.json`, `witness_family.json`, `refine_results.json`, `certificate_input.json`, `certificate_result.json`, `matlab_run.txt`, `source_manifest.json`. Temporary logs and Python caches were removed. No STATE/ledger changes, commit or push were made; no files outside lower/n3 were written.
