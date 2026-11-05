# Fixed-location rational weight certificates

Candidate A is `candidate_05.json`: relative to R002 candidate002, atom 51 has mass 4/5000. Candidate B is `transfer_34_27.json`: from A transfer 1/10^6 mass from atom 34 to atom 27. Original atoms are one-indexed. The two original area algorithms and their actual output logs are included; both independently return the same exact fraction for each input. For B their fraction is

4629593106418826487624326959763271900003966948610305702577 / 2718886251002898428308996140604808231510448537600000000000.

The first-core OS PID for B was not captured; the other recorded process information appears in `run-receipt.json` and `*-core2-run.json`. Runs use Dr=2^40 and the same finite reach as R002; these give strict inner-radius certificates, not exact full halos. Descriptive floating metadata inherited from the original seed is ignored; atom coordinates and masses are the input.

`python3 -I round004/weights/transfer_finite.py` and `wolframscript -file round004/weights/transfer_scalar.wl` reproduce arithmetic checks. The original and relocated finite-window arithmetic stdout are byte-identical and deduplicated as `transfer_finite.log`; the relocated symbolic output is `portable-transfer-scalar.log`. Their original copies remain local. The N=10^12 finite-window lower bound checked by these scripts exceeds 1.7027535096. The periodic-window and smoothing transfer is a separate proof, not supplied by these scalar calculations. No local/global optimality or weak-type upper bound is claimed.
