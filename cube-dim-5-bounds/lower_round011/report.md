# Lower round 011: corrected staggered-row batch

Round 010's claimed 2.341885634028755 is withdrawn. Its source and lock used the wrong second-coordinate extrema; they remain unchanged and hash-identified in `../lower_round010/withdrawal.json`. The three original candidates were registered again against corrected source before any corrected execution. No new grid or model was added.

The independent audit supplied the N=4 block `(j,k,l)=(0,0,(1,4,1,1,1))`, r=u²/2. The old box contains the interior receiver `(u²/4,s,s/4,s/4,s/4)`, but only one of the four selected atoms lies in its cube. Row 3 has y-distance 2s>r. The corrected intersection `[3s-r,r]` is empty. `regression_row_guard.py` actually verifies this example with exact arithmetic.

| tau | Exact value, decimal display | Seconds | Sweep nodes | Retained boxes |
|---|---:|---:|---:|---:|
| 1/8 | 1.767237955729165538391124630583537927257736648741578393 | 11.738660416 | 722 | 160 |
| 1/4 | 1.644984489473170924403645484959616455675896286004086154 | 11.515311416 | 579 | 149 |
| 1/2 | 1.501818584711004131365039638214149168391899747471351916 | 10.132086334 | 250 | 96 |

All three exact ratio intervals lie strictly below the standing complete C. Full algebraic coefficients in Q(u,v), rational isolating intervals, difference intervals and actual timings are in the three `corrected_tau_*_results.json` files. The 120-second/200,000-node limit was respected separately for each candidate. No candidate lock or formal bound promotion is warranted.

The source checks distinct N and enumerates every selected atom of every retained witness; both receiver endpoints are checked against every atom coordinate. The corresponding coordinate-interval check totals are 49180, 46700, 33605 (each entails two inequalities). Two exact union engines agree in every run. Maximum absolute period indices among retained witnesses are (1,1,2,2,2).

The independent audit reconstructs tau=1/8 without producer imports. Its exact coefficient dictionary equals this corrected result; the comparison was executed, not inferred from the displayed decimals. Audit result SHA256: `1ad913f408e9d7adf918e370b974f55d84806a51cecd9274873fca5858ba3932`. Audit source SHA256: `9072a7f83977829a8a1c97b99a6cf92b4357ac41ac3629792ad4834d55e6e9d2`. The minimal counterexample is preserved in this round; the full audit is external to this round and owned by the audit agent.

`mechanism_transfer.md` gives an explicit width-10 period collar, exact finite mass 2(M+20)^5, reflected receiver volume 8 M^5 V_epsilon, and the common-radius choice r'=r-epsilon/2 after inward shrink epsilon. Unit bumps with support delta<epsilon/4 give smooth nonnegative L1 inputs at strict level 1. Taking M to infinity then epsilon to zero proves the finite-union ratio is a valid general-L1 all-scale lower bound. It does not make this below-C batch a new record.

Reproduce from the repository root:

```sh
python3 cube-dim-5-bounds/lower_round011/regression_row_guard.py
python3 cube-dim-5-bounds/lower_round011/corrected_two_row.py 1/8
python3 cube-dim-5-bounds/lower_round011/corrected_two_row.py 1/4
python3 cube-dim-5-bounds/lower_round011/corrected_two_row.py 1/2
```

Corrected source SHA256: `38e61c858c466723a04b3d05e597363d6d8d338bc00605fb251e335f38d8800c`. Result hashes:

- tau=1/8: `4622799d91dc83af85d3ebaae23216a9506aae32d20072b7cb9eebe2e27badc2`
- tau=1/4: `ac21060d8db15cf10f45ccd8e85feb9421182779a653a2490831bae54b049536`
- tau=1/2: `1e4294c870529e3ef9a4ef22850c695bb0a6747649bb0eecaa2e84dd30ecf47f`

Registration, source, baseline, result, regression and transfer hashes are collected in `completion_receipt.json`. Actual local Python was used throughout; no MATLAB was used. The registered batch is complete. This agent made no commits, did not alter the formal ledger, and did not change the original round-010 source/result/lock. No unrelated staged or untracked files were modified.

Open question: better point configurations or a richer valid witness family may improve C, but this fixed three-parameter batch gives no evidence of improvement. This round closes without extending a completed search.
