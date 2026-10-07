# Frozen tensor soft future-cap diagnostics

Four selected receivers reuse the complete old tensor measure, old lambda and old A-stream seeds. They are selected from hard-band/nonconcentrated short-shell rows before looking at soft responses. No population frequency or dimensional order is inferred.

| n / saved pair | round | all L nodes | visible raw sigma | K/q at finite candidate | K/q scalar bounds* |
|---|---:|---:|---:|---:|---:|
| 128 / pair0 | 1 | 529 | 0.691758176 | 0.432004187 | [0.428778367, 0.435254099] |
| 128 / pair1 | 1 | 524 | 0.453240246 | 0.273001523 | [0.271010986, 0.275006595] |
| 512 / pair0 | 1 | 2043 | 0.161192959 | 0.00168794606 | [0.00163940341, 0.00173792692] |
| 512 / pair1 | 1 | 2040 | 0.241129316 | 0.00303698941 | [0.00294911893, 0.00312747783] |
| 128 / pair0 | 2 | 609 | 0.691768508 | 0.432021797 | [0.431213055, 0.432832045] |
| 128 / pair1 | 2 | 604 | 0.453247745 | 0.273012522 | [0.272513508, 0.273512445] |
| 512 / pair0 | 2 | 2123 | 0.161196506 | 0.00168819949 | [0.00167592873, 0.00170056014] |
| 512 / pair1 | 2 | 2120 | 0.241134269 | 0.003037464 | [0.0030152504, 0.00305984124] |
| 128 / pair0 | 3 | 705 | 0.691771222 | 0.432026382 | [0.431824051, 0.432228807] |
| 128 / pair1 | 3 | 699 | 0.453249621 | 0.273015266 | [0.272890425, 0.273140164] |
| 512 / pair0 | 3 | 2219 | 0.161197408 | 0.00168826332 | [0.00168518706, 0.0016913452] |
| 512 / pair1 | 3 | 2216 | 0.241135526 | 0.00303758416 | [0.00303201518, 0.00304316337] |

*The scalar quadrature remainder is analytic, but the floating function evaluation and guard are not outward certified. These are numerical bounds, not formal machine intervals.

All original coordinate arrival radii remain closed-cube cuts. Every jump has a saved left/right response trace. Targeted rounds retain prior cuts and add up to64 geometric L midpoints per point, besides the nested logL and sigma grids. The raw bisection root is only the last visible crossing of the finite physical-node maximum, not a true last-crossing boundary.

Possible full last-crossing sigma corridors (round3):

| n / pair | last clear finite-node violation lower | full continuous numerical cap start upper | finite K / theory epsilon_n | root-L hard jump |
|---|---:|---:|---:|---:|
| 128 / pair0 | 0.6875 | 0.816771222 | 242.519 | 0.131518 |
| 128 / pair1 | 0.453125 | 0.515749621 | 153.258 | 0.173458 |
| 512 / pair0 | 0.15625 | 0.176822408 | 2.56329 | 0.453629 |
| 512 / pair1 | 0.234375 | 0.248948026 | 4.61196 | 0.228626 |

The final numerical corridors lie above1/n and S0=1/(16sqrt n). This supplies middle-softness candidates for the original shared-input necessary cap layer, without proving actual FIRST/GOOD/score/CP/GP/continuation qualification. The numerical fullcap upper is at a sigma greater than the finite visible root; it is not equality P=q.

K is the original weighted future average with alpha=ceil sqrt n. Positive Bernstein integration eliminates v quadrature: the field is a polynomial of degree n, and weights are alpha/(n+alpha) product_{j=1}^{alpha-1}(k+j)/(n+j). Independent Gauss integration at a polynomial-exact order differs by at most4.94e-15 on the final four candidate pairs.

At the four finite candidate pairs, scalar-lower K/q exceeds the prechosen theory threshold epsilon_n=log(n+2)^(-4). The n128 values also exceed1/16; n512 values do not. This is not evidence of a remaining low-future-average branch. K over the entire possible sigma/L winner box remains UNKNOWN: actual Ls and actual FIRST are not certified, so no actual original output is declared deleted.

Finite-grid future maxima alone only exclude starts when a node is above q; all nodes below q cannot establish continuous cap. The separate cell oracle uses a positive coordinate endpoint envelope in sigma and (u/l)^n P_{sigma,u} in L, retains every hard boundary, and subdivides unresolved sigma cells. Failure or budget remainder stays UNKNOWN. Completed numerical envelopes cover every sigma/L cell, but machine rounding prevents promoting them to continuous certificates.

No old hard energy simulation was rerun. All old source/round/NPZ files were read-only. Three soft rounds, amendments, source coordinates/seeds, complete response grids and jump traces are saved. There is no new weak-type fee or fitted exponent.
