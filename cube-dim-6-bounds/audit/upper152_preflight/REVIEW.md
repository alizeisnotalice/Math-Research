# C152 independent static gate — PASS, concrete Terra packet pending

Reviewed source directory: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round152/laplace_kernel_first_case`. This review reads `PROOF.md`, `GENERAL_2D_INTERFACE.md`, `prepare_for_terra.py`, `single_case.py`, and `HANDOFF_TO_TERRA.md`. No mathematical program, prepare/selection, Arb/FLINT/SymPy call, fixture, integration or preflight was run. Byte-hash/JSON metadata handling is the only executable inventory operation. Actual library compatibility and numerical success remain untested by this reviewer.

## Kernel and closed transform

For the normalized rate-one Gamma(1/2) density, multiplying by `exp(a)` is essential: the untruncated expectation is `exp(-a)` and the full kernel tends to one as h tends to infinity. At finite h the stated kernel is positive and bounded by `p_h=erf(sqrt(h))`; its Gaussian envelope makes the transform finite for every real lambda.

Positive Tonelli interchange is valid. Completing the square in a gives `sqrt(pi*y)*exp((lambda-1)^2*y)*erfc((lambda-1)*sqrt(y))`; multiplication by the original Gamma density cancels `sqrt(pi*y)` and leaves exactly `exp(lambda*(lambda-2)*y)*erfc((lambda-1)*sqrt(y))`. Integration of its derivative gives the displayed signed numerator divided by `lambda*(lambda-2)`. The origin boundary is one, and the derivative's `y^(-1/2)` term is integrable.

The removable branches are correct:

- `I0=h*(1+p_h)-p_h/2+sqrt(h)*exp(-h)/sqrt(pi)`;
- `I2=h*(1-p_h)+p_h/2-sqrt(h)*exp(-h)/sqrt(pi)`.

They imply `I0+I2=2h`; the lambda-one identity is `1-exp(-h)`. The current release evaluates only lambda30. No numerical derivative/degenerate test is implied by this static check.

## Genuine m2 contribution and scale

The source-freeze code requires the exact audited C147 prototype and saved audit manifests, cone zero, pair (16,17), h=1/128 and alpha=1/171992678400. It checks the actual cone's two coefficient rows. The two surviving Exp(1) inputs give `A=E5/60+E6/30`, `B=E5/60`; the four zero columns integrate with probability mass one. The inverse is `E5=60B`, `E6=30(A-B)`, with absolute inverse Jacobian 1800 and support `0<B<A`. Thus the density is `1800*exp(-30A-30B)` and has mass `1800/(2*30^2)=1`.

The function `exp(-30(A+B))*k_h(A)*k_h(B)` is symmetric in its two physical coordinates. Its two complementary triangles have equal integrals; the diagonal has zero measure. Therefore this correlated law has expectation `900*I30^2`. This is a symmetry identity, not independence of A and B.

The complete group-reduced C130 positive expression contributes `(46080/pi^3)*alpha0*900*I30^2` for this one cone and one unordered pair. The factor 46080 is inherited from the already qualified full symmetry normalization. There is no extra 496, second triangle, pattern multiplicity, or column-scale factor. This is a positive part of the full beta2 sum. Its upper endpoint cannot serve as an upper endpoint for the complete sum, and it yields no new M/C6 certificate.

## Accurate enclosure and one-run gate

The lambda30 branch forms the entire signed numerator in Arb512 and divides by 840. It does not clip summands. It then forms 900 times its square and the positive common pi factor and alpha in the same outward interval arithmetic. Exact rational endpoints are extracted for final comparisons. Shared intermediate interval dependencies may widen the enclosure but cannot make it invalid.

The strict tests require positive I and pair/contribution, `I_U<p_L/30`, `E_U<p_L^2`, and widths of I/E/partial contribution at most 10^-80. The first two upper comparisons are conservative sufficient checks using the already proved kernel bound and complete density mass one. There is exactly one released lambda30 transform, one sqrt(h), erf, exp, erfc and pi constant; the final counter equality enforces these intended calls. The code contains no quadrature, pointwise kernel loop, density compiler, MC/native kernel or old functional replay.

`single_case.py` requires the exact source/frozen-input/registration hashes and both review approvals. It validates every frozen source before creating an exclusive started marker, uses nice at least ten, checks the four one-thread environment variables, fixes Arb512/version requirements and stops on failure. Exclusive creation of the marker prevents a second evaluator in the same directory even if two callers race. No precision increase or automatic retry exists. An external Terra observer must supply complete original-PID exit/reap/wall/CPU evidence; the code's internal timing is explicitly before its final receipt. If a process is interrupted, unobserved scalar-call progress is unknown rather than zero. This static review does not create a release or authorize an independent second run.

Terra's source-only preparation must still freeze the concrete input constants, source references and registration and return their exact hashes. Those constants must preserve pair_factor=900, alpha0, h/lambda and common numerator=46080. The final concrete-hash dual gate remains necessary before the unique execution.

## General two-dimensional interface

The future interface correctly retains the ray Jacobian, exact binomial expansion and both Gamma normalizers `p!/d1^(p+1)` and `q!/d2^(q+1)` when both ray rates are positive. Signed coefficients remain signed; error accumulation needs their absolute weights. Gamma rectangle monotonicity and the complete residual probability give valid lower/upper boxes without renormalizing truncated mass.

For a nonpositive ray rate, the stated finite-interval exponential-polynomial recurrence is valid, including its separate exact d=0 branch. Only the complete signed density sum can be intersected with the known positive probability range. The proposed Chernoff tail uses the original nonnegative Exp-column law and an admissible theta; either omitted ray coordinate forces the displayed A+B threshold, so the union bound pays the missing region. The common-pi/weighted pattern-width accounting is valid for complete positive pattern expectations. None of this prepays actual decay, cancellation, coverage, width or computational budgets for all 9961 patterns. The suggested thresholds are future conditional budgets, not results of this single-case packet.

No source149/151 algorithm or live output, lower-bound task, Git index or remote was modified. The final static whitelist contains this report and the source-hash inventory only. Concrete Terra prepared hashes and runtime evidence remain outstanding.
