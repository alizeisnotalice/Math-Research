# B024: one released polynomial construction and exact isolation batch

The main agent's explicit release (after the earlier prepare-only assignment):

> F023已final PASS，report SHA2851b03ea55e4ce3e95e7bbbae0b97daa9724a01f0de051aadc586255427544e、manifestc8b3f27cb4a4cada0993e50e89366b7a5f9a4997946f3986a70a31093cf9d70a。主已读B/F完整报告，独立核相位/fold/最高pole/G积分推导，接受stage1。现在明确release B024：先完成本阶段registration和code/inputhash，再唯一有限N/D/P合成+exact正根隔离批，原20min/degree512/≤4positivebrackets/0新K/H/native/失败即停不变。请把此release与F锁hash写入，不必再等主单独release文件；从实际batch start计时。N/D/P code若发现数学接口新不确定先停报。主当前你的/round024/kernel_polynomial目录尚未落盘，故未审核任何未写代码，证明责任由明确验收/checker与后续独审承担，不能声称已主审运行源码。

The main agent also expressly accepted the exact Descartes completeness
shortcut below. Code itself is not claimed to have had a prior main review.
All writes belong only to this directory; B023 and F023 are read-only.

## Fixed inputs and mathematical scope

Original body32 frame and the C022 point x0=(50,50,50,45,27,0)/50. Use only
B023's frozen real principal table, true417-degree factor list, and exact
Gamma_G. F023's independent full-symbolic audit is a required pinned source.
No new point, positive-time kernel, K/H/native call or full symbolic G pass.

Build R(z)=N(z)/D(z), with D=prod_(s>0)(z+s)^m monic and degree417. Combine
local rational terms by a balanced exact tree, not floating interpolation.
Require gcd(N,D)=1, N(0)=0 and degree(N)<=414. Check every D factor positive
for z>0. Then set c=1/4100 and

`ell=3141592653/10^9 < pi < u=3141592654/10^9`,
`P=N-64*ell^6*c*D`.

One Arb512 pi enclosure must prove the displayed strict outer bounds. One
Arb512 G enclosure from the OLD exact Gamma_G/(128*pi^5) is permitted for a
large-time tail bound. This is saved-coefficient postprocessing, not a new G
full computation. P/(64*ell^6 D) is a conservative density-minus-c comparison;
its roots are not asserted to be roots of the actual p-c involving pi.

Clear denominators with a POSITIVE rational scale to a primitive integer P.
Preserve its sign. Exact squarefree factorization must reconstruct it and
give all multiplicities. P(0)<0 and its leading coefficient is negative.

## Exact positive-root completeness and finite gates

Let S be the squarefree part. Require S(0)!=0. Descartes positive coefficient
sign variations V are only an UPPER bound. If V<=4, exact VAS isolation may
produce at most4 positive brackets. When exactly V disjoint positive brackets
have exact endpoint sign changes (or exact rational roots), their lower count
equals V, proving completeness and uniqueness in each bracket. Otherwise,
including V>4, obtain exact `Poly.count_roots(0,oo)` (local source uses Sturm)
before accepting completeness. If V>4, this count precedes any isolation, and
count>4 immediately stops. A count mismatch, endpoint-root ambiguity, absent
exact sign evidence or any other interface problem stops the whole batch.
No approximate root finder or approximate root count is allowed.

VAS uses exact rational `eps=2^-80`, `inf=0`, squarefree S. If a left endpoint
is0, use half the exact reciprocal-Cauchy lower bound for all nonzero roots;
verify the new endpoint signs. Distinct positive brackets may touch only at
a nonroot endpoint. Each bracket's original multiplicity is the exponent of
the unique squarefree factor having its exact root/sign change there.

Map each bracket to a rational positive T bracket by exact integer-square-root
rounding. Try the finite bit list32,64,128,256,512 only; prove lower_T^2<=lower_z
and upper_z<=upper_T^2 and disjoint T brackets. A rational perfect-square
root may remain an exact singleton. No kernel is evaluated in these brackets.

Degree gate<=512, at most4 DISTINCT positive brackets (multiplicities stored),
and1200 seconds TOTAL from actual batch start. A parent Python process kills
and waits its sole worker at that deadline. Stage checkpoints precede and
follow constants, synthesis, integer/squarefree construction, root count/
isolation, and saved-data verification. Sample resident memory at checkpoints
and stop above768MiB; any MemoryError, timeout, domain error, source mismatch,
degree/root gate failure stops with the last completed artifacts. No retries,
precision/gate expansion or second batch.

The minimal checker reconstructs D, checks the N/D/P identity and positive
scaling, squarefree factors, exact root-count evidence, individual endpoint
signs/multiplicities and T-square bounds. It is run inside the same deadline.
If no positive P roots exist, P(0)<0 yields P<0 on z>0, so p_T<c and K_T/T<c
for every T>0 after the source/pi/sign audit. Even then this is a fixed-model
obstruction, not a C6 bound. With roots, this stage does not evaluate K there
or assert an all-time bound.

The actual local API inspection found python-flint0.9.0 `numer`, `denom`,
`gcd`, `factor_squarefree`, and SymPy1.14.0 exact Sturm `count_roots` plus
VAS `intervals`. No target polynomial was constructed during that inspection.
No MATLAB, subagents, git stage, commit or push is authorized.
