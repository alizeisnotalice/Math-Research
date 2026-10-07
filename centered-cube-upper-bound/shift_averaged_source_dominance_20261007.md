# 平移平均来源支配：无平移数费用、packet 容纳概率与尺度权衡

2026-10-07。新增本稿，不改已有稿，不新增数值。以下是直接可测性、Tonelli 与覆盖推导，不调用新的外部定理。一般短壳主目标仍未闭合。

## 1. 冻结完整输入和随机网格

原有限正 Borel 来源为 `μ`，质量 `W>0`，原可测赢家／半径 `R(x)∈[a,b]`，
\[
m(x)=\mu(Q(x,R(x))),\quad u(x)=m(x)/R(x)^n,
\quad E=\{2\lambda<u(x)\le4\lambda\}.
\tag{1}
\]
保持原完整输入、同一 `λ,E,R`；下面的覆盖费不要求 `R` 最大，故适用于有限或连续原允许尺度集。

先固定格边长 `\ell>0`。平移 `zeta∈[0,ell)^n` 按归一化 Lebesgue 概率 `dπ(zeta)=ell^{-n}dzeta` 选取，半开格为
\[
C_{\zeta,k}=\zeta+\ell(k+[0,1)^n),\quad k\in\mathbb Z^n.
\tag{2}
\]
对每个平移它们完整分割原来源；无论平移数或格数，`Σ_k μ(C_zeta,k)=W`。平移不是实际 firstjump 的随机标签，也没有改变原来源概率 `μ/W`。

固定 `η∈(0,1]`，定义
\[
D_\zeta=\{x\in E:\exists k,\ 
\mu(C_{\zeta,k}\cap Q(x,R(x)))>\eta m(x)\},
\qquad p(x)=\int\mathbf1_{D_\zeta}(x)d\pi(\zeta).
\tag{3}
\]
这里分母必须是完整 `m`，不是单格质量或已覆盖子源质量。采用 `≥` 也可，但必须一致分配等号；本稿沿用已有搜索的严格支配、等号留残余规范。

可测性逐项成立：`(zeta,x,y)↦1_{C_zeta,k}(y)1_{Q(x,R(x))}(y)` 为 Borel，积分有限固定 `μ(dy)` 得联合可测 captured mass；可数并集给 `1_D(zeta,x)` 联合可测，再积分 `π` 给 `p` 可测。网格边界用半开规范处理，来源可以有原子；不要求每个平移的来源边界质量为零。

## 2. 平移平均只收 `1/kappa`，不收平移数

令
\[
L_\ell=(1+\ell/a)^n.
\tag{4}
\]
对固定平移，原覆盖证明将 `a/n` 换成 `ell`：支配格总质量 `M_zeta,k` 满足 `M_zeta,k>2ηλR^n`，输出距格中心至多 `R/2+ell/2≤(1+ell/a)R/2`。体积求和得到
\[
\lambda|D_\zeta|\le\frac{L_\ell}{2\eta}W,
\qquad\int_{D_\zeta}u(x)dx\le\frac{2L_\ell}{\eta}W.
\tag{5}
\]
Tonelli 对非负积分严格给
\[
\lambda\int_Ep(x)dx\le\frac{L_\ell}{2\eta}W,
\qquad\int_Ep(x)u(x)dx\le\frac{2L_\ell}{\eta}W.
\tag{6}
\]
对 `κ∈(0,1]`，取真实输出事件
`E_κ={x∈E:p(x)≥κ}`。点态 `κ1_Eκ≤p`，所以
\[
\boxed{\lambda|E_\kappa|\le\frac{L_\ell}{2\eta\kappa}W,
\qquad\int_{E_\kappa}u(x)dx\le\frac{2L_\ell}{\eta\kappa}W.}
\tag{7}
\]
在 `ell=a/n` 时 `L_ell≤e`。这证明“至少 `κ` 比例的预先固定平移使某格支配”的一般输出费只付 `1/κ`，**不付平移数**。无限平移本身没有产生额外费用。

对任何保留 actual 门、来源 `ν≤μ` 的当前 residual 子交通，(7) 仍支付其在 `E_κ` 的完整交通，因为其行响应不超过原 `u`。实际种子可同时保留，再 Tonelli；不能将 `p(x)` 当作其它随机门的独立事件。原 `E_κ` 与已经支付输出取残余交集，不重复收费。

## 3. `Γ_1` 能量：同样只付一次 `1/kappa`

令 `r_t=ae^{t/n}`，`β=nlog(R/a)`，`B=nlog(b/a)`。对任意输出权 `w(x)∈[0,1]`，保持完整来源 `μ/W` 定义
\[
\Gamma_w(t)=W^{-1}\int d\mu(y)\int d\varsigma_n(\omega)
\mathbf1_E(x)w(x)\mathbf1_{\{t\le\beta(x)\le t+1\}},
\quad x=y+\tfrac12r_t\omega.
\tag{8}
\]
cone Jacobian 给 `∫_0^B Γ_w≤eW^{-1}∫_Ewu`。因此
\[
\int_0^B\Gamma_p^2\le\int_0^B\Gamma_p
\le\frac{2eL_\ell}{\eta},
\]
\[
\boxed{\int_0^B\Gamma_{\mathbf1_{E_\kappa}}(t)^2dt
\le\frac{2eL_\ell}{\eta\kappa}.}
\tag{9}
\]
在 `ell=a/n` 时为 `2e²/(ηκ)`。不能只从 `Γ_Eκ≤κ^{-1}Γ_p` 平方再得到不必要的 `κ^{-2}`；正确证明使用 `Γ_Eκ≤1` 及 (7) 的一次 L¹ 费。

另由非负 Tonelli 有精确 `Γ_p=∫Γ_1Dζ dπ`；Jensen 也给第一条能量界。两种证明都始终以原 `W` 归一化，未按支配平移、来源格或覆盖质量重新归一。

## 4. 小包完整容纳概率：乘积而非维数无关常数

固定来源包 `μ_k` 的最小轴包络宽度为 `d_j`，假设 `0≤d_j≤ell`。一个随机平移的网格完整容纳这个包于某一格，等价于每个坐标的格线未穿过该包络的内部；端点遇格线对连续平移概率为零。各平移坐标独立，故完整容纳的概率严格为
\[
\boxed{\kappa_{\rm contain}
=\prod_{j=1}^n(1-d_j/\ell).}
\tag{10}
\]
宽度 `>ell` 时概率为零。对集中测度可用各坐标 essential 上下确界定义宽度；只有零质量端点的开闭差异不影响上述概率。

若仅知包 `side≤c a/n`、网格 `ell=a/n`，则 (10) 给下界 `(1−c)^n`；若每个坐标恰好宽 `c a/n`，就是等号。`c∈(0,1)` 固定时这是指数小的保证，不是常数 `κ`。退化少数坐标宽或 `Σd_j/ell` 很小则可能更好，不能按最大宽度自动认定所有包同样差。

若输出有这份固定包捕获 `>ηm`，在完整容纳事件上该格捕获至少同样质量，所以 `p(x)≥κ_contain`。这是充分条件；**实际格支配概率可以更高**，因为包未完全容纳时某个片仍可能支配。由 (10) 的指数小概率不能断言一般输出必须付指数费用，只能指出“用完整容纳认证平均支配”这条路线的损失。

## 5. 严格 packet 障碍：格支配本身也可能指数稀少

这里给一个仅用于检验平均网格代理的合法 hardband 输入，不研究特殊族的端点。取 `d=c a/n`、`c∈(0,1)` 固定，完整来源为轴盒 `[-d/2,d/2]^n` 上的均匀 `L¹` 密度、质量 `W`。在接收盒
`E_0=[-(a−d)/2,(a−d)/2]^n`，所有 `R≥a` 都捕获全部输入，因此连续 `[a,b]`（或含 `a` 的原尺度集）唯一赢家是 `a`。取 `λ=W/(3a^n)`，就有 `u=W/a^n=3λ`，同一 hardband 成立。

对 `ell=a/n` 的随机网格，某坐标不切包的概率为 `1−c`；切包的概率为 `c`，条件切点在包内均匀。该坐标的最大格内质量份额为
\[
Z_j=\begin{cases}1,&\text{没有切分},\\
\max(U_j,1-U_j),&\text{切分},\quad U_j\sim\mathrm{Unif}[0,1].
\end{cases}
\]
完整均匀 product 包的最大格份额正好是 `Z=∏Z_j`。所以
\[
\mathbb EZ=(1-c/4)^n,\qquad
\boxed{p(x)=\Pr(Z>\eta)
\le\eta^{-1}(1-c/4)^n\quad(x\in E_0).}
\tag{11}
\]
这是解析 Markov 界，不是数值推断。它证明即便整个来源是一个合法 `side≤a/n` 的 packet、完整赢家与 hardband 都成立，固定细网格的平均支配概率仍可能指数小。预先选这个 packet 作为来源分配则由上一稿的直接盒覆盖付常数费；随机格代理不一定保留这份容易的支配。

该例不反驳 (7)、(9)，不反驳一般 `Γ_1` 能量候选，也不认证实际 softFIRST；它只否定“所有小 packet 都有维数无关平均网格 dominance”这种额外猜测。

## 6. 放大格尺度的真实权衡

取 `ell=L a/n`、`L>c`，包每坐标宽 `c a/n`。覆盖膨胀与完整容纳概率分别为
\[
L_\ell=(1+L/n)^n,\qquad
\kappa_{\rm contain}=(1-c/L)^n.
\tag{12}
\]
若只通过完整容纳支付 (7)，其损失是
`L_ell/κ_contain`，对中间尺度近似 `exp(L+cn/L)`。`ell=a/√n` 即 `L=√n`，膨胀本身约 `exp(√n)`，容纳下界约 `exp(−c√n)`；仍不适合多项式或 polylog 主账。

该收费因子还可精确优化。令 `α=c/n`、`δ=ell/a>α`，则
\[
F(\delta)=\frac{L_\ell}{\kappa_{\rm contain}}
=\left[\frac{\delta(1+\delta)}{\delta-\alpha}\right]^n,
\qquad
\frac1n\frac{d}{d\delta}\log F
=\frac{\delta^2-2\alpha\delta-\alpha}
{\delta(1+\delta)(\delta-\alpha)}.
\tag{12a}
\]
分母为正，分子在合法域仅有一个零点；两端 `F` 均趋于无穷，故唯一极小为
\[
\delta_* =\alpha+\sqrt{\alpha^2+\alpha},\qquad
\min F=(\sqrt{1+\alpha}+\sqrt\alpha)^{2n}
=\exp\bigl(2n\operatorname{arsinh}\sqrt\alpha\bigr).
\tag{12b}
\]
代入或使用 `δ_*²−2αδ_*−α=0` 即核验该极小值。对固定 `c>0`，其指数为 `2√(cn)+O(n^{-1/2})`；所以即使允许所有格尺度，**完整容纳收费**的最优结果仍为 `exp(Θ(√n))`。这不推断真实平均支配概率的最优费用。

可给严格范围下界而不只作渐近：若 `c/n<δ=ell/a≤1`，则
\[
\log\frac{L_\ell}{\kappa_{\rm contain}}
=n\log(1+\delta)-n\log(1-c/(n\delta))
\ge\frac{n\delta}{2}+\frac c\delta
\ge\sqrt{2cn}.
\tag{13}
\]
若 `δ≥1`，覆盖膨胀已至少 `2^n`。这只是**完整容纳充分条件所产生的费用因子**的障碍；不把它误称为最优平均格支配的必要费。

## 7. 可用组合与多尺度弱预算

一个一般平移混合引理是：预先固定可数尺度／来源分配方案 `j`，概率权 `τ_j≥0`、`Στ_j=1`，各方案有自己的 `ell_j,η_j` 及平均支配权 `p_j(x)`。令
\[
p_*(x)=\sum_j\tau_jp_j(x),\qquad
A_* =\sum_j\tau_j\frac{(1+\ell_j/a)^n}{\eta_j}<\infty.
\tag{14}
\]
先混合再阈值 `p_*≥κ`，同一 Tonelli 证明给
\[
\int_{\{p_*\ge\kappa\}}u\le\frac{2A_*}{\kappa}W,
\qquad
\int_0^B\Gamma_{\{p_*\ge\kappa\}}^2\le\frac{2eA_*}{\kappa}.
\tag{15}
\]
每份方案都使用同一完整输入；`τ_j` 使多方案重复查看原来源的费用成为概率平均。若改取“任意方案成功”的并集，不再有这份免费平均，而需支付各方案费用之和或另证覆盖重数。

更直接的来源分配多尺度工具通常更强：固定 `Σ_k μ_k≤μ`，第 `k` 包的各坐标轴盒宽为 `d_k,j`，允许不同包尺度，输出支配阈值为 `η_k`。同一中心覆盖逐坐标给
\[
D_k=\prod_{j=1}^n(1+d_{k,j}/a),\qquad
\lambda|E_{\rm dom}|\le\frac12\sum_k\frac{D_kM_k}{\eta_k},
\]
\[
\boxed{\int_{E_{\rm dom}}u\le2\sum_k\frac{D_kM_k}{\eta_k},
\qquad
\int_0^B\Gamma_{\rm dom}^2
\le\frac{2e}{W}\sum_k\frac{D_kM_k}{\eta_k}.}
\tag{16}
\]
这里支配仍是 `μ_k(Q)>η_k μ(Q)`。证明使用 `R≥a` 逐坐标将 `R+d_k,j` 放宽为 `(1+d_k,j/a)R`，再按包质量求和；没有要求包盒互不交或包数受限。

因为 `D_k≤exp(Σ_j d_k,j/a)`，(16) 把可支付性的条件准确化为来源质量加权的 `ΣD_kM_k/η_k` 预算。包宽每坐标 `≤c a/n` 时直接 `D_k≤e^c`，比先随机容纳的小概率工具可靠。对于更宽包，若其总质量相应较小，也可支付；没有这种质量／几何合同则不能免费推广。

以上都是有费用的覆盖支。未达到平均阈值的输出，或每个固定包份额都低于对应阈值的输出，继续保留原完整来源、hardband、赢家和全部 actual 门；无一般结论说这些残余为空。

**状态。** 平移平均 dominance 可严格只付 `1/κ`，`Γ_1` 能量也只有一次该因子；可测性与完整来源账成立。高维 packet 的完整容纳概率是坐标宽度乘积，并且存在实际 hardband 与唯一赢家下的指数稀少格支配例。可用弱进展是 (15) 的预固定混合与 (16) 的来源质量加权多尺度预算；它们没有闭合一般协同余项。
