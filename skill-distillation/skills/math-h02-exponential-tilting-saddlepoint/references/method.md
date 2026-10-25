# 指数倾斜与鞍点：推导、定理条件与迁移边界

## 1. 精确指数倾斜恒等式

设 `K(s)=log E exp(sS)<∞` 于含 `0,ŝ` 的开区间，且 `K′(ŝ)=x`、`K″(ŝ)>0`。令
\[
\frac{d\mathbb P_{\hat s}}{d\mathbb P}=e^{\hat s S-K(\hat s)}.
\]
对上尾 `ŝ>0`，直接换测得
\[
\mathbb P(S\ge x)=e^{K(\hat s)-\hat s x}\,
\mathbb E_{\hat s}\!\left[e^{-\hat s(S-x)}\mathbf1_{\{S\ge x\}}\right].
\]
这是恒等式；从它到密度在鞍点附近为 Gaussian 还需额外局部极限/正则性。下尾应对 `−S` 使用同一上尾式。把负 `ŝ` 直接代入上尾的 `1/ŝ` 会给出错误符号。

## 2. 已读的条件 Lugannani–Rice 定理（原文原条件）

来源为 Niu, Ray Choudhury & Katsevich, arXiv:2407.08915v3（身份及完整阅读范围见 `references/papers/HK-P-d13d1a9028d29a0f.json`）。令 `W_{1n},…,W_{nn}` 条件于 `F_n` 独立且条件均值为 0，
\[
K_{in}(s)=\log\mathbb E(e^{sW_{in}}\mid\mathcal F_n),\qquad
K_n(s)=n^{-1}\sum_iK_{in}(s).
\]
原文 Theorem 1 的假设是以下两项之一加方差条件：

- CSE：存在固定 `β>0`，使条件尾由 `θ_n e^(−βt)` 控制，`θ_n` 为 `F_n`-可测、有限 a.s. 且 `θ_n=O_P(1)`；或
- CCS：`|W_{in}|≤ν_{in}<∞` a.s.，其中 `ν_{in}` 为 `F_n`-可测，且 `n^(−1)Σ_iν_{in}^4=O_P(1)`；并且
- `n^(−1)Σ_i E(W_{in}²|F_n)=Ω_P(1)`，阈值 `w_n` 为 `F_n`-可测且 `w_n=o_P(1)`。

在概率趋于 1 的事件上，`K_n′(ŝ_n)=w_n` 在原文固定邻域内有唯一根。定义总 CGF `Λ_n(s)=nK_n(s)`、总和阈值 `x_n=nw_n`，并记
\[
\lambda_n=\hat s_n\sqrt{\Lambda_n''(\hat s_n)}
 =\hat s_n\sqrt{nK_n''(\hat s_n)},\qquad
r_n=\operatorname{sgn}(\hat s_n)\sqrt{2\{\hat s_nx_n-\Lambda_n(\hat s_n)\}}.
\]
根号内为非负时按此式定义；原文对其余分支及 `w_n=0` 给出符号/连续约定。记 `Q(r)=1−Φ(r)`。原文 Theorem 1 结论为
\[
\mathbb P\!\left(n^{-1}\sum_iW_{in}\ge w_n\mid\mathcal F_n\right)
=\left[Q(r_n)+\phi(r_n)\left(\frac1{\lambda_n}-\frac1{r_n}\right)\right](1+o_P(1)).
\]
括号内是完整 LR 近似；误差是相对 `o_P(1)`，没有给出确定误差常数或收敛速率。原文的条件独立数组定理不要求额外非格点条件，但只处理 `w_n=o_P(1)`；它不覆盖固定非零 cutoff 大偏差、一般依赖和、studentized 统计量或置换依赖。对应限定见原文 §5。

原文 Appendix H 给出 `λ_n/r_n→_P1`（在正上尾分支并满足其定理条件时）。这项比例结论不能脱离定理假设单独调用。

## 3. 从 LR 到简化前因子：仅在正且发散的 `r_n` 区间

在定理适用且 `w_n>0` 以概率趋于 1 时，`ŝ_n>0,r_n>0`。令 `M(r)=rQ(r)/φ(r)`。精确代数恒等式为
\[
\frac{Q(r_n)+\phi(r_n)(\lambda_n^{-1}-r_n^{-1})}
     {\phi(r_n)/\lambda_n}
=1+\frac{\lambda_n}{r_n}\{M(r_n)-1\}.
\]
对 `r>0`，Mills 界给出
\[
\frac{r^2}{r^2+1}<M(r)<1,
\]
故 LR 与 `φ(r_n)/λ_n` 的比值偏低于 1，且其相对亏损至多 `λ_n/[r_n(r_n²+1)]`。若进一步 `r_n→_P+∞`，并用原文的 `λ_n/r_n→_P1`，则 LR 与 `φ(r_n)/λ_n` 相对等价。又因
\[
\frac{\phi(r_n)}{\lambda_n}
=\frac{\exp\{\Lambda_n(\hat s_n)-\hat s_nx_n\}}
 {\hat s_n\sqrt{2\pi\Lambda_n''(\hat s_n)}},
\]
这给出 Niu 定理模型内的简化上尾前因子。它是附加 `r_n→∞` 后的渐近式，仍没有定量总误差率。对 `r_n=O_P(1)`，该化简无效；中心例子中原 LR 修正项是必要的。

注意这与经典固定非零偏移的大偏差 saddlepoint 展开不是同一套条件。Niu 结果的 cutoff 趋于零；固定 `ŝ≠0` 时，格点修正、非格点密度、局部极限定理和误差条件须由对应的大偏差来源另行提供。该迁移不得只凭同一指数速率完成。

## 4. Gaussian Mills 界（独立初等证明）

令 `Q(u)=∫_u^∞φ(t)dt`，`u>0`。因为 `t/u≥1`，
\[
Q(u)\le u^{-1}\int_u^\infty t\phi(t)\,dt=\frac{\phi(u)}u.
\]
分部积分给出
\[
Q(u)=\frac{\phi(u)}u-\int_u^\infty\frac{\phi(t)}{t^2}\,dt,
\qquad
\int_u^\infty\frac{\phi(t)}{t^2}\,dt\le\frac{Q(u)}{u^2},
\]
从而
\[
\frac{\phi(u)u}{u^2+1}<Q(u)<\frac{\phi(u)}u.
\]
所以简化 Mills 首项 `A(u)=φ(u)/u` 是**上界**；相对首项的亏损满足 `0<(A−Q)/A<1/(u²+1)`，而相对真值的超额满足 `0<A/Q−1<1/u²`。仅当 `u→∞` 时 `Q(u)∼A(u)`。这是标准 Gaussian 的精确界，不是任意倾斜变量的有限样本界。

## 5. 记录中的数值反例/示例

`audit_current/hi/run_h02_i01_cases.py` 有标准正态中心 `r=1` 与 moderate-tail `r≈31.62` 两例，以及一个 Rademacher 格点数组实例。中心例子中简化项比精确尾大约 52.5%；大 `r` Gaussian 例子的相对误差约 0.10%。格点例子的偏差约 3.1%，显示有限样本格点效应可以明显大于连续 Mills 尾项。它们都是输入级检查，不用于证明第 2–4 节命题。

Package-relative `audit_current/...` locators in this file are resolved with the package-root resolver described in [portable audit access](portable-audit-access.md); they are not Skill-local file links.
