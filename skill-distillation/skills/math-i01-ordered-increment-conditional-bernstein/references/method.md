# 有序增量条件 Bernstein 矩：条件、自包含证明与来源边界

## 假设

在滤过 `(F_k)_{k=0}^N` 上，`X_k` 为 `F_k`-可测鞅差，`c>0` 为确定常数，`v_k≥0` 为 `F_{k−1}`-可测。对全部整数 `p≥2`，假定
\[
\mathbb E(X_k\mid\mathcal F_{k-1})=0,\qquad
\mathbb E(|X_k|^p\mid\mathcal F_{k-1})\le(p!/2)v_kc^{p-2}.
\]
并要求 `V_N=Σ_{k=1}^N v_k≤v` a.s.，其中 `v` 是确定常数。可将有限时域换成有证明的停时构造，但不可直接对随机 `V_N` 当作常数优化。

## 条件 MGF 与尾界证明

对任意 `0<θ<1/c`，用条件中心化和矩级数（非负参数下绝对可积由假设保证）得
\[
\begin{aligned}
\mathbb E(e^{\theta X_k}\mid\mathcal F_{k-1})
&=1+\sum_{p\ge2}\frac{\theta^p}{p!}\mathbb E(X_k^p\mid\mathcal F_{k-1})\\
&\le1+\sum_{p\ge2}\frac{\theta^p}{p!}\mathbb E(|X_k|^p\mid\mathcal F_{k-1})\\
&\le1+\frac{\theta^2v_k}{2}\sum_{j\ge0}(c\theta)^j\\
&\le\exp\!\left(\frac{\theta^2v_k}{2(1-c\theta)}\right).
\end{aligned}
\]
所以
\[
M_j=\exp\!\left(\theta\sum_{k=1}^jX_k-
\frac{\theta^2}{2(1-c\theta)}\sum_{k=1}^jv_k\right)
\]
是非负超鞅且 `E M_N≤1`。在 `Σv_k≤v` 下，若 `ΣX_k≥t`，则
\[
M_N\ge\exp\!\left(\theta t-\frac{\theta^2v}{2(1-c\theta)}\right),
\]
故
\[
\mathbb P\!\left(\sum X_k\ge t\right)
\le\exp\!\left(-\theta t+\frac{\theta^2v}{2(1-c\theta)}\right).
\]
取合法参数 `θ=t/(v+ct)`，有 `1−cθ=v/(v+ct)>0`，代入即
\[
\mathbb P\!\left(\sum X_k\ge t\right)
\le\exp\!\left(-\frac{t^2}{2(v+ct)}\right).
\]
这是条件下独立推导，无需增量之间无条件独立。若 `v=0`，则 `v_k=0` a.s.，二阶条件矩为零并使 `X_k=0` a.s.；正阈值尾概率因此为零。

令 `t=√(2vx)+2cx`。展开可得
\[
t^2-2x(v+ct)=2cx\sqrt{2vx}\ge0,
\]
所以尾界至多 `e^(−x)`。这是松弛的常见阈值形式，不是前一分式指数界的等号点。

## Grey (2010) 的实际适用范围

独立完整阅读卡为 [`evidence-card-Grey2010-P-073066f98d171070.md`](evidence-card-Grey2010-P-073066f98d171070.md)。论文 §1 Assumption 1（PDF p.3）假定存在 `θ>0`，使 `q=lim_n E exp(−θS_n)∈(0,∞)`；其伴随随机游走构造用双端 Assumption 2（PDF p.3；所有 `k,B∈F_{−k,k}`，`m,n→∞` 联合极限），Wald 型鞅用单端 Assumption 2*（PDF p.5；所有 `k,B∈F_k`，`n→∞` 极限）。论文没有给出上述逐增量条件中心化或 `p!v_kc^{p−2}/2` 矩控制。

§2.2 Gaussian 例（PDF pp.8–9）还要求 `Σ_{r≥1} r|ρ_r|<∞`，令 `R=Σ_{r≥1}ρ_r`、`S=Σ_{r≥1}rρ_r`，且排除 `1+2R=0`；在这些条件下
\[
\theta=\frac{2\mu}{\sigma^2(1+2R)},\qquad
q=\exp\!\left[-\frac{4\mu^2S}{\sigma^2(1+2R)^2}\right].
\]
这是相关随机游走的特定换测计算，不可迁移成有序增量 Bernstein 结论。Grey 2001 的收敛定理只是其参考文献，此处没有独立重证。

## 当前证明状态

上面的条件矩到尾界链条在列明假设下已逐步推导；这不是把 Grey 当作定理来源，也不声称某个应用已满足输入假设。脚本 `audit_current/hi/run_h02_i01_cases.py` 中一个条件中心化满足的 Rademacher 例与一个共享符号、条件中心化失败的反例均已执行；它们只作实例检查。
6. 条件绝对阶乘矩假设的相邻接口：Fan–Grama–Liu 2012 年论文 P-53739dd293f0de1c Proposition 8.1(III), PDF pp.23–24，证明其 Bernstein 条件与 `E(|η_i|^k|F_{i−1})≤(1/2)k!ρ^(k−2)E(η_i²|F_{i−1})` 等价。该式使用实际条件二阶矩；只有当 I01 取 `v_i=E(X_i²|F_{i−1})` 时才是同一个接口。若 `v_i` 是更大的可预测上界，则 I01 仍成立，但来源等价命题不直接给出它。条件 MGF 级数求和、超鞅构造、确定总预算 Chernoff 界均是本目录独立推导。2015 年 P-048bfe95a33fe7e7 是另一篇论文，其 Theorem 2.6 的条件 MGF 最大值界只作邻近对照，不支持该矩到 MGF 推导。

Package-relative `audit_current/...` locators in this file are resolved with the package-root resolver described in [portable audit access](portable-audit-access.md); they are not Skill-local file links.
