---
name: math-i01-ordered-increment-conditional-bernstein
description: "在有序滤过下核验条件 Bernstein 矩条件，并推导可选指数超鞅和尾界；不把无条件矩或相邻文献当成条件假设。"
---

# I01 · 有序增量条件 Bernstein 矩

当前证据状态：**假设条件下的自包含推导已核；外部文献只作相邻背景**。此结论不表示某一课题已满足这些假设。Grey (2010) 已由本轮独立逐页重读，但只支撑特定平稳增量过程中的换测鞅背景，不能支撑下面的 Bernstein 矩条件或尾界。来源与阅读状态见[来源与阅读状态](references/provenance.md)。旧入口存档于 `references/history/SKILL.pre-evidence-revision-20261007.md`。

## 输入与结论

给定有限序列与滤过 `(F_k)_{k=0}^N`。对每个 `k`，要求 `X_k` 为 `F_k`-可测、`c>0` 为确定常数，`v_k≥0` 为 `F_{k−1}`-可测，并且
\[
\mathbb E(X_k\mid\mathcal F_{k-1})=0,\qquad
\mathbb E(|X_k|^p\mid\mathcal F_{k-1})\le\frac{p!}{2}v_kc^{p-2}
\quad(p=2,3,\ldots).
\]
若 `Σ_{k=1}^N v_k≤v` a.s.，其中 `v>0` 为确定常数，则对 `t>0`
\[
\mathbb P\!\left(\sum_{k=1}^NX_k\ge t\right)
\le\exp\!\left(-\frac{t^2}{2(v+ct)}\right).
\]
同一形式用于负尾：对 `−X_k` 应用假设并检查条件。若 `v=0`，假设强制所有增量为零 a.s.，正阈值尾概率为 0。随机方差和不能冒充确定 `v`；只有经有效停时/可预测控制后才能代入相应界。

## 执行步骤

1. 逐项检查条件中心化、条件绝对矩、可测性、可预测 `v_k`、共同确定尺度 `c`，以及方差和的控制方式。无条件矩检查不能替代条件检查。
2. 按[方法说明](references/method.md)从指数级数推出 `0<θ<1/c` 时的条件 MGF 界，构造非负超鞅，并应用 Chernoff。采用 `θ=t/(v+ct)` 可得上述显式界；它是方便的合法选择，不声称为对原始分式指数目标的精确最优解。
3. 若以 `x>0` 参数化，阈值 `t=√(2vx)+2cx` 满足 `t²≥2x(v+ct)`，因此尾概率不超过 `e^(−x)`。需要 `c/3` 常数、随机方差界、向量/矩阵版本或停止时间版本时，另行核验相应假设和定理。

Fan–Grama–Liu 2012 年论文只给本项条件矩接口的受限对照：其 Proposition 8.1(III) 的右端使用实际条件二阶矩 `E(X_k²|F_{k−1})`。若本项的可预测上界 `v_k` 严格大于该矩，二者不是同一个假设；当前 MGF 级数求和、超鞅与 Chernoff 链条完全按上述 `v_k` 自行推导。另一个 P-048bfe95a33fe7e7 是 2015 年论文，其 Theorem 2.6 可作条件 MGF 最大值界的邻近对照，但不证明当前的矩到 MGF 步骤。

## 失败处理

条件中心化缺失时不得构造上述超鞅；只知道无条件矩时不得通过；`Σv_k` 无确定/停时控制时不得代入确定 `v`。用[适用案例](examples/positive.md)和[条件缺失案例](examples/negative.md)检查执行边界。测试案例只检查实例，不证明一般命题。

H–I 操作步骤与来源定理的逐项证据、证明状态及未覆盖接口见[H–I 操作证据与共享审计文件定位说明](references/portable-audit-access.md)。
