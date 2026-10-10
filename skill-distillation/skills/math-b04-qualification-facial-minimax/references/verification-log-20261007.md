# B04 当前复核记录（2026-10-07）

本记录用于收敛当前入口的修订状态；`verification-log-20261006.md` 保留历史复核，不删除。

- P-552143febc99a746 的 Proposition 2（PDF pp.3–4）与 Theorem 10（PDF pp.14–15）已按其源条件核对。R+×R+ 示例 `A=I,b=c=0` 的原问题可行但无 primal-PPS，而对偶侧满足 dual-PPS；该例仅表明普通可行性不等于 PPS，并不推翻命题。
- 对共享证据源 P-5ae169908ed70fca，PDF p.22 的原式 `sqrt(λ*/σ_K(A*y))` 与推导中的齐次性不符。若 `x*∈λ*∂σ_K(A*y)` 且 `σ_K(A*y)>0`，使用 `α=λ*/σ_K(A*y)` 有 `σ_K(αA*y)∂σ_K(αA*y)=λ*∂σ_K(A*y)`；还需源中共享法向和 Fenchel–Rockafellar 稳定性条件。平方根缩放是来源错误，局部修复未证明该命题在任意更宽条件下成立。
- 精确标量例：`K=[−1,1], A=1, b=x*=2, λ*=2, y=3`，支持值 `σ_K(A*y)=3`。比值缩放给 `y*=2`，使 `J(q)=½q²−2q` 达到唯一最小；印刷根号给 `q=√6`，目标值严格更大。详见 `audit_current/independent-ab/cases/case-results-initial.json` 的 `B01.p5ae_proposition_4_1_rescaling`。
- 当前能力限于明确列出的有限维闭凸锥、PPS、有限值和对偶拓扑假设；无穷维/一般 minimax 扩展仍待另外证明。全库相关性与其余跨主题 source SHA 筛查仍进行中。
