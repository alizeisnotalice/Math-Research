# J03 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-3ca42b9ac0084b5b eq16–18 lines625–759以有限jump-range Ω和interaction volume ΩI处理外部落点；非local boundary-only。P-fa9183276c325867是Brownian连续出口邻近基线，不证明jump balayage。

## 证明骨架与已检查推导

定义τD与actual Xτ→exit measure在D^c→harmonic/Dynkin pairing需domain与可积性；totalmass还扣nonexit/killed事件。

本次实际检查：0状态吸收到−2概率1/3、2概率2/3，外部g(−2)=1,g(2)=4，u0=3；压到±1且g(±1)=1给错值1。

## 正反例与失败处理

正例：有限状态吸收链的离散退出测度由吸收概率给出，可与解线性调和方程相核对。

条件缺失例：对跳过程把所有退出质量压到最近边界点，或忽略正概率 τ_D=∞ 的质量。

## 中心立方体迁移推断

中心立方体迁移接口：若以跳跃方式搜索 cube-indexed 参数域，定义退出集为整个参数域补集并跟踪越界落点；对空间立方体上的非局部算子，须在整个域外定义边界数据，不能压回几何边界。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

volume constraint source不等于一般Riesz balayage；cube的jump退出不能强压∂D。

非局部跳跃会越过几何边界；仅给 ∂D 数据可能不足以定义非局部 Dirichlet 问题。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
