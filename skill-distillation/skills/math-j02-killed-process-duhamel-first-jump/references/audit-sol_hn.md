# J02 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-c580aef7365f8786 Thm3.4 lines531–549在(A1–A4)、weighted Sobolev/boundedκ下给weaksolution和a.e.x Feynman–Kac；Cor5.5 lines1050–1081仅finiteλ域外occupation penalty。

## 证明骨架与已检查推导

合法killed evolution→firstjump互斥分解；Duhamel从operator difference积分导出，并保domain、strongcontinuity/可积性。

本次实际检查：常κ2,t1 survival e^{−2}=.135335，墓地mass1−e^{−2}=.864665，合1。

## 正反例与失败处理

正例：常数杀死率独立于状态时，生存半群乘 e^(−λt)，首次杀死时密度可直接算。

条件缺失例：令跳出域即杀死却继续用完整概率核积分到域内，并不减去墓地质量。

## 中心立方体迁移推断

中心立方体迁移接口：若按首次越过半径、阈值或预算来终止 cube-state 路径，先定义杀死率/墓地状态与真实首次跳标签，再推导 Duhamel 分解；几何立方体边界不等同于跳过程杀死边界。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

finiteλ→∞不自动等first exit：boundary/path occupancy与limit exchange需另证，不能从作者discussion认证。

杀死、越界跳跃和普通转移可能是不同机制；遗漏墓地质量或把边界杀死当一次存活跳跃会破坏质量守恒。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
