# I01 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

本项conditional Bernstein矩前提是对每个p≥2、每个过去滤过满足(p!/2)vk c^{p−2}控制；Fan source的epsilon normalization与Freedman c/3不可互换。

## 证明骨架与已检查推导

条件power series与1−cθ geo度量 sum→nonnegative exponential supermartingale→合法终端/停时→Chernoff；t/(v+ct)是可行参数，不称精确minimizer。

本次实际检查：取v=c=t=1给θ=.5，尾bounde^{−1/4}=.77880；actualRademacher尾.5。√(2vx)+2cx验证t²≥2x(v+ct)。

## 正反例与失败处理

正例：令 X=±1 等概率，c=v=1；其 mgf cosh(θ)≤exp(θ²/[2(1−θ)]) 对 0<θ<1 成立，可直接核对级数界的尺度。

条件缺失例：令一行只有 X_{n,1}=n(E_n−1)，E_n∼Exp(1)，mgf=e^(−nθ)/(1−nθ) 仅在 θ<1/n 有限；逐行存在 mgf 邻域不能推出共同 Bernstein 尺度。若除以标准差n，变量变成E_n−1，仍非随n正态化的和，故统一尺度与小单项条件必须分别检查。

## 中心立方体迁移推断

中心立方体迁移接口：若按半径/尺度排序增量，定义对应滤过 F_r 和 X_Q−X_parent(Q)；须逐项验证条件 Bernstein 矩与路径/停止线上的总预测方差上界，不能只用无条件方差。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

无条件moments或随机variance总和没有determinate上界不能给同一tail；v0须独立degenerate处理。

只核验无条件矩会漏掉条件 Bernstein 假设；随机方差和不能当确定 v。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
