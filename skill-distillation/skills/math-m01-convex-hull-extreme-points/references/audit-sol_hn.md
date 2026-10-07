# M01 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-3d2c7d539b055402 lines138–146/177–182是locallyconvex紧metrizable集的Krein–Milman–KyFan exposed/affine结果；Carathéodory有限维与Choquet integral需各自版本。

## 证明骨架与已检查推导

R^d有限凸组合多于d+1时以affine依赖adjustweights直到一项0并迭代；infinitedimcompact只给闭凸包，closure不可省。

本次实际检查：三角形(1/3,1/3)由三个顶点各1/3；[0,1]Lebesgue概率不可能finiteDiraccombo（有限原子处Lebmass0）。

## 正反例与失败处理

正例：三角形 conv{(0,0),(1,0),(0,1)} 中，(1/3,1/3) 是三个顶点以 1/3 权重组成的凸组合；R^2 中 Carathéodory 至多需要 3 个顶点。

条件缺失例：[0,1] 上概率测度的弱星紧凸集，其极点是 Dirac 测度。Lebesgue 测度属于这些极点的闭凸包，但不是有限个 Dirac 测度的凸组合。

## 中心立方体迁移推断

只有在可行凸集具备相应紧性和拓扑条件时，才用极点化简凸松弛或矩集合。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

拓扑、紧性、metrizability不自动属于目标可行域；名字Krein space不证明convex theorem。

无限维 Krein--Milman 一般只给闭凸包，不保证有限凸组合。Krein 空间是带不定内积的空间，不能凭同名认定为 Krein--Milman 结果。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
