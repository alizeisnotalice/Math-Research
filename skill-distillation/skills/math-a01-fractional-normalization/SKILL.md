---
name: math-a01-fractional-normalization
description: "用于齐次归一化与分式规划的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> 外部文献、claims/cases 与 PDF 属于完整证据包，不随单个 Skill 安装。请使用[便携证据访问指南](references/evidence-guide.md)，按本 Skill 的[来源索引](references/handoff-evidence-index.csv)通过统一 resolver 定位。registry/显式包根路由均不依赖当前工作目录；SHA 检查只确认当前字节，不等于全文阅读或数学验收。


# A01 · 齐次归一化与分式规划

当前证据状态：**部分可用**。支持清单中明确假设下的单比值CCT、射线放松、SOS层级及题型专属审计；不自动覆盖零/变号分母、不同不确定性结构或课题结论。P-fdba 弦/切线界须用已修正的覆盖域和方向；其余引用范围见逐篇卡。当前全库相关性重筛仍在进行。

## 输入与产出

可行集 X、分子 f、分母 g、优化方向、齐次/凸结构声明，以及中心立方体最大算子的函数与阈值接口。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 先定义实际单比值f/g、可行集K与优化方向；逐可行点检g>0，零分母与趋零面另处理，不由fractional名称推断算法。
2. 精确CCT设s=x/g(x),t=1/g(x)>0，s/t∈K且tg(s/t)=1，目标tf(s/t)，逆映射x=s/t逐约束核对。线性与一般多项式perspective表示不可混用。
3. 全文来源P-33450c5f81cf2615把tg=1松成tg≥1时还要求f≥0，沿射线τ≥1目标按τ缩放才不改善最小值；若f可负，该松弛可能无界。值相等不等于最优点集合相同。
4. 其凸多项式达到定理需要K非空、f≥0、g>0、f/−g/h_i凸且sup_Kg<∞；CCT强对偶另检原约束Slater，不能从正分母或凸性单独推。
5. 其SOS渐近层级还需有界K、Slater及g_min>0，保L>域半径、C>2/g_min的罚常数与每层下界方向；无通用速率。首层精确仅SOS-convex版本，普通凸不够。
6. 矩恢复需来源具体rank-one与非零y₀条件或SOS-convex Jensen；浮点近rank-one不认证可行/最优。固定λ residual f−λg的方向保留，不对鲁棒sup/inf任意交换。
7. 若实际问题是binary多比值、鲁棒比值、proximal迭代、projective凸包或chance constraint，先读[分模型来源与定位](references/frozen-source-locators.md)，按题型打开相应全文阅读卡；该表明确标注P-334/P-fdba之外的五篇来源。保共同binary单项式、相关不确定性和正分母；P-fdba03cbf10513db的目标曲率/尾分位缺口及弦包络收紧约束的最大化下界方向另列，不沿用其卡中upper推论。
8. 齐次归一化先检缩放保全部约束；中心立方体接口固定M_c f(x)=sup_{Q:center(Q)=x}|Q|⁻¹∫_Q|f|，维数、阈值、尺度和同输入耦合保留，CCT不证明其解析界。

## 证据与失败处理

单比值正分母的精确CCT、f≥0的射线放松、sup g有限的达到、Slater强对偶与SOS-convex首层精确是不同层次。外引凸多项式达到/SOS扰动/Jensen依赖仍未独立全证；鲁棒量词、零分母和私有耦合缺口不自动消失。 多比值不能各项独立凸化后宣称共同可实现；proximal驻点不默认全局最优，Gaussian chance reformulation的凸性和近1分位尾域仍有来源缺口。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。