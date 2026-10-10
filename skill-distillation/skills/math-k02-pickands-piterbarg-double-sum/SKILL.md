---
name: math-k02-pickands-piterbarg-double-sum
description: "用于Pickands–Piterbarg 双重求和的数学研究：核对输入与假设，组织方法和证据；不凭名称补造内部接口。"
---

# K02 · Pickands–Piterbarg 双重求和

> 审核边界：本流程整理给定对象的定义、前提、推导与证据；它本身不构成专题定理证明。数值或玩具例只核对对应实例。来源结论须在账本列出的条件内使用，未独立核验的迁移必须标注“未证迁移”。本轮逐条断言、原文定位与实际案例记录见工作包 `audit_current/jk/`。
> 文献定位见 `references/handoff-evidence-index.csv`；完整原文应以工作包 `evidence/papers/<paper_id>/paper.pdf` 及其 SHA-256 为准，自动转换文本只作检索辅助。

## 输入与产出

连续高斯过程/场、方差函数、局部相关展开、索引域与高阈值。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 验证样本路径可测/连续、方差上界与最大方差点结构；记录局部相关 1−r(t,t+h) 的幂指数 α、慢变项和常数。
2. 按相关长度与阈值选择小块，使用 Pickands/Piterbarg 单块超越渐近并保留对应常数 H_α/Piterbarg 常数。
3. 执行double-sum时先固定block长度参数T并令阈值u→∞，再令T→∞；分别给相邻边带、中距和固定紧域远距交集预算。对有限块并集，除非联合律另给独立性，必须以实际交集概率处理 Bonferroni；两边际相同或各自标准正态不蕴含独立。只有在已核验块事件相互独立时才可写 P(E₁∪E₂)=2p−p²。P-6c5bcc98b64321b6 Lemma5只写固定t₀的sufficiently-large-u，不能直接代入随u增长的中距索引；该uniform pair界尚待补。
4. 核对域维数、网格稠密度与连续上确界版本；常数、指数项和多重峰贡献分别列出。

## 证据与失败处理

从独立安装位置访问文献 PDF 或共享审计记录时，先按 [便携证据访问指南](references/portable-audit-access.md) 定位证据包并调用共享解析器。

常数取决于局部相关、方差峰和索引维度；单块渐近不能在无联合界时直接乘块数。 Michna首个完整double-sum来源的局部条件极限需uniform entropy/Borell domination，不能省固定域依赖；远距相关隙应在紧区间[ε/4,p]取min，不能用全半轴inf。Theorem1的Hα显式下界只有外引证明；作者整体证明的中距uniform量词尚未独立补齐。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
