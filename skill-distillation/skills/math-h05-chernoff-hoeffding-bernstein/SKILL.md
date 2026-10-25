---
name: math-h05-chernoff-hoeffding-bernstein
description: "用于Chernoff、Hoeffding、Bernstein 不等式的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

# H05 · Chernoff、Hoeffding、Bernstein 不等式

> 便携证据定位见[本Skill证据索引](references/handoff-evidence-index.csv)。读取外部文献时，将 `EVIDENCE_ROOT` 设为单独提供的数据包根目录，再按索引的 `resolved_portable_path` 定位文件。文献文件属于额外数据输入，不是本Skill目录内的必需文件。

当前证据状态：**部分可用：基础推导已检查，专题证据仍待完整验收**。全文转换、论文阅读和证明核验是三个独立状态。

## 输入与产出

随机变量/鞅差、独立或条件结构、范围界、方差和矩母函数上界。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 明确单侧/双侧事件和依赖结构；先证明相应矩母函数乘积或条件矩母函数上界。
2. Chernoff：对任意 θ>0 用 P(S≥t)≤e^(−θt)E e^(θS)，再在已证明的定义域中优化 θ。
3. 有界独立变量用 Hoeffding 引理得到范围平方和尺度；中心化且 Bernstein 矩条件成立时用方差和线性尺度得到 Bernstein 界。
4. 将最终阈值、概率水平和常数回代原式；区分独立、条件独立和仅不相关。
5. Linial–Luria P-64dfefc5434002d4允许依赖indicator，但必须验证所有大小k子集的共同成功概率：P(ΣXi≥βn)≤Σ_{|S|=k}P(全1)/binom(βn,k)，k是0<k<βn的整数。KL后果的k=(β−α)n/(1−α)要可行；仅边际均值或pairwise相关不够，任意取整需另证。
6. Hertz P-14bb456bd53bdaeb对EX=0,X∈[a,b],s>0的mgf可取v=|a|b（|a|≥b）或(b−a)²/4（|a|≤b）。独立和Sn正确尾界为exp(−t²/[2Σv_i])；负尾要对−Xi交换区间另取v。原Eq28–36的n归一化错误已用n=2,Rademacher,t=2反例核，禁止照抄。

## 证据与失败处理

界的适用结构和范围/矩母函数假设不能省略；优化参数超出 mgf 定义域无效。 Har-Peled P-c918dd6c7cd3df44只自含证明6M独立faircoins的P(≤M)≤2^−M，QuickSort3/4条件收缩仅remark未细证。Hertz有效部分是单向mgf和Eq23，打印的归一化尾公式不能直接用。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。

H–I 操作步骤与来源定理的逐项证据、证明状态及未覆盖接口见[H–I 操作证据与共享审计文件定位说明](references/portable-audit-access.md)。
