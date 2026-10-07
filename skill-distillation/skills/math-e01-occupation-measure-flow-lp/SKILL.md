---
name: math-e01-occupation-measure-flow-lp
description: "用于占用测度 LP 与流守恒的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# E01 · 占用测度 LP 与流守恒

当前证据状态：**待证据验收：可执行方法工作流，不是已证明的课题结论**。全文转换、论文阅读和证明核验是三个独立状态。（此行声明已被下方【修订】标注替代：该 skill 经外部审核提升为「部分可用」，详见对应修订行与 references/verification-log-20261006.md）
  【修订 2026-10-06（外部审核 R26）】既有验证记录支撑档位提升：VP 有效域文献确认、Cantelli 紧性完整证明、VP 极值取等、CVaR/Thm5.5 双反例复算（4 条子断言）。证据状态提升为：部分可用：基础推导已检查，专题证据仍待完整验收。源定理证明本身未重证（守卫保持）。

## 输入与产出

状态空间 S、动作集 A、初始分布 μ₀、转移核 P（或生成元）、折扣因子 γ/时域、奖励或成本 r(s,a)。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 先选有限γ折扣MDP、定时域确定性控制、扩散/跳控制、停止风险、均场相互作用或变分域模型；时间/折扣/停止occupation的质量及界方向不同，均场二次目标不默认纯LP。
2. 有限MDPγ∈[0,1)、μ₀概率：x(s,a)=EΣγᵗ1_(sₜ,aₜ)≥0，检Σ_ax(s,a)−γΣ_(s′,a)P(s|s′,a)x(s′,a)=μ₀(s)，总质量1/(1−γ)。正占用处恢复π=x/Σx，零占用状态另定。
3. 确定性OCP对C¹测试q保初末Liouville等式。扩散P-3b44b82c06eccd85需C^{1,2}及Dynkin可积性，Aq=∂_tq+f·∇q+(1/2)tr(ggᵀ∇²q)；终端质量1、运行occupation质量T。跳过程改完整跳生成元，不能漏跳项。
4. 最小成本弱measure LP可能严格低于原轨迹控制最小值；多项式紧域球约束下moment-SOS层级渐近趋该弱LP，不认证任意有限阶或真实控制恢复。
5. 均场来源P-078f65ff156d77ea与P-22e8303b48a27276用同一时刻人口边缘的二次交互；凸性检矩阵核对所有有符号测度PSD及交互权符号。冻结梯度后的LMO线性分离不证明原目标凸；exact LMO、可测经典最优选择、状态约束内点和smoothness/曲率条件才调用FW速率。详见[分模型来源与定位](references/frozen-source-locators.md)。
6. 变分源P-b8164dc1ce4dcf5d的无gap限连通分片C¹域、矩/紧支撑及体内/边界联合凸性或OC1–OC3；一般仅M_aff≤M_r≤M。reachability体积、危险区期望停留时间和最小成本的目标/界方向分开；有限矩匹配不能恢复真实控制轨迹，详见补充审计。
7. 局部分片来源按适配分区保反对称signed界面通量并验证聚合抵消；扩散空间界面要连续，确定性单向放宽不适用扩散。作者有限阶更紧比较另需闭片的严格更丰富多项式表示，非爆炸/停时及通量可定义性另核。
8. P-610d8ca4e35b8a23的peak最大化弱measure/SOC给上界：Cantelli r=√(1/ε−1)，VP另需单峰且0<ε≤1/6；方差c²+a²≤b保统计矩约束。Problem3.1/5.1 的t*明确为确定性标量；(29)未约束共同terminal time，checked-local双分支cutoff ODE在ε=3/4给固定时刻ES最优2/3、measure LP最优1，故Theorem5.5的literal固定时刻no-gap为假；Theorem5.4的upperbound仍可用。
9. 其ES用0<ε≤1及概率ψ=p#μτ，εν+νhat=ψ、ν质量1编码可拆分分位原子的标准CVaR；打印整原子事件公式在V1/0各概率1/2、ε1/4给2而标准式/LP给1，不沿用错误等价。
10. 作者tail-SOC强对偶附录用max p控制|E p|有符号缺口；可由a²≤b≤Π₂将trace粗界修为6(1+2Π₂+Π₂²)，仍保抽象对偶/停止表示/层级外引依赖。
11. 闭环控制、稀疏分区实验、数值重构与论文理论分开；无限离散跳态的半代数外包、所有共同输入和中心立方体流接口另证。
12. P-0297391ab7ac5c7b finite-mode切换：μ_j为u_j加权投影；逆向以Σμ_j及δ_ej合并，weak LP等价不等于原switching无gap。用v=t检Σ质量T，time-only tests检Σπ_tμ_j=dt；joint reconstruction须保共同time marginal并真实积分验证。Archimedean polynomial hierarchy仅渐近给weak LP lower bounds；source p24乘性dynamic与linear chattering计算矛盾，Cor2 support无权积分暂不作认证。

## 证据与失败处理

occupation质量与界方向取决于模型；弱LP未必无gap或可恢复轨迹，有限阶数值亦不认证。随机局部分片需合法通量/界面条件；peak来源的确定性共同终止时刻问题与未限制共同终止时刻的停止测度LP不等价、ES分位原子错误及tail-SOC符号界缺口保留，不能当一般无条件occupation定理。 人口交互一般为二次测度泛函，核PSD或权重符号缺失不认证凸性；精确oracle速率不覆盖Adam/网格近似，无gap变分定理不泛化到任意动态控制。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
