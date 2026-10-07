# 独立复核记录（append-only，2026-10-06）

本文件由外部审核批次写入，不修改本 skill 任何既有内容。

## 重要更正（R14，方法错误第 19 次）

**撤回早前评价**：R12 曾称 I01 的 sources.json「信息量严格多于其余 63 项的单卡格式」。
该判断只看了 schema 字段丰富度，**未检查 read_scope 内容**。实际核对：

- 40 条 attachment_entries 的 read_scope **全部**为「Title and first-page abstract/
  introduction lines7–38 only; initial relevance screen」——即**全部 40 篇文献只做了
  题名+摘要初筛，无一篇全文深读**。
- 因此 I01 的「待证据」档位**完全正确且必要**：它比其余 63 项更缺全文证据。
- 单卡格式（references/papers/*.json）虽然字段少，但其余 63 项中大量卡记录了
  全文阅读范围（见 pending_source_read_scan.json 的 full 列统计）。

## 验证状态

- I01 自身步骤含条件 Bernstein（中心化）与有序增量结构——案例层未执行（40 卡均初筛，
  无可用全文模型）。
- 待证据：全部 40 篇文献的全文获取与深读；H05 依赖条件与有序增量 Bernstein 的
  具体化（源卡 review 中已把相关文献标为 adjacent）。


---

## 更正条目（R27）：缺口性质修正——PDF 全部本地可得

核查 evidence 池：**40/40 篇文献的 PDF 全部在本地**（evidence/papers/<paper_id>/paper.pdf），
零缺失。因此本项缺口是「全文深读未执行」而非「来源不可得」——补读在原则上是可行的，
只受工作量约束。

**身份抽样核对（3 篇，人工比对首页）**：
- E0824 P-fe2dd431ab2300be：标题 "Processes with block-associated increments" (Jakubowski) 与文件名吻合 ✓
- E0825 P-4d166c976407386e："Inference for parameters identified by conditional moment restrictions using a generalized Bierens..." 吻合 ✓
- E0826 P-81afa2f95543677b："Weighted Gaussian Approximations for Increments of the Uniform Empirical and Quantile Processes" 吻合 ✓
（自动前缀匹配报 False 是因文件名带序号/作者/年份前缀——方法假象，非身份问题。）

**direct 论文定位**：40 篇中仅 1 篇判 direct——E0847 P-073066f98d171070
（Grey 2010, The associated random walk and martingales in random walks with non-iid increments；
摘要将 Wald 鞅/associated random walk 从 iid 扩至 stationary ergodic 增量，p1 定义
E exp(−θX)=1 与指数换测——与 I01/I02 指数鞅结构直接相关）。其余 18 adjacent / 21 irrelevant。

**状态：保持待证据**。全文深读（至少 direct 的 Grey 2010 + 18 篇 adjacent 中的
条件 Bernstein 相关者）仍未执行；这是本项提升为部分可用的前置条件。


---

## 全文深读记录（R28）：direct 论文 Grey 2010 完成 13/13 页

**证据卡**：`references/evidence-card-Grey2010-P-073066f98d171070.md`（完整：
定义、Assumption 1/2/2*、主构造逐步、三应用、开放问题、可迁移方法、未核验环节）。

**独立验算（4 处，全部通过）**：
1. M/M/1 θ=μ−λ 精确（(λ+θ)=μ、(μ−θ)=λ 抵消，1.000000000000）；
2. Gaussian θ=2μ/(σ²(1+2R)) 二次方程正根；
3. q=exp(−4μ²S/(σ²(1+2R)²)) 恒等式（指数=−σ²Sθ²，1e−14 精确）；
4. 预约系统 (μ/λ)e^{1−μ/λ}<1 ⟺ log u<u−1，θ 更大 ⟹ 尾更薄——方向正确。

**关键端点捕获**：Gaussian 情形 1+2R=0 须排除（存在反例 X_n=μ+Z_n−Z_{n−1}）——
迁移时的"方差增长退化"端点警示；对偶需 mixing（强于遍历），一般情形开放。

**未核验环节**（如实）：Grey 2001 收敛定理引用未重证；Markov p*_ij 归一化 routine 未逐步；
Gaussian 2* 验证论文写"similar calculations"未全文给出。

## 档位提升依据（R28）

待证据 → **部分可用**。依据：
1. 唯一 direct 文献全文深读完成 + 证据卡（含全部假设/端点/未核验环节）；
2. 鞅构造（RN 导数链 V_k=E(V_{k+1}|F_k)）经阅读核验为正确（3 行证明完整）；
3. 40/40 PDF 本地可得，补读可行；3 篇抽样身份人工确认；
4. 换测/鞅接口的非 iid 推广路径（本文核心）已进入 skill 可引用的证据层。

仍待证据（如实保留）：18 篇 adjacent 文献未深读（条件 Bernstein 的具体化依赖其中
若干篇）；Grey 2001 收敛定理未重证；专题迁移主张保持守卫。
