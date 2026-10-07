---
name: math-f05-l2-spectral-parameter-schur-test
description: "用于L² 谱参数与 Schur 估计的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **交接准备状态（结构性）**：本副本只整理技能结构与引用可移植性；未进行全文深读或数学验证。原评级、可用状态、假设、证明限制和未决条件均沿用原件。结构检查不构成科学验证。
> 外部文献证据由第二证据包提供，是独立数据输入，不随单个 Skill 安装。合并交接包默认将 `EVIDENCE_ROOT` 设为包内 `evidence/`；单独安装时由调用方传入 `EVIDENCE_ROOT`，按本地[证据索引](references/handoff-evidence-index.csv)中的 `portable_path` 查找。无需全局安装，也不要把外部 PDF 当作 Skill 内文件。


# F05 · L² 谱参数与 Schur 估计

当前证据状态：**待证据验收：可执行方法工作流，不是已证明的课题结论**。全文转换、论文阅读和证明核验是三个独立状态。（此行声明已被下方【修订】标注替代：该 skill 经外部审核提升为「部分可用」，详见对应修订行与 references/verification-log-20261006.md）
  【修订 2026-10-06（外部审核 R20）】步骤 3 的 Schur test 推导经独立完整重证（Cauchy–Schwarz 拆分+交换非负和+列界 ⟹ ‖A‖≤√(RC)），行列 Schur 与逐项乘子的区分、双谱积分条件均与源卡一致；证据状态提升为：部分可用：基础推导已检查，专题证据仍待完整验收。双谱积分/Peller 外引证明未另核；具体 A(λ) 实例化待做。

## 输入与产出

Hilbert 基或谱分解、参数矩阵/核 A(λ)、输入输出测度、所声称的行和与列和估计。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 在有限支撑稠密子空间定义A(λ)，写明两Hilbert空间、底测度及λ的含义；矩阵行列Schur、逐项Schur乘子和谱分解是不同命题。
2. 取正权p_i,q_j，验证Σ_j|a_ij(λ)|q_j≤R p_i及Σ_i|a_ij(λ)|p_i≤C q_j，两侧同时成立且R,C对所需λ范围统一。
3. 已检查Cauchy–Schwarz构件：把每行Σ_j|a_ij||x_j|拆为√(|a_ij|q_j)与√(|a_ij|/q_j)|x_j|；平方后求i、交换非负和并用列界，得||Ax||²≤RC||x||²。
4. 积分核同理先检可测性、正权、Tonelli与有限支撑估计，再密度延拓；权可依赖λ但界及定义域必须可控。
5. P-3da1b59567d35c36的Luna全文卡现覆盖99页；第3节行列和界为√(RC)，第4节正PSD逐项乘子另为||H∘A||≤sup h_ii||A||，不能相互替换。全文覆盖不等于全部外引定理独立认证。
6. 若实际问题是逐项乘子且H=LM*，给D_L=sup_iΣ|l_ir|²与D_M=sup_jΣ|m_jr|²，才能用√(D_LD_M)界；Bennett逆向、Peller与Birman–Solomyak双谱积分是作者转述，外部原始证明未另核。
7. survey双谱积分入口另需正交谱测度E,F和核h(mu,lambda)=∫m(mu,x)l(lambda,x)dx，两因子L2 norm的esssup受控；仅矩阵行列界或有界h都不能替代该因子证书。Birman–Solomyak/Peller等原始外引证明仍未另读认证。
8. 参数端点、极点、零权及正交谱块分别检查；常数1的有限Weyl谱/奇异值幂和比较也不提供具体A(λ)一致界。范数统一界不证明核逐点收敛、Hilbert矩阵最优权或任何私人谱接口。

## 证据与失败处理

两侧正权行列界缺一不能调用Schur；survey全文卡已覆盖99页，但不含私有参数λ一致定理，且Bennett/Peller/Birman–Solomyak为外引依赖。独立Cauchy–Schwarz只认证该弱构件，不能代替真实核、谱重叠、端点或private接口。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
