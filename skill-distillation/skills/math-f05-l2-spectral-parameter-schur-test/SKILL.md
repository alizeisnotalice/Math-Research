---
name: math-f05-l2-spectral-parameter-schur-test
description: "用于L² 谱参数与 Schur 估计的数学研究：核对输入与假设，组织方法和证据，说明中心立方体极大算子的迁移条件；不凭名称补造内部接口。"
---

> **生产修订（2026-10-07）**：当前断言与证据等级记录位于 `audit_current/f/production_audit_20261007.md`（包内路径）；独立安装时按[包级审计访问说明](references/portable-audit-access.md)解析。此前交接状态和追加修订历史位于 `audit_current/f/history/pre-audit-20261007/`。示例运行只检验具体输入，不作为一般定理证明。

来源索引路径及便携解析方式见[来源访问说明](references/source-access.md)；解析器会校验单项来源 SHA，不能替代阅读与数学核验。


# F05 · L² 谱参数与 Schur 估计

当前断言与等级以原子账本为准；本条 Skill 不因某一基础引理通过就整体标记为已验收。

## 输入与产出

Hilbert 基或谱分解、参数矩阵/核 A(λ)、输入输出测度、所声称的行和与列和估计。

产出应包含适用性判断、使用的假设、数学步骤、证书或界的方向、常数依赖，以及尚未解决的缺口。

## 执行步骤

1. 在有限支撑稠密子空间定义A(λ)，写明两Hilbert空间、底测度及λ的含义；矩阵行列Schur、逐项Schur乘子和谱分解是不同命题。
2. 取正权p_i,q_j，验证Σ_j|a_ij(λ)|q_j≤R p_i及Σ_i|a_ij(λ)|p_i≤C q_j，两侧同时成立且R,C对所需λ范围统一。
3. 已检查Cauchy–Schwarz构件：把每行Σ_j|a_ij||x_j|拆为√(|a_ij|q_j)与√(|a_ij|/q_j)|x_j|；平方后求i、交换非负和并用列界，得||Ax||²≤RC||x||²。
4. 积分核同理先检可测性、正权、Tonelli与有限支撑估计，再密度延拓；权可依赖λ但界及定义域必须可控。
5. P-3da1b59567d35c36第3节的两侧行列和给√(RC)界；第4节的正PSD逐项乘子是另一命题。有限正半定 H 可作 Gram 分解 H=LL*；由 (H∘A)_{ij}=H_{ij}A_{ij} 的因子化可独立得 ||H∘A||_{2→2}≤(max_i h_ii)||A||_{2→2}。这不是行列和检验，也不提供参数族的一致界。
6. 对有限同尺寸矩阵，若 `H=LM*`、`h_ij=Σ_r l_ir conjugate(m_jr)`，且 `D_L=sup_iΣ_r|l_ir|²<∞`、`D_M=sup_jΣ_r|m_jr|²<∞`，则 **F05-LMSTAR-005** 的独立有限支撑证明给出 `||H∘A||_(2→2)≤√(D_LD_M)||A||_(2→2)`：对有限支撑测试向量展开乘积，按因子指标 `r` 求和，以 `||A||` 控制每项，再对 `r` 用 Cauchy–Schwarz 和两个行 `ℓ²` 上界。此界不要求 `H⪰0`，也不等于两侧行列和 Schur 检验。来源 §4 (4.6)–(4.10) 的正向版本及可选积分/非方阵扩展另作来源定位；当前原子证明只认证有限矩阵。**F05-BENNETT-006** 单独登记 Bennett 逆向：作者在 (4.11) 转述其定理，有限/无限矩阵上若 Schur 乘子对所有同尺寸 `A` 一致有界为 `D`，其中 `D_L=sup_p Σ_r|l_pr|²`、`D_M=sup_q Σ_r|m_qr|²`，则每个 `ε>0` 存在 `H=LM*` 因子化且 `√(D_LD_M)<D+ε`。它保持 `citation_only`，不得据此认证 Bennett/Peller 的外部证明。
7. survey 的 Stieltjes 双算子积分模型要求 Λ、M 为可测空间，E(dλ)、F(dμ) 是同一可分 Hilbert 空间 H 上的完备正交谱测度（E(Λ)=F(M)=I）；双算子积分按文中所述意义存在。若可测核有因子化 h(μ,λ)=∫_X m(μ,x)l(λ,x)dx，其中 X 带非负 σ 有限测度 dx，且 C_m=esssup_μ∫|m(μ,x)|²dx、C_l=esssup_λ∫|l(λ,x)|²dx 均有限，则 survey 转述的估计分别为 ‖T_h‖_{B(H)→B(H)}≤√(C_mC_l) 与 ‖T_h‖_{S₁(H)→S₁(H)}≤√(C_mC_l)。这是作者转述 Birman–Solomyak 结果，未独立证明；Peller 的反向因子化刻画也只是转述。有界 h、矩阵行列界或 PSD Schur 界均不能替代上述谱测度和因子证书。
8. 参数端点、极点、零权及正交谱块分别检查；常数1的有限Weyl谱/奇异值幂和比较也不提供具体A(λ)一致界。范数统一界不证明核逐点收敛、Hilbert矩阵最优权或任何私人谱接口。

## 证据与失败处理

两侧正权行列界缺一不能调用Schur；survey全文卡已覆盖99页，但不含私有参数λ一致定理，且Bennett/Peller/Birman–Solomyak为外引依赖。独立Cauchy–Schwarz只认证该弱构件，不能代替真实核、谱重叠、端点或private接口。

先读 [来源与阅读状态](references/provenance.md)。只有实际阅读并核对的定理才能作为文献结论引用；候选题录及自动转换不证明其适用。
需要完整数学步骤时读 [方法说明](references/method.md)；涉及当前课题时读 [课题接口](references/cube-interface.md)。
遇到定义缺失、条件不满足或来源未核验，指出具体缺项，保留可证明的弱结论。不得将迁移推断写成原文定理。
依赖表提供方法路线，不表示所有工具之间存在无条件定理蕴含。

用 [适用案例](examples/positive.md) 与 [条件缺失案例](examples/negative.md) 检查适用边界。
