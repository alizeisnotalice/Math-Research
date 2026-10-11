# L² 谱参数与 Schur 估计：方法与验收边界

1. 在有限支撑稠密子空间定义A(λ)，写明两Hilbert空间、底测度及λ的含义；矩阵行列Schur、逐项Schur乘子和谱分解是不同命题。

2. 取正权p_i,q_j，验证Σ_j|a_ij(λ)|q_j≤R p_i及Σ_i|a_ij(λ)|p_i≤C q_j，两侧同时成立且R,C对所需λ范围统一。

3. 已检查Cauchy–Schwarz构件：把每行Σ_j|a_ij||x_j|拆为√(|a_ij|q_j)与√(|a_ij|/q_j)|x_j|；平方后求i、交换非负和并用列界，得||Ax||²≤RC||x||²。

4. 积分核同理先检可测性、正权、Tonelli与有限支撑估计，再密度延拓；权可依赖λ但界及定义域必须可控。

5. P-3da1b59567d35c36的Luna全文卡现覆盖99页；第3节行列和界为√(RC)，第4节正PSD逐项乘子另为||H∘A||≤sup h_ii||A||，不能相互替换。全文覆盖不等于全部外引定理独立认证。

6. **一般 Schur 因子界的独立有限证明（F05-LMSTAR-005）**：令有限行指标 `I`、列指标 `J`，公共因子指标为 `r`，并设 `h_ij=Σ_r l_ir conjugate(m_jr)`、`D_L=max_iΣ_r|l_ir|²<∞`、`D_M=max_jΣ_r|m_jr|²<∞`。对 `x∈ℓ²(J)`、`y∈ℓ²(I)`，有限支撑求和给
   `⟨(H∘A)x,y⟩ = Σ_r ⟨ A(x·overline{m_r}), y·overline{l_r}⟩`，
   其中矩阵内积约定为 `⟨u,v⟩=Σ_i u_i overline{v_i}`（对第一变量线性），`(m_r)_j=m_jr`、`(l_r)_i=l_ir`。对 `r` 求和绝对收敛，因为 `I,J` 有限且每个因子行属于 `ℓ²(r)`。因而其绝对值不超过
   `||A|| [Σ_r||x·overline{m_r}||²]^(1/2)[Σ_r||y·overline{l_r}||²]^(1/2) ≤ ||A||√(D_M D_L)||x||||y||`，最后一步按 Tonelli/单调收敛交换非负级数并应用两项行 `ℓ²` 上界。取单位向量上确界即得 `||H∘A||_(2→2)≤√(D_LD_M)||A||_(2→2)`。该独立证明限有限矩阵/有限支撑求和及可行矩阵乘法；没有把来源的更一般积分因子、非方阵或无限指标情形冒充为本原子证明。
   Dym–Katsnelson survey 的 §4 PDF pp.22–23 equations (4.6)–(4.10) 报告上述正向估计，并允许非方阵及适当的 `L²(X,dx)` 因子。其 §4 PDF pp.23–24 equation (4.11) 另转述 Bennett 逆向定理：若有限或无限 `H` 满足对所有同尺寸 `A` 的 Schur 乘子范数界 `D`，其中 `D_L=sup_p Σ_r|l_pr|²`、`D_M=sup_q Σ_r|m_qr|²`，则每个 `ε>0` 有 `H=LM*` 且 `√(D_LD_M)<D+ε`。此逆向是 **F05-BENNETT-006 citation-only**，原文指向 Bennett Theorem 6.4，未独立阅读/证明；不与正向独立证明合并。Peller 与 Birman–Solomyak 的其他转述仍保持原有外部证明缺口。

7. survey 的 Stieltjes 双算子积分模型要求 Λ、M 为可测空间，E(dλ)、F(dμ) 是同一可分 Hilbert 空间 H 上的完备正交谱测度（E(Λ)=F(M)=I），且双算子积分按文中所述意义存在。核需可测并分解为 h(μ,λ)=∫_X m(μ,x)l(λ,x)dx，其中 X 带非负 σ 有限测度 dx，C_m=esssup_μ∫|m(μ,x)|²dx、C_l=esssup_λ∫|l(λ,x)|²dx 均有限。survey 转述分别给出 ‖T_h‖_{B(H)→B(H)}≤√(C_mC_l) 和 ‖T_h‖_{S₁(H)→S₁(H)}≤√(C_mC_l)；这是 Birman–Solomyak 结果的作者转述，未独立证明，Peller 反向因子化也未独立证明。不能仅凭有界 h、矩阵行列界或 PSD Schur 界调用该 DOI 结论。

8. 参数端点、极点、零权及正交谱块分别检查；常数1的有限Weyl谱/奇异值幂和比较也不提供具体A(λ)一致界。范数统一界不证明核逐点收敛、Hilbert矩阵最优权或任何私人谱接口。

### 三种不同的 ℓ² 估计

- **行列和 Schur 检验**作用于矩阵/核的绝对行和、列和，需两边同时有界，结论为算子范数 ≤√(RC)。
- **正半定 Schur 乘子**是有限矩阵逐项乘积。若 H⪰0，取 Gram 向量 u_i 使 h_ij=⟨u_j,u_i⟩，并令 V e_i=e_i⊗u_i，则 H∘A=V*(A⊗I)V，‖V‖²=max_i h_ii；所以 ‖H∘A‖≤(max_i h_ii)‖A‖。这是 Gram 因子化的独立推导，不能替代前一项的行列和条件。
- **Stieltjes 双算子积分**是谱测度模型。按来源 survey 的 pp.25–26 (4.12)–(4.20)，Λ、M 是可测空间，E(dλ),F(dμ) 是同一可分 Hilbert 空间 H 上满足 E(Λ)=F(M)=I 的完备正交谱测度，且 DOI 按文中所述意义存在；可测核需有 h(μ,λ)=∫_X m(μ,x)l(λ,x)dx 因子化，X 带非负 σ 有限测度 dx，并满足 C_m=esssup_μ∫|m(μ,x)|²dx<∞、C_l=esssup_λ∫|l(λ,x)|²dx<∞。survey 分别报告 B(H)→B(H) 与 S₁(H)→S₁(H) 两个算子范数界，常数均为 √(C_mC_l)；这是 Birman–Solomyak 结果的作者转述，不是本组独立证明，Peller 反向因子化也仅按转述记录。具体定位和缺口见 [P-3da1b59567d35c36 来源卡](papers/EG-P-3da1b59567d35c36.json)。

## 不可省略的限制

两侧正权行列界缺一不能调用Schur；survey全文卡已覆盖99页，但不含私有参数λ一致定理，且Bennett/Peller/Birman–Solomyak为外引依赖。独立Cauchy–Schwarz只认证该弱构件，不能代替真实核、谱重叠、端点或private接口。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。

局部证据补充见 [ROOT记录](audit-root.md)；其审读者与范围以记录为准。
