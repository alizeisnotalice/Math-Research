# 独立阅读卡：Grey 2010（I01 adjacent/background，全文重读 13/13 页）

**paper_id**: P-073066f98d171070 ｜ arXiv:1006.4465v1 [math.PR] 23 Jun 2010
**题名**: The associated random walk and martingales in random walks with stationary increments
**作者**: D. R. Grey（Sheffield）｜ 13 页 ｜ 本轮独立阅读范围：PDF p1–13（2026-10-07）；历史 R28 声明单独保留，不作为本轮阅读依据
**SHA256（sources.json 登记）**: 073066f98d171070514958d8eff6d59ac50ee69fabf330c051d34099073a2ca5

## 定义与假设

- 背景：iid 增量、正漂移随机游走 S_n；Cramér 根 θ>0 满足 F̂(θ)=E(e^{−θX₁})=1；
  伴随游走 dF*=e^{−θx}dF（对偶：F̂*(−θ)=1）；Wald 鞅 V_n=e^{−θS_n}。
- **Assumption 1**：∃θ>0 使 q := lim_{n→∞} E(e^{−θS_n}) 存在且 ∈(0,∞)。
  θ 唯一（S_n→∞ a.s. 时：φ<θ ⟹ E(e^{−φS_n})→0；φ>θ ⟹ →∞，论文给了完整论证）。
- **Assumption 2**（双端版本，用于伴随游走）：对 k≥1、B∈F_{−k,k}，
  q(B) := lim_{m,n→∞} E(e^{−θS_{−m,n}}; B) 存在。
- **Assumption 2\***（单端版本，用于鞅）：对 B∈F_k，r(B) := lim_{n→∞} E(e^{−θS_n}; B) 存在。
- 增量过程：平稳遍历（stationary ergodic）；扩为双端平稳序列（Breiman Prop 6.5）。

## 主构造（证明结构，每步已核）

1. P*_{m,n}(B) = E(e^{−θS_{−m,n}}; B)/E(e^{−θS_{−m,n}}) → q(B)/q；由逐点收敛定理
   （Grey 2001，Vitali–Hahn–Saks 特例）得 P* := q(B)/q 为概率测度；
   k 的一致性 + Carathéodory 扩张到 F_{−∞,∞}。**伴随游走由此定义**；
   平稳性保持（双端极限所需）；遍历性与对偶性一般**未解决**（§3）。
2. 鞅：r 是 F_k 上总质量 q 的测度，P≪r 的逆——r≪P（e^{−θS_n}dP 在零测集上为 0）；
   V_k := dr/dP|_{F_k}，链式 V_k = E(V_{k+1}|F_k) ⟹ **鞅**。
   注：V_k = lim E(e^{−θS_n}|F_k) a.s. 在一般情形"难以作为定义"（论文原话）。
   iid 特例 V_k = e^{−θS_k} ✓ 复原 Wald。

## 三个应用（关键公式全部独立验算）

### 2.1 平稳 Markov 链增量（可数态空间、不可约非周期）
- 正则条件：Q=(p_ij e^{−θj}) 有 PF 特征值 1，左右特征向量 v^T、c，v^Tc=1，Q^n→cv^T。
- Assumption 1: E(e^{−θS_n}) = π^TQ^n1 → π^Tc·v^T1 ✓（矩阵极限）。
- 伴随链：**p*_ij = p_ij·e^{−θj}c_j/c_i，π*_i = c_i v_i**（三因子分解+µ^Tc=π^Tc 恒等式，
  论文 p7 逐步推导，已核）。
- 鞅：**V_k = c_{X_k} e^{−θS_k}**（取 v^T1=1；同 Lu 1991）。

### 2.2 平稳 Gaussian 增量（μ>0, σ²>0）
- 正则条件：Σ_{r≥1} r|ρ_r| < ∞；R := Σρ_r，S := Σrρ_r（有限）。
- 必要 1+2R ≥ 0（var S_n = σ²(n(1+2R)−2S)+o(1)）；**极端情形 1+2R=0 须排除**
  （例 X_n := μ+Z_n−Z_{n−1}, Z_n iid N(0,σ²/2)——存在非平凡核）。**端点条件，重要**。
- **θ = 2μ/(σ²(1+2R))**（由 −μθ+½σ²(1+2R)θ²=0；已独立验算）。
- **q = exp(−4μ²S/(σ²(1+2R)²))**（指数化简为 −σ²Sθ²；已独立验算到 1e−14）。
- 伴随游走：协方差结构与原过程相同、漂移 −μ（向下）。
- 鞅：V_k = exp(γ^TX+δ)（指数仿射形式）。

### 2.3 排队应用（G/GI/1）
- Lindley (1952)/Kingman (1964)：W_n 平衡分布 = 无约束随机游走全时最小值的负值；
  增量 X_n := T_{−n}−U_{−n}。
- M/M/1：E(e^{−θX}) = λμ/{(λ+θ)(μ−θ)}，θ = μ−λ（**已独立验算精确=1**）。
- 预约系统（误差 iid，Laplace ψ）：θ 满足 **φ(−θ)exp(−λ^{−1}θ) = 1**（与 ε_n 分布无关）；
  与随机到达比较：(μ/λ)e^{1−μ/λ}<1 ⟹ 预约系统 θ 更大 ⟹ 等待尾更薄（**已验算**，
  log u < u−1）。

## §3 对偶与渐近独立性（开放问题）

对偶一般**不自动成立**：需 S_{−r,−m−1}、S_{n+1,s} 与 B 的近似独立性——
一种 **mixing 条件，强于遍历性**（Bradley 2005）。两个特例（Markov、Gaussian）中对偶成立。
**论文诚实标注**：一般平稳遍历情形伴随过程的遍历性、对偶性均开放。

## 可迁移方法（对 I01/I02 接口）

1. **指数换测的非 iid 推广**：θ 由渐近斜率确定（Assumption 1 的唯一性论证可复用）；
2. **RN 导数构造鞅**（V_k = dr/dP|_{F_k}）——不依赖增量独立性，只依赖 r(B) 极限存在；
3. **端点警示**：Gaussian 情形 1+2R=0 必须排除——迁移到中心立方体接口时，
   类似的"方差增长退化"端点需单独检查；
4. **遍历 ≠ 充分**：对偶需 mixing——凡引用"伴随/换测对偶"的迁移推断必须额外假设混合条件。

## 未核验环节（如实登记）

- Grey (2001) 的逐点收敛定理（VHS 特例）**引用未重证**（原文给"straightforward proof"，
  主体依赖 Dunford–Schwartz III.7.2–4）；
- Markov 情形 p*_ij 构成转移概率、π* 为平衡分布——论文称 routine，未逐步展开（结构可信，
  行/列和归一未独立逐项验算）；
- Gaussian 情形 Assumption 2* 的"similar calculations"未在论文中完整给出；
- 遍历性/对偶性开放问题（论文自declare）。

## 与 I01 skill 的关系（本轮修正）

本文与有序滤过中的换测鞅有概念邻近性，但不能支撑 I01 的条件 Bernstein 矩输入或尾界。Assumption 1/2/2* 是构造伴随随机游走和 RN 导数鞅的渐近极限条件；文中没有要求逐个增量条件中心化，也没有 `E(|X_k|^p|F_{k-1})≤(p!/2)v_k c^(p-2)` 假设。因此 E0847 当前归类为 **adjacent/background**，不是 I01 的 direct theorem source。

本轮独立逐页重读全文 lines 1–571、PDF p.1–13，视觉核验 pp.3、5、8–9 的假设和公式。Assumption 1 在 PDF p.3；双端 Assumption 2 在 p.3；单端 Assumption 2* 在 p.5；Gaussian 例横跨 pp.8–9。I01 当前条件 Bernstein 推导的证明见同目录 `method.md`，该证明仅以 I01 明列的条件矩假设为前提，不借用 Grey。

Grey (2001) VHS 收敛定理仅在本文参考文献中被引用，本轮未读未证；不能把本文 full-read 状态传递给该参考文献。
