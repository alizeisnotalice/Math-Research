# 独立复核记录（append-only，2026-10-06）

## 引用忠实性深查 + 核心推导完整证明（R20）

### 步骤 3 的 Schur test 推导 = A 级（完整证明）
Σ_j|a_ij||x_j| = Σ_j √(|a_ij|q_j)·√(|a_ij|/q_j)|x_j|
  ≤ √(Σ_j|a_ij|q_j)·√(Σ_j|a_ij||x_j|²/q_j)   [Cauchy–Schwarz]
  ≤ √(R p_i)·√(Σ_j|a_ij||x_j|²/q_j)
⟹ Σ_i|(Ax)_i|² ≤ R Σ_i p_i Σ_j|a_ij||x_j|²/q_j ≤ R·C·Σ_j|x_j|²   [列界 Σ_i|a_ij|p_i≤Cq_j]
⟹ ‖Ax‖² ≤ RC‖x‖²，即 ‖A‖≤√(RC)。与源卡 C_A≤√(ζ(A)κ(A)) 一致。□

### 引用忠实性（对照 4/4 全文已读源卡）
- 步骤 5：行列和界 √(RC) 与 PSD 逐项乘子 ‖H∘A‖≤sup h_ii‖A‖ 的区分 —— 与源卡两条
  theorem_cards 的明确区分一致 ✓
- 步骤 6：D_L、D_M 因子分解与 √(D_LD_M) 界 —— Schur 乘子标准形式 ✓
- 步骤 7：双谱积分需正交谱测度 E,F 与 h(μ,λ)=∫m·l，两因子 L² esssup —— Birman–Solomyak
  标准形式 ✓；skill 如实声明 Peller/Birman–Solomyak 原始证明未另核
- 步骤 8 边界（离散 HL 极大、维数无关）与 EG-S 系列卡一致 ✓

## 档位提升依据（R20）
待证据 → **部分可用**（核心推导完整证明 + 引用忠实性通过）。
仍待证据：双谱积分/Peller 外引证明；具体 A(λ) 的参数实例化。


---

## 更正条目（R41）：卡↔PDF 保真性——F05 的 Schur 源通过

- **加权 Schur test**（Dym–Katsnelson 综述 PDF p17）：ζ_r/κ_r 定义与
  **C_A ≤ √(ζ_r(A)κ_r(A))** (3.15) 与源卡逐字一致；(3.16) 插值形式与 R20 的
  独立证明结构对应（R=ζ, C=κ，两路径互证）；Toeplitz/Hankel/Hilbert 例确认
  加权版必要性。
- 卡↔PDF 累计 12/13 卡逐字通过、1 处转录偏差（E03，已修正）。
