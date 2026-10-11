# 适用案例

取非零频率 ξ≠0 的单频乘子模式 f(x)=e^{iξx}。因为 ∂ₜPₜf=−|ξ|e^{-t|ξ|}f，故逐模式积分 ∫₀∞t|∂ₜPₜf|²dt=|f|²/4。plane wave 仅用于乘子核验，并非 L²(R) 输入；一般 L²(R) 版本由 Plancherel 成立，非负自伴 L 的谱公式为 (1/4)||(I−Π_{ker L})f||₂²。

已执行检查：`audit_current/f/cases.json` 中的 `F04-POS`；逐项输出见 `audit_current/f/example_execution_20261007.txt`。这些有限计算只核验该例，不证明一般定理。
