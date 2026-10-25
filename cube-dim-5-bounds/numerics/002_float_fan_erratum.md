# 002 勘误：旧浮点 Fourier 锚点输出无效

`body_fan_check.py` 的旧 `heat_density` 将 `matrices[c,m,j]` 的物理坐标轴与 ray 轴混淆，使用了错误的 `einsum` 下标。故 `fan_results.json` 中的浮点 anchor 值和由同一实现得到的齐次回归全部无效；它们不用于任何结果比较。

7 锥精确覆盖、方向符号检查、精确 basis 审计和 `body_fan_arb_anchor.py` 的逐 ray 点积均不使用该函数，因而不受影响。修复后的浮点实现必须在新的输出中与 Arb 点值交叉核验。
