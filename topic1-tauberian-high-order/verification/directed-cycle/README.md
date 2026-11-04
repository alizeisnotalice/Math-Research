# 第16轮：必要性尺度的最小符号核验

本机 Mathematica/Wolfram Language 15.0.0 ARM 实际运行 `scale.wl`（单内核PID62439），7项精确检查全部True、FAILURES0：临界指数两种写法、平衡恒等式、q≥3时β>6的全实域违反集为空、q2回归，以及坏量下界/可删除预算的两个幂次比。实际输出为 `scale.log`。

复现：`wolframscript -file scale.wl > scale.log`。无私有绝对路径，无自动时间限。

它只检验代数。精确外包邻孔枚举、可行性迫使e=o(u)、实际潜在坏块下界和有限凸截断的排除由 `proofs/016-directed-cycle-exclusion.md` 的独立短证明承担。本轮不归档不再用于主攻目标的更晚项系数展开、重复逐例结果或临时缓存；这些历史探索保留在仓库外。
