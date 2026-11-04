# R19 局部窗口的标量核验

`local_window.wl` 在 Mathematica 15.0.0 中实际执行，输出 `local_window.log`。这只是对局部铺满归纳中若干标量步骤的精确核验，不证明几何归纳或一般输入到 LT 的桥。

在 R≥1、r>0、c≥0、c≤R+r/2 且 c>max(1,r) 下，`Reduce` 检查 r≥2R 和 c≥2R 的违例集合，均返回False。若原单位半盒已覆盖候选中心，使它到某单位中心的距离≤1/2，则此距离不可能同时满足 >max(1,r)；该违例集合也为False。

覆盖归纳得到半径 ρ_m=1/2+m/4。只需 m≥8R+8，已有 ρ_m−(2R+1)≥3/2；因此 m=ceil(8R+8) 符合要求。代码将m按任意实数处理，比所需整数条件更宽。

严格分离不可省掉：若只取非严格条件，R=1、r=c=2满足接触假设且r=2R、c=2R。若物理模型采用其他开闭边界，需另行处理这些等号情形。代码无自动时限，复现命令为 `wolframscript -file local_window.wl`。

## Hölder 指数

`holder_exponents.wl`/`.log` 实际核验一般实q≥2：令ζ=q/(2q−1)，有

  q+6ζ−6 = (2q−3)(q−2)/(2q−1)≥0，  1−ζ≥1/3。

当q≥3时前一差值至少3/5，q=3取等；三个违反集均由Reduce返回False。因此 0<v,μ≤1 时，Hölder产生的 v^(q+6ζ)μ^(1−ζ) 不超过 v⁶μ^(1/3)。此外在线性预算 Σε=O(v⁵S)、F=O(v⁶S)、S_b=O(v⁵S) 下，v²Σε、v²F、v⁴S_b分别为 O(v⁷S)、O(v⁸S)、O(v⁹S)。这些指数核验不替代相应几何计费条件。

本包的全部标量已实跑；完整一般输入桥及Hölder的几何计费条件现由项目 proofs/019-input-tiling.md 给出，独立审核见 reviews/019-input-tiling.md。保留日志原始PID/版本。两个脚本不依赖私有路径，均可用wolframscript -file 文件名实际复现。
