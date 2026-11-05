# 非中心强极大算子端点常数

研究对象为 Lebesgue 测度下轴平行非中心强矩形极大算子，统一使用
\(\Phi_n(t)=t(1+\log^+t)^{n-1}\)、自然对数及阈值 1。
本项目研究最优比值 \(C_n=\sup_f |\{M_Sf>1\}|/\int\Phi_n(f)\)。

## 当前状态（2026-11-06，第 3 轮）

- 已复核：指数率 \(C_n^{1/n}\to B\)，以及 \(\limsup C_n/B^n\le\sqrt B\) 的证明路线；B 满足 \(\log B=1+B^{-1}\)。详细推导见第 3 轮。
- 已实现：等参数单层见证的 Gamma–二项混合、数值优化和 Wolfram 高精度交叉复算。
- 待补证明：附件所称整个双曲示性函数族 H_n 的全局三项展开，特别是全空间 halo 与有限 L 的统一比较。
- 已完成：有限多层轮廓的严格极限归约与固定网格凸振幅优化；n=4…256 的候选均优于同维等参数单层数值对照。
- **第3轮严格下界**：n≥1024 时，log C_n ≥ n mu − sigma sqrt(n(log n−2loglog n)) − log(2n)。新增正的 sigma loglog n sqrt(n/log n) 项；首个负项系数仍为 sigma。完整证明见 rounds/r003/report.md。
- 主攻：有限阶梯定向认证、完整凸轮廓的渐近障碍及一般上界近阈值核心。
- 一般输入的负 \(\sqrt{n\log n}\) 次项仍未证明。

## 复现

MATLAB R2026a（Optimization Toolbox）：
```matlab
addpath('code');
sm_checks;
rows=sm_baseline;
```
Wolfram 15：`wolframscript -file code/crosscheck.wl`。
固定种子不适用于本轮：所有算法均为确定性。
结果摘要保存在 `rounds/r001/`；数值精度检查不等于区间认证。

## 研究约定

每轮分别记录严格定理、受限模型结论、数值候选与证明缺口。
每轮结束仅归档必要报告、代码和摘要；一轮一个提交。
下一轮接续点见最新轮次报告，不以低维改进宣称高维次项突破。

第2轮复现：MATLAB `rows=sm_multilayer;`，Wolfram `wolframscript -file code/multicheck.wl`。

第3轮复现：MATLAB `rows=sm_exponential;`；Wolfram `code/softproof_checks.wl`、`code/exponential_check.wl`；Python `code/soft_constants.py`。
