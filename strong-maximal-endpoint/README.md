# 非中心强极大算子端点常数

研究对象为 Lebesgue 测度下轴平行非中心强矩形极大算子，统一使用
\(\Phi_n(t)=t(1+\log^+t)^{n-1}\)、自然对数及阈值 1。
本项目研究最优比值 \(C_n=\sup_f |\{M_Sf>1\}|/\int\Phi_n(f)\)。

## 当前状态（2026-11-06，第 14 轮）

- 已复核：指数率 \(C_n^{1/n}\to B\)，以及 \(\limsup C_n/B^n\le\sqrt B\) 的证明路线；B 满足 \(\log B=1+B^{-1}\)。详细推导见 rounds/r001/report.md。
- 已实现：等参数单层见证的 Gamma–二项混合、数值优化和 Wolfram 高精度交叉复算。
- 待补证明：附件所称整个双曲示性函数族 H_n 的全局三项展开，特别是全空间 halo 与有限 L 的统一比较。
- 已完成：有限多层轮廓的严格极限归约与固定网格凸振幅优化；n=4…256 的候选均优于同维等参数单层数值对照。
- **第3轮严格下界**：n≥1024 时，log C_n ≥ n mu − sigma sqrt(n(log n−2loglog n)) − log(2n)。新增正的 sigma loglog n sqrt(n/log n) 项；首个负项系数仍为 sigma。完整证明见 rounds/r003/report.md。
- 第4轮：n=4…256 的七个多层候选已通过定向算术认证；完整端点和输入hash见 rounds/r004/certificates.json。
- **第5轮机制障碍**：任意非减轮廓、任意坐标比例向量的单一双曲接收见证，其最优 Q_n 满足 (n mu−log Q_n)/sqrt(n log n)→sigma。这是受限机制定理，不能用作 C_n 的一般上界。
- 第6轮：非负加权深层单调轮廓逐项被 Q_n 支配；分块张量示性函数的产品见证逐项被普通单层示性见证支配。均有明确固定维数深层极限限制。
- **第7轮几何障碍**：共同轴向多中心/对角缩放双曲并集满足 V(Kt)≤K V(t)。统一振幅下，仅合并各分量原有双曲见证不能超过单分量深层基准；不控制新增跨分量矩形。
- **第8轮完整机制展开**：log Q_n = n mu − sigma sqrt(n log n) + sigma[loglog n+log(sigma sqrt(2pi))]sqrt(n/log n) + o(sqrt(n/log n))，并作为一般 C_n 的严格下界。
- **第9轮重排障碍**：允许非单调标量轮廓、任意接收层及逐点选择比例后，锚定族的模块弱常数 W_n 仍精确等于 Q_n。此结论不覆盖非锚定真实矩形。
- 第10轮：闭合一维非负缺口及端点四项账本；证明仅振幅上界不能强迫纯Lp稳定性损失，必须保留Markov和模数效率。另给真实强矩形/坐标复合、真实区间/锚定区间的精确严格反例。
- 第11轮：一般输入严格归约到可预测反馈的一维上界模型 mathcalR_n，满足 Q_n≤C_n≤mathcalR_n。Bellman基准已计算至256步，反馈严格优于任意固定两步比例。
- **第12轮一般定理**：log C_n = n mu − sigma sqrt(n log n) + O(sqrt n)。证明一般首个缺口系数为sigma，并推出 C_n/B^n→0。完整一般上界见 rounds/r012/report.md。
- **第13轮一般三项定理**：第8轮公式现对一般C_n及反馈上界mathcalR_n成立；证明截停累计均值损失，再逐项核验FGL条件鞅相对中偏差。详见 rounds/r013/report.md。
- **第14轮隐式根定理**：令barPhi(z_n)=sigma z_n/sqrt n、d_n=sigma sqrt n z_n，则nmu−d_n−O(log n)≤log C_n≤nmu−d_n+O(n^(3/8)(log n)^(1/8))；因而确定任意固定阶逆对数展开。见 rounds/r014/report.md。
- 主攻：分析固定比例优化和可预测反馈在更小尺度的差别。
- 一般前三项现已闭合；乘法渐近等价、后续锐余项和极值结构仍待研究。

## 复现

MATLAB R2026a（Optimization Toolbox）：
```matlab
addpath('code');
sm_checks;
rows=sm_baseline;
```
Wolfram 15：`wolframscript -file code/crosscheck.wl`。
基准卷积、优化和Bellman递推均为确定性；第7、9轮几何诊断使用明确记录的固定种子。
结果摘要保存在 `rounds/r001/`；数值精度检查不等于区间认证。

## 研究约定

每轮分别记录严格定理、受限模型结论、数值候选与证明缺口。
每轮结束仅归档必要报告、代码和摘要；一轮一个提交。
下一轮接续点见最新轮次报告，不以低维改进宣称高维次项突破。

第2轮复现：MATLAB `rows=sm_multilayer;`，Wolfram `wolframscript -file code/multicheck.wl`。

第3轮复现：MATLAB `rows=sm_exponential;`；Wolfram `code/softproof_checks.wl`、`code/exponential_check.wl`；Python `code/soft_constants.py`。

第4轮复现：`python3 code/certify_staircase.py`；MATLAB `rows=sm_dense;`。

第5轮：完整首系数障碍及有界核心归约见 rounds/r005/report.md；单盒数值核对 `rows=sm_box;`。

第6轮：加权与分块精确支配见 rounds/r006/report.md；复现 `rows=sm_weighted;`、`code/weighted_checks.wl`。

第7轮：精确多中心膨胀障碍见 rounds/r007/report.md；MATLAB `result=sm_union_checks;`。

第8轮：完整证明与一般上界的显式维数障碍见 rounds/r008/report.md；MATLAB `rows=sm_large_soft;`，Wolfram `code/large_soft_check.wl`。

第9轮：非单调轮廓重排见 rounds/r009/report.md；MATLAB `result=sm_rearrangement_checks;`。

第10轮：缺口恒等式、缩放障碍与算子反例见 rounds/r010/report.md；Wolfram `code/deficit_checks.wl`。

第11轮：一般反馈上界归约见 rounds/r011/report.md；MATLAB `sm_bellman(h,256,outputPath)`，Wolfram `code/bellman_checks.wl`。

第12轮：一般首个次指数项与似然比鞅证明见 rounds/r012/report.md；MATLAB `result=sm_tilt_checks;`，Wolfram `code/tilt_checks.wl`。

第13轮：一般完整三项展开见 rounds/r013/report.md；MATLAB `sm_controlled_gaussian_checks`，Wolfram `code/stopping_checks.wl`。

第14轮：隐式正态根和任意固定阶逆对数展开见 rounds/r014/report.md；MATLAB `sm_root_checks`，Wolfram `code/root_checks.wl`。
