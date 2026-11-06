# 非中心强极大算子端点常数

研究对象为 Lebesgue 测度下轴平行非中心强矩形极大算子，统一使用
\(\Phi_n(t)=t(1+\log^+t)^{n-1}\)、自然对数及阈值 1。
本项目研究最优比值 \(C_n=\sup_f |\{M_Sf>1\}|/\int\Phi_n(f)\)。

## 当前状态（2026-11-06，第 37 轮）

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
- **第15轮精确指数控制**：负Laplace目标的可预测最优值为exp(−n lambda c_lambda)，由唯一固定动作达到，并有逐步非负损失恒等式。此障碍仅限指数目标，不覆盖尾概率或端点；见 rounds/r015/report.md。
- **第16轮真实分位数定理**：固定c=mu−(sigma/2)sqrt(log n/n)相对c=mu的临界下分位数改善为(mu−1)log n/8+o(log n)。证明使用精确CDF参数导数与自包含Fourier密度估计；仅限该测试族，不是C_n的新对数级项。见 rounds/r016/report.md。
- **第17轮实际下界族改善**：在固定斜率1/n的指数软轮廓中，指定比例偏移使真实见证评分的对数提升(mu−1)log n/8+o(log n)；已结算严格接收与I_n分母。此为家族内改善，未缩小一般C_n上界。见 rounds/r017/report.md。
- **第18轮绝对下界项**：固定mu指数斜坡的软根已精确到o(1)，从而得到该实际下界族自身的乘法渐近等价；结合第17轮偏移，一般log C_n≥nmu−d_n−(mu+3)log n/24+o(log n)。这仍不是一般C_n的乘法等价。见 rounds/r018/report.md。
- **第19轮窄平台改善**：宽w=log n/2+loglog n的高度1平台使实际分母因子趋于1，只将软根降低1/2+o(1)；一般log C_n≥nmu−d_n+(9−mu)log n/24+o(log n)，对数级系数转为正。见 rounds/r019/report.md。
- **第20轮一般上界**：有限时域Gaussian离散屏障给log C_n≤nmu−d_n+O((log n)²)。一般上下界对数差距现为O((log n)²)，真正log n系数及一般乘法等价仍未闭合。完整证明见 rounds/r020/report.md。
- **第21轮一般对数级上界**：精确非线性Gaussian移位、C³截断与热半群积分给log C_n≤nmu−d_n+O(log n)。结合第19轮下界，一般上下界对数差距现为O(log n)，最优系数及乘法等价仍待研究。见 rounds/r021/report.md。
- **第22轮Q_n精确对数项**：任意确定性独立比例向量和任意单调轮廓的log Q_n=nmu−d_n+(9−mu)log n/24+o(log n)；机制上界匹配第19轮下界。此结论未扩展至一般C_n或反馈模型。见 rounds/r022/report.md。
- **第23轮实际相关输入反例**：一个缺角两层输入严格优于相同振幅上限下的全部常幅矩形；振幅上限1.01时认证比值≥1.020472，而常幅矩形≤1.019754。固定二维近单位收益至少gamma epsilon+O(epsilon²)，gamma=(6log(4/3)−1)/8>0；不改变一般C_n渐近结论。见 rounds/r023/report.md。
- 第24轮：高度1+1/n的满立方体近常值输入类最优比值趋于I_0(2sqrt2)，并证明三个小振幅尺度上的满立方体支撑障碍。见 rounds/r024/report.md。
- 第25轮：宏观角分块的完整锚定weak-halo仍被Q_n控制；辅助反馈热方程的精确线性响应含g log t，但未据此声称离散反馈渐近。见 rounds/r025/report.md。
- 第26轮：固定块数、次指数高度及任意确定比例的单锥接收评分≤C_m sqrt(n/∏b_j) Q_n；平衡块损失(m−1)log n/2。见 rounds/r026/report.md。
- 第27轮：共同参数的增量卷积密度在四阶起对数凹，三阶有精确反例；完成曲线接收体积及尾截断，给出精确反馈核参数导数和局部telescope。见 rounds/r027/report.md。
- **第28轮全输入固定比例障碍**：双射线全局重排使任意非负输入的固定比例反射锚定族端点精确等于对应标量Q（含有限高度）；K组固定锚点/比例的逐点最大至多KQ_n。另严格量化中等下尾带的真实一步贪婪收益p x² f/(8k²)，累计占据量仍待证明。见 rounds/r028/report.md。
- **第29轮反馈模型分离**：精确Poisson校正子、全状态比较及有界软轮廓证明liminf(log mathcalR_n−log Q_n)/log n≥(mu−1)/8；从而log mathcalR_n≥nmu−d_n+(3+mu)log n/12+o(log n)。这仅是上界模型的下界和机制障碍，未提高真实C_n下界。见 rounds/r029/report.md。
- **第30轮连续预算与熵障碍**：有限真实接收点预算有达到的精确强对偶及盒/尾认证接口；单步核移动跳点使KL随动作偏移线性，排除二次Fisher路径比较。独立比例交替仅在离散族内改善，仍低于单块参考；固定旧输入的跳跃分段积分把原网格平均1.00061修正为约0.99809，未登记连续下界。见 rounds/r030/report.md。
- **第31轮块正则性**：b≥2块核具有二次局部KL；有效界为O(bh²)。严格证明I_b(mu)→p_mu²/2及精确矩下界，同时给出b h⁴漂移障碍，禁止误用全块长统一O(h²)界。旧接收输入连续数值修复后平均约1.0000304，但仍低于单块参考。见 rounds/r031/report.md。
- **第32轮任意输入接口**：对所有可预测策略证明实际Orlicz预算的Fenchel高层尾吸收，将任意单调轮廓归约到指数/线性单终端，误差O(log n/n)。精确独立桥公式给出g log n源质量及全局对数包络反例；有限正比例参数无一阶随机序。局部接收比例小收益被统一比例对照超过，未新增C_n认证下界。见 rounds/r032/report.md。
- **第33轮正确终端结构**：证明payoff Bellman全局唯一动作根与前沿恒等式、非负running-cost等价、带中心余项的局部全动作收益；给出固定短时Gamma平移初始化和正指数外侧屏障。完整MAX递推至128仅为插值数值，global-c校准后一个独立坐标正向探针未增益，sharp全局校正仍未闭合。见 rounds/r033/report.md。
- **第34轮全控制外侧接口**：构造精确仿射/Gamma终端修补，给出自然时间管两侧的全策略边界及固定短时起点的参考桥外带源求和；纯Exp占据反例与加权大跳相消排除两种拼接捷径，结合增长屏障证明有缓冲内部的全动作跨界界。小坐标扰动获益仍被统一比例细化超过，未新增C_n认证下界。见 rounds/r034/report.md。
- **第35轮局部源与全局反例**：匹配参考局部全动作收益保留bQ/n并有可求和误差；证明固定参考多项式校正直接与两侧屏障取min在正大偏差区必失败。给出真实单分量halo外的跨分量接收方盒，以及n≥8同步二维积见证的有限分母障碍。一般锐系数仍未闭合。见 rounds/r035/report.md。
- **第36轮非线性局部闭合**：参考payoff对数凹性给出全动作Doob指数尾，控制指数化校正的Jensen缺陷；精确矩消除闭合自然时间管内MAX势，临界附加误差仅O(f_n)。全状态延拓仍未闭合。曲线/source联合变形按有限锥体积验收后无收益，三种接收点依赖方向也未胜过同输入对照。更正第34轮接收试验维数为n=64（256为求积节点）。见 rounds/r036/report.md。
- **第37轮晚时段全空间比较**：匹配仿射余量使局部非线性势在 t≥N^gamma 可与全局外屏障严格拼接；初始正缺口面积为 O(log t)，但常值偏移给出任意面积误差传播的反例，尚缺带控制损失的初始化比较。Mathematica 30组双侧变换诊断及5组真实 mixed source/receiver 均未登记新认证下界。见 rounds/r037/report.md。
- 主攻：带漂移损失的初始化误差传播；下界继续联动增长模板数、source、receiver 与实际并集体积。
- 一般前三项现已闭合；乘法渐近等价、后续锐余项和极值结构仍待研究。

## 复现

计算优先使用WolframScript。MATLAB按需短批运行，保存结果后立即退出，不保留各子代理的常驻引擎。

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

第15轮：精确指数控制和局部动作偏移见 rounds/r015/report.md；MATLAB `sm_laplace_checks`，Wolfram `code/laplace_checks.wl`。

第16轮：真实下分位数的对数级改善见 rounds/r016/report.md；MATLAB `sm_quantile_derivative_checks`，Wolfram `code/quantile_derivative_checks.wl`。

第17轮：实际指数软轮廓下界的对数级改善见 rounds/r017/report.md；MATLAB `sm_soft_shift_checks`，Wolfram `code/soft_shift_checks.wl`。

第18轮：固定比例斜坡的乘法渐近及一般显式下界见 rounds/r018/report.md；MATLAB `sm_edgeworth_ramp_checks`，Wolfram `code/edgeworth_ramp_checks.wl`。

第19轮：窄平台实际下界见 rounds/r019/report.md；MATLAB `sm_plateau_checks`，Wolfram `code/plateau_checks.wl`。

第20轮：一般上界的对数平方余项见 rounds/r020/report.md；MATLAB `sm_gaussian_barrier_checks`，Wolfram `code/gaussian_barrier_checks.wl`。

第21轮：一般对数级上界见 rounds/r021/report.md；MATLAB `sm_shift_barrier_checks`，Wolfram `code/shift_barrier_checks.wl`。

第22轮：独立比例机制的精确对数项及初始化障碍见 rounds/r022/report.md；MATLAB `sm_independent_log_checks`，Wolfram `code/independent_log_checks.wl`。

第23轮：实际缺角输入与近单位振幅见 rounds/r023/report.md；Python `code/cert_capped_corner.py`，MATLAB `sm_capped_corner_checks`，Wolfram `code/capped_corner_checks.wl`。

第24轮：近常值满立方体类的精确极限与三尺度支撑障碍见 rounds/r024/report.md；MATLAB `sm_nearflat_ceiling_checks`，Wolfram `code/nearflat_ceiling_checks.wl`。

第25轮：宏观角分块的锚定障碍与反馈辅助热方程见 rounds/r025/report.md；MATLAB `sm_feedback_probe`，Wolfram `code/angular_feedback_checks.wl`。

第26轮：微观块凸预算与单锥接收的对数损失见 rounds/r026/report.md；MATLAB `sm_microblock_dual_checks`，Wolfram `code/microblock_checks.wl`。

第27轮：共同参数卷积的对数凹证书、曲线接收边界及精确反馈核见 rounds/r027/report.md；Python `code/cert_block_logconcavity.py --verify`，MATLAB `sm_microblock_halo`，Wolfram `code/block_logconcavity_checks.wl`。

第28轮：固定比例反射锚定族的全输入精确端点与带内反馈收益见 rounds/r028/report.md；MATLAB `sm_greedy_band_checks`，Wolfram `code/fixed_anchor_checks.wl`。
