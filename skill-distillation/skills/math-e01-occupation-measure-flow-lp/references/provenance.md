# E01 来源与实际处理状态

附件归属 42 条；本专题关联不同PDF 42 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E1970|26_Han_2018_Controller Synthesis for Discrete-Time Polynomial Systems vi.pdf|P-aac215965bebc368|direct|离散多项式控制系统以受控 Liouville 方程和占用测度 LP 表示，直接展示流守恒建模。|
|E1971|36_Borkar_2011_Asymptotics of the Invariant Measure in Mean Field Models wi.pdf|P-6fe71a3ceded8417|pending|摘要讨论耦合 Markov 链不变测度与大偏差；尚未定位该文主定理，暂不判作占用 LP 证据。 主结论未在已读范围内定位，暂列待核。|
|E1972|05_Mangold_2024_SCAFFLSA_ Taming Heterogeneity in Federated Linear Stochasti.pdf|P-5df1d37d0abf87e8|irrelevant|联邦线性随机逼近与 TD 学习的样本/通信复杂度，与占用测度流 LP 无关。|
|E1973|27_Josz_2018_Transient Stability Analysis of Power Systems via Occupation.pdf|P-b3dff2f748c1f95d|direct|全文已读：以占用测度和 Liouville 方程 LP 估计 ROA 体积，给出矩/SOS 外逼近；可作占用流模型范例，但迁移至立方体极大算子需另证。|
|E1974|17_Miller_2023_Peak Estimation of Time Delay Systems using Occupation Measu.pdf|P-a97f4407517c8a48|pending|以占用测度 LP 上界估计时滞微分方程轨迹峰值。 主结论未在已读范围内定位，暂列待核。|
|E1975|08_Korda_2022_The gap between a variational problem and its occupation mea.pdf|P-1377aa658f9785aa|pending|研究变分问题与占用测度松弛的等价性及凸性条件。 主结论未在已读范围内定位，暂列待核。|
|E1976|34_Krivulin_2012_Using max-algebra linear models in the representation of que.pdf|P-3d4a2e81b8be7fbf|irrelevant|max-plus 线性代数用于排队系统表示，与测度 LP 的流约束不同。|
|E1977|31_Claeys_2014_Reconstructing trajectories from the moments of occupation m.pdf|P-c4534358e3e04fe8|adjacent|全文已读：研究由 occupation measure 有限阶矩恢复近似轨迹/控制节点，是流 LP 的逆问题后处理；重构不保证可行，需局部直接法验证。|
|E1978|17_Ahmad_2025_Impact of occupancy patterns on buildings performances gaps.pdf|P-e8c31b56f76bacb8|irrelevant|医院术后床位占用预测是调度应用，不是数学占用测度。|
|E1979|15_Alexander_2021_Late levels of nested occupancy scheme in random environment.pdf|P-b9e72f8116a5d4d3|irrelevant|随机环境下 balls-in-boxes 的嵌套占用方案，术语相同但对象不同。|
|E1980|29_Kuntz_2018_Approximations of countably-infinite linear programs over bo.pdf|P-92b542665234ed99|adjacent|给出有界测度空间上的可数无限 LP 有限维近似与误差界，支持离散化思路但非特定控制流式。|
|E1981|24_Chen_2019_Safety Verification of Nonlinear Autonomous System via Occup.pdf|P-38eb110ed131a199|direct|全文含附录已读：occupation-measure 流 LP 精确编码固定 ODE 与初始分布下危险区暴露时间，并证明 LP 强对偶和矩/SOS 上界收敛；几何极大函数迁移仍需新建模型。|
|E1982|40_Henrion_2008_Nonlinear optimal control synthesis via occupation measures.pdf|P-0a6262177bc3b3c4|direct|全文17页已读：连续时间多项式OCP的三测度流守恒、矩/SOS下界与HJB次解直接展示占用LP操作；控制构造只给算法和例子，不是一般最优性定理。|
|E1983|32_Claeys_2014_Modal occupation measures and LMI relaxations for nonlinear.pdf|P-0297391ab7ac5c7b|direct|32页全文给出finite-mode occupation projection/merge、weak LP与HJB dual、Archimedean moment/SOS hierarchy及joint time-marginal reconstruction；原switching无gap条件和source example矛盾分开保留。|
|E1984|37_Zinchenko_2010_Shrink-Wrapping trajectories for Linear Programming.pdf|P-30419b83bb1d60e6|irrelevant|超双曲规划中的 shrink-wrapping 轨迹几何，与占用测度流无关。|
|E1985|18_Lahiri_2025_Linear models of dynamic optimization with linear constraint.pdf|P-0f3c0580a1eee3c4|irrelevant|无限时域线性动态规划与 Euler/横截条件，非测度占用 LP。|
|E1986|19_Jaffuel_2020_Occupation measures arising in finite stochastic games.pdf|P-a90820e5d37395da|adjacent|研究有限随机博弈的折扣极限占用测度，可参考 MDP 流解释但目标是渐近博弈。|
|E1987|03_Henrion_2023_Occupation measure relaxations in variational problems_ the.pdf|P-b8164dc1ce4dcf5d|direct|正文研究 calculus-of-variations 的 occupation-measure 线性松弛、弱 Liouville 约束与原问题无间隙；核心覆盖凸情形、非凸 affine 松弛与 convex-envelope 等价，并推广到最优控制。PDF 首页题名、作者与 arXiv:2303.02434v1 已核对。|
|E1988|28_Spratt_2018_Reducing post-surgery recovery bed occupancy with a probabil.pdf|P-15010e8e8618b009|irrelevant|术后恢复床位随机排程中的普通 occupancy 指标。|
|E1989|07_Yu_2026_An Occupation-Measure and Frank-Wolfe Framework for Heteroge.pdf|P-22e8303b48a27276|direct|论文直接将异质多人口均值场控制写成占用测度优化，以线性弱 Liouville 约束编码动力学，并证明对称化核凸性条件、人口分离 FW oracle 与有限曲率率；PDF 首页题名/作者/arXiv:2607.09907v2 已核对。|
|E1990|25_Mania_2019_Certainty Equivalence is Efficient for Linear Quadratic Cont.pdf|P-726e57aeb50cc57f|irrelevant|未知线性系统上的 certainty-equivalent 控制性能界，不是占用 LP。|
|E1991|20_Ekbatani_2021_Circuit imbalance measures and linear programming.pdf|P-63696ad045602440|irrelevant|线性规划 circuit imbalance 与 polyhedral 算法性质。|
|E1992|22_Suttle_2022_Occupancy Information Ratio_ Infinite-Horizon, Information-D.pdf|P-ff629c9346640b96|adjacent|MDP 策略的 state occupancy entropy ratio 与占用分布有关，但非流守恒 LP 主结果。|
|E1993|15_Hinz_2025_On fractal minimizers and potentials of occupation measures.pdf|P-d5edae09546b81fc|adjacent|研究高斯场占用测度的势及其连续性，属于 potential of occupation measures，不是占用 LP。|
|E1994|11_Miller_2023_Peak Value-at-Risk Estimation of Stochastic Processes using.pdf|P-610d8ca4e35b8a23|direct|用 occupation-measure Dynkin LP 构造随机过程沿轨迹的峰值 VaR/ES 上下界与 moment-SOS SDP；是风险峰值分析对占用测度方法的直接应用。|
|E1995|35_Dai_2011_A Gel'fand-type spectral radius formula and stability of lin.pdf|P-3c9937944470e5da|irrelevant|线性切换系统谱半径与稳定性公式。|
|E1996|12_Feng_2020_Dynamics Compensation in Observation of Abstract Linear Syst.pdf|P-e209a35d71d8eb8d|irrelevant|抽象线性系统传感器动态补偿和观测器稳定性。|
|E1997|04_Holtorf_2022_Stochastic Optimal Control via Local Occupation Measures.pdf|P-3b44b82c06eccd85|direct|扩散和跳过程随机最优控制的局部占用测度、弱 Dynkin 传输关系及分片 moment-SOS 松弛是论文核心主题。|
|E1998|21_Yu_2026_Occupation-Measure Mean-Field Control_ Optimization over Mea.pdf|P-078f65ff156d77ea|direct|本文直接研究 mean-field control 的 occupation-measure 优化和线性 Liouville 流约束，给出测度目标凸性/存在、FW oracle classical-trajectory 实现及 O(1/k) 收敛结果。PDF 首页题名、作者及 arXiv:2603.16094v1 已核对。|
|E1999|06_Piunovskiy_2023_Extreme occupation measures in Markov decision processes wit.pdf|P-31afba7e98fa0ddc|direct|Borel MDP 有限占用测度的极点刻画及确定性平稳策略生成。|
|E2000|23_Shen_2026_Distributed Variational Quantum Linear Solver.pdf|P-9a687ff59723f6a5|irrelevant|分布式变分量子线性求解器与 LP 占用流无关。|
|E2001|33_Deshpande_2014_On the role of Föllmer-Schweizer minimal martingale measure.pdf|P-1d85e482cfab9fa2|irrelevant|风险敏感资产管理中的最小鞅测度和控制博弈。|
|E2002|10_Dereziński_2026_Towards Universal Convergence of Backward Error in Linear Sy.pdf|P-9b72dfbb79eaa866|irrelevant|线性方程迭代求解器 backward error 收敛分析。|
|E2003|39_Malesevic_2008_Differential transcendency in the theory of linear different.pdf|P-1906c827cdb5a436|irrelevant|常系数线性微分系统与微分超越性。|
|E2004|14_Fantuzzi_2022_Sharpness and non-sharpness of occupation measure bounds for.pdf|P-f8e1ed1dfb1932d9|direct|比较占用测度与点态对偶松弛，研究变分问题的尖锐性与非尖锐例子。|
|E2005|01_Liang_2025_Mean Field Game with Reflected Jump Diffusion Dynamics_ A Li.pdf|P-84922ac336b62d22|direct|反射跳扩散均值场博弈的占用测度 LP 与弱松弛控制等价。|
|E2006|13_Feng_2020_Actuator Dynamics Compensation in Stabilization of Abstract.pdf|P-42be6b04b99d8d48|irrelevant|执行器动力学补偿的抽象线性系统稳定化。|
|E2007|30_Chermakani_2015_Optimal Aggregation of Blocks into Subproblems in Linear-Pro.pdf|P-aeac4bf34c3b2121|irrelevant|按计算平台聚合分块 LP 子问题，研究复杂度而非占用测度。|
|E2008|02_Piunovskiy_2020_Aggregated occupation measures and linear programming approa.pdf|P-f4bae8c8bfb9d2df|pending|受约束冲击控制以占用测度和聚合占用测度构造等价 LP。 主结论未在已读范围内定位，暂列待核。|
|E2009|38_Morio_2009_On the Computation of $π$-Flat Outputs for Linear Time-Delay.pdf|P-e71e3e6254b8ca98|irrelevant|线性时滞系统 π-flatness 与多项式矩阵分解。|
|E2010|09_Alfaro_2020_Computing sandpile configurations using integer linear progr.pdf|P-2d44b1b98bb30f36|irrelevant|沙堆群的 recurrent configuration 与整数规划。|
|E2011|16_Piunovskiy_2023_On the continuity of the projection mapping from strategic m.pdf|P-b8fd135108e28ef5|adjacent|刻画吸收 MDP strategic measure 到 occupation measure 的弱拓扑连续性与一致吸收条件。|

## 实际阅读卡

- [P-0297391ab7ac5c7b：full_read](papers/TEAM_EF-P-0297391ab7ac5c7b.json)
- [P-078f65ff156d77ea：full_read](papers/EG-P-078f65ff156d77ea.json)
- [P-0a6262177bc3b3c4：full_read](papers/EG-P-0a6262177bc3b3c4.json)
- [P-22e8303b48a27276：full_read](papers/EG-P-22e8303b48a27276.json)
- [P-38eb110ed131a199：full_read](papers/EG-P-38eb110ed131a199.json)
- [P-3b44b82c06eccd85：full_read](papers/EG-P-3b44b82c06eccd85.json)
- [P-610d8ca4e35b8a23：full_read](papers/EG-P-610d8ca4e35b8a23.json)
- [P-aac215965bebc368：full_read](papers/EG-P-aac215965bebc368.json)
- [P-b3dff2f748c1f95d：full_read](papers/EG-P-b3dff2f748c1f95d.json)
- [P-b8164dc1ce4dcf5d：full_read](papers/EG-P-b8164dc1ce4dcf5d.json)
- [P-c4534358e3e04fe8：full_read](papers/EG-P-c4534358e3e04fe8.json)

## 本地原文定位

原附件PDF：工作包根目录的 `evidence/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`evidence/supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。
