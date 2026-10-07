# C04 来源与实际处理状态

附件归属 42 条；本专题关联不同PDF 41 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E1592|04_Bauß_2023_Adaptive Improvements of Multi-Objective Branch and Bound.pdf|P-83bd4ae20d216ff0|direct|摘要明确提出多目标优化B&B算法并研究上界集、下界集、节点顺序与改进策略；直接涉及多目标节点界和全局搜索。|
|E1593|33_Christophel_2017_Shaping and Trimming Branch-and-bound Trees.pdf|P-f5acc87285f8c7a3|direct|提出MILP的offshoot搜索树，讨论改变分支变量顺序(shaping)和删去不必要分支变量(trimming)，且普通B&B是其特例。|
|E1594|05_Dym_2020_Quasi Branch and Bound for Smooth Global Optimization.pdf|P-922909d5e8e47a30|direct|C04 直接相关：全文给出平滑立方体全局优化的 quasi-BnB 规则、二阶/三阶界及广度优先有限终止证明。关键范围 lines 252–279、309–444、471–603、992–1140；限制在连续 f、已知 Lipschitz 界、规定盒分割等条件，不能据此宣称中心立方体算子的解析估计。|
|E1595|07_Bauß_2023_Augmenting Bi-objective Branch and Bound by Scalarization-Ba.pdf|P-c1ec8eed4a071c5b|direct|改进双目标B&B，用标量化信息和超体积间隙选择节点并强化搜索；直接讨论上下Pareto界与节点策略。|
|E1596|17_Wu_2024_A Branch and Bound Algorithm for Multiobjective Optimization.pdf|P-4734890e4114b7ab|direct|提出带一般序锥偏好的多目标B&B，利用锥支配做节点剔除并近似Pareto有效解集。|
|E1597|11_Muñoz_2022_Compressing Branch-and-Bound Trees.pdf|P-517b37b72548ef9b|direct|研究可认证整数规划对偶界的B&B树压缩问题，分析后处理改进与算法树规模，直接关系搜索证明结构。|
|E1598|13_Woodman_2025_Modern Computational Methods in Reinsurance Optimization_ Fr.pdf|P-ef9138d2fd5ce0fd|adjacent|应用论文比较再保险优化的模拟退火与量子可处理B&B，含实践模型和实验；可迁移搜索方案，但领域特定且非通用界/收敛理论。|
|E1599|37_Hager_2009_An ellipsoidal branch and bound algorithm for global optimiz.pdf|P-aec71cadd9f56fda|direct|构造椭球分割的连续全局优化B&B，以仿射下估计提供节点下界，并求解凸松弛。|
|E1600|38_Pataki_2009_Basis Reduction, and the Complexity of Branch-and-Bound.pdf|P-d3ca9e7e3e081918|direct|研究短而近正交约化基对整数可行问题B&B最坏/平均复杂度的影响，给出复杂度保证。|
|E1601|15_Wu_2023_Reference-Point-Based Branch and Bound Algorithm for Multiob.pdf|P-f47854600849e70a|direct|基于参考点偏好设计多目标B&B丢弃测试，证明感兴趣Pareto区域内ε有效解和有限迭代界。|
|E1602|27_Dey_2020_Branch-and-Bound Solves Random Binary IPs in Polytime.pdf|P-dc62681cfc628ced|direct|理论分析标准变量分支B&B在随机二元整数规划上的多项式时间表现；直接提供复杂度条件。|
|E1603|39_Goldsztejn_2008_Revisiting the upper bounding process in a safe Branch and B.pdf|P-ca0d63c86e40ac6c|direct|连续安全B&B中改进可行点和节点上界的认证计算策略，直接对应安全剪枝所需的可行incumbent。|
|E1604|32_Liu_2017_A branch and bound algorithm for the robust parall machine s.pdf|P-ca132dfd25f64256|direct|为带随机加工时间和顺序相关准备时间的并行机调度构造B&B，给出节点上下界、分支规则与支配规则。|
|E1605|03_Bauß_2023_Adapting Branching and Queuing for Multi-objective Branch an.pdf|P-7775bcfc546bbbc1|direct|研究多目标B&B的分支变量与节点队列策略，关注上下界Pareto集的规模及分布。|
|E1606|14_Byun_2024_Branch-and-bound algorithm for efficient reliability analysi.pdf|P-1ff127abaef3ffee|adjacent|在相干系统可靠性分析中用B&B构造系统失效/存活事件的可复用表示；属特定概率可靠性应用，节点剪枝思想可迁移。|
|E1607|09_Naguib_2023_On Statistical Learning of Branch and Bound for Vehicle Rout.pdf|P-b1fd9bcac1b78e51|adjacent|用图神经网络模仿CVRP的强分支决策，是B&B分支策略的机器学习启发式；不提供正确性界或全局最优证明的新定理。|
|E1608|02_Aldana-López_2023_Designing controllers with predefined convergence-time bound.pdf|P-974c0c803e3a7dd2|irrelevant|首页和摘要研究预设时间收敛控制器及有界时变增益；“bound”是时间界，未研究branch-and-bound、节点界或搜索树。|
|E1609|13_Gupta_2022_Branch-and-Bound Performance Estimation Programming_ A Unifi.pdf|P-6e660372e6874d0a|direct|提出BnB-PEP，以定制B&B全局求解非凸QCQP来自动设计最优一阶优化方法；直接研究B&B构模与全局求解。|
|E1610|23_Knaeble_2024_Branch and Bound to Assess Stability of Regression Coefficie.pdf|P-6d4cc057dc81beda|direct|用B&B搜索离散正则化回归模型空间，计算不确定模型下回归系数的最大/最小值并提供数学支撑。|
|E1611|12_Lindstrom_2020_Error bounds, facial residual functions and applications to.pdf|P-4bf3149d7443892c|irrelevant|摘要处理指数锥可行问题误差界和facial residual functions；与B&B搜索树/节点剪枝无关，属于锥可行性误差分析。|
|E1612|24_Dey_2021_Lower bound on size of branch-and-bound trees for solving lo.pdf|P-40446b9407b8d603|direct|构造lot-sizing实例族，证明即使允许一般split disjunction，任何B&B证明树都需指数节点；直接是树大小下界。|
|E1613|01_Sciandra_2024_Supplementary Materials to Graph Convolutional Branch and Bo.pdf|P-3701734af53a6778|adjacent|这是图卷积增强B&B的补充材料，包含学习节点优先级/最优性评分的模型信息；启发式可迁移，但其补充性质及理论界需与主文一并核对。|
|E1614|40_Chen_2008_Parallel Branch and Bound Algorithm for Computing Maximal St.pdf|P-a171e89f663f3d98|direct|提出并行B&B计算structured singular value上界/最大值，不要求逐频率紧界，属于有保证的连续鲁棒控制全局搜索。|
|E1615|30_Kazachkov_2021_An Abstract Model for Branch and Cut.pdf|P-642917e1219e52ae|adjacent|抽象建模branch-and-cut（B&B搜索加切平面），扩展到切平面策略权衡；与标准B&B相邻，不能不加区分地当作纯B&B算法结论。|
|E1616|29_Huang_2021_Branch and Bound in Mixed Integer Linear Programming Problem.pdf|P-2d230e96b5235f02|adjacent|综述MILP B&B四个组件（分支、节点选择、剪枝和切平面）及学习加速；是方法地图而非原始定理证据。|
|E1617|25_Huang_2025_Non-linear Multi-objective Optimization with Probabilistic B.pdf|P-07a73c50d2215d3f|direct|提出非线性多目标概率B&B，处理概率可行性和全局Pareto优化；直接对应B&B节点和概率界机制。|
|E1618|31_Luo_2019_A Geometric Branch and Bound Method for a Class of Robust Ma.pdf|P-e2c7deeb84592370|direct|针对有限不确定集下凸函数鲁棒最大化提出几何B&B，以分段线性近似生成节点上下界并证明ε有限收敛。|
|E1619|10_Gupta_2022_Branch-and-Bound Performance Estimation Programming_ A Unifi.pdf|P-6e660372e6874d0a|direct|与E1609同一paper_id和相同BnB-PEP全文，仅归档文件名不同；重复目录项，不重复阅读。|
|E1620|08_Häner_2024_Solving QUBOs with a quantum-amenable branch and bound metho.pdf|P-ecab669061ed516b|direct|构造精确QUBO B&B并利用低成本节点界，报告实验；直接研究QUBO搜索与可认证剪枝。|
|E1621|26_Wu_2024_Obtaining properly Pareto optimal solutions of multiobjectiv.pdf|P-a101e5ef3c9891f8|direct|用B&B为多目标问题输出properly Pareto-optimal解，减少完整Pareto前沿带来的决策负担。|
|E1622|20_Gurvits_2020_Capacity Lower Bounds via Productization.pdf|P-c3c8be58dbbbecf3|irrelevant|摘要研究实稳定齐次多项式capacity界和scaling算法；未涉及branch-and-bound/搜索树，分类词面误命中。|
|E1623|35_Mohammadi_2016_Combining SOS and Moment Relaxations with Branch and Bound t.pdf|P-deb72842d1067b4f|direct|将SOS/矩松弛接入多项式盒B&B；原8页PDF完整，旧两栏次序导致的截断判断已撤回。Lemma3严格域a≤c<d≤b已局部独立核验；完整SOS层/定理1及外引证明链仍待核，本补正新增全文阅读计数为0。|
|E1624|22_Gläser_2023_Sub-Exponential Lower Bounds for Branch-and-Bound with Gener.pdf|P-dd61b6dcf6e0fb40|direct|针对一般整数分支不等式B&B树，证明特定紧凑整数规划实例存在次指数规模下界；直接分析证明树复杂度。|
|E1625|36_Kitahara_2011_A Bound for the Number of Different Basic Solutions Generate.pdf|P-b5b02b5d8e1be50f|irrelevant|摘要研究单纯形法产生的不同基本可行解数量和迭代界；“basic solutions”并非B&B树，主题误命中。|
|E1626|19_Qu_2022_Yordle_ An Efficient Imitation Learning for Branch and Bound.pdf|P-cae1d2a50e7a1ed0|adjacent|提出YORDLE模仿学习以替换组合优化B&B中的启发式决策，主旨是数据驱动分支加速而非可认证正确性。|
|E1627|21_Legat_2020_Abstraction-based branch and bound approach to Q-learning fo.pdf|P-6a8a1d69a6a50866|adjacent|把抽象与B&B用于混合最优控制/近似动态规划；提供控制特定框架，迁移相关但与通用整数规划节点界不同。|
|E1628|16_Scavuzzo_2024_Machine Learning Augmented Branch and Bound for Mixed Intege.pdf|P-9e98fef276e7de34|adjacent|为MILP B&B加入学习方法辅助决策，系统总结与实验学习策略；正确性仍由底层精确求解器承担。|
|E1629|18_Paulus_2023_Learning To Dive In Branch And Bound.pdf|P-92a05ebbdab4c481|adjacent|学习B&B diving primal heuristics以寻找可行整数解；服务incumbent构造但不产生独立节点上界证书。|
|E1630|34_Adelgren_2017_Branch-and-bound for biobjective mixed-integer linear progra.pdf|P-d06622775790e3db|direct|构造双目标MILP的通用B&B求解Pareto解，给出节点对偶界、fathoming和gap等算法组件。|
|E1631|28_Chakrabarti_2022_Universal Quantum Speedup for Branch-and-Bound, Branch-and-C.pdf|P-36f66cd67e257dcd|adjacent|分析量子加速B&B/B&C/tree-search的算法复杂度；能迁移到树搜索工作量估计，但量子模型和目标不同于经典认证界。|
|E1632|06_Wu_2022_The Hybridization of Branch and Bound with Metaheuristics fo.pdf|P-235b40772166d034|direct|将多目标进化算法的上下界与B&B结合处理非凸多目标优化，并以B&B保证全局收敛。|
|E1633|10_Dey_2021_Lower Bounds on the Size of General Branch-and-Bound Trees.pdf|P-12964b7ec923eda9|direct|为general branching rules的B&B树构造packing、set covering和TSP下界实例，证明指数规模；直接分析B&B证明复杂度。|

## 实际阅读卡

- [P-922909d5e8e47a30：full_read](papers/AD-P-922909d5e8e47a30.json)
- [P-deb72842d1067b4f：scoped correction，整条SOS链pending](papers/AD-P-deb72842d1067b4f.json)

## 本地原文定位

原附件PDF：工作包根目录的 `evidence/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`evidence/supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。

## 有效状态与历史

本表与当前纸卡采用[来源完整性更正](source-completeness-correction.md)的有效状态；[原生成卡历史](history/AD-P-deb72842d1067b4f-before-completeness-correction.json)保留旧字节与旧判断，不能继续把旧partial_read的物理截断理由当作当前事实。原AD卡未修改，未新增full-read计数。
