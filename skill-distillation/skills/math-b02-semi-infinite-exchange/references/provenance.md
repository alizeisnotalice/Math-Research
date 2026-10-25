# B02 来源基线清单（旧标签待复核）

> 当前筛查说明：下方标题、作者/版本、direct/adjacent/irrelevant及理由均是基线记录或旧人工备注，不代表本轮已复核；本轮状态以交付根 `audit_current/pdf_screening_and_reading_queue.csv` 和共享读卡为准。列表中的专题附件数是路径/条目口径，不可直接当SHA去重的论文数。题录TXT不是全文证据。


附件归属 40 条；本专题关联不同PDF 40 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E1445|37_Wang_2013_Semidefinite relaxations for semi-infinite polynomial progra.pdf|P-ecb91aed2e290a23|direct|研究半无限多项式规划，并用交换/约束生成算法结合SDP求解；紧索引集与非紧集齐次化都直接覆盖本专题。|
|E1446|03_Huertas_2025_Constraint Programming Models For Serial Batch Scheduling Wi.pdf|P-a2dad4c87754714b|irrelevant|题名和摘要是批次排程的constraint-programming模型；此处constraint programming指约束规划建模，不是半无限约束生成。|
|E1447|13_Wehbeh_2024_Semi-Infinite Programs for Robust Control and Optimization.pdf|P-e82723d1821371e3|direct|鲁棒最优控制被明确表述为连续不确定集上的半无限规划，并扩展local-reduction方法处理存在约束；属于半无限算法应用。|
|E1448|14_Guo_2021_An SDP method for Fractional Semi-infinite Programming Probl.pdf|P-3370cc64f9cf50d4|direct|研究SOS凸分式半无限多项式规划及其收敛SDP层级，直接相关。|
|E1449|29_Zhang_2017_Sparse Semidefinite Programs with Guaranteed Near-Linear Tim.pdf|P-4582d054b75b8d94|adjacent|研究稀疏SDP的dualized-clique-tree转换与复杂度；是SDP求解器技术，可辅助松弛实现，但不研究半无限约束/生成。|
|E1450|32_Basu_2015_Strong duality and sensitivity analysis in semi-infinite lin.pdf|P-21cdac95bef06357|direct|研究半无限线性规划的强对偶、灵敏度与约束空间扩张；可补充半无限问题对偶基础。|
|E1451|34_Basu_2013_Projection_ A Unified Approach to Semi-Infinite Linear Progr.pdf|P-d0dfb02b32ca64a5|direct|把Fourier–Motzkin投影扩展到半无限线性规划，并给出投影/对偶算法框架。|
|E1452|23_Dolgopolik_2019_A New Constraint Qualification and Sharp Optimality Conditio.pdf|P-0adbe9aa88a63b24|adjacent|研究一般非光滑数学规划的约束资格和最优性条件；与B04交叉，但题目/摘要未见半无限索引集。|
|E1453|25_Okuno_2018_Primal-dual path following method for nonlinear semi-infinit.pdf|P-ed4991b823121218|direct|研究带半定约束、无限凸不等式的非线性半无限规划，并给出primal-dual路径跟踪方法。|
|E1454|11_Wehbeh_2026_Generalized Semi-Infinite Programming for Robust Optimal Con.pdf|P-25ca932784db7e58|direct|研究决策相关不确定性下的广义半无限鲁棒控制及GSIP求解，直接覆盖不确定集依赖决策的半无限结构。|
|E1455|38_Mordukhovich_2011_Constraint Qualifications and Optimality Conditions for Nonc.pdf|P-87abf351b97093dc|direct|研究非凸半无限/无限规划的约束资格与最优性条件，是本专题的理论基础来源。|
|E1456|08_Hu_2023_Polynomial Optimization Relaxations for Generalized Semi-Inf.pdf|P-fb59ad9aea84dab7|direct|研究多项式广义半无限规划及其多项式优化松弛层级，直接覆盖GSIP算法。|
|E1457|20_Harwood_2019_A note on generalized semi-infinite program bounding methods.pdf|P-646526af87253fe6|direct|专门研究广义半无限规划的全局下界与bounding方法，直接相关。|
|E1458|26_Okuno_2018_An interior point sequential quadratic programming-type meth.pdf|P-ac31185ed67f3b33|direct|研究带log-det和无限凸不等式的半无限规划，并提出内点SQP型方法。|
|E1459|28_Chow_2018_On Dynamic Programming Principle for Stochastic Control unde.pdf|P-935d169aba15d33c|irrelevant|研究随机控制中的逐时刻期望约束与动态规划原理；约束是概率/期望型而非对连续索引集合的无限不等式。|
|E1460|30_Guo_2015_Semidefinite programming relaxations for linear semi-infinit.pdf|P-4033503c69070497|direct|研究线性半无限多项式规划及其SDP松弛层级，直接相关。|
|E1461|18_Arima_2024_Exact SDP relaxations for a class of quadratic programs with.pdf|P-01257b41bc5a72eb|direct|研究具有有限与无限二次约束的非凸二次规划之精确SDP松弛，属于半无限二次约束求解。|
|E1462|21_Mordukhovich_2019_New extremal principles with applications to stochastic and.pdf|P-ff08a2dba86a5855|adjacent|变分分析极值原理用于随机与半无限规划；可迁移的理论工具，主贡献是极值原理而非SIP算法。|
|E1463|09_Das_2024_Data-driven distributionally robust MPC for systems with mul.pdf|P-25edfbd4a946fb56|direct|明确提出分布鲁棒MPC的半无限半定规划模型与算法，属半无限鲁棒优化应用。|
|E1464|19_Harwood_2019_A note on semi-infinite program bounding methods.pdf|P-3580cf604bdcbfd2|direct|研究半无限规划的bounding方法，可直接迁移到约束生成过程中的上下界控制。|
|E1465|35_Basu_2013_On the sufficiency of finite support duals in semi-infinite.pdf|P-b7e6d327ca15e5eb|direct|研究半无限线性规划有限支撑对偶何时充分，直接涉及有限活跃约束/对偶结构。|
|E1466|10_Huang_2024_Piecewise SOS-Convex Moment Optimization and Applications vi.pdf|P-ec78cd78268be31b|adjacent|研究SOS凸分段函数的矩优化与精确SDP；与多项式/SIP松弛邻近，但摘要不涉及半无限约束。|
|E1467|40_Gouveia_2009_A new semidefinite programming hierarchy for cycles in binar.pdf|P-1610c4a501d1f46d|irrelevant|研究二元拟阵环与图割的theta-body半定层级；“semi-definite”只是SDP术语，不是semi-infinite规划。|
|E1468|16_Wehbeh_2025_State-Dependent Uncertainty Modeling in Robust Optimal Contr.pdf|P-ef89e86fccc03703|direct|用广义半无限规划建模状态相关不确定性的鲁棒最优控制，直接相关。|
|E1469|15_Gimelfarb_2024_Constraint-Generation Policy Optimization (CGPO)_ Nonlinear.pdf|P-fae2e5521424b24c|adjacent|研究混合MDP策略优化中的非线性规划constraint-generation；可借鉴切约束实现，但模型不是半无限规划。|
|E1470|17_Hoa_2020_Asymptotic behavior of Integer Programming and the stability.pdf|P-54cd3c43f77bc518|irrelevant|研究整数规划与Castelnuovo–Mumford regularity的渐近关系；无半无限约束或切生成。|
|E1471|07_Oustry_2023_Minimal-time nonlinear control via semi-infinite programming.pdf|P-a4c1aa39cd6f2829|direct|以半无限规划层级逼近非线性最短时间控制，属于直接应用及算法建模。|
|E1472|36_Canovas_2013_Calmness modulus of linear semi-infinite programs.pdf|P-28bf059385cfa185|direct|研究线性半无限规划最优解映射的calmness及稳定性，直接属于SIP性质。|
|E1473|22_Xu_2019_On Solving a Class of Linear Semi-Infinite Programs by the T.pdf|P-fb822fc3ab40cf73|direct|利用三角矩问题/半定规划求解一类线性半无限规划，直接相关。|
|E1474|05_Guo_2020_On Solving a Class of Fractional Semi-infinite Polynomial Pr.pdf|P-9b209dd1268468de|direct|研究分式半无限多项式规划、对偶锥重构和SOS/SDP松弛，直接相关。|
|E1475|02_Terrien_2026_An iterative Constraint Programming approach to integrate ma.pdf|P-c5002886017f314b|irrelevant|研究带最大工作量约束的抢占式jobshop排程constraint-programming；constraint programming不表示半无限规划/约束生成。|
|E1476|31_Papp_2015_Semi-infinite programming using high-degree polynomial inter.pdf|P-eaa51628b8b3a5a3|direct|以高次多项式插值和SDP表示半无限约束，直接对应约束半无限化处理。|
|E1477|06_Oustry_2023_Convex semi-infinite programming algorithms with inexact sep.pdf|P-6823eca2545a7dbf|direct|研究凸半无限规划在不精确分离预言机下的算法与收敛性，直接命中约束生成/分离步骤。|
|E1478|04_Huertas_2024_Parallel Batch Scheduling With Incompatible Job Families Via.pdf|P-1802c52b23a7dc81|irrelevant|研究批次排程的constraint-programming模型；与半无限规划的约束生成算法无关。|
|E1479|27_Wei_2018_An Inexact Primal-Dual Algorithm for Semi-Infinite Programmi.pdf|P-9a632a388e489afb|direct|提出带误差界和非负测度近端更新的半无限规划inexact primal-dual算法，直接相关。|
|E1480|01_Pantuso_2025_The L-Shaped Method for Stochastic Programs with Decision-De.pdf|P-c87123ba8b3ccbda|adjacent|研究决策相关不确定性的两阶段随机规划L-shaped cuts；cut generation相邻，但不是半无限模型。|
|E1481|12_Aravind_2024_Distributed alternating gradient descent for convex semi-inf.pdf|P-636829975d80f85d|direct|研究凸半无限规划的分布式交替梯度算法并给出收敛保证。|
|E1482|33_Mehrotra_2013_A cutting surface algorithm for semi-infinite convex program.pdf|P-d6b57dfc8a330d19|direct|分析凸半无限规划的central cutting-surface算法并用于矩鲁棒优化，直接相关。|
|E1483|24_Wei_2018_Corporative Stochastic Approximation with Random Constraint.pdf|P-bb1c8c906cbcfa9c|direct|针对半无限规划cut-generation子问题使用随机采样和随机近似，给出误差界与收敛率。|
|E1484|39_CÁnovas_2011_Quantitative Stability and Optimality Conditions in Convex S.pdf|P-6d1aed887aa3e5be|direct|研究凸半无限/无限约束系统的定量稳定性与最优性条件，直接覆盖SIP理论。|

## 实际阅读卡

- [P-3580cf604bdcbfd2：full_read](papers/TEAM_AB-P-3580cf604bdcbfd2.json)
- [P-646526af87253fe6：full_read](papers/TEAM_AB-P-646526af87253fe6.json)
- [P-ecb91aed2e290a23：full_read](papers/AD-P-ecb91aed2e290a23.json)

## 本地原文定位

原附件PDF：工作包根目录的 `evidence/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`evidence/supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。
