# C01 来源与实际处理状态

附件归属 41 条；本专题关联不同PDF 41 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E0486|27_Louca_2014_Nondegeneracy and Inexactness of Semidefinite Relaxations of.pdf|P-1d50f53dd9e0c24a|direct|OPF非凸QCQP的SDP松弛精确性与rank-one解存在条件；提供明确的tightness/非精确判据。|
|E0487|33_Wang_2013_Semidefinite relaxations for semi-infinite polynomial progra.pdf|P-ecb91aed2e290a23|adjacent|半无限多项式规划用SDP松弛与交换算法，针对无限约束索引集而非一般lifting接口。|
|E0488|23_Shah_2016_Biconvex Relaxation for Semidefinite Programming in Computer.pdf|P-8acd1d26c6add036|adjacent|把大型SDP重写为低维双凸近似并交替求解；重点是数值算法，且结果为近似。|
|E0489|08_Burer_2023_A Slightly Lifted Convex Relaxation for Nonconvex Quadratic.pdf|P-783dff58d54a2ee6|direct|非凸球交QCQP的辅助变量β lifting与单矩阵Beta松弛；对两球给出精确性定理，已全文读。|
|E0490|17_Eltved_2018_On the Robustness and Scalability of Semidefinite Relaxation.pdf|P-4ff5b4866621bdf0|adjacent|大规模OPF SDP松弛求解器的可扩展性、稳健性和数值gap比较；主要是数值实现。|
|E0491|03_Cory-Wright_2026_Compact Lifted Relaxations for Low-Rank Optimization.pdf|P-d6da1bb3a2238f8d|direct|任意结构低秩二次优化的紧凑半定lifting与冗余矩阵块消除；直接研究低秩二次优化的紧凑提升凸松弛。|
|E0492|10_Holmes_2023_On Semidefinite Relaxations for Matrix-Weighted State-Estima.pdf|P-d427c4b0a7c06e3e|adjacent|矩阵加权定位估计的SDP松弛tightness与误差/噪声限制；应用专用。|
|E0493|34_Xia_2013_A new semidefinite relaxation for $_ell_{1}$-constrained qua.pdf|P-e45748c9e7fcaa7b|direct|ℓ1球上的非凸二次优化变量拆分与SDP松弛，并比较包含关系及扩展。|
|E0494|16_Khoo_2019_Semidefinite relaxation of multi-marginal optimal transport.pdf|P-f8a87989c82417ca|adjacent|多边际最优输运量子模型的SDP凸松弛和弱Slater下对偶达到；领域专用但可用作锥松弛例。|
|E0495|30_Lee_2014_Lower bounds on the size of semidefinite programming relaxat.pdf|P-df6ea7e760c17f43|direct|研究多项式/可扩展SDP lift的最小维度下界及extension complexity，直接对应提升规模成本。|
|E0496|18_Heimendahl_2022_A semidefinite program for least distortion embeddings of fl.pdf|P-9b689f486c804a6f|adjacent|平坦环面嵌入问题的无限维SDP；属于半定表示但对象和无限维分析特殊。|
|E0497|37_Paparella_2012_A note on the Lovasz-Schrijver Semidefinite Programming Rela.pdf|P-f9aa8ab78f03a701|adjacent|二元整数规划的Lovasz–Schrijver半定lift-and-project松弛；专门层级应用。|
|E0498|35_Nie_2013_Semidefinite Relaxations for Best Rank-1 Tensor Approximatio.pdf|P-430da7e9916c90b1|direct|秩一张量近似转为球面多项式优化并构造SOS/SDP松弛，研究近似保证。|
|E0499|09_Guedes-Ayala_2024_Sparse Sub-gaussian Random Projections for Semidefinite Prog.pdf|P-76fc5e5ff59f9005|adjacent|随机投影压缩SDP变量并给近似误差界；重点是投影算法，不是精确矩阵提升。|
|E0500|01_Xu_2023_New semidefinite relaxations for a class of complex quadrati.pdf|P-fc57ada29665bec4|direct|复数非凸QCQP中构造新的lift矩阵有效约束，比较松弛强度；直接支持复数lift扩展。|
|E0501|18_Filho_2018_On the integrality gap of the maximum-cut semidefinite progr.pdf|P-6559b62bf4af5e0b|adjacent|最大割SDP松弛的integrality gap factor-revealing优化；研究松弛性能界而非构造接口。|
|E0502|36_Lasserre_2012_A Lagrangian relaxation view of linear and semidefinite hier.pdf|P-6d53b70110481e4c|direct|将拉格朗日放松、LP/RLT与SOS/SDP层级关联，给出多项式优化凸松弛构造。|
|E0503|32_Fawzi_2013_Equivariant semidefinite lifts and sum-of-squares hierarchie.pdf|P-83fc85ecba810ce3|direct|多面体的等变PSD lifts和SOS层级，定义lift为PSD锥仿射切片投影，分析lift size。|
|E0504|12_Kanoh_2020_Centering ADMM for the Semidefinite Relaxation of the QAP.pdf|P-11ef2b6538707606|adjacent|QAP的SDP松弛上提出Centering ADMM和收敛分析；关注数值算法。|
|E0505|26_Fawzi_2014_Equivariant semidefinite lifts of regular polygons.pdf|P-d282095d5c0595a9|direct|规则多边形的等变PSD lifts与Lasserre/SOS阶数，直接分析半定扩展表示和规模。|
|E0506|24_Fawzi_2015_Lieb's concavity theorem, matrix geometric means, and semide.pdf|P-cb8e06d9cb6a933c|adjacent|Lieb矩阵函数的显式SDP表示及矩阵几何均值；提供锥表示技术，非原变量非凸QCQP典型lift。|
|E0507|40_Klerk_2009_On semidefinite programming relaxations of the traveling sal.pdf|P-f0a49be48d514d20|adjacent|TSP的QAP导出SDP松弛和界强度比较；组合优化特定。|
|E0508|28_Sevcovic_2014_Solution to the Inverse Wulff Problem by Means of the Enhanc.pdf|P-7f19c1a7fb3e0b92|adjacent|逆Wulff各向异性函数问题经增强SDP松弛与二次线性约束求解；是非凸松弛的应用案例。|
|E0509|11_Hu_2024_Affine Facial Reduction for Semidefinite Relaxations of Bina.pdf|P-d689fd295109b74a|direct|二元/混合二元SDP松弛的affine facial reduction预处理，处理Slater失败和lift大矩阵。|
|E0510|39_Bachoc_2010_Invariant semidefinite programs.pdf|P-65e42f19d2b84512|direct|对称不变SDP的结构约化和表示；一般对称性降维lift方法。|
|E0511|22_Paredes_2017_A Low-Rank Rounding Heuristic for Semidefinite Relaxation of.pdf|P-0e4758bbe823504f|adjacent|水电机组调度SDP松弛的low-rank rounding启发式；回构候选而非证明松弛精确。|
|E0512|04_Roux_2023_Instance-specific linear relaxations of semidefinite optimiz.pdf|P-6dc834a2467f2aa6|direct|基于目标/约束矩阵可交换性构造带保证的实例专属LP近似SDP，并分析两者同值条件。|
|E0513|07_Chen_2024_Exactness Conditions for Semidefinite Relaxations of the Qua.pdf|P-e2fe503f77e0d383|direct|QAP多种SDP松弛的输入矩阵精确性充分条件；明确界定可无gap使用的场景。|
|E0514|19_Gutekunst_2017_The Unbounded Integrality Gap of a Semidefinite Relaxation o.pdf|P-f2e3423894f2cf30|adjacent|TSP/k-cycle-cover SDP松弛的无界integrality gap；是重要反例/失效证据。|
|E0515|38_Helton_2011_Semidefinite programming in matrix unknowns which are dimens.pdf|P-5d6dc4f622d558cb|adjacent|矩阵未知量和非交换多项式的dimension-free SDP综述；抽象对象与有限维实变量提升不同。|
|E0516|25_Sevcovic_2014_Application of the Enhanced Semidefinite Relaxation Method t.pdf|P-8afb30223b637d67|adjacent|用增强SDP松弛求最优各向异性函数并给充分精确条件；专门应用。|
|E0517|05_Yang_2025_A New Semidefinite Relaxation for Linear and Piecewise-Affin.pdf|P-a73a127fa16f8f55|adjacent|时标变量使线性/分段仿射最优控制非凸，构造二阶SDP松弛；面向控制场景。|
|E0518|31_Andersen_2013_Reduced-Complexity Semidefinite Relaxations of Optimal Power.pdf|P-fe54da374477782f|adjacent|通过弦图转换删约束得到更便宜的OPF SDP松弛；关注稀疏计算和精度。|
|E0519|13_Yuan_2021_Semidefinite Relaxations of Products of Nonnegative Forms on.pdf|P-f3b991cc0787e7d4|direct|非负二次型几何平均最大化用紧凑SDP/SOS relaxation，并给常数因子保证。|
|E0520|15_Zhang_2019_Calculating Entanglement Eigenvalues for Non-Symmetric Quant.pdf|P-39f841345a095ad5|direct|复张量主特征值转为等式约束多项式优化并用Jacobian SDP层级松弛。|
|E0521|29_Huang_2014_Scalable Semidefinite Relaxation for Maximum A Posterior Est.pdf|P-3eb1975583042c90|adjacent|离散MRF MAP估计的可扩展SDP松弛；组合推理应用，不主张通用tightness。|
|E0522|20_Cifuentes_2017_On the local stability of semidefinite relaxations.pdf|P-5600afe0464247f4|direct|研究QCQP在名义点精确后对参数扰动的局部tightness稳定性及定量界。|
|E0523|14_Chevalier_2026_Low-Rank and Lifted Semidefinite Programming for Mixed-Integ.pdf|P-b0204c99397d2551|direct|混合整数多项式电网优化采用低秩SDP与Lasserre矩lifting迭代tightening。|
|E0524|21_Paredes_2017_Inexactness of the Hydro-Thermal Coordination Semidefinite R.pdf|P-b30716d325a37d73|adjacent|水火电协调问题特定Shor松弛的非精确性条件及修正方案；提供失效案例。|
|E0525|02_Au_2020_On Connections Between Association Schemes and Analyses of P.pdf|P-e3411dc0f93e9a88|direct|Association schemes与Lovasz–Schrijver半定lift-and-project层级，研究rank和图优化界。|
|E0526|06_Sinjorgo_2024_Cuts and semidefinite liftings for the complex cut polytope.pdf|P-9e5c272e6deccbe0|direct|复cut多面体的有效切与更小型第二半定lifting，并证明与同阶复Lasserre lift等价。|

## 实际阅读卡

- [P-783dff58d54a2ee6：full_read](papers/AD-P-783dff58d54a2ee6.json)
- [P-ecb91aed2e290a23：full_read](papers/AD-P-ecb91aed2e290a23.json)

## 本地原文定位

原附件PDF：工作包根目录的 `evidence/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`evidence/supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。
