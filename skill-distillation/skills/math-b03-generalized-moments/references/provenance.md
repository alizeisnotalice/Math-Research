# B03 来源基线清单（旧标签待复核）

> 当前筛查说明：下方标题、作者/版本、direct/adjacent/irrelevant及理由均是基线记录或旧人工备注，不代表本轮已复核；本轮状态以交付根 `audit_current/pdf_screening_and_reading_queue.csv` 和共享读卡为准。列表中的专题附件数是路径/条目口径，不可直接当SHA去重的论文数。题录TXT不是全文证据。


附件归属 41 条；本专题关联不同PDF 38 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E0983|05_Gamertsfelder_2025_The Effective Countable Generalized Moment Problem.pdf|P-858d0c6799128b32|direct|全文摘要明确为带可数矩约束的广义矩问题(GMP)，研究Moment-SoS松弛收敛率及唯一优化测度收敛；是广义矩优化的直接数值理论。|
|E0984|21_Wu_2022_A Non-Gaussian Bayesian Filter Using Power and Generalized L.pdf|P-05eb64cf961db17d|adjacent|研究非高斯Bayesian filter，以power/generalized-logarithmic moments构造密度近似并证明参数映射性质；使用矩统计但不是在矩歧义集上求最坏分布界。|
|E0985|39_Bayraktar_2007_On the One-Dimensional Optimal Switching Problem.pdf|P-c2c9b1b15ecc5674|irrelevant|研究一维扩散过程的最优切换/停止、动态规划与价值函数；标题/摘要没有广义矩或矩约束分布优化，切换意义的“最优”不属于矩问题。|
|E0986|04_Huang_2022_Generalized truncated moment problems with unbounded sets.pdf|P-75a1e84b124307ee|direct|研究非紧支撑上的广义截断矩问题、矩锥及其SOS近似，并给出可行测度/有限原子表示；与测度矩约束直接相关。|
|E0987|01_Narayanan_2020_Moment-Based Ensemble Control.pdf|P-71440609e6e30463|adjacent|把Hausdorff矩序列用于集合体控制并研究moment-system的可控性；矩问题是表征控制系统的工具，不是分布歧义下的锐界优化。|
|E0988|23_Zhu_2020_Worst-Case Risk Quantification under Distributional Ambiguit.pdf|P-21b418417150d856|direct|明确以kernel mean embedding定义矩约束歧义集，用广义矩优化最大化风险/约束违背概率；直接对应最坏分布概率界。|
|E0989|40_Yannick_2022_Option-Implied Moments_ A Generalized Moment Problem Approac【题录】.txt|R-9c75a57dbe1ff058|direct|题名即“Option-Implied Moments: A Generalized Moment Problem Approach”，可确认属于广义矩方法；但仅有无摘要题录且全文未取得，结果与假设仍待核。|
|E0990|25_Dio_2018_The multidimensional truncated Moment Problem_ The Moment Co.pdf|P-dbf5b46fca0494a4|direct|研究多维截断矩锥的面结构、Carathéodory原子数及最大质量问题，提供矩可行性和有限原子归约的直接理论基础。|
|E0991|19_Wu_2023_General Discrete-Time Fokker-Planck Control by Power Moments.pdf|P-3c1dffa1064f978f|adjacent|用power moments建模离散时间Fokker–Planck分布控制并求控制实现；目标是分布动态控制，不是对满足矩约束的所有分布取最坏值。|
|E0992|14_Kirschner_2021_Convergence rates of RLT and Lasserre-type hierarchies for t.pdf|P-9a1b62ec2ad3c819|direct|直接研究simplex/sphere上的GMP之RLT与Lasserre层级，给出明确O(1/r)、O(1/r²)收敛率；属于广义矩界的计算理论。|
|E0993|10_Mula_2022_Moment-SoS Methods for Optimal Transport Problems.pdf|P-0fe0e3690c911ee7|direct|把多边最优传输写成GMP，对边际矩做Moment-SoS层级并证明收敛；目标是耦合测度的线性泛函，方法可迁移到分布矩优化。|
|E0994|30_Hieu_2016_Parallel hybrid methods for generalized equilibrium problems.pdf|P-f0bf2e545e1e005f|irrelevant|“generalized”修饰均衡问题；摘要讨论伪收缩映射及广义均衡算法，没有概率矩、测度矩或矩锥。|
|E0995|17_Hu_2024_Positivstellensätze and Moment problems with Universal Quant.pdf|P-af82c7b6d22f9a8e|direct|研究带全称量词的Positivstellensatz与moment问题，核心是矩约束/非负证书之间的对偶；对构造广义矩上下界的证书直接相关。|
|E0996|34_Karlsson_2015_The Multidimensional Moment Problem with Complexity Constrai.pdf|P-cf70ceb9362e7615|direct|研究带复杂度约束的多维矩问题、正多项式对偶锥及可行表示；属多维矩锥的直接理论。|
|E0997|06_Guo_2022_A Unified Framework for Generalized Moment Problems_ a Novel.pdf|P-75068ca1d23232fc|direct|摘要和题名均以Generalized Moment Problems为主题，提出适用一般可测矩函数的统一primal-dual框架和最优性条件，属于核心理论来源。|
|E0998|26_Klerk_2018_A survey of semidefinite programming approaches to the gener.pdf|P-9e49e5c1d29f474f|direct|综述GMP的半定规划方法，覆盖测度问题、矩松弛和SOS对偶证书；直接提供方法综述。|
|E0999|11_Henrion_2023_Algebraic certificates for the truncated moment problem.pdf|P-8ae2d5ff87c54133|direct|研究半代数支撑下截断矩可表示性与代数证书；与矩锥/正多项式分离直接相关。|
|E1000|09_Korda_2022_Exploiting ideal-sparsity in the generalized moment problem.pdf|P-d2ee235cf7790ae9|direct|针对GMP提出ideal-sparsity分解与层级方法，优化变量仍为支撑受约束的测度；直接相关。|
|E1001|12_Huang_2024_Moment-SOS relaxations for moment and tensor recovery proble.pdf|P-c7080e7b60af9ef8|adjacent|Moment-SOS用于矩/tensor recovery并研究张量分解；可迁移稀疏矩松弛，但目标是恢复张量而不是概率分布的最坏界。|
|E1002|15_Wang_2024_On stochastic control problems with higher-order moments.pdf|P-969c8d817643b4fa|adjacent|研究含高阶中心矩目标的时间不一致随机控制及均衡控制；矩出现于控制目标/动态方程，并非广义矩歧义优化。|
|E1003|25_Henrion_2023_Algebraic certificates for the truncated moment problem.pdf|P-8ae2d5ff87c54133|direct|与E0999是同一paper_id的重复条目；同篇研究半代数截断矩问题的代数证书，只需后续深读一次。|
|E1004|03_Sultana_2023_Variational Reformulation of Generalized Nash Equilibrium Pr.pdf|P-2af7eb18048dea36|irrelevant|研究非序偏好下的广义纳什均衡变分重构；“generalized”指均衡博弈，非广义矩问题。|
|E1005|27_Zhu_2018_On the Uniqueness Result of Theorem 6 in _Relative Entropy a.pdf|P-343e568df3f7547a|adjacent|以广义矩约束谱估计为背景，分析相对熵/矩映射唯一性并给出反例；可提醒对偶正则性与唯一性风险，但不直接求概率上/下界。|
|E1006|18_Henrion_2026_Extreme points and faces in the moment problem.pdf|P-b77447cac7b06cc4|direct|研究广义矩问题可行测度集合的极点与面，给出有限原子极点及将积分下确界限制到极点的结果；直接对应极端分布/有限支撑方法。|
|E1007|29_Georgiou_2016_Likelihood Analysis of Power Spectra and Generalized Moment.pdf|P-ca0893e2ca4a464c|adjacent|研究功率谱似然与广义矩谱估计之间的映射/参数化；核心为谱密度估计及结构唯一性，和矩约束最坏分布界仅相邻。|
|E1008|02_Huang_2024_Piecewise SOS-Convex Moment Optimization and Applications vi.pdf|P-ec78cd78268be31b|direct|研究分段SOS凸矩优化与应用，目标是利用矩/SOS证书求全局最优；与广义矩优化的计算性方法直接相关。|
|E1009|07_Halaseh_2026_Duality attainment and strict feasibility of the generalized.pdf|P-c5b0ca086449f4b8|direct|研究GMP及其松弛的对偶可达性和严格可行条件，直接支持何时矩对偶界能达到/强对偶成立。|
|E1010|33_Mohammadi_2016_Combining SOS and Moment Relaxations with Branch and Bound t.pdf|P-deb72842d1067b4f|adjacent|把SOS与moment松弛结合分支定界求多项式全局优化；使用矩层级但不是分布歧义集上的锐概率/期望界。|
|E1011|31_Molzahn_2016_Moment Relaxations of Optimal Power Flow Problems_ Beyond th.pdf|P-ab3e7532bb324b7b|adjacent|moment relaxations用于非凸optimal power flow及网络稀疏化；是特定工程应用，不是一般矩界理论。|
|E1012|13_Eekelen_2023_A generalized moment approach to sharp bounds for conditiona.pdf|P-91433be844e70d18|direct|提出以GMP与锥对偶得到已知条件事件下conditional expectation的sharp bounds，并扩展到带side information的DRO；本组已对该文完整通读。|
|E1013|36_Molzahn_2014_Sparsity-Exploiting Moment-Based Relaxations of the Optimal.pdf|P-67b88482d416115f|adjacent|稀疏moment-based relaxations用于OPF问题；属于矩优化应用，未给一般概率分布尾界框架。|
|E1014|28_Chen_2017_The discrete moment problem with nonconvex shape constraints.pdf|P-6f0d6890eb9f2625|direct|研究含非凸形状约束的离散矩问题，表征最坏/极值分布并给出有限结构算法；和用矩与分布形状锐化界直接相关。|
|E1015|38_Bandeira_2013_Approximating the Little Grothendieck Problem over the Ortho.pdf|P-be84558b32a777c1|irrelevant|Little Grothendieck问题在正交群上的近似与舍入；未出现广义概率矩优化，和“moment”专题无关。|
|E1016|24_Zhu_2023_A Weaker Regularity Condition for the Multidimensional $ν$-M.pdf|P-06dbd1e1b3e57478|adjacent|研究多维ν-moment谱估计问题的正则性条件，目标在谱参数化/似然映射；不是测度矩锥上的最坏界。|
|E1017|22_Sklyar_2020_Truncated Hausdorff Moment Problem and Analytic Solution of.pdf|P-c8906ecb0fe43090|adjacent|确实研究截断Hausdorff矩问题，但其应用主线是时间最优控制的解析解；可迁移有限矩可行性思路，非分布鲁棒概率界。|
|E1018|16_Henrion_2024_Solving moment and polynomial optimization problems on Sobol.pdf|P-947cbaf82e4e3078|direct|研究Sobolev球上的矩与多项式优化，扩展有限维矩/SOS优化到函数空间；测度矩优化仍是主题。|
|E1019|20_Zhang_2020_Safe Screening Rules for Generalized Double Sparsity Learnin【题录】.txt|R-cba2c6567fda5b48|irrelevant|标题/摘要中的generalized修饰double-sparsity learning；变量选择和混合正则化问题，没有概率/测度矩约束。|
|E1020|37_Mehrotra_2013_A cutting surface algorithm for semi-infinite convex program.pdf|P-d6b57dfc8a330d19|direct|研究半无限凸规划cutting-surface算法并直接应用于给定任意阶/非多项式矩的DRO；是GMP概率界的优化算法来源。|
|E1021|08_Huang_2025_An ideal-sparse generalized moment problem reformulation for.pdf|P-2e7ec0a57bbbceca|direct|将完全正张量分解精确重构为ideal-sparse generalized moment problem并分析Moment hierarchy收敛；虽以张量为应用，核心确为GMP建模。|
|E1022|32_Molzahn_2016_Computational Analysis of Sparsity-Exploiting Moment Relaxat.pdf|P-8390636c956069de|adjacent|Moment-SOS稀疏层级用于电力流问题，主题在OPF计算与全局解；属于矩优化工程应用。|
|E1023|35_Molzahn_2015_Solution of Optimal Power Flow Problems using Moment Relaxat.pdf|P-8fe67628583bc153|adjacent|用Lasserre moment relaxations与目标惩罚求解OPF；不是一般分布最坏界，作为应用文献邻近。|

## 实际阅读卡

- [P-91433be844e70d18：full_read](papers/AD-P-91433be844e70d18.json)
- [P-deb72842d1067b4f：partial_read](papers/AD-P-deb72842d1067b4f.json)

## 本地原文定位

原附件PDF：工作包根目录的 `evidence/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`evidence/supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。


## 本轮复核状态（2026-10-07）

- **E0997 / P-75068ca1d23232fc**：PDF 身份/版本按首页核验为 Guo, He, Jiang, Wang，arXiv:2201.01445v3；本地 PDF SHA-256 与 56 页逐页卡一致，主文、e-companion 附录和参考文献均记录覆盖。具体相关性为 B03 direct。仅将 §2 Theorem 1 与三步候选—验证法登记为条件性证据备注；`SKILL.md` 默认操作流程尚未启用此路线，定理的独立复核仍待完成。卡片中具体定理条件及边界见 [E0997 逐页卡](papers/AD-P-75068ca1d23232fc.json)。外部出版元数据未独立核验。
