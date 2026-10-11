# J01 来源与实际处理状态

附件归属 39 条；本专题关联不同PDF 37 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

> **E0785 来源警告（P-86e3b72f；精确源 PDF SHA-256 `86e3b72f3b605def55817230693c22b778a410721cab800e95065b6a4c27617f`）** 该文 Theorem 1 声称的 weighted-f / total-variation 几何遍历结论在印刷假设下被独立反例否定：取 λ₀=1、μ(λ)=1−λ、B=0、ν=Unif(0,1)、f(λ)=1+λ²，则假设与 Af≤2−f 成立，但 λ>1 时有限时刻转移 Dirac 分布不趋于 δ₁。反例只否定该几何遍历结论；不否定显示的生成元公式、弱收敛或例中 δ₁ 的不变性/唯一性。当前 J01 SKILL、method、claims 未把该结论用作操作依赖。证据：共享全文卡 `audit_current/shared-reading/revise_e/86e3b72f3b605def55817230693c22b778a410721cab800e95065b6a4c27617f.json`（当前 SHA `e7eef8e3f63d5d7f8503c561b95fbf1f3a2cc54c3a64b024905071ba503f52b4`；其窄审回执绑定前一卡SHA `120f199b5415bdfa8c431dcf62f671e1af249c43a6ecf99bbb96489948a5f836`），独立JK窄审回执 `audit_current/independent-jk/source_review_P-86e3b72f_theorem1_geometric_ergodicity_20261008_r2.json`（SHA `e4ed7cb3542db1243dbfd2f6888c59ecede64283dff517aabd9f154d6c999528`；只复核关键页与反例，不是第二次全文阅读）。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E0785|13_Eyjolfsson_2021_Multivariate self-exciting jump processes with applications.pdf|P-86e3b72f3b605def|adjacent（22页全文；Theorem 1警告）|仅邻近：模型特定extended generator不构成一般J01生成元定理。**来源警告**：独立JK回执发现Theorem 1按印刷Assumption 1推出的weighted-f/TV几何遍历结论有反例（缺少可达/φ-不可约性）；只否定该遍历结论，不否定显示生成元、弱收敛或例中δ1不变分布。当前SKILL/method/claims均未引用该遍历结论。完整卡SHA `e7eef8e3f63d5d7f8503c561b95fbf1f3a2cc54c3a64b024905071ba503f52b4`；窄审回执SHA `e4ed7cb3542db1243dbfd2f6888c59ecede64283dff517aabd9f154d6c999528`。|
|E0786|04_Feinberg_2021_Kolmogorov's Equations for Jump Markov Processes and their A.pdf|P-ea2d2ca36be09621|direct（全文复核）|18页全文已核；PDF p.4 Assumption 2.3，pp.6–7 Theorem 3.5/Example 3.2：前向式只用于 (q,s)-bounded B 和 a.e. 终时，regular 才唯一，B=X 反例右端未定义。通用定理外引 [9] 未独立证明，见 atomic claims。|
|E0787|39_Mingang_2020_Quantized $$H_{_infty }$$ Filtering for Continuous-Time Nonh【题录】.txt|R-e1cbb50828e5d3df|pending|题录已读，标题涉及H∞ filtering/state estimation非齐次Markov jump systems；摘要和全文缺失，不仅凭题名排除或认证。|
|E0788|28_Ratanov_2013_Damped jump-telegraph processes.pdf|P-5ca861a4e5706dbc|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0789|02_Macrina_2021_Captive Jump Processes.pdf|P-61e8f1cc14624bef|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0790|29_Zeng_2013_How Vertex reinforced jump process arises naturally.pdf|P-9dca3131d0950261|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0791|40_Yufeng_2022_New Approach to $$H_{_infty }$$ State Estimation for Continu【题录】.txt|R-1b295980cfa76267|pending|题录已读，标题涉及H∞ filtering/state estimation非齐次Markov jump systems；摘要和全文缺失，不仅凭题名排除或认证。|
|E0792|08_Sabot_2021_The _-Vertex-Reinforced Jump Process.pdf|P-350a6f3d1a4ab464|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0793|36_Burdzy_2010_Stationary distributions for jump processes with memory.pdf|P-14e61c3e62f5ead6|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0794|35_Mischaikow_2010_Topology-guided sampling of nonhomogeneous random processes.pdf|P-b1a1cc0fe6978c26|irrelevant|随机函数nodal/excursion组件离散采样，非跳跃生成元或Q-function构造。|
|E0795|31_Bogdan_2012_Boundary Harnack inequality for Markov processes with jumps.pdf|P-5c064f0f0892f54d|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0796|33_Mukhamedov_2011_On $L_1$-Weak Ergodicity of nonhomogeneous discrete Markov p.pdf|P-0b9e2f2a28245e42|irrelevant|离散时间Markov链L1弱遍历，非连续时间跳率/生成元。|
|E0797|18_Miles_2017_Jump Locations of Jump-Diffusion Processes with State-Depend.pdf|P-ee636b8a2c2329ad|direct|全文22页：状态相关hazard的 killed-between-jump方程、跳位置resolvent迭代；原文索引/符号/λ=0/推前与附录矩公式缺口逐项保留。没有冻结坐标补偿契约。|
|E0798|24_Yang_2015_Multifractality of jump diffusion processes.pdf|P-8e60655fe6cbdf80|irrelevant|跳扩散多分形路径正则性；非所需时变Q-function或首跳构造。|
|E0799|15_Rabehasaina_2019_Multitype branching process with nonhomogeneous Poisson and.pdf|P-2c58bb2b0a01fae5|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0800|09_Sikorski_2020_The Augmented Jump Chain -- a sparse representation of time-.pdf|P-674e36bed18cb66c|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0801|26_Burch_2014_The exit-time problem for a Markov jump process.pdf|P-3ca42b9ac0084b5b|adjacent|复用HK已读P-3ca42b9ac0084b5b：非局部有限跳程退出模型的生成元/forward-kernel构件；不提供一般时变Q-function。|
|E0802|16_Guo_2019_Periodic solutions of hybrid jump diffusion processes.pdf|P-01fcd4a98b04dfba|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0803|03_Stettner_2025_Long run control of nonhomogeneous Markov processes.pdf|P-b67984eba11757b2|irrelevant|离散时间受控Markov转移与long-run Bellman方程，非连续时间跳强度。|
|E0804|38_Bayraktar_2007_Analysis of the optimal exercise boundary of American option.pdf|P-c7e3f4012a2fbb57|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0805|37_Barral_2009_A pure jump Markov process with a random singularity spectru.pdf|P-008e733184819c44|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0806|25_Oertel_2015_On Jump Measures of Optional Processes with Regulated Trajec.pdf|P-2f81500e5dfa70c6|irrelevant|首页摘要和全文核验显示本文构造跳测度，不推出非齐次跳过程的真实生成元；无时间/状态依赖生成元定理。|
|E0807|23_Sarantsev_2015_Explicit Rates of Exponential Convergence for Reflected Jump.pdf|P-cd94492f5c41c046|irrelevant|反射跳扩散的平稳分布指数收敛率，非目标非齐次生成元构造。|
|E0808|30_Barczyk_2012_Scaling limits of coupled continuous time random walks and r.pdf|P-548cf76ea89ebd8f|irrelevant|耦合CTRW的缩放极限与residual order statistics，非Q-function/补偿接口。|
|E0809|32_Mikulevicius_2011_On the rate of convergence of simple and jump-adapted weak E.pdf|P-8b935a3036f27915|irrelevant|Lévy驱动SDE的weak Euler数值误差，非时变跳强度构造。|
|E0810|22_Zhang_2015_On the nonexplosion and explosion for nonhomogeneous Markov.pdf|P-500a71307a38c036|direct-but-scoped|精确SHA全读卡复用，28/28页含附录和参考文献；只支持稳定保守非齐次 Q 模型的非爆炸/漂移判据及对应条件，不支持任意测试函数生成元域或一般Kolmogorov前向式。|
|E0811|07_Kohatsu-Higa_2021_Density estimates for jump diffusion processes.pdf|P-3556faeb8af5a8b6|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0812|05_Liu_2021_On the comparison between jump processes and subordinated di.pdf|P-abe34d350b944d66|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0813|27_Feinberg_2013_On solutions of Kolmogorov's equations for jump Markov proce.pdf|P-64a7ecbfc2a44f93|direct|正文对可测、真实时变 Q-function 构造多跳 Markov 过程及补偿随机测度，证明首跳递归、backward/forward 方程、极小性；符合非齐次 jump generator 主题。|
|E0814|12_Larsson_2024_Inverting the Markovian projection for pure jump processes.pdf|P-8c2298dab4e50832|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0815|01_Disertori_2023_Transience of vertex-reinforced jump processes with long-ran.pdf|P-0baa94015b19fb47|irrelevant|顶点强化跳过程的transience/random environment，非此项确定时变Q-function。|
|E0816|11_Bandini_2020_Stochastic filtering of a pure jump process with predictable.pdf|P-feb30426b68c574f|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0817|06_Gérard_2023_A multi-dimensional version of Lamperti's relation and the M.pdf|P-cf022f495a45ff81|irrelevant|相互作用Brownian与Bessel桥的Lamperti/Matsumoto–Yor关系，非跳跃生成元。|
|E0818|20_Budiana_2016_Flowgraph Models and Analysis for Markov Jump Processes.pdf|P-fd12b6570711afdd|irrelevant|齐次有限链flowgraph等待时间mgf计算；本J01要求时变生成元，此项不提供。|
|E0819|21_Martin_2016_Efficient posterior inference on the volatility of a jump di.pdf|P-9da14ac3fc6ae10a|irrelevant|基于离散样本的volatility Bayesian后验推断，非跳生成元/补偿。|
|E0820|10_Sarantsev_2020_Convergence rate to equilibrium in Wasserstein distance for.pdf|P-33e40ef1d76ecd7a|irrelevant|反射跳扩散Wasserstein平稳收敛率，非时变Q-function构造。|
|E0821|17_Chen_2018_Martingale solutions for the three-dimensional stochastic no.pdf|P-020260a5fca4f4df|irrelevant|三维随机非齐次Navier–Stokes鞅解，“非齐次”指流体密度，非目标跳过程生成元。|
|E0822|34_Whitehead_2011_Occupation Times for Jump Processes.pdf|P-f0f533eb54a0dcde|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|
|E0823|14_Kagan_2025_Averaging principle for jump processes depending on fast erg.pdf|P-caf84402ee14b035|pending|摘要与状态相关跳强度/生成元、非爆炸、非局部过程或首跳构造邻近；需要主定理及全文核验才能确定可迁移范围。|

## 实际阅读卡

- [P-2f81500e5dfa70c6：full_read](papers/HK-P-2f81500e5dfa70c6.json)
- [P-3ca42b9ac0084b5b：full_read](papers/HK-P-3ca42b9ac0084b5b.json)
- [P-64a7ecbfc2a44f93：full_read](papers/HK-P-64a7ecbfc2a44f93.json)
- [P-ea2d2ca36be09621：full_read](papers/HK-P-ea2d2ca36be09621.json)
- [P-ee636b8a2c2329ad：full_read](papers/TEAM_JK-P-ee636b8a2c2329ad.json)

## 本地原文定位

原附件PDF：工作包根目录的 `corpus/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。


## E0810：Zhang exact-SHA full-read integration

- arXiv:1511.05011v2; PDF SHA `500a71307a38c036d97bca3c6b958c710c22abd9ca1479c201a1ed17713b424f`; 28/28 pages, including Appendix A pp.23–26 and references pp.27–28. Exact full-read card: [independent audit card](../../audit_current/shared-reading/audit_jk/500a71307a38c036d97bca3c6b958c710c22abd9ca1479c201a1ed17713b424f.json), card SHA `cd12ef9aff0a82aac45d1dc503cce2177f92b52ec95fd229fe4f6d0195ce68d2`. Visual pages p.2–4, 8, 10, 12, 15–17 were checked in the review.
- Direct J01 scope: conservative/stable time-inhomogeneous Q-function construction; Proposition 2.2 forward equation only on the source’s q-bounded set; Theorem 3.1 nonexplosion equivalence; Theorem 3.2 drift sufficiency and regularity-bounded necessity; Corollary 3.2 and Theorem 4.1 drift/Dynkin criteria under their stated model.
- The paper does not establish an unrestricted generator domain for arbitrary test functions or a general Kolmogorov forward equation. Theorem 3.1/3.2 proofs rely in part on cited stability-theory inputs that remain external.
