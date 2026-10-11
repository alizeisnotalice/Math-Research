> **冻结入库快照（2026-10-06）**：以下附件行与当时的相关性/阅读状态是历史证据，不是本轮最终覆盖表。历史原文另存于 `../../../audit_current/f/history/pre-audit-20261007/`。当前逐条数学断言请看 `../../../audit_current/f/claims.json`，实际运行请看 `../../../audit_current/f/cases.json`；本轮候选池身份筛查和逐页阅读账以 `../../../audit_current/independent-f/` 的审核报告为准。不得把本文中的 `pending` 自动解释为无关，也不得把旧卡的 `full_read` 标签视为本轮已复核。

便携包或已安装 Skill 的本地来源访问方法见[来源访问说明](source-access.md)。单条索引查找会验证来源 SHA；解析成功不等于已阅读或已证明。

# F02 来源与实际处理状态

附件归属 41 条；本专题关联不同PDF 40 份。
题名按附件文件名原样保留，未自动认定身份。下表的人工判定与读取范围来源于真实审读记录；pending不能当作无关。

|条目|文件名|SHA-256 / 论文ID|人工判定|理由|
|---|---|---|---|---|
|E0001|38_Jaming_2008_On the Fourier transform of the symmetric decreasing rearran.pdf|P-7d16bea48a9dfeca|irrelevant|首页摘要实际为symmetric decreasing rearrangement后的Fourier L1/L2局部范数比较与Sobolev regularity，不给tensor kernel正定性/非负Fourier或softmax平方上确界；本组不迁移其结论。|
|E0002|24_Stolyarov_2016_Functions whose Fourier transform vanishes on a surface.pdf|P-3e8f1ef605bea15a|irrelevant|首页为Fourier transform在codimension-one surface上消失的Lp闭包/密度问题，没有tensor正定核或softmax操作。|
|E0003|15_Bogo_2020_Modular forms, deformation of punctured spheres, and extensi.pdf|P-2490753fc3a62bdd|irrelevant|首页为modular forms与punctured-sphere deformation/tensor representations，tensor是表示论，不是张量卷积核Fourier正性。|
|E0004|25_Choi_2015_Averaged decay estimates for Fourier transforms of measures.pdf|P-b68257475d2ab063|irrelevant|首页为曲线上measure Fourier averaged decay/energy，不是正定性刻画或tensor soft-supremum。|
|E0005|26_Guella_2015_Strictly positive definite kernels on a product of circles.pdf|P-8205499d1bddf331|direct|全文刻画圆直积上的连续实各向同性正定核；严格正定等价于双频谱正支集穿透所有平移矩形子格，直接给出张量频率正性的公开构件。|
|E0006|28_Tephnadze_2014_A note on the Fourier coefficients and partial sums of Vilen.pdf|P-6d2d8a7efe41cb96|irrelevant|首页为bounded Vilenkin group上的Fourier coefficients、Paley/Hardy-Littlewood与partial-sum strong convergence，未提供F02所需核正定性或softmax。|
|E0007|12_Zelent_2025_Time frequency localization in the Fourier Symmetric Sobolev.pdf|P-71ef66b90e5cba7e|irrelevant|首页为Fourier symmetric Sobolev space concentration spectrum/RKHS，核心是eigenvalue plunge；非tensor Fourier-positive kernel或softmax平方上确界。|
|E0008|23_Madiman_2016_The norm of the Fourier transform on compact or discrete abe.pdf|P-b5117a459bd01a66|irrelevant|首页为compact/discrete infinite Abelian group上Fourier Lp→Lq operator norm/entropy uncertainty，没有kernel positivity或softmax预算。|
|E0009|02_Beltran_2022_$L^p-L^q$ local smoothing estimates for the wave equation vi.pdf|P-440bf8c97a0e0230|irrelevant|首页为wave local smoothing的k-broad restriction/broad-narrow argument，不是tensor positivity或softmax核入口。|
|E0010|18_Nicola_2019_A note on the HRT conjecture and a new uncertainty principle.pdf|P-339239bfb93d5350|irrelevant|HRT/Gabor time-frequency shifts与STFT uncertainty，不给tensor positive-definite kernel或softmax。|
|E0011|35_Grafakos_2011_On Fourier transforms of radial functions and distributions.pdf|P-35ab79d8223c1786|adjacent|径向Fourier变换的n→n+2 recurrence可迁移为kernel Fourier sign/normalization检查，并可用于F03旧径向幂次回源；需要全文核假设和Bessel公式。|
|E0012|13_Rozendaal_2021_Rough pseudodifferential operators on Hardy spaces for Fouri.pdf|P-186eaf5447f261f2|irrelevant|rough pseudodifferential/FIO Hardy-Sobolev mapping，不给张量正定核或softmax预算。|
|E0013|31_Csordas_2013_Fourier transforms of positive definite kernels and the Riem.pdf|P-333c28bbdea9a7a7|direct|首页明确把positive-definite canonical kernels与Fourier/Laguerre inequalities连接，可提取正性定义与scope；Riemann ξ应用需与一般tensor/softmax区分，排全文。|
|E0014|04_Liu_2024_On the positivity of Fourier transform of the stretched Gauß.pdf|P-2ebc46987d4f4c81|adjacent|全文证明 0<s≤2 的径向 stretched Gaussian Fourier 变换处处严格正；提供完全单调函数的 Gaussian 正混合方法，但非一般张量核或软最大定理。|
|E0015|03_Orponen_2023_On the Fourier decay of multiplicative convolutions.pdf|P-52342022c2d36a84|irrelevant|multiplicative convolution Fourier decay由energy hypotheses推出，非kernel Fourier positivity/softmax。|
|E0016|11_Langowski_2020_On derivatives, Riesz transforms and Sobolev spaces for Four.pdf|P-9d07ad98b6473ec5|irrelevant|Fourier-Bessel derivative/Riesz/Sobolev boundary问题，不给F02张量核正性或softmax。|
|E0017|29_Tephnadze_2014_A note on the Vilenkin-Fourier coefficients.pdf|P-2ed90a941d5fb895|irrelevant|bounded Vilenkin Fourier coefficient estimates，非tensor核正定性。|
|E0018|22_D'Agnolo_2017_Topological computation of some Stokes phenomena on the affi.pdf|P-79868a353a7b0a93|irrelevant|holonomic D-module的Fourier-Laplace/Fourier-Sato及Stokes，非harmonic-analysis tensor kernel。|
|E0019|17_Song_2024_Positive definiteness of fourth order three dimensional symm.pdf|P-ceb4efdd564d74cc|irrelevant|4阶3维±1 symmetric tensor的quartic polynomial正定性，与positive-definite convolution kernel/tensor product Fourier正性对象不同。|
|E0020|32_Hitzer_2013_The Orthogonal 2D Planes Split of Quaternions and Steerable.pdf|P-9b959d658de77160|irrelevant|quaternion Fourier transform及orthogonal-plane split/FFT实现，不是Abelian tensor正定核。|
|E0021|14_Neerven_2026_The paper _On the constant in a transference inequality for.pdf|P-2c45c32048568b99|adjacent|periodized squared-sinc powers的global minimum是Fourier transfer常数的explicit kernel比较；可迁移为F02 sinc正性/平方化与维数常数检查，排全文；不宣称softmax。|
|E0022|37_Kozma_2010_Singular distributions, dimension of support, and symmetry o.pdf|P-537bd19a19c0a929|irrelevant|circle singular distributions的support dimension与one-sided Fourier symmetry，非tensor正定核或softmax。|
|E0023|19_Castro_2018_Characterization of Strict Positive Definiteness on products.pdf|P-ac20188cfcb81698|direct|complex spheres product上positive-definite/disc-polynomial coefficients严格正性刻画，直接相关F02乘积结构，须保p/q=1/infinity分情形与离散频谱条件，排全文。|
|E0024|21_Gát_2017_Almost everywhere convergence of Fejér means of two-dimensio.pdf|P-28ab949dfb84247f|irrelevant|triangular Walsh-Fejer summability/a.e. convergence，不是F02张量核正性或softmax。|
|E0025|27_Xu_2015_Positivity and Fourier integrals over regular hexagon.pdf|P-49e834dfe88ef038|direct|hexagonal norm Fourier Riesz means的kernel positivity门槛及H-radial positive-definite functions，直接相关F02核正性/几何依赖，排全文。|
|E0026|01_Jaye_2025_A High-Frequency Uncertainty Principle for the Fourier-Besse.pdf|P-d327d8a324d06269|irrelevant|Fourier-Bessel high-frequency PLS uncertainty和damped wave observability，非tensor核Fourier正性。|
|E0027|20_Goolish_2018_Sparse Averages of Partial Sums of Fourier Series【题录】.txt|R-1771a0444e4f5604|irrelevant|完整题录摘要仅为稀疏Fourier partial-sum averages的uniform convergence，不是F02核正定性或softmax。无本地正文，身份和证明未核。|
|E0028|05_Fraser_2021_On the Fourier dimension of $(d,k)$-sets and Kakeya sets wit.pdf|P-a83271b921ca77ce|irrelevant|restricted Kakeya/(d,k)-set的Fourier dimension，非tensor positive-definite kernel。|
|E0029|30_Carneiro_2013_Entire approximations for a class of truncated and odd funct.pdf|P-5ef7684716744d0c|irrelevant|truncated/odd functions的entire L1 extremal approximation，没有F02指定tensor正性/softmax操作；本组不迁移其指数类型majorant。|
|E0030|16_Leclerc_2021_Julia sets of hyperbolic rational maps have positive Fourier.pdf|P-08d95ba9b15ccf53|irrelevant|hyperbolic rational Julia sets上的Fourier dimension/polynomial decay，positive修饰dimension，不是kernel positive-definite。|
|E0031|40_Dhaouadi_2007_Functions of q-positive type.pdf|P-2594d4a54906e45b|adjacent|q-Bessel Fourier transform的q-positive type/Bochner analogue给正性检验的明确不同底空间版本；可迁移到定义对齐审查，但不直接证明经典tensor核，排全文。|
|E0032|39_Li_2008_Discrete Fourier analysis on fundamental domain of $A_d$ lat.pdf|P-203094ba50b4455f|irrelevant|A_d lattice fundamental-domain Fourier/interpolation/cubature和Chebyshev nodes，未给F02正定tensor核或softmax。|
|E0033|08_Caputo_2025_Fourier transform of vector-valued graph signals.pdf|P-a404dc2b29e56805|irrelevant|Banach-valued graph signal spectral Fourier/operator norms，不是LCA群卷积tensor正定核。|
|E0034|06_Leclerc_2025_Fourier decay in parabolic $C^{1+α}$ systems with overlaps.pdf|P-cc140dcd0b2bb495|irrelevant|parabolic IFS equilibrium measure的Fourier decay，positive Fourier dimension不等于Fourier-positive kernel。|
|E0035|10_Khare_2025_The entrywise calculus and dimension-free positivity preserv.pdf|P-60830612beeb212e|direct|dimension-free entrywise positivity preservers综述连接positive-definite functions/Schur calculus；可用于F02 kernel entrywise transform与tensor coupling正性审查，需全文保每域假设和外引状态。|
|E0036|34_Talvila_2011_Fourier series with the continuous primitive integral.pdf|P-157cd2b54d703584|irrelevant|continuous primitive distributions/Alexiewicz norm的Fourier series与summability，非kernel positivity/softmax。|
|E0037|36_Bourgain_2010_Sur les séries de Fourier des fonctions continues unimodulai.pdf|P-73709ca192e5f48b|irrelevant|连续unimodular circle functions的topological degree与weighted Fourier coefficients，非tensor正定核。|
|E0038|33_Bie_2011_The class of Clifford-Fourier transforms.pdf|P-9be9c3d5ba89bc37|irrelevant|Clifford Fourier transform PDE class/eigenvalues/inverse，非Abelian tensor核Bochner正性。|
|E0039|27_Huang_2025_On the supremum of random cusp forms.pdf|P-11894fc52f6dd9c0|irrelevant|random cusp forms的expected supremum/concentration，未定义F02 softmax平方核或tensor Fourier positivity。|
|E0040|09_Bulj_2022_Multi-parameter maximal Fourier restriction.pdf|P-69448beefbc02932|irrelevant|multi-parameter Fourier restriction/maximal and Menshov-Paley-Zygmund theorem，不给F02 tensor positive-definite kernel/softmax。|
|E0041|07_Oganesyan_2022_Two-dimensional Hardy-Littlewood theorem for functions with.pdf|P-bb8dd5522407a1d5|irrelevant|2D general monotone Fourier coefficients的Hardy-Littlewood norm relation允许signed coefficients，非核正定性或softmax。|

## 实际阅读卡

- [P-2ebc46987d4f4c81：full_read](papers/EG-P-2ebc46987d4f4c81.json)
- [P-8205499d1bddf331：full_read](papers/EG-P-8205499d1bddf331.json)
- [S-01a4d818f5fb3fd6：full_read](papers/ROOT-S-01a4d818f5fb3fd6.json)
- [S-10977fcb18554549：full_read](papers/EG-S-10977fcb18554549.json)
- [S-3ba98fcad8583c4d：full_read](papers/EG-S-3ba98fcad8583c4d.json)
- [S-cc47322e4d57f42a：full_read](papers/ROOT-S-cc47322e4d57f42a.json)

## 本地原文定位

原附件PDF：工作包根目录的 `evidence/papers/<paper_id>/paper.pdf`；全文 `paper.md`；页码映射 `pages.json`。
补充PDF：`evidence/supplementary/<paper_id>/`。完整原文不放入轻量Skill ZIP；可用SHA-256与附件成员名定位原文件。
逐篇完整处理状态见工作包 `catalog/entries-reviewed.json`、`catalog/reading-queue.json`。
