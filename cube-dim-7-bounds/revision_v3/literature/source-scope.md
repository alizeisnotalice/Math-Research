# 正文引用范围（主稿可直接采用）

- `Aldaz2017Covering`：J. M. Aldaz, *The Stein–Strömberg covering theorem in metric spaces*, J. Math. Anal. Appl. 449(2) (2017), 1741–1753, DOI10.1016/j.jmaa.2016.11.070。实际数学阅读版本为arXiv1605.05596v2（2017-01-07），Definition2.2 p2、Theorem5.2及证明pp14–15。v1同一定理5.2已比对，不存在系数2至1的版本变更。精确量词为任意固定ε>0，n≥N(ε)时，对Lebesgue任意范数球的全部尺度中心算子有弱常数≤(2+3ε)nlogn。不要直接取ε随n变化；不要将稀疏半径Theorem5.1常数6用于全尺度。正式版题录核实，未另取正式版PDF逐字对照。
- `Melas2003`：Antonios D. Melas, *The best constant for the centered Hardy–Littlewood maximal inequality*, Ann. of Math.157(2) (2003),647–688, DOI10.4007/annals.2003.157.647。原文Corollary1（printed p649，PDF第3页）精确给C1=(11+sqrt61)/12且尖锐；definition(1.1)为全L1、所有h>0的中心区间Lebesgue平均。数字2不适合当一维最优常数。
- `Tomar2026Squares`：Tejasvi Singh Tomar, *Bounds for the weak type (1,1) constant of the centred maximal operator over squares*, Zenodo preprint v4（2026-10-22），DOI10.5281/zenodo.23554628。实际下载官方12页PDF并读pp1–5、8–10、12；原文Theorem1.2为C2≤3.59123，Theorem1.3为C2≥1.7027237335，Theorem2.2及pp4–5延拓覆盖全L1与全正尺度。正文标“公开预印本”或“预印本及其公开证书仓库”，不称正式期刊论文。固定代码bca281a706bbf6375b04d601ad19ba4ad0bf28b6。作者报告两核验器及Lean证明，本轮未重跑，不写成本轮新核验。
- `Aldaz2020Kissing`：v2 Theorem1.2+Corollary3.5支持测度类统一尖锐2^n，Lebesgue特例≤2^n；不支持Lebesgue尖锐性、给定维数最好公开上界或“七维最好128”。

本轮有限12查询及原始来源跟进未核到同一全尺度Lebesgue中心方体模型的公开系数1或更好渐近结果。仅据此可说“相对于Aldaz上述明确估计改善前导系数”，不能说“世界首次/最新最好”。全部内部上下界仍为既有接受证据的只读汇总；精确端点在internal_baselines.json。本轮未启动下界或上界计算。
