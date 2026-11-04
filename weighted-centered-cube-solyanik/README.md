# 项目2：加权中心方体 Solyanik 估计

当前归档为第30轮。研究内禀权平均、外层加权光晕及共享中心/半径的多集合公共交集归一化。

- `paper/paper.tex`：当前中文主稿入口。包含原有内禀端点、水位传递、一般无原子测度临界与留数，第18轮可比齐次权的外层精确首项和两档二阶范围，第19轮等价双测度的弱型极限、临界传递和外层显式界，第20轮分离奇点混合权的全阈值精确公式，第21轮低阶背景保持、重合奇点分组与角向合并，第22轮径向匹配的尖锐小亏损稳定性，第23轮三次矩全局相对L1最优二次界，第24轮全p精确角向稳定包络及反演，第25轮小τ的单原子弱型测试和临界下界，第26轮全0<τ<1的双原子严格改进，以及第27轮完整小齐次次数渐近，并由第28轮统一到角向A1分类、无界轮廓及对数余项；第29轮以集合浓度推广近单位首项，允许无界或趋零角向轮廓；第30轮闭合负角幂临界次数的精确对数首项及非尖锐相对余项。
- `paper/homogeneous_extension.tex`、`homogeneous_applications.tex`：统一外层证明、几何匹配准则和可计算算例，取代此前重复的三段专章。
- `paper/mixture_extension.tex`：分离齐次奇点正和的精确最大值公式、非齐次端点与奇点合并反例。
- `proofs/`：一页独立量化平移证书、编译入口及旧专章的历史提交索引。旧轮报告中的路径以当时提交为准。
- `src/`、`experiments/`：必要的Python、Wolfram、MATLAB程序与最小记录；第18轮复现说明在`src/homogeneous/README.md`，第19轮在`src/outer_critical/README.md`，第20轮在`src/separated_mixture/README.md`，第21轮在`src/low_order_background/README.md`和`src/angular_grouping/README.md`，第22轮在`src/angular_stability/README.md`，第23轮在`src/cubic_l1/README.md`，第24轮在`src/angular_envelope/README.md`，第25轮在`src/atomic_weak/README.md`。第27轮在`src/small_homogeneity/README.md`，第28轮在`src/angular_a1/README.md`，第29轮在`src/round29/`，第30轮在`src/critical_angular/`。有限实验不替代一般证明。
- `rounds/`：逐轮报告与独立审阅；`STATE.md`标记当前有效结论。
- `sources/`：原项目输入、版本与文献核对范围。

在`paper`目录用XeLaTeX编译`paper.tex`两遍；本机字体为宋体、黑体与Times New Roman。`proofs/certificates.tex`在其目录单独编译。PDF、OCR、渲染、完整本地检查点和ZIP保留于本机，Git只保留精简归档。

第1–5轮在逐轮推送要求之前于本地完成，未伪造对应Git记录；第6轮建立基线，此后每轮独立提交并推送。第18轮清理的重复TeX仍可从第15–17轮历史及本机交付包读取，未删去相应数学结论。
