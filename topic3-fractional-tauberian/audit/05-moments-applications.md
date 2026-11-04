# 第 1 轮核查：高矩、弱型与势函数应用

## 范围与判定口径

原稿第 11 节与第 14--16 节分别位于原 PDF 第 30--36 页和第 46--53 页。逐项核查时固定区分：Lebesgue 弱型与任意底测度弱型、固定输入的渐近与对全体集合取上确界的渐近、势函数给出的下界与 PDE 正则性所需的反向估计。

## 第 11 节：弱型与高矩

| 原稿 | 页 | 复核结论 | 新位置 | 状态与改动 |
|---|---:|---|---|---|
| 引理 11.1 | 31 | 有限中心方体贪心选择正确。须排除有限个坐标超平面后按象限证明重叠至多为 \(2^n\)。 | Lemma `mom:selection` | 已证明；补足零测集处理。 |
| 引理 11.2 | 31--32 | 固定半径平均的平移连续性足以产生稳定见证；紧性才允许有限化。 | Lemma `mom:stable-witness` | 已证明；明确只用于紧子集。 |
| 定理 11.3 | 32--33 | 常数 \(2^{n/q_\beta}=2^{n-\beta}\) 正确，前提是 Lebesgue 测度、\(f\in L^1\)。 | Theorem `mom:centered-weak` | 已证明；以紧内正则性补足从紧集到开水平集的步骤。 |
| 式 (63) 与注释 | 33--34 | \(C_{\beta,\mathcal Q_c}(\alpha)\le2^n\alpha^{-q_\beta}\) 仅为统一上界，不能称端点 sharp。 | `mom:centered-tauberian-bound`, Remark `mom:weak-scope` | 已修正措辞。 |
| 引理 11.4 | 34--35 | 方体的 \(3^n\) 覆盖与强矩形的 \(L^2\) 上界只作为低阈值支配使用。 | Proposition `mom:low-level` | 已补全紧集有限严格见证、不交子族及 \(\lambda^{-q_\beta}\sum s_j^{q_\beta}\le\lambda^{-q_\beta}|E|^{q_\beta}\)；强矩形部分注明逐坐标迭代或 de Guzmán 的 \(L^2\) 有界性。 |
| 推论 11.5 | 35 | 固定 \(Q_0\) 的精确高矩首项由第 4 节已证明的固定模型切向极限和全域 \((1-\alpha)^n\) 上界推出。 | Theorem `mom:fixed-cube-moment` | 已证明；未升级为全体极值。 |
| 定理 11.6 | 36 | 非中心的 \(1+O(p^{-1})\)、中心的 \(O(p^{-n})\)，以及两条非中心 \(\liminf\) 都由统一 Tauberian 上界和固定方体/矩形测试得到；固定 \(Q_0\) 给中心 \(\liminf\)。 | Theorem `mom:restricted-growth` | 已补回 \(2n/(n-\beta)\) 与 \(2n^2/(n-\beta)\) 的下界。原稿已明确固定 \(Q_0\) 首项和全体 \(\liminf\) 的区别；本次保留该区别。 |
| 注 11.7 | 36 | “光晕端点阶的泛函后果”正确，但必须有全体集合上界，单一测试集不足。 | Remark `mom:global-gap` | 已明确。 |

## 第 14 节：势函数与 PDE 接口

| 原稿 | 页 | 复核结论 | 新位置 | 状态与改动 |
|---|---:|---|---|---|
| 定理 14.1 | 46--47 | Riesz 核比较方向、常数 \(n^{-(n-\beta)/2}\) 和 \((2/\sqrt n)^{n-\beta}\) 正确；集中近点质量给出 sharpness。 | Theorem `app:riesz` | 已证明。 |
| 注 14.2 | 47 | 两个相同的 \(2^{n-\beta}\) 数值因子来自不同对象，不能推断两条不等式等价。 | Remark `app:factor` | 已保留并压缩。 |
| 定理 14.3 | 47--48 | Wolff 势尾积分给出常数 \(\gamma^{-1}n^{-\gamma/2}\)，中心情形再乘 \(2^\gamma\)。Dirac 质量给出 sharpness。 | Theorem `app:wolff` | 已证明；明确 \(\gamma=(n-p)/(p-1)\)。 |
| 推论 14.4 | 48 | Riesz 势下界把已构造的方体光晕包含送入势的超水平集，方向正确。 | Corollary `app:fractional-laplace` | 已证明；仅限方体测试资料。 |
| 推论 14.5 | 48--49 | 点态势上界只能排除资料过度集中，不能给出解正则性。 | Corollary `app:concentration` | 已证明。 |
| 推论 14.6 | 49 | 结论依赖已知的下势估计和该点 Wolff 势可用；对任意有限测度不能无条件宣称。 | Corollary `app:plaplace` | 已补充全部前提。 |
| 14.4 与问题 14.7 | 49--50 | 单向势比较不足以推出 good--\(\lambda\) 或 Schauder 型结论。 | Remark `app:boundary` 及本节末问题 | 已降级为开放接口。 |
| 第 15--16 节 | 50--51 | 全局中心 sharp 首项、一般凸体极值、双参数极限及一般多输入常数仍未闭合。 | `app:section` 的末段 | 原稿当时列为开放；第9轮gl:main闭合中心方体全体首项，其余范围继续保留。 |

## 原始文献与使用范围

- P. A. Hagelstein and I. Parissis, *Solyanik estimates in harmonic analysis*, arXiv:1310.3771：术语与经典光晕背景。
- J. M. Aldaz, *Kissing numbers and the centered maximal operator*, arXiv:1811.10478：中心极大函数、Besicovitch 交叠的背景；不作为本节 Lebesgue 选择引理的替代证明。
- T.-N. Nguyen and M.-P. Tran, *Level-set inequalities on fractional maximal distribution functions and applications to regularity theory*, J. Funct. Anal. 280 (2021), 108797, arXiv:2004.06394：说明需要额外反向分布机制的 PDE 语境。
- M. Kwaśnicki, *Ten equivalent definitions of the fractional Laplace operator*, Fract. Calc. Appl. Anal. 20 (2017), 7--51, arXiv:1507.07356：分数 Laplace 与 Riesz 核归一化的来源。
- T. Kilpeläinen and J. Malý, *The Wiener test and potential estimates for quasilinear elliptic equations*, Acta Math. 172 (1994), 137--161；T. Kuusi and G. Mingione, *Guide to nonlinear potential estimates*, Bull. Math. Sci. 4 (2014), 1--82：Wolff 势下势估计使用时应从具体 PDE 定理核实全部前提。

## 可复现实验与非证明性检查

本轮不以浮点计算承担任何定理。高矩公式的数值实现应在主项目的验证脚本中使用层析恒等式、Beta 函数和第 4 节给出的精确 \(|S_{n,c}|\) 交叉检查。数值拟合只能检验固定输入候选，不能替代统一 Tauberian 上界。

## 第9轮状态更新

保持原稿区分固定资料与全体liminf的历史核查结论；新增gl:main的全体匹配上界现将mom:restricted-growth的中心结论提高为n!L p^(−n)+O_J(p^(−n−1))。证明使用逐输入统一层析上界、固定Q0下界及低阈值弱型支配，不交换sup与积分。势函数单向比较及PDE限制没有扩大。

## 第18轮次阶状态更新

二维半阶pb:volume的逐输入全体上包络经层析给mom:pair-bracket：p³(F−4p⁻²)的liminf/limsup夹在[0,8]。上界才使用C的全体上包络；下界固定Q0，不把矩形上确界R当作一份输入的全部水平集。下一系数存在性、精确值及PDE反向传递仍未解决。
