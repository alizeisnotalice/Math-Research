# 课题1：直接相关文献核验

检索/核验日期：2026-10-29。范围限于指定的中心方体 Tauberian 常数、Besicovitch/weak 型常数，以及用户点名的 Solyanik 相关原文。以 arXiv 论文全文和出版社/期刊题录为主；未开展穷尽性综述。

## 主要结论（对象与已知阶数）

记
\[
 C^{\mathrm c}_{\ell_\infty,n}(\alpha)=\sup_{0<|E|<\infty}
 \frac{|\{x:M_{\mathrm c,\ell_\infty}\mathbf1_E(x)>\alpha\}|}{|E|},
\]
即 Lebesgue 测度下、以同中心轴平行方体为平均集的线性 Tauberian 常数。它与多线性共享半径且以公共交集归一化的课题1对象不同。

- **已有同一线性中心方体对象的一阶量级结果。** Solyanik 1993 的结果由 Hagelstein–Parissis 2014 明确转述为：\(C^{\mathrm c}_{\mathrm{HL},c}(\alpha)\to1\)，且
  \[
  C^{\mathrm c}_{\mathrm{HL},c}(\alpha)-1\asymp_n 1/\alpha-1\asymp_n 1-\alpha.
  \]
  这里 \(\asymp_n\) 是上下界至多相差维数依赖常数；它**不**给出精确一阶系数 \(2^n\)，也不提供扣除线性项后的二阶/高阶渐近展开。
- **Aldaz 的 \(2^n\) 不是上面固定 Lebesgue Tauberian 常数的系数。** 其 Theorem 1.2 识别的是对所有局部有限、\(\tau\)-可加 Borel 测度取上确界的中心 maximal weak \((1,1)\) 常数；Corollary 3.5 给出 \(\ell_\infty\) 范数下该统一弱型常数恰为 \(2^n\)，并可达到。此量词、范数和对象均须与 Lebesgue Tauberian halo 常数区分。
- **Hagelstein–Parissis 2015 的高阶/正则性结果不能直接当作中心方体渐近。** arXiv:1409.3811 的主设定是 homothecy-invariant 的凸集基所定义的几何 maximal operator；论文具体处理 strong maximal operator 及未中心的球/方体 Hardy–Littlewood maximal operator，并推出 halo 函数的局部 Hölder 正则性。其方体基是“点在所选方体中”的未中心对象，不是每个 halo 点必须作中心的 centered cube 对象；也没有给出 \(C^{\mathrm c}_{\ell_\infty,n}\) 的 \(2^n\) 首项系数或更高阶项。
- 本次核对**未发现**所列原文证明固定 Lebesgue 中心方体 Tauberian 常数的一阶精确系数为 \(2^n\)，也未发现其真实二阶/高阶展开；这只是对下列有限、针对性文献的阅读结论，**不能据此宣称新颖性或“首次”**。已有文献已经建立相同线性对象的精确阶数 \(\asymp_n(1-\alpha)\)。

## 可引用来源

### 1. 原始 Solyanik 论文（俄文原版及正式英文译本）

A. A. Solyanik, “Bordering functions for differential bases”, *Matematicheskie Zametki* **54**(6) (1993), 82–89; 英译 “Bordering functions for differential bases,” *Mathematical Notes* **54**(5–6) (1993), 1241–1245. DOI（英译）：[10.1007/BF01209086](https://doi.org/10.1007/BF01209086). 题录由 [MathNet 原刊记录](https://www.mathnet.ru/eng/mzm2453) 与 [Springer 英译记录](https://link.springer.com/article/10.1007/BF01209086) 核对。MathNet 原刊页面可读，全文 PDF 在本次 web 读取时超时；因此对具体定理范围采用下条论文对 Solyanik 定理的逐项转述，不假称已从俄文全文独立复核。

**范围/对象：** Hagelstein–Parissis 将该文列为中心 HL 算子与 strong maximal 算子的 Solyanik 极限结果来源；centered cube 是直接相关对象。与课题1的线性单指标版本同类，但不是共享半径多线性公共交集版本。

### 2. Hagelstein–Parissis，Solyanik 估计会议论文

Paul Hagelstein and Ioannis Parissis, “Solyanik Estimates in Harmonic Analysis,” in C. Georgakis, A. M. Stokolos, and W. Urbina (eds.), *Special Functions, Partial Differential Equations, and Harmonic Analysis: In Honor of Calixto P. Calderón*, *Springer Proceedings in Mathematics & Statistics* **108**, Springer, Cham (2014), 87–103. DOI: [10.1007/978-3-319-10545-1_9](https://doi.org/10.1007/978-3-319-10545-1_9). 正式版书章与题录见 [Springer 书章页](https://link.springer.com/chapter/10.1007/978-3-319-10545-1_9)；开放预印本为 [arXiv:1310.3771v2](https://arxiv.org/abs/1310.3771v2)（2014-05-26 版本）。

**范围/对象：** 预印本 Theorem 1（PDF 第4页）明确区分未中心与中心算子，陈述中心方体和中心球 Tauberian 常数趋于1，并给出 \(C^{\mathrm c}_{\mathrm{HL},c}(\alpha)-1\asymp_n1/\alpha-1\)。同定中心方体的 Lebesgue 线性 Tauberian 对象，直接支持“已有线性阶量级”；没有精确系数或余项展开。其对未中心 balls/strong maximal 的结论是不同对象，不可移植为 centered cube 结论。

### 3. Hagelstein–Parissis，halo 函数 Hölder 正则性

Paul A. Hagelstein and Ioannis Parissis, “Solyanik estimates and local Hölder continuity of halo functions of geometric maximal operators,” *Advances in Mathematics* **285** (2015), 434–453. DOI: [10.1016/j.aim.2015.07.032](https://doi.org/10.1016/j.aim.2015.07.032). 正式题录见 [Elsevier](https://www.sciencedirect.com/science/article/pii/S0001870815002984)；开放预印本 [arXiv:1409.3811v2](https://arxiv.org/abs/1409.3811v2)（2015-07-10 版本）。

**范围/对象：** Theorem 1 将强 maximal 及未中心的方体/球 maximal 的已知 Solyanik 阶数作为前提/输入；Corollary 4 给出方体基 halo 函数 \(C^{1/n}\) 正则性。中心方体算子因其中心依赖，不属于文中所定义的 homothecy-invariant set-basis maximal operator；故不是课题1中心方体高阶渐近的直接结果，也不能证明或否定 \(2^n\) 精确首项。

### 4. Aldaz，中心 maximal weak 型常数与 \(\ell_\infty\) 方体

J. M. Aldaz, “Kissing Numbers and the Centered Maximal Operator,” *Journal of Geometric Analysis* **31** (2021), 10194–10214. DOI: [10.1007/s12220-021-00640-1](https://doi.org/10.1007/s12220-021-00640-1). 正式版题录见 [Springer](https://link.springer.com/article/10.1007/s12220-021-00640-1)；指定全文为 [arXiv:1811.10478v2](https://arxiv.org/abs/1811.10478v2)（2020-03-10 版本）。

**范围/对象：** Theorem 1.2：Besicovitch 常数等于对所有局部有限、\(\tau\)-可加 Borel 测度取上确界的 centered maximal weak \((1,1)\) 算子范数。Theorem 2.10：反例可离散有限化。Theorem 3.2：有限维任意范数下该几何常数等于 strict Hadwiger number。Theorem 3.3：存在达到常数的离散测度。Corollary 3.5：\(\ell_\infty\) 范数（轴平行中心方体）下统一 weak \((1,1)\) 常数恰为且达到 \(2^n\)。**它研究所有 Borel 测度上的函数 weak 型常数，不是固定 Lebesgue 测度的 Tauberian 常数随 \(\alpha\uparrow1\) 的展开。**

**题录修正：** 此文 DOI 核为 `10.1007/s12220-021-00640-1`；此前本地长版正文参考文献出现的 `10.1007/s12220-020-00584-w` 与该正式题录不符，应由 Bib 核验者修正。

## 来源与阅读记录

- 读取了 Aldaz v2 的 arXiv PDF 全文，定位并核对 Theorems 1.2, 2.10, 3.2, 3.3 和 Corollary 3.5；并以 Springer 正式题录核验出版年、卷页和 DOI。
- 读取了 Hagelstein–Parissis 1310.3771v2 的 arXiv HTML/PDF 中引言及 Theorem 1；以 Springer 正式书章页核验书目信息及其对 Solyanik 1993 的引用。
- 读取了 Hagelstein–Parissis 1409.3811v2 的 arXiv HTML，核对摘要、Theorem 1 和 Corollary 4 的对象；以 Elsevier 正式期刊记录核验卷期页和 DOI。
- Solyanik 原文由 MathNet 与 Springer 元数据页核验题录；原文全文未能通过 MathNet 下载读取。本文关于其中心方体定理的数学范围依赖 Hagelstein–Parissis 对原定理的明确陈述。
