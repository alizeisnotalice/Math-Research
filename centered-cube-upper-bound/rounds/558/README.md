# 第558轮：完整核族接收几何与全行动预算

**一般中心立方体极大算子的 $n^{1/2+o(1)}$ 一致弱型上界仍未证明。** 本轮取得改变构造决策的完整核族估计，以及原完整概率场的一项真实几何系数改进。局部接口、核正性、有限数值及条件性归约均不算一般闭合。

## 完整核族结论

对第556—557轮顶点生成元、$Q=[-1,1]^n$，完整归一切换核 $K_{n,k}=(P_k+F*P_k)/5$ 满足
\[
\operatorname*{ess\,inf}_QK_{n,k}
\le C_s2^{-n}(\log n/n)^s,\qquad n\ge2,\quad 1\le k<\infty,
\]
其中任意固定 $s\in(1/2,1)$，常数独立于 n 和 k。证明使用完整共同跳幅的双余弦负矩、四阶矩与条件浓缩。实际接收集 $A_n\subset Q$ 的 Lebesgue 体积至少为 $3\cdot2^n/4$，且与 k 无关。因此还排除跨有限 k 的概率混合达到固定对数幂根阶。**这不否定一般极大算子目标。**

两个真实 Cauchy 坐标支付全部大时间尾，进一步证明
\[
R_\lambda^A\delta_0(A_n)\le\frac{C_s}{1+\lambda}(\log n/n)^s,\qquad\lambda\ge0.
\]
结合全过程寿命下界，对全部 $G=R_{\lambda+h1_{D^c}}\delta_0$，$D\supset Q$、$\lambda>0$、有限 $h\ge0$，有
\[
\frac{\int G}{2^n\operatorname*{ess\,inf}_QG}
\ge\frac{c_s}{1+\log(8\log(4n))}(n/\log n)^s.
\]
全部参数可依赖 n，包括 $\lambda\downarrow0$。未将自由过程的条件独立性移植到被杀过程。此结论不涵盖 Q 内非恒定杀死或其他生成元。

另一全维生成元为
\[
\nu_\alpha(dz)=\frac{\alpha}{n2^n}\|z\|_\infty^{-n-\alpha}dz,\qquad0<\alpha\le1/2.
\]
正跳构造支付 $Cn^{1+\alpha}/\alpha$ 的全 Q 预算；完整多跳与共同时间分析又给
\[
2^n\lambda\operatorname*{ess\,inf}_QR_\lambda\delta_0
\le\frac{\alpha}{4n}+(2047/2048)^n.
\]
其全部 resolvent 及正确定性时间混合的归一成本至少 $n/1025$。取 $\alpha_n=(2+\log n)^{-1}$ 时，最佳 resolvent 成本为 $\Theta(n\log n)$。仅调整跳数或共同时间不能达标；空间/路径依赖权重仍在定理范围之外。

## 原完整场的改进及剩余几何

保留全部父对 Gram 交叉，$E=\|A_\pi\|^2$。逐父更新严格给
\[
E_{\rm new}\le E-\int\delta_B^2/e_B\,dB.
\]
保端点原子的有限微步构造还保证每宏阶段能量减半，或得到真正全行动缺口 $\mathcal K\le E/16$。原子、选择误差、阶段转换及终端已计账；近平价分支仍可能出现，不是无条件收缩。

在 $\mathcal K\le\theta E$ 下，完整探索律 $p_G\,dg/g+p_G\delta_G$，$p_G=(1+\log G)^{-1}$，严格给
\[
(1-\theta)E\le\frac{C_*(G)S A_n}{1+\log G}W,\qquad
A_n=1+n\log(b_{\rm geom}/a_{\rm geom}).
\]
原场立方体参数是边长，真实中心列积分没有遗漏 $2^n$，全部相位与跨父交叉保留。该律在共同、来源无关膨胀律的径向包络方法中最优，仍不支付多对数预算。放大 G 会增加完整母账中的 $G^2J_{\rm ph}$ 和 $C_*(G)$。

真正未证的几何蕴含是
\[
\mathcal K\le E/16\Longrightarrow
E\le\xi Q_{\rm bulk}/G^2+L(n)W,\qquad
\xi<\frac{(1-\beta)^2(1-3h_{\rm win})^2}{J_{\rm ph}},
\quad L(n)\lesssim\log^B(n+2).
\]
原入口、宽距离、低强度、初始、删除跳跃、未来被删来源、极端尺度比、全尺度及终端费用全部保留。内部高杀死另有合法检查点，但外区实际续程接收仍未证，形式上抵消选择概率不算几何改善。

## 实际验证及复现

- 指定本机 Python3.12.14、SymPy1.14.0、python-flint0.9.0/Arb256。MATLAB 本轮未运行。
- $s=3/4,q=40$：SciPy 与50位 mpmath 独立方法分别得到 $I_{1/2}\approx-4.731777548642816$、$-4.731777548604387$，差值在所报浮点估误内。这只是数值探路。
- $q=512$：解析分段覆盖原点、全部余弦零点与无穷尾；有理/Arb严格证书给 $I_{1/2},I_\infty\le-243474392372/70233453593<0$、$\mathbb Ew(R)\le18729711000/70233453593<1$。主代理实际独立复跑通过。
- 早期 n=2 三档截断试验发现未截断 cell 面邻接权发散，粗截断显著改变 cell 均值，因而停止扩大模型。cell 均值不是 essinf；δ 源及连续误差未被假定已付。
- 源分布、负矩、真实接收体积、寿命、角点 essinf、探索列界与回代分别经过独审。职责、范围及原稿哈希见包内 AUDIT.md。

总报告 REPORT.tex、详细 LaTeX、最小复现代码、实际运行摘要和逐成员 SHA-256 均保留。重复审稿全文、旧轮解压副本和可重生数据留在仓库外，不编译 PDF。最新协作按任务启停，不为满员派工；没有逐模型实际扣额数据，不估算成本比例。

## 归档

研究包26个成员，92,314字节，SHA-256：02b0954c48e12282a09eb42278e5528dc1ac18d6dbf20ed313072745d2ddb80b。ZIP CRC、逐成员哈希及主代理复跑通过。

仓库仅新增本轮 README、research.zip，并更新项目总览。已核对 [GitHub 官方规模说明](https://docs.github.com/en/repositories/creating-and-managing-repositories/repository-limits)：本轮远低于单对象100 MB及单次推送2 GB；网页diff的300文件显示上限不是仓库总文件数上限。
