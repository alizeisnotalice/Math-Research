# 固定单位方体有限尺度归约的独立复核

## 结论

该推导的主结论成立，但只能作为固定输入
\(Q=[-1/2,1/2]^2\) 的高水平光晕精确公式；它不产生全体集合
Tauberian 常数的等式。四个面积公式和列出的 \(t^2,t^3\) 系数均正确。

有两项必须在进入正文前补足：

1. “原点邻域”及相应 \(\alpha _0(\beta)\) 必须依赖于固定的
   \(\beta\)。该邻域不能在 \(\beta\to1/2,1,3/2\) 时取一致值。例如
   \(\beta=1.49\) 时，若 \(s\) 大于约 \(0.01\)，即使 \(r-s\) 很小，
   候选 \(D\) 也可以超过 \(C\)；只有缩小到随 \(2-\beta-1/2\)
   变化的邻域后，\(C\) 才统一占优。
2. 候选比较的草稿论证应改成有量化余量的局部引理，特别是
   \(B/C\) 和 \(C/D\) 的比较。仅说“导数在原点为正/负”不足以控制
   \(h=r-s\) 与 \(s\) 同时趋零的全部路径。

## 有限尺度归约

令 \(d_j=|x_j|<1/2\)。边长为 \(\ell\) 的、以 \(x\) 为中心的方体
与 \(Q\) 的交长度为
\[
 q_j(\ell)=\min\left\{1,\ell,\frac{1+\ell}{2}-d_j\right\}.
\]
因此
\[
 F_x(\ell)=\ell^{\beta-n}\prod_{j=1}^nq_j(\ell).
\]
在由有限集合 \(\{1-2d_j,1+2d_j\}\) 切出的每个对数区间内，
\(\log F_x(e^u)\) 是仿射函数与若干
\(\log(e^u+a)\)（\(a>0\)）之和，故为凸函数。每个闭区间上的
最大值在端点取得；两个无界端的极限分别为零，因为
\(F_x(\ell)\sim\ell^\beta\)（\(\ell\downarrow0\)）及
\(F_x(\ell)\sim\ell^{\beta-n}\)（\(\ell\uparrow\infty\)）。故
\[
 M_{\beta,\mathcal Q_c}\mathbf 1_Q(x)
 =\max_{1\leq j\leq n,\,\varepsilon\in\{-1,1\}}
 F_x(1+2\varepsilon d_j).
\]
这一步对全部 \(x\in\operatorname{int}Q\) 成立。原稿定理 6.2 的近重合
中心方体论证提供高水平光晕 \(x=O(1-\alpha)\) 的局部化，故在
\(\alpha\) 充分接近一时可以合法地只用该内部公式。

## 二维候选与比较

取 \(r=d_1\geq s=d_2\geq0\)。四个接触尺度产生
\[
 A=(1-2r)^\beta,\quad
 B=(1-2s)^{\beta-1}(1-s-r),
\]
\[
 C=(1+2s)^{\beta-2}(1+s-r),\quad
 D=(1+2r)^{\beta-2}.
\]
代入交长度可直接验证这些式子。

令 \(h=r-s\)、\(y=1-2s\) 及 \(z=h/y\)。则
\[
 \log\frac BA=\log(1-z)-\beta\log(1-2z).
\]
而 \(\log(1-z)/\log(1-2z)<1/2\) 对 \(0<z<1/2\) 成立，极限为
\(1/2\)。这给出 \(\beta\geq1/2\) 时 \(B\geq A\)，并给出
\(\beta<1/2\) 时在一个 \(\beta\)-依赖邻域内 \(A\geq B\)。

写 \(c=2-\beta\)。有
\[
 \frac DC=\frac{(1+2h/(1+2s))^{-c}}{1-h}.
\]
若 \(c\leq1/2\)，利用
\((1+2h/(1+2s))^{-c}\geq(1+2h)^{-1/2}\geq1-h\)
可得 \(D\geq C\)。若 \(c>1/2\)，则
\[
 \partial_h\log(D/C)|_{h=0}=1-\frac{2c}{1+2s}.
\]
固定 \(c>1/2\) 后选择 \(s\) 的邻域，使此量有严格负上界，再以连续性
得到 \(C\geq D\)。这就是上述第一项限制的来源。

最后
\[
 \log(C/B)=(\beta-2)\log(1+2s)-(\beta-1)\log(1-2s)
 +\log(1-h)-\log(1-2s-h).
\]
在 \(\beta=1\)，\(C-B=2s(r+s)/(1+2s)\geq0\)；对 \(\beta\)
严格递增，故 \(\beta\geq1\) 时 \(C\geq B\)。对每个固定
\(\beta<1\)，上式关于 \(s\) 的导数在 \((s,h)=(0,0)\) 为
\(4(\beta-1)<0\)，从而缩小联合邻域后有 \(B\geq C\)。

综上，局部最大候选依次为
\[
 A\;(0<\beta<1/2),\qquad B\;(1/2\leq\beta<1),
\]
\[
 C\;(1\leq\beta<3/2),\qquad D\;(3/2\leq\beta<2).
\]
临界点的相等只发生在低维边界集合或较高阶比较中；题稿选取右侧分支给出了
正确的光晕面积。

## 面积与级数核验

令 \(V_\beta(\alpha)=|H_{\beta,\mathcal Q_c}(\alpha,Q)|\)。在
\(\alpha\in(\alpha_0(\beta),1)\)，直接积分有
\[
 V_\beta(\alpha)=
 \begin{cases}
 (1-\alpha^{1/\beta})^2,&0<\beta<1/2,\\
 2-\dfrac{4\alpha}{2-\beta}
 +\dfrac{2\beta\alpha^{2/\beta}}{2-\beta},&1/2\leq\beta<1,\\
 4\left(\dfrac{(2-\beta)\alpha^{-1/(2-\beta)}}{3-\beta}
 -1+\dfrac\alpha{3-\beta}\right),&1\leq\beta<3/2,\\
 (\alpha^{-1/(2-\beta)}-1)^2,&3/2\leq\beta<2.
 \end{cases}
\]
第二、三行分别是八个有序第一象限扇形的面积；积分上下限为
\((1-\alpha^{1/\beta})/2\) 和
\((\alpha^{-1/(2-\beta)}-1)/2\)。

令 \(t=1-\alpha\)。符号级数的独立计算给出
\[
 V_\beta(1-t)=L_\beta t^2+T_\beta t^3+O_\beta(t^4),
\]
其中
\[
 (L_\beta,T_\beta)=
 \begin{cases}
 (\beta^{-2},-(1-\beta)\beta^{-3}),&0<\beta<1/2,\\
 (2/\beta,-4(1-\beta)/(3\beta^2)),&1/2\leq\beta<1,\\
 (2/c,2(1+2c)/(3c^2)),&1\leq\beta<3/2,\ c=2-\beta,\\
 (c^{-2},(1+c)c^{-3}),&3/2\leq\beta<2,\ c=2-\beta.
 \end{cases}
\]
Wolfram 15 和 SymPy 1.14 分别核验了候选代入及全部四个泰勒系数。另作
20 万个局部随机样本的候选比较，只作为反例搜索：在足够小的固定邻域内未见
反例；在 \(\beta=1.49\) 而邻域未经缩小时确实出现 \(D>C\)，与上述
非一致性警告相符。

## 发表建议

建议将该结果作为“固定单位方体的二维局部精细渐近”命题写入正文或附录，并在
定理中明示 \(\alpha_0=\alpha_0(\beta)\)。不得把它表述为二维全体
Tauberian 常数的首项或二阶项。
