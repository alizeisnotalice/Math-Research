# 第 10 节外层权模型审计

审计范围：原稿第 23--30 页，原编号 (46)--(59)。本节中的权均为外层测度；光晕定义中的平均和阈值均保持 Lebesgue 版本。除特别说明外，\(0<\beta<n\)、\(\theta=\beta/n\)。

|原编号与页码|原结论及必要前提|本次核查依据与处理|状态|新标签|
|---|---|---|---|---|
|(46), p.23|定义 \(C^w_{\beta,\mathcal B}\)，要求 \(0<|E|<\infty\)、\(0<w(E)<\infty\)。|定义保留，并明确不是内禀 \(w\,dx\) 平均。|proved/definition|`wgt:definition-concentration` 前的定义|
|(47), p.23|定义 \(\Omega_{\mathcal D,w}\)，\(\mathcal D\in\{\mathcal Q,\mathcal R\}\)，分母仅在 \(w(D)>0\) 时使用。|保留；补足 \(0<u<1\) 和正分母约定。|proved/definition|`wgt:definition-concentration`|
|定理10.1, (48)--(50), pp.23--24|集中模量的正向控制；使用单一锚定严格见证、先前的双密度误差和见证聚拢。|逐步复核：\(B_0^\ast\setminus B_0\)、\(B_0\setminus E\)、\(B_0^\ast\setminus E\) 三个相对 Lebesgue 误差均为 \(O(\rho_\theta)\)；非中心情形得到 \(1+\omega/(1-\omega)^2\)，\(\omega\le1/2\) 时至多 \(1+4\omega\)。中心情形使用同心小方体，其体积比例为 \(O(\rho_\theta^n)\)。|proved|`wgt:theorem-forward`; `wgt:Q-forward`; `wgt:R-forward`; `wgt:Qc-forward`|
|定理10.2, (51), pp.24--25|反向集中估计；须有 \(w>0\) 几乎处处，且基是非中心 \(\mathcal Q\) 或 \(\mathcal R\)。|保留。\(E=D\setminus F\) 时，同一个 \(D\) 为其每一点见证；严格阈值由 \((1-\delta)^{1-\theta}>1-u\) 给出。中心基不能照搬，因为固定 \(D\) 未必以其中每点为中心。|proved|`wgt:theorem-reverse`; `wgt:reverse`|
|推论10.3, (52)--(53), p.25|外层加权 Solyanik 极限分别刻画 \(A_\infty\)、\(A_\infty^\ast\)。|正文定义 \(A_\infty=\bigcup_{p>1}A_p(\mathcal Q)\)、\(A_\infty^\ast=\bigcup_{p>1}A_p^\ast\)，并引用 \(\Omega\to0\) 的标准集中刻画，见 \cite{HagelsteinParissis2016,HagelsteinParissis2018}。正反集中估计完成本稿中的传递；未推广至内禀加权平均。|proved under stated hypotheses|`wgt:corollary-Ainfty`|
|推论10.4, p.25|非中心/强矩形的最佳幂指数等于集中指数，且与固定 \(\beta\) 无关。|幂型集中界是对应 \(A_\infty\)/强 \(A_\infty^\ast\) 的标准定量形式（同上两篇文献）；正反比较和 \(\rho_\theta(\alpha)\asymp1-\alpha\) 给出指数恒等。结论只涵盖 \(\mathcal Q,\mathcal R\)。|proved|`wgt:exponent-preservation`|
|推论10.5, (54), pp.25--26|幂型集中条件给出两个非中心基的 \(\delta\) 阶与中心基的 \(n\delta\) 上界。|已拆分前提：方体集中只推出 \(\mathcal Q\) 与 \(\mathcal Q_c\) 结论；矩形集中才推出 \(\mathcal R\) 结论，并自动蕴含方体前提。中心指数一般不应称为 sharp，特别是各向异性奇性权。|proved as upper bound only|`wgt:corollary-power-conc`; `wgt:power-concentration-cube`|
|引理10.6, p.26|正有限 Lebesgue 点处，缩放集合的加权比趋于 Lebesgue 比。|原证明正确；重写时明确 \(E,F\) 有界，并以共同包络方体控制平均振荡。|proved|`wgt:proposition-localization`; `wgt:localization`|
|引理10.7, (55), p.26|径向核 \(|x|^{-a}\) 的集合积分上界，\(0<a<d\)。|层饼、Tonelli 和精确分段积分复核无误；常数为 \(d(d-a)^{-1}\omega_d^{a/d}\)。|proved|`wgt:lemma-singular-integral`; `wgt:singular-integral`|
|命题10.8, (56), pp.27--28|一坐标权 \(|x_1|^{\gamma-1}\) 对 \(\mathcal Q,\mathcal R\) 的精确阶，\(0<\gamma<1\)。|一维集中估计经远离/邻近奇点分解成立；层截与 Jensen 将其提升至任意矩形。平移薄条 \(E_h,F_h\) 给出两个基的匹配下界；已补 \(E_h\subset H\)，故 \(E_h\cup F_h\subset H\) 才能导出 \(C^w-1\) 的下界。|proved|`wgt:proposition-one-coordinate`; `wgt:one-coordinate`; `wgt:one-dimensional-concentration`|
|命题10.9, (57)--(59), pp.28--30|横向径向权 \(|y|^{-a}\) 对中心方体的 \((1-\alpha)^{n-a}\) 阶，\(0<a<k\)。|严格保留两类输入：(i) 方体集中 \((57)\) 只用于证明 \(w(B_0\setminus E)=o(w(B_0))\)；(ii) 同心缩放 \((59)\) 才产生 \(n-a\)。因此没有把一般集中指数 \(1-a/k\) 代入通用 \(n\delta\) 上界并误称为 \(n-a\)。固定 \(Q_0\) 的内部方体给出匹配下界。|proved|`wgt:proposition-radial`; `wgt:radial-concentration`; `wgt:radial-dilation`; `wgt:radial`|

## 实质边界

- 正向权估计需要的是外层权的集中模量；若将平均改成 \(w\)-平均，原证据链不适用。
- 反向估计和 \(A_\infty\) 刻画都依赖 \(w>0\) 几乎处处；允许权在正 Lebesgue 测度集上为零时，补集测试无法保证 \(w(E)>0\)。
- 中心通用式给出的 \((1-\alpha)^{n\delta}\) 仅为上界。\(w_{k,a}\) 的 \((1-\alpha)^{n-a}\) 由专门的同心缩放估计取得。
- 带 `wgt:` 的所有公式现由主文的章节自动编号；不再使用手工 `(W1)` 等标签，避免与全篇编号体系冲突。
