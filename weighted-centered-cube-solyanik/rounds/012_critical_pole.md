# 第十二轮：超临界临界线处的一次极点

## 结论与范围

令 \(\lambda_c=P_{\boldsymbol a}(1)\)，并令 \(\rho_c\) 解
\[
\sum_i(1-a_i/\rho_c)_+=1.
\]
在 \(s=1\) 的左邻域，令 \(I_c=\{i:a_i<\rho_c\}\)、
\(q=|I_c|\)、\(A=\sum_{I_c}a_i\)。严格解析水位计算给出
\[
P_{\boldsymbol a}(1-\beta)=\lambda_c
\left(1+\frac{\beta}{q-1}\right)^A,
\qquad \rho_c=\frac A{q-1}.
\]
这在 \(a_i=\rho_c\) 的活动集切换点也有效，因为该坐标从左仍有零缺陷。
故
\[
P'_{\boldsymbol a}(1-)=-\rho_c\lambda_c,
\qquad
\beta(\lambda)\sim\frac{\lambda-\lambda_c}{\rho_c\lambda_c}.
\]

更强地，令 \(W_\mu=\|M_{c,\mu}\|_{L^1\to L^{1,\infty}}\)。对每个非零
无原子局部有限 Radon 测度，有限 Boolean 原子上的确定性抽样证明
\[
\lim_{\beta\downarrow0}\beta C^{\rm int}_{c,\mu}(\beta)=W_\mu.
\]
它把有界紧支 \(f\) 的每个原子 $B$ 选为质量
\(\eta\int_Bf\,d\mu\) 的子集；零质量边界参考方体和以 \(\eta\) 为尺度的
中心扰动余量保证该结论对每个充分小 \(\eta\)，而非仅子列成立。顺序上先固定
\(\eta\)，逐点选参考方体及满足该 \(\eta\) 误差预算的中心邻域，再由紧性取
有限子覆盖、固定 Boolean 原子，最后才选择 \(E_\eta\)；因而没有把随
\(\eta\) 缩小的邻域误当成固定覆盖。接入已有的
超临界传递后，固定 \(\mu\) 的多线性常数满足精确临界首项系数公式
\[
\lim_{\lambda\downarrow\lambda_c}(\lambda-\lambda_c)
C^{\boldsymbol a}_{c,\mu}(\lambda)=\rho_c\lambda_cW_\mu.
\]
原先的 \(\asymp_{\mu,\boldsymbol a}(\lambda-\lambda_c)^{-1}\) 是直接推论。
无原子性不可删：\(\mu=\delta_0\) 时指标函数比值给出
\(\beta C_{c,\mu}^{\rm int}(\beta)\to0\)，而弱型范数为 1。

对 \(d\mu=w\,dx\) 的全权族，只有在已有方体全权族精确公式下才得到
精确临界首项系数：
\[
\lim_{\lambda\downarrow\lambda_c}(\lambda-\lambda_c)
\sup_w C^{\boldsymbol a}_{c,w}(\lambda)=2^n\rho_c\lambda_c.
\]
当 \(a_i=1\) 时，这就是 \(2^n(1-1/m)^{m-1}\)。

## 粗阶下界的独立核验

先固定 \(\mu\)，由内正则性选定一次且永不随 \(\beta\) 改变的有界紧
\(K\) 使 \(0<\mu(K)<\infty\)。再选固定 \(R,D\)，满足
\(K\subset Q(x,R)\subset D\) 对所有 \(x\in K\)。当
\(2\beta\mu(D)<\mu(K)\) 时，非原子性选出 \(E\subset K\)、
\(\mu(E)=2\beta\mu(D)\)。于是每个 \(x\in K\) 都有严格平均
大于 \(\beta\)，故比值至少为 \(\mu(K)/(2\beta\mu(D))\)。这既不使用
第十一轮的临界以下结论，也不把依赖 \(\mu\) 的下界常数宣称为一致常数。

## 实际计算核验

`critical_pole_certificate.py` 是无随机性的 Python 有理算术证书；它检查
六组整数指数、末段活动集、\(P(1)\)、局部精确公式和线性系数
\(\rho_c\lambda_c\)，并特别检查 \((1,1,2)\) 的切换
\(a_3=\rho_c\)。`critical_pole_check.wl` 由 Wolfram Language 做符号
导数、反演级数和全权族临界首项系数极限交叉核验。计算只核验标量代数，
不替代证明中的测度论步骤。

2026-11-05 实际运行的 Python 证书六例全部通过：
\((1,1),(1,1,2),(1,2,3),(1,1,3,5),(2,2,2,7),(1,1,1)\)。例如
\((1,1,2)\) 给出 \(\rho_c=2\)、\(\lambda_c=1/4\)、缺陷
\((1/2,1/2,0)\) 及线性系数 \(1/2\)，从而直接覆盖临界水位碰到指数的
情形。Wolfram Language 15.0.0 实际返回：左导数恒等式为 `True`；反演
级数的一阶项为
\((q-1)(\lambda-\lambda_c)/(A\lambda_c)\)；将全权族表达式相乘后的
临界首项系数为 \(c_nA\lambda_c/(q-1)=c_n\rho_c\lambda_c\)（方体族取
\(c_n=2^n\)）。

一维 Lebesgue 情形的 \(W_{dx}=(11+\sqrt{61})/12\) 采用
Melas, \emph{The best constant for the centered Hardy--Littlewood maximal inequality},
Annals of Mathematics 157 (2003), 647--688,
[DOI: 10.4007/annals.2003.157.647](https://doi.org/10.4007/annals.2003.157.647)；
TeX 片段使用 `\cite{Melas2003}`，由主稿维护该参考文献。

## 成稿检查

`critical_pole_extension.tex` 在 `extension.tex` 之后实际以 XeLaTeX 编译。
在临时主稿补入 `Melas2003` 题录后，`latexmk` 成功完成，日志无未定义引用、
公式溢出或重新运行提示；渲染末两页后检查了长证明的分页、公式和引文。主稿
只需加入同一 `Melas2003` 题录即可输入该片段。

`critical_pole_extension.tex` 可在 `extension.tex` 之后输入；它不依赖
第十一轮仍在核验的临界以下发散主张。


## 主代理验收与量词核对

主代理重跑Python六组精确水位实例和Wolfram15，均通过；另执行melas_coefficient_check.wl，验证已出版弱型常数的二次方程及m=2,3,4,8的临界首项系数。m=2时固定一维Lebesgue临界首项系数为(11+√61)/24≈0.7837604031627773，而全权族为1。Melas公式来自官方期刊摘要的既有定理，不是本轮计算证明。

最终主文将小阈值引理放在临界渐近定理之前，并明确先固定η、逐点选取有界误差邻域、再紧性取有限覆盖、再分割Eη。几何与Boolean家族可依赖η，故下界对全部充分小阈值成立，不只是子列。新增全权族分段公式由第11轮及既有超临界传递直接推出。

最终终审确认：已显式令Eη为各原子所选部分之并并在D外取空；整合稿12页，XeLaTeX双遍通过，无未定义引用或Overfull。新增章节与参考文献已渲染检查。现有特殊构造保留为历史复现附件，本文不重复使用其较弱结论。
