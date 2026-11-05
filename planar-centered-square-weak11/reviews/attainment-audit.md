# 达到性与覆盖对偶模块：冻结稿独立终审

2026-11-05。**裁决：冻结稿全部明示结构结论已独立重建，无必要数学修补。** 本报告不判定全局 \(L^1\) 极值存在或不存在，也不计算锐 \(C_2\)。结论由以下证明核验支持，不以其他审核标签为前提。

## 0. 冻结对象与实读证据

已完整实读
/Users/zhengzhihao/Desktop/T/planar-centered-square-weak11/scratch/attainment-dual/attainment_and_cover_dual.md，
345 行、18645 bytes、无 CR；实际 SHA256 为
94fd2370609f8d3ec321ac7a89158c8e28dc4207f13bc0fe69cede82f0b11348。
已完整读 source/research-request.md。作者冻结前，本代理从定义独立写出 interface_reconstruction.md，未读动态稿；冻结后逐节核对了整份主稿，包括新增的联合半连续、移动阈值和强紧性等价表述。

实际读了作者两组配对计算源码与 stdout/JSON：finite_checks.wl/log 及 two_blobs.py/log/json。另有本代理独立实际运行的 finite_interfaces.wl/log 与 two_blob_lower.py/json/log。代码、输入几何与证据边界在第 9 节列明。只改本审查目录，作者稿、其他代理文件和仓库均未改。

核验范围为固定函数的阈值结构、强 \(L^1\) 传递、固定有限见证的两个 LP、一般中心权的条件覆盖证书与互补条件。尚未供应的全局锐常数、强预紧性、无穷见证强对偶或最优覆盖权，不因这些局部结构被解除。

## 1. 固定 Lebesgue 背景下 \(A=9\) 的自足证明

稿 §1 保留分母 \(4r^2\)、接收面积 Lebesgue、输入有限正测度三者的区别。任一严格水平集紧子集 \(K\)，每点有一枚原严格见证；方体内部含其中心，故可从这些内部取有限覆盖。

按半径降序，仅保留与此前保留闭方体都不相交者。每个被删方体与某个较大保留方体相交，若半径为 \(r_i\le r_j\)，中心距 \(\le r_i+r_j\)，故其整个闭方体落在 \(3Q_j\)。因此三倍保留族覆盖原有限族的整个并，进而覆盖 \(K\)；不是仅把原有限中心被覆盖误称为 \(K\) 被覆盖。

保留闭方体真正不交，所以即使原子位于边界也有
\[
\lambda|K|\le9\lambda\sum_j4r_j^2
<9\sum_j\mu(Q_j)\le9\mu(\mathbb R^2).
\]
对有界紧集内逼近得到严格弱型，非严格版再以阈值从下逼近。输入函数取 \(|f|\,dy\) 即得固定 Lebesgue 版本。没有调用任意 Borel 背景模型的常数 4。

一般输入测度时，固定半径平移质量只保证 Borel，不能假称连续。稿末已另证 strict halo 开：若 \(\mu(Q(x,r))>4\lambda r^2\)，小增半径至 \(r+h\) 仍使原质量 \(>4\lambda(r+h)^2\)，近中心的增大方体含原方体。可数有理半径从上逼近与有限测度从上连续性又给 Borel supremum。这同时使内正则及正阈值有限面积合法。

## 2. 低端质量、高端零与输入类别

稿定理 1 的量词严格限固定非零 \(f\in L^1\)。

对 \(g\ge0\) 支撑于 \(Q(0,R)\)、质量 \(m_0>0\)，设 \(s_\lambda=\sqrt{m_0/(4\lambda)}\)。当 \(\|x\|_\infty<s_\lambda-R\)，半径 \(\|x\|_\infty+R\) 捕获全部质量并给严格平均。另一方面在盒外，任何碰到支持的半径至少 \(\|x\|_\infty-R\)，平均至多 \(m_0/[4(\|x\|_\infty-R)^2]\)。两类水平集面积夹在 \(4(s_\lambda-R)_+^2\) 与 \(4(s_\lambda+R)^2\)，乘 \(\lambda\) 后趋 \(m_0\)。左侧使用开端严格条件；稿件没有把临界闭边的等号误当严格见证。

一般 \(f\) 先取绝对值，用空间截断 \(g_R\) 与质量尾 \(\varepsilon_R\)。阈值分裂在 strict 与 nonstrict 各自成立：
\[
F_f^\diamond(\lambda)\le
\frac{F_{g_R}^\diamond((1-\eta)\lambda)}{1-\eta}
+\frac{A\varepsilon_R}{\eta}.
\]
先 \(\lambda\downarrow0\)，再 \(R\to\infty\)，最后固定 \(\eta\downarrow0\)；下侧由单调支配给质量。因此两类低端都为 \(\|f\|_1\)，参数顺序没有交换不受控极限。

高端取 \(g_B=\min(|f|,B)\)、\(h_B=(|f|-B)_+\)。因 \(Mg_B\le B\)，非严格水平集包含于 \(\{Mh_B\ge\lambda-B\}\)；弱型给 \(\lambda>2B\) 时 \(F_f^\ge(\lambda)\le2A\|h_B\|_1\)。先 \(\lambda\to\infty\)，再 \(B\to\infty\)，得零。strict 是其子集。

单原子 \(m\delta_0\) 的每个严格接收面积恰为 \(m/\lambda\)，高端仍是 \(m\)，故不能把该高端零照搬到有限原子或一般 Borel 输入。冻结稿明确作了这一类别区分。

## 3. 阈值半连续及有限正阈值的达到条件

正阈值有限面积使所有从上集合连续性合法。对从相应侧趋近的阈值，稿表 (9) 四个极限均正确：strict 右极限为 strict，左极限为 nonstrict；nonstrict 左极限为 nonstrict，右极限为 strict。乘以阈值仍保这些极限。

因此 \(F^>\) 是右连续、下半连续；\(F^\ge\) 是左连续、上半连续。\(f=\mathbf1_Q\) 在 \(Mf=1\) 的正面积 plateau 上确给 strict 版非上半连续例子。

以 \(s\uparrow\lambda\) 得每个固定函数
\[
\sup_\lambda F_f^>(\lambda)=\sup_\lambda F_f^\ge(\lambda).
\]
该恒等式不声称固定阈值集合相同。若 \(N(f)>\|f\|_1\)，第 2 节端点把所有充分近最大阈值限制于一个正紧区间；非严格上半连续性给 finite argmax。strict finite max 的充要条件正是某个 nonstrict argmax 的 plateau 零测：若 strict 值为 \(N\) 而 plateau 正面积，同阈值 nonstrict 值会超过 \(N\)，与两 supremum 相同矛盾。

稿 (13)–(15) 的真实 \(L^1\) 双 blob 提供独立 \(C_2>1\)，不依赖外部数值基线。每 blob 半宽 \(1/100\)、质量 4，单见证 \(99/100\)、双见证 \(139/100\) 的平均分别为 \(10000/9801\)、\(20000/19321\)，都严格大于 1。两个单接收方形的并面积 \(9457/1250\)，双接收矩形只在上下条补 \(86/125\)，总面积 \(10317/1250\)、质量 8，故比值 \(10317/10000>1\)。双矩形的水平投影在单方形投影并内，故新增面积确只计上下条，没有重复面积。

所以任何真的 \(L^1\) 极值函数必有有限正**非严格**等号阈值。这是必要结构；不是该函数存在性的证明。

## 4. \(\sqrt{N}\) 连续模量与强收敛传递

对几乎处处两边有限的点，
\(|Mf-Mg|\le M(f-g)\)；有限性由弱型给出。对 \(0<\eta<1\)，strict 阈值分裂得稿 (17)：
\[
N(f)\le\frac{N(g)}{1-\eta}+\frac{A\|f-g\|_1}{\eta}.
\]
设 \(U=N(g),V=A\|f-g\|_1\)。精确恒等式
\[
\frac U{1-\eta}+\frac V\eta-(\sqrt U+\sqrt V)^2
=\frac{[\sqrt U\,\eta-\sqrt V(1-\eta)]^2}{\eta(1-\eta)}
\]
给其最小值；若 \(U,V\) 为零，以端点极限处理。交换 \(f,g\)，得到
\[
|\sqrt{N(f)}-\sqrt{N(g)}|\le\sqrt{A\|f-g\|_1}.
\]
没有把 weak \(L^1\) 准范数当作常数一的三角范数。

质量 1 的极值序列若实际强 \(L^1\) 收敛，极限质量仍为 1，以上模量给 \(N(f)=C_2\)。只弱测度收敛、tightness 或质量不逃逸都没有被升级为强收敛。

移动阈值 (19) 由 \(\lambda\pm2\epsilon\) 包含及
\[
|\{|Mf_j-Mf|>\epsilon\}|\le A\|f_j-f\|_1/\epsilon
\]
得到。先 \(j\to\infty\)，再 \(\epsilon\downarrow0\)，给 strict 集的 liminf 下界与 nonstrict 集的 limsup 上界。正阈值面积有限，乘 \(\lambda_j\to\lambda>0\) 后得到强 \(L^1\times(0,\infty)\) 的联合 nonstrict 上半连续、strict 下半连续。

若强 \(L^1\) 收敛的非零输入同时带近极值阈值，阈值不能有趋零或无穷子列：固定极限 \(f\) 与误差 \(f_j-f\) 再作上述 split，端点子列分别只有质量及零的 limsup；\(C_2>1\) 排除它们。抽取有限正极限阈值后，联合上界与全局 nonstrict 常数才给等号，plateau 未检验前不能换成 strict。

## 5. 新增的强紧性等价表述

稿 §4.2 是合法条件等价，不是全局存在性断言：

- 若质量 1、阈值 1 的 strict 近极值序列有强 \(L^1\) 子列，则 \(N(f_j)\ge F_{f_j}^>(1)\to C_2\) 且 \(N(f_j)\le C_2\)，由第 4 节得到极值极限。
- 若极值 \(f\) 存在，设非严格最大阈值为 \(\lambda_*\)，取 \(s_j\uparrow\lambda_*\)。阈值左极限给 \(s_j|\{Mf>s_j\}|\to C_2m\)。归一化
  \(g_j(x)=s_j^{-1}f(\sqrt{m/s_j}\,x)\)
  确有质量 1、阈值 1 接收面积趋 \(C_2\)。其振幅和空间缩放参数趋有限正值，\(L^1\) 膨胀连续性给强收敛到同样归一化的 \(f\)。

这里所需膨胀连续性可先对 \(C_c\) 函数证明，再以 \(L^1\) 逼近；归一化算子对误差的 \(L^1\) 范数恰为原误差除以 \(m\)，故逼近在参数上统一。尚未供应任意近极值序列的此类强预紧性。

## 6. 固定有限见证 LP 的完备性与强对偶

闭 \(Q_i\)、\(r_i>0\)、\(b_i=4r_i^2>0\) 固定。非空 incidence cell \(A_J\) 包含所有确切属于 \(J\) 中方体、且不属于其它方体的点，Borel 且有界。输入并集外的质量不帮助约束，可先删除。

measure 版允许全部非空 cells；以任一点承载 \(m_J\) 精确还原其全部 closed-square memberships。\(L^1\) 版只能用正面积 cells，以 \(m_J|A_J|^{-1}\mathbf1_{A_J}\) 精确还原。零面积边、点不能借密度承载；也不能把 cell 的闭包当成其原 membership 后删掉约束。

两版原问题均可行且质量可限于 \(\sum b_i\)。有限非负变量的该次水平集紧，所以 nonstrict minimum 存在。每行至少有一个允许列：measure 因 \(Q_i\ne\varnothing\)，density 因 \(Q_i\) 正面积且由有限 cells 分割。故 dual 的每个 \(z_i\le1\)，其可行集紧。

强对偶的分离论证符号已逐项核对。有限生成锥
\((Am-s,\mathbf1^\top m+t)\) 闭；对 \(p<P_*\) 分离得到
\(\alpha\le0,\beta\ge0,A^\top\alpha+\beta\mathbf1\ge0\)，且 \(\alpha\cdot b+\beta p<0\)。
若 \(\beta=0\)，与原可行性矛盾；归一化 \(z=-\alpha/\beta\) 便得 dual 可行且 \(b\cdot z>p\)。令 \(p\uparrow P_*\)，配合 weak dual 与 dual 紧性即得 (25) 两边达到且相等。无数值求解器假设。

固定输入支持的矩阵版本也准确：每行至少一个 1 时可行、强对偶达；若某行全零，primal 不可行而 dual 沿该行无界。只在候选输入支持核对负载不能给所有输入证书。

finite strict 约束的 inf 与 nonstrict min 相同，缩放 \(1+\epsilon\) 即给一侧；任何 strict 可行输入有有限正比例 slack，可稍缩整体质量，故 strict inf 不达。此处的优化变量是固定有限见证的输入质量，不能把其不达移植成全局 \(N(f)/\|f\|_1\) 不达。

稿 (27) 的相邻面积 1 闭方形给 measure opt 1、density opt 2。我的独立核使用边长 2 的缩放版本，得到 4 与 8，二者完全相容。公共零面积边上的原子能同时满足两约束；density 不能。对偶 \(z_1=z_2=1\) 只违反共享边上的 pointwise 负载，故合法 density 证书不是合法任意 measure 证书。

## 7. 中心权测度、所有 \(y\) 与闭边 atoms

源 §6 允许正 \(\sigma\)-有限 Borel 中心权 \(\nu\)，比有限权略广；该扩展在此证明中合法。Borel \(E,r\) 给联合 Borel 核，输入 \(\mu\) 有限且 \(\nu\) 为 \(\sigma\)-有限，Tonelli 正确：
\[
\int_E\mu(Q_x)\,d\nu(x)=\int L_\nu(y)\,d\mu(y).
\]
原约束必须对 \(\nu\)-a.e. 中心合法；Lebesgue-a.e. 中心合法只对 \(\nu\ll dx\) 自动传递，不能任意改用 atomic \(\nu\)。

输入 measure 版的统一负载要对所有 \(y\) 控制；固定 \(\mu\) 可只在 \(\mu\)-a.e. 控制；\(L^1\) 版可只在 Lebesgue-a.e. 控制。这三者在闭边输入 atoms 上有严格差别。稿件已明确，没有用零面积免责。

中心权不必有密度。有限 LP 的 dual 恰是 labelled centers 上的原子中心权；若同一个坐标有多种半径，稿件也正确提升到标签空间或 \((x,r)\) 空间，而未把它们强行拼成单值 \(r(x)\)。

在可行负载和原约束下，(30) 强制 \(B_\nu\le\mu(\mathbb R^2)<\infty\)，即使 \(\nu\) 仅 \(\sigma\)-有限也没有无限积分相减。slack 恒等式
\[
\mu(\mathbb R^2)-B_\nu
=\int(1-L_\nu)\,d\mu+\int_E[\mu(Q_x)-b(x)]\,d\nu
\]
由有限的 Tonelli 积分直接展开；两项 a.e. 非负，等号恰对应两个 a.e. 饱和条件。若 \(\nu\ne0\) 且原见证 \(\nu\)-a.e. strict，第二项确实严格正。

但取得几何目标 \(B_\nu\ge|E|/K\) 是另一供货；稿件没有从有限 LP 强对偶自动推出无穷见证强对偶、最优 \(\nu\)、统一 sharp cover 或全局极值不可达。它还指出阈值 supremum 可从下逼近 nonstrict plateau，故单个 finite strict 证书的 slack 不自动给统一间隙。这些量词隔离正确。

## 8. Melas 引文范围独立回源

已直接打开 [arXiv 原 PDF](https://arxiv.org/pdf/math/0311452) 与 [Annals 题录页](https://annals.math.princeton.edu/2003/157-2/p08)，并将本地同文 PDF 的第 4 页实际渲染为 melas650.png 后查看。印刷页 650 的 Theorem 3 限于正 Dirac 线性组合及每个阈值的严格不等式；同文印刷页 649 将该用语定义为有限 \(n\) 个正原子。Theorem 3 前一段把 \(L^1\) 或更一般有限测度的等号问题另行提出。因此冻结稿没有把一维 finite-atom 结论当作二维 \(L^1\) 或一般有限 Borel 输入不可达结论。

本报告只核该引用的准确范围，不以 Melas 结论证明本模块的新定理。Annals 页码、年份与作者也与原 PDF 一致。

## 9. 实际计算、可复现文件与最终状态

作者配对代码与原始输出已实读：

- finite_checks.wl/log，Wolfram 15 ARM，PID 96736：面积 1 接触方形的 measure/density primal-dual 为 1/2，平方优化残差 0。
- two_blobs.py/log/json，Python 3.15 Fraction，PID 96760：真实内接收矩形、strict 平均及 \(10317/10000\) 精确回归。

本代理新实际运行并保留完整输入／输出：

- finite_interfaces.wl/log，Wolfram 15 ARM，PID 97283，exit 0：缩放接触例 4/8、阈值 split 的平方恒等式残差 0、紧支撑两端面积极限均为 \(m\)。
- two_blob_lower.py/json/log，Python 3.15 Fraction，PID 98542，exit 0：独立同一真实 \(L^1\) 双 blob 的几何不等式和面积归一化，所有断言通过。

这些有限输出不证明 weak 覆盖、所有阈值、强收敛传递、LP incidence 完备性或连续 cover 目标；对应推导已在第 1–7 节独立给出。原主稿冻结哈希在审核后重新核验未变，文件清单另见 freeze_check.json。

最终：没有发现必要数学修补。可将冻结 94fd2370... 的上述结构命题作为已重建供货；全局极值存在／不存在、锐值、无穷见证证书存在性和任何强预紧性仍保持未证。
