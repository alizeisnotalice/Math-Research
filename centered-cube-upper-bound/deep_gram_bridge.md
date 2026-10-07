# 实际深重父组的有符号平方：Gram 结构与不能删去的负项

2026-10-07。解析审计；未运行脚本、求解器或数值试验。仅新增本文件。

**本次结果。** 已建立同一来源参考上的六个真正 Gram 核及其精确有符号分解。非负、输出依赖的节点权本身不破坏 Gram；但父减子、来源对门以及删掉负 discrepancy 项这三步不能混为一谈。任何非零实际 LCA 交通核在固定种子的来源空间上都不是 PSD；若仍保留原远对门，种子平均核也不是 PSD。这是原实际核的条件结构命题，不是构造了违反平方根目标的新输入。未得到深重交通的新平方根费用。

## 1. 已读来源与冻结接口

使用 math-g03-directional-capture-tree-quadratic-audit，已读 SKILL、provenance、method、cube-interface。没有启用其未定义 TC-A14 私有接口，也没有调用 Bloch/Bergman 来源定理。下述论证完全基于实际来源映射和原森林。

已完整读取 occupation identity 与 heavy_mass_fee_review_20261006.md，包括浅分支已计算的同源对偶证书；读取概率接口阶段证明.tex 2690–2762 行的固定种子森林、半开格、唯一 LCA、两端落于不同孩子、全部原资格及剩余系数；沿用已读 10 月 6 日 geom 提示词的实际归一化、条件 LCA 捕获与早第一跳定义。

更新主账与新增条件 Bessel 按“提供的接口前提”使用；没有补造缺失的 10 月 5/6 日上游证明。保留浅分支已审费用
\[
R_{\le D}\le J(1+D)W,\qquad D=\sqrt n/J.
\]
本稿只研究精确补支 \(R_{>D}\)。质量深度是审计文件中
\[
d_P=\min\left\{n\log(b/a),\log\frac{h_P}{2\tau\lambda a^n}\right\},
\qquad h_P>2\tau\lambda a^n,
\]
不是 martingale 时间、信息量或来源熵。

## 2. 同一来源上的实际软硬响应

设 \(W>0\)、\(\nu=\mu/W\)，\(w=d\mu_{\rm hi}/d\mu\in[0,1]\)。冻结原输入产生的 \(E,x,\sigma(x),L_s(x),R_h(x),q,u(x)\)、出生 \(B\)、seed \(\vartheta\) 和所有门。软硬完整行归一化保留原值，\(u=0\) 行交通为零。

在共同参考 \(\nu\) 上定义
\[
\alpha_x(y)=\frac WqP_{\sigma(x),L_s(x)}(x-y)w(y),\qquad
\beta_x(y)=\frac W{u(x)}h_{R_h(x)}(x-y)w(y)\mathbf1_B(y).
\tag{1}
\]
于是 \(d\bar p_x=\alpha_xd\nu\)、\(d\bar r_x=\beta_xd\nu\)。对父格 \(P\) 及其孩子 \(C\)，
\[
a_{x,P}=\alpha_x\mathbf1_P,\quad b_{x,P}=\beta_x\mathbf1_P,\qquad
s_P=\int a_{x,P}d\nu,\quad t_P=\int b_{x,P}d\nu.
\tag{2}
\]
这些都是同一份完整来源的函数；没有把高子源再归一化。

无额外资格的原 LCA 块核为
\[
K_{x,P}(y,z)=a_{x,P}(y)b_{x,P}(z)
       -\sum_{C\in\operatorname{ch}(P)}a_{x,C}(y)b_{x,C}(z).
\tag{3}
\]
它恰等于 \(\alpha_x(y)\beta_x(z)\) 乘“两端在 \(P\) 的不同孩子”的指标，逐点非负。其在真实来源上的积分是
\[
c_P=s_Pt_P-\sum_Cs_Ct_C\ge0.
\tag{4}
\]

**实际早第一跳与余项门不可漏掉。** 原目标在 (3) 上还乘早期积分子核与完整软核的比值、原源对分数、删除、GOOD、score、owntrace、出生、严格 CP/GP、FIRST/全部未来约束等实际资格。记其父组归一化行交通为 \(\kappa_P\)。门不一致时只有
\[
0\le\kappa_P\le c_P,
\tag{5}
\]
以及相应树子块支配；不能写 \(\kappa_P=c_P\)。早期比值指已经积分第一跳时间、落点与全部 continuation 后的 \(Q_x^{\rm early}/P_{\sigma,L_s}\)，不把固定未积分历史密度比当作概率。

## 3. 六个真 Gram 核，与精确有符号分解

对固定 \(x,P\)，令 \(e_{x,P}=a_{x,P}-b_{x,P}\)，与质量深度 \(d_P\) 区分。代数对称化给
\[
\operatorname{Sym}(a_P\otimes b_P)
=\tfrac12(a_P\otimes a_P+b_P\otimes b_P-e_P\otimes e_P).
\tag{6}
\]
因此
\[
\begin{split}
\operatorname{Sym}K_{x,P}
=\tfrac12\big[&
a_P\otimes a_P-\sum_Ca_C\otimes a_C\\
&+b_P\otimes b_P-\sum_Cb_C\otimes b_C\\
&-e_P\otimes e_P+\sum_Ce_C\otimes e_C\big].
\end{split}
\tag{7}
\]
这里六种正求和分别是 Gram 核；**(7) 是这些 Gram 的有符号差，不是一个 Gram 核。** 对称化只处理代数，不交换原软硬资格。

设可测节点权 \(g(x,\vartheta,P)\ge0\)，暂只依赖 \(x,\vartheta,P\)，不依赖当前 \(y,z\)。取真实索引测度
\[
d\rho=\mathbf1_E(x)\frac{u(x)}Wdx\,d\Pr(\vartheta)\,d\#P.
\tag{8}
\]
例如第一种 Gram 有特征 \(\sqrt g\,a_{x,P}\)，孩子 Gram 有特征 \(\sqrt{g(x,\vartheta,P)}a_{x,C}\)。对任何使积分有限的真实测试 \(\varphi\)，
\[
\begin{split}
&\iint \varphi(y)\varphi(z)
 \int g a_{x,P}(y)a_{x,P}(z)d\rho\,d\nu(y)d\nu(z)\\
&\hspace{2em}=\int g\left|\int\varphi a_{x,P}d\nu\right|^2d\rho\ge0.
\end{split}
\tag{9}
\]
其它五种同理。输出依赖不改变此身份，参考始终是 \(\nu\)。可先在有限空间盒、有限格索引与有限软度截断上证明，再对非负平方项分别取单调极限；有符号相减仍须验证各项可积。

所以“receiver-dependent 门必然破坏 PSD”是错误表述。正确区别是：

- 非负 scalar 节点/行权可放进同一 Gram 的平方根特征；
- 同一来源上的因子 \(v(y)v(z)\) 可以同时乘特征，保 PSD；
- 一般来源对门 \(A(x,y,z)\) 对 Gram 核逐点乘 \(A\) 没有 PSD 保证；
- 父减子及其它带符号线性组合，即使每个成员真为 Gram，组合仍需单独认证。

若 pair 门本身是可测 PSD Gram 核，Schur 乘积可通过张量特征保 PSD；本题不能免费假定这一点。原“不同孩子”LCA 门已经不满足该充分条件。

## 4. discrepancy 差不是可以删掉的负平方

对任意可加节点数 \(v_P=\sum_Cv_C\)，记
\[
\mathcal D_P(v)=v_P^2-\sum_Cv_C^2.
\]
标量 (4) 准确满足
\[
c_P=\tfrac12\big[\mathcal D_P(s)+\mathcal D_P(t)
                         -\mathcal D_P(s-t)\big].
\tag{10}
\]
\(s,t\) 为正测度的质量，前两项非负。但 \(s-t\) 是有符号测度，\(\mathcal D_P(s-t)\) 可为负；最后一项因而可增加交通，不能以“减平方”删掉。

精确二孩子算术：取
\[
(s_1,s_2)=(9/10,1/10),\qquad(t_1,t_2)=(1/10,9/10).
\]
则 \(c_P=41/50\)，\(\mathcal D_P(s)=\mathcal D_P(t)=9/50\)，而 \(\mathcal D_P(s-t)=-32/25\)。故
\[
41/50\not\le\tfrac12(9/50+9/50)=9/50.
\tag{11}
\]
这仅是量词/符号的精确代数检查，不声称四个值已由完整 FIRST/早第一跳实例实现。下一节给原实际核本身的更强条件障碍。

## 5. 原实际非零 LCA 核的 PSD 障碍

### 命题 A：固定 seed 的实际来源核

冻结 seed \(\vartheta\)。令 \(C_{f,\vartheta}^{>D}(y,z)\ge0\) 是原实际深重源对核，归一化使
\[
R_{>D,\vartheta}/W=\iint C_{f,\vartheta}^{>D}d\nu d\nu.
\]
它仍只允许两端在该有限森林某父格的不同孩子，并有唯一 LCA；所有 FIRST、早第一跳、续传播与其它门保留。令 \(C_j\) 为该 seed 的最细半开格。于是
\[
C_{f,\vartheta}^{>D}=0\quad\nu\otimes\nu\text{-a.e. on }C_j\times C_j
\quad\text{for every }j.
\tag{12}
\]
若 \(R_{>D,\vartheta}>0\) 且核可积，就有两格 \(C_i,C_j\)、\(i\ne j\)，满足
\[
I_{ij}=\iint_{C_i\times C_j}
\operatorname{Sym}C_{f,\vartheta}^{>D}(y,z)d\nu(y)d\nu(z)>0.
\]
对**同一真实来源**取 \(\varphi=\mathbf1_{C_i}-\mathbf1_{C_j}\in L^2(\nu)\)。由 (12)，
\[
\iint\varphi(y)\varphi(z)\operatorname{Sym}C_{f,\vartheta}^{>D}(y,z)
d\nu(y)d\nu(z)=-2I_{ij}<0.
\tag{13}
\]
因此该实际对称核不是 PSD。没有采用抽象树样本或假设来源坐标独立；只是对任何已有正实际深交通的原核提取它自己的两块来源测试。若目标交通恰为零，则命题不声称存在违反点。

### 命题 B：种子平均，须保持确定远对门

不能把命题 A 中随 seed 变的测试直接拿去证明平均核不 PSD。原森林段还明确保留远对资格
\[
\|y-z\|_1>H_{\rm mark}>a.
\tag{14}
\]
若当前深补支继承 (14)，取一个**确定、不随 seed 变化**的边长 \(a/n\) 网格 \(D_j\)。同格来源对距离 \(\le a\)，所以平均核
\[
C_f^{>D}=\mathbb E_\vartheta C_{f,\vartheta}^{>D}
\]
在每个 \(D_j\times D_j\) 上为零。只要 \(R_{>D}>0\)，重复 (13) 得同一固定 \(\varphi=\mathbf1_{D_i}-\mathbf1_{D_j}\) 的负二次型。这个版本真正作用于 seed 平均后的同源核。

**资格守卫。** (14) 已在读取的旧实际森林段中出现；若更新余项没有继承它，必须先核验映射，此时只能报告命题 A，不能报告平均版本。最新资料缺失不能靠同名余项补齐此门。

命题 A/B 不否定所求预算。源核逐点非负，仍可采用 \(N+\mathrm{Gram}\) 对偶证书；只说明“把实际交通核自己变成正平方”这一步不可能。有效证书可以增加乘子、来源势、非负松弛及共同约束，使完整证书落入充分锥；证书的 Gram 不需要等于交通核。尤其不能把 (13) 的负方向误说成对偶证书必须填补的负方向：证书含的是 \(-C\)，该测试方向上的 \(-C\) 反而为正。

## 6. 带实际节点门的 telescoping：确切边界项

在任意有限根树上，对节点可加数 \(v\) 与节点权 \(g_P\)，逐节点收集平方系数得到
\[
\begin{split}
\sum_{P\ {\rm internal}}g_P\mathcal D_P(v)
={}&\sum_{Q\ {\rm root}}g_Qv_Q^2\\
&+\sum_{P\ {\rm nonroot\ internal}}(g_P-g_{\operatorname{par}(P)})v_P^2
-\sum_{C\ {\rm leaf}}g_{\operatorname{par}(C)}v_C^2.
\end{split}
\tag{15}
\]
标量乘积版本是
\[
\begin{split}
\sum_Pg_Pc_P
={}&\sum_Qg_Qs_Qt_Q\\
&+\sum_{P\ {\rm nonroot\ internal}}(g_P-g_{\operatorname{par}(P)})s_Pt_P
-\sum_{C\ {\rm leaf}}g_{\operatorname{par}(C)}s_Ct_C.
\end{split}
\tag{16}
\]
(15) 对 \(s,t,s-t\) 分别成立，再由 (10) 复原 (16)。这不是省略适应选择的普通 telescoping。

若只保留重捕获、质量深度门
\[
g_P=\mathbf1_{\{s_P>\tau,t_P>\tau,d_P>D\}},
\tag{17}
\]
则在共同阈值和固定父组质量的同一森林中该门向祖先保持：\(s,t,h\) 随祖先不减，故 \(g_P\le g_{\operatorname{par}(P)}\)。此时 (16) 非根项与叶项非正，确有
\[
\sum_Pg_P\kappa_P\le\sum_Pg_Pc_P
\le\sum_{Q\ {\rm top}}g_Qs_Qt_Q\le\sum_{Q\ {\rm top}}g_Qt_Q.
\tag{18}
\]
(15) 对 discrepancy 仍不是 PSD：符号测试 (13) 和算术 (11) 阻止把这一步称为 Gram 校准。若额外实际节点门不向祖先保持，(15)/(16) 的正变化项还必须另付；若原门依赖源对，只能先用 (5) 正扩大到 (17)，不能对 actual \(\kappa_P\) 宣称精确 telescoping。

## 7. 层数可避免的硬端粗上界；它不是新 Gram 费用

(18) 可纠正父组逐层收费的 \(J\) 因子。对固定 seed 的 top 格 \(Q\)，令
\[
h_Q=\mu_{\rm hi}(Q\cap B),\qquad
H_Q(x)=h_{R_h(x)}*(\mu_{\rm hi}|_{Q\cap B})(x).
\]
在 top 重门上 \(t_Q>\tau\)，沿审计文件的实际硬幅度 \(u>2\lambda\)，原赢家满足
\[
R_h(x)\le c_Q=\min\{b,(h_Q/(2\tau\lambda))^{1/n}\}.
\]
用审计中已证的硬包络积分 \(1+d_Q\) 与 (18)，
\[
R_{>D}\le\mathbb E_\vartheta\sum_{Q\ {\rm top,\ active}}h_Q(1+d_Q)
\le[1+n\log(b/a)]W.
\tag{19}
\]
top 格不相交，质量只收一次，没有 \(J\)。全部 actual 门只是正支配中被丢掉，不声称扩大核继承资格。

相比把所有深父组各用 \(h_P(1+d_P)\) 支付，(19) 去掉人工层数损失。不过它也可由整个实际 hard field 的固定窗口强包络直接得到，是**已有粗天花板**，不能登记为新 Gram 或平方根分支进展。倍尺度窗中仍是 \(O(n)W\)，不满足目标 \(\sqrt n\,n^{o(1)}W\)。本稿不会用 (19) 将最新主账标为闭合。

## 8. 最小未付转换与下一证书对象

本稿完成的是明确源函数、真假 Gram、带门边界与实际 PSD 障碍。未得到一条通过完整 FIRST/早第一跳门的新深费用式。

可继续寻找的对象是：在同一 \(\nu\) 上，将 (7) 的六 Gram **有符号差**与实际源对门、完整未来行约束、早出口 occupation 一起校准，使
\[
\begin{split}
&\frac{\eta(y)+\eta(z)}2-\operatorname{Sym}C_f^{>D}(y,z)\\
&+\frac12\int[F_\xi(y)\Lambda_\xi(z)+F_\xi(z)\Lambda_\xi(y)]d\rho(\xi)
=N(y,z)+G(y,z),
\end{split}
\]
其中 \(N\ge0\)、\(G\) 真为同源 Gram，且**实际计算的**势/乘子成本达到平方根阶或可吸收小行费用。命题 (13) 排除了把交通核本身称作 Gram 的捷径，并没有给出证书成本，也没有否定完整对偶证书。

occupation 文件中的 \(R_{\rm heavy}=\mathcal B_0+\mathcal D_f\) 可作为保留后续传播的参考差额，但这个恒等式同样不是费用改善，不能用它支付 (15) 的负 discrepancy 或未知列负载。

**最终判定：** 原 scalar LCA 极性分解成立；输出依赖 scalar 门下六个同源 Gram 可合法构造；父减子 discrepancy 不是 PSD，实际非零 cross-child 交通核也不是 PSD。一般源对资格门不能免费 Schur 到 Gram。深补支的平方根成本仍未证，主账不修改。
