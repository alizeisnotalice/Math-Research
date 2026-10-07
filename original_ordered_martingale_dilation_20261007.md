# 原 ordered 核的非齐次 martingale dilation：同源停止部分与投影重叠余项

2026-10-07。固定原物理尺度及 c>0，完整同一来源。不改主账或旧文件。本轮证明了原 K_r 的准确非齐次对称 Markov dilation，并将 Ω 外选尺度的输出精确分为质量 ≤W 的辅助停止部分与非负 projection-overlap 余项。余项尚未支付，故**没有证明新的 weak(1,1) 或 actual geom 合同**。特别没有把条件期望的 weak 范数当作整体 weak 范数。

## 1. 查重与方法范围

已读 [ordered_propagation_endpoint](ordered_propagation_endpoint_20261007.md) §3：旧表示使用固定 Q_c=G_c^{1/2} 的二步投影以及中间 commuting-projection semigroup；已指出外层条件期望不保 weak L¹。本轮的新增表示是随 r 有真实 Markov 增量的 P_r=K_r^{1/2}，不是重复旧固定 Q_c 表示。

应用并读取 [E04](/Users/zhengzhihao/.codex/skills/math-e04-snell-envelope-optimal-stopping/SKILL.md) 与 [I04](/Users/zhengzhihao/.codex/skills/math-i04-random-carleson-occupancy-llogl/SKILL.md) 的 method、provenance、cube-interface。只用有限时间可积 martingale 的条件期望/停止质量身份；随机停止允许未停收益为 0，未假定停止总质量为 1。没有引用这些 skill 中尚未认证的 LlogL 嵌入或连续无限时域 Snell 定理。查到 Rota 原论文的出版链接，但 AMS 原 PDF 获取返回 403，未将未读定理作为证据；下面所需有限网表示直接证明。

## 2. 原 P_r 及真正的非齐次增量

原一维轴核为 G_c=(I+cB)^{-1}，B 的 Fourier symbol 为

\[
B(v)=\frac v{\log(1+v)}-1
=\int_0^1[(1+v)^t-1]dt,\qquad v\ge0.
\tag{1}
\]

右式说明 B 是 Bernstein 函数：0<t<1 时 t(1+v)^{t−1} 完全单调；也可用其正 Lévy 表示建立 Brownian 从属概率核。固定一个坐标的 q=cB(ξ_i²)，原符号为

\[
\widehat K_r(\xi)=\prod_i\frac{1+(1-r)q_i}{1+q_i},\qquad
\widehat P_r(\xi)=\prod_i\sqrt{\frac{1+(1-r)q_i}{1+q_i}}.
\tag{2}
\]

对 0≤r<s≤1，记真实增量 V_{s,r}=P_s/P_r。其每坐标指数为

\[
\psi_{s,r}(q)=\tfrac12\{\log(1+(1-r)q)-\log(1+(1-s)q)\}.
\tag{3}
\]

置 A=1−r、D=1−s（D=0 取极限），有正 Lévy 表示

\[
\psi_{s,r}(q)=\frac12\int_0^\infty
(1-e^{-qt})\frac{e^{-t/A}-e^{-t/D}}t\,dt.
\tag{4}
\]

因此 ψ_{s,r}(cB(ξ_i²)) 由原 B 的 Markov 从属产生对称概率卷积增量；tensor 乘积同样保正、保质量。没有仅由乘子非负来断言 Markov。增量满足

\[
V_{t,s}V_{s,r}=V_{t,r},\quad P_s=V_{s,r}P_r,\quad
P_0=I,\quad P_r^2=K_r.
\tag{5}
\]

对 r<1，其真实生成元为

\[
\mathscr A^P_r
=\frac12\sum_i\frac{cB_i}{I+(1-r)cB_i}
=\frac1{2(1-r)}\sum_i(I-G_{c(1-r),i}),
\qquad \partial_rP_r=-\mathscr A^P_rP_r.
\tag{6}
\]

K_r 的生成元为 2𝒜^P_r。这正是原 softness 改变的非齐次轴过程，未换成 reset 或固定 G_c 过程。r=1 的概率核由 (2)–(5) 存在；不在该端点假定有 bounded generator。以下先用有限 r 网，包含端点也可直接使用 Markov 增量，无需无限时间可选抽样。

## 3. Lebesgue σ 有限路径测度与 reverse martingale

取 0=r_0<r_1<⋯<r_J≤1。构造路径测度

\[
d\mathfrak m(x_0,\ldots,x_J)
=dx_0\prod_{j=1}^JV_{r_j,r_{j-1}}(x_{j-1},dx_j).
\tag{7}
\]

每个 X_j 的边缘都是 Lebesgue 测度。这是 σ 有限测度，不是“均匀概率分布于 Rⁿ”。**不能**先截断 X_0 到一个有限盒并归一化、再沿用 stationarity 或相同的 reverse 条件核；截断会改变后验。本稿直接在 (7) 下使用条件期望和可积函数，全部积分由保质量/Tonelli 合法化。

令 𝔽_j=σ(X_j,…,X_J)，为递减滤过。对完整 f≥0、f∈L¹，W=∫f，定义

\[
M_j=\mathbb E_{\mathfrak m}[f(X_0)\mid\mathcal F_j]
=P_{r_j}f(X_j),\qquad
\mathbb E_{\mathfrak m}[M_j\mid X_0]=K_{r_j}f(X_0).
\tag{8}
\]

第一式由对称性、Markov 性与 (5) 得到；第二式是两次 P_{r_j}，不是一次。M_j 是非负 reverse martingale，∫M_j d𝔪=W。逆序 j=J,…,0 后是有限 forward martingale，其终端为 f(X_0)。

若 Z=max_j M_j，则有限首次越阈停止给 α𝔪{Z>α}≤W。证明只需越阈停止事件与对应滤过可测、以及停止收益积分等于该事件内终端 f 的积分；在 σ 有限测度下同样成立。这个路径 weak 界不等于接收点 K_*f 的 weak 界，因为 (8) 还含外层 conditional projection。

## 4. Ω 外 selector 的精确 covariance

现专取原饱和障碍来源 f=ν_b，支撑于 Ω，满足

\[
Su=\nu_b-\mu_b,\quad u\ge0,\quad\Omega=\{u>0\},\quad
0\le\mu_b\le\kappa,\quad W=\int\nu_b=\int\mu_b,
\quad S=\sum_i(I-G_{c,i}).
\tag{9}
\]

完整输入的 good 来源若存在仍须单独保留；本节不把总 ν 当 ν_b。κ>0 固定。取 finite-grid 外域超阈集合

\[
A=\{x\in\Omega^c:\max_jK_{r_j}\nu_b(x)>\kappa\}.
\tag{10}
\]

将 A 按原最大响应的最小索引 tie rule 分为可测不交 A_j；r_0=0 在外域没有贡献。∫K_rν_b=W 给 |A|≤JW/κ，只为证明有限性，不作为目标费。a_j=1_{A_j}(X_0) 已知于 reverse 过程的**终端**，一般不是合法的 reverse 停止标签。其条件投影为

\[
\pi_j=\mathbb E[a_j\mid\mathcal F_j]
=P_{r_j}1_{A_j}(X_j),\quad 0\le\pi_j\le1.
\tag{11}
\]

真实选中输出精确等于

\[
T_A=\sum_j\int_{A_j}K_{r_j}\nu_b,dx
=\mathbb E_{\mathfrak m}\sum_ja_jM_j
=\mathbb E_{\mathfrak m}\left[f(X_0)\sum_j\pi_j\right].
\tag{12}
\]

由于 f 支撑 Ω 而 A_j⊂Ωᶜ，a_j f(X_0)=0 几乎处处，所以其条件期望也为零。于是

\[
\operatorname{Cov}_{\mathfrak m}(a_j,f(X_0)\mid\mathcal F_j)
=-\pi_jM_j,
\quad
T_A=-\sum_j\mathbb E_{\mathfrak m}
\operatorname{Cov}(a_j,f(X_0)\mid\mathcal F_j).
\tag{13}
\]

这不是把源、接收点或后验当独立。外域支持恰使原 terminal product 为零；输出全由条件投影中的负 covariance 产生。原分区 Σa_j≤1 不推出 Σπ_j≤1，因为各标签在不同滤过上投影。

## 5. 合法停止源 ≤ν_b，加上未付 overlap

按降序 j=J,…,1 观察 𝔽_j，并在每个阶段揭示一枚额外独立 Uniform(0,1) coin；尚未停止时，以概率 π_j 停止。给定原路径的停止权重为

\[
b_j=\pi_j\prod_{k>j}(1-\pi_k),\qquad
\sum_jb_j=1-\prod_j(1-\pi_j)\le1.
\tag{14}
\]

各 π_k、k>j 在 𝔽_j 可测，因此 b_j 为适应权；加 coin 后确为 reverse 顺序下的随机停止规则。这里不将一个 future 事件当作 predictable 指示；离散阶段的 coin 在已观察当前 π_j 后才使用。允许不停止，未停止收益为 0。有限条件期望身份给

\[
T_{\rm stop}=\sum_j\mathbb E[b_jM_j]
=\mathbb E\left[f(X_0)\left(1-\prod_j(1-\pi_j)\right)\right]
\le W.
\tag{15}
\]

对应完整来源空间中的正停止源是

\[
d\nu_{\rm stop}(y)=f(y)
\mathbb E\left[\sum_jb_j\mid X_0=y\right]dy\le d\nu_b(y).
\tag{16}
\]

这是一份合法 source-once 的**辅助 dilation 停止源**，不是重新定义原 actual FIRST。真正选中输出尚有精确非负余项

\[
\boxed{T_A=T_{\rm stop}+D_A,\quad
D_A=\mathbb E\left[f(X_0)
\left(\sum_j\pi_j-1+\prod_j(1-\pi_j)\right)\right]\ge0.}
\tag{17}
\]

点态把 π_j 视为独立 Bernoulli 概率的数值，可得 elementary union bound/Bonferroni 身份

\[
0\le\sum_j\pi_j-1+\prod_j(1-\pi_j)
\le\sum_{i<j}\pi_i\pi_j,
\quad D_A\le\sum_{i<j}\Gamma_{ij},\quad
\Gamma_{ij}=\mathbb E[f(X_0)\pi_i\pi_j].
\tag{18}
\]

这里 Bernoulli 仅用于合法独立辅助 coins 的代数，没有假定原捕获后 mask 独立。若能对原 selector 证明 D_A≤C_nW，则由 κ|A|≤T_A 得整体 finite-grid weak 常数 ≤1+C_n。若常数与网无关，原 K_r 是有限个 G_Aν_b 的 polynomial：选共同满测集上全部系数有限的版本，r↦K_rν_b(x) 连续，故 dense-grid 穷尽可给全 r family。**当前没有 D_A 的所需预算。**

## 6. Pair overlap 的原 P 转移核空间表示

写 P_i=P_{r_i}、V_{j,i}=P_j/P_i，为 (3)–(5) 的真实 Markov 增量，i<j。Markov 性精确给

\[
\mathbb E[\pi_i\pi_j\mid X_0=y]
=P_i\big[(P_i1_{A_i})\,(V_{j,i}P_j1_{A_j})\big](y)
=P_i\big[(P_i1_{A_i})\,(P_iV_{j,i}^{\,2}1_{A_j})\big](y).
\tag{19}
\]

因此需要支付的是具体的同源三因子空间积分

\[
\boxed{\Gamma_{ij}
=\int_{\mathbb R^n}(P_i\nu_b)(x)
 (P_i1_{A_i})(x)
 (P_iV_{j,i}^{\,2}1_{A_j})(x)dx.}
\tag{20}
\]

V_{j,i}²=K_{r_j}/K_{r_i} 是原完整 ordered 增量，不是固定 G_c 或自由 Gaussian 替代。式 (20) 保留同一原 ν_b、共同 pivot x 与原接收分区，没有将 hard 来源替换为 soft 来源。

A_j 不交只推出 Σ_j1_{A_j}≤1；各 j 的转移 V_{j,i}² 不同，不能免费断言 Σ_jV_{j,i}²1_{A_j}≤1。P_iν_b 也不由 μ_b≤κ 自动封顶。式 (20) 是三因子 occupation，不是固定 Gram 平方，不能直接以 Bessel 或谱正交支付。

精确 D_A 还含更高重叠的取消：对索引 j_1<⋯<j_k，其 Γ_{j_1,…,j_k}=E[f(X_0)∏_lπ_{j_l}] 由 (7) 的同一源路径和真实 V 链表示，且

\[
D_A=\sum_{k\ge2}(-1)^k
\sum_{j_1<\cdots<j_k}\Gamma_{j_1,\ldots,j_k}.
\tag{21}
\]

这在有限网是准确有限和。用 pair 上界是一个明确的正放大，会丢掉 (21) 的高阶取消；本稿未登记它为已可支付的新强合同。

## 7. 饱和障碍与原停止路径目前提供什么

令

\[
C_A^{\rm def}(y)=\mathbb E\left[
\sum_j\pi_j-1+\prod_j(1-\pi_j)\mid X_0=y\right].
\]

0≤C_A^{def}≤Σ_jK_{r_j}1_{A_j}，所以它有界、可积，且 ∫C_A^{def}≤|A|。用原固定 S 的对称性及 (9)，得到

\[
D_A=\int\mu_b C_A^{\rm def}
       +\int u\,S C_A^{\rm def}
\le\kappa|A|+\int u\,S C_A^{\rm def}.
\tag{22}
\]

仅代入 μ_b≤κ 给系数 1 的 κ|A|，原不等式 κ|A|≤W+D_A 因而不能吸收；S C_A^{def} 的符号/来源预算也未证。这明确定位了障碍仍需控制的 selector potential/commutator 项，没有通过自由 potential energy cap 关闭。

也可对 P 或 K 的真实 forward Markov 过程在 Ω 外首次退出。其正退出源总质量 ≤W；对 K 的增量 V² 使用一次退出分解，正好回到旧稿尚未支付的原时间依赖 continuation maximal。对 P 的一次退出后再投影 P_r，则还包含 P_r 对未退出 live 来源的外域输出。这些重写不解决 (17)，不能用质量 ≤W 就继承完整传播的整体 weak 界。

标准截断 martingale 论证仅给 LlogL：分出 f≤t/2 后，𝔪{Z>t}≤(2/t)∫f1_{f>t/2}。因 K_*f≤E[Z|X_0]，对外域 (10) 可得

\[
\kappa|A|\le4\int f\log_+(4f/\kappa).
\tag{23}
\]

它不支付任意 L¹ 的 W 端点。饱和来源也没有免费 entropy/W cap；这是已有尖峰边界，本轮仅引用其适用限制，不登记为新反例突破，不重复实验。式 (17)–(20) 比 (23) 更具体：真正的新待估项是原转移核下的 source-weighted selector overlap。

## 8. 与截帽 reference 吸收的关系

[capped_reference_absorption](capped_reference_absorption_20261007.md) 只需一个同源整体非负 B 的独立 weak 常数，不要求 ∫B 强范数。若将来 (17) 的原 D_A 预算对所有阈值与网成立，得到整体原 K_*ν_b 的 weak 界，则在**另外核准**实际交通 T≤B+e、T≤cλ 及剩余误差同源预算后，可用该吸收引理。没有先把各路径或各参数的 conditional weak 界混合成 B 的整体 weak 界。

目前路径 Doob weak 界尚有 (17) 的投影损失；原 actual 交通另含 moving hard 物理尺度、soft winner、来源对和 FIRST 门，仍需保持其完整桥接。辅助停止源 (16) 不更换这些原合同。

## 9. 三轮有理代数守卫与范围

执行前保存同前缀 `_registration_20261007.json`。新 guard 一次执行，5956 项 Fraction 精确检查全部通过；不导入或重跑旧 spike、entropy、TV 脚本。

|round|原 symbol tensor n|有限代数滤过层数 J|terminal atoms|检查范围|
|---:|---:|---:|---:|---|
|1|4|3|16|BF 增量导数符号、tensor 微分、tower/covariance、停止薄化和 overlap|
|2|16|6|128|同上|
|3|64|9|1024|同上|

原 symbol 部分对登记的非负有理 q、r 网检查 (2)–(6) 的精确代数；不以有限 q 检查代替 (4) 的全部 Bernstein/Markov 证明。有限滤过部分是固定完整 terminal 函数与 exterior selector 的嵌套 partition 模型，只核条件期望和 (13)–(18) 的普适代数，**不是原 c=1 空间核样本**，也不是饱和障碍、FIRST/geom 或维数阶数压力。它没有替换原 Gi 来验证任何空间弱端点。

保存结果含完整 terminal f、receiver labels、M_j、π_j、b_j 和精确 stopped/deficit/pair 数值，可以从 Fraction 数据重构；无浮点判正、容差、Monte Carlo 或置信区间主张。

本轮已付的是 (15) 的辅助同源停止部分；(17) 的原 projection-overlap、(20) 的同源三因子占用及 actual 回代仍未闭合，主账不变。
