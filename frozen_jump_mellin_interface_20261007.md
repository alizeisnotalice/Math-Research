# 固定 jump 的导数 Mellin 接口与高频缺口

2026-10-07。只新增本稿及同前缀 registration/guard/results。不修改已核 Abel、obstacle、总稿或 cube 主账；未重跑旧数据。

**终态：低频接口成立；冻结轴向核的统一高频 \(L^1\) 费用尚未证明。** 已证明的障碍与 Abel 工具给出低频 \(\log^2(e\Lambda)\) 最大 \(L^2\) 界。若补足本文 (11) 的高尾费用，单个冻结 semigroup 可得 \(O(\log(e\Lambda))\) 弱型界。逐 Poisson 跳层取绝对值再积分会发散，不能用于补足该缺口；此发散是方法失败，不是原高尾或 actual cube 的反例。

## 1. 阅读、查重及不变的对象

已读 [Spector–Stockdale, arXiv:2609.05377v1, §4](https://arxiv.org/html/2609.05377v1)，重点其 (4.1)–(4.6) 的高 Mellin 计算，及 §4 后段/§5 的低频转换。高频依赖 Euclidean heat kernel 的显式空间 Mellin 核及障碍 potential 的域外积分，不能由一般 fixed semigroup 或 subordination 的名称继承。本文不复述原文完整证明。

旧总 tex 1748 附近的高频是 cube 核的空间 Fourier 参数审计，5916 后是热流初始平均线性损失，均不是这里的 jump logtime 高频预算。已有 public_search_scope 仅记录阅读；旧 nonlocal_abel_absorption_review 与 fixed_jump_obstacle_interface 分别证明 fixed Abel 及单个 frozen 生成元的 density obstacle。本稿不把这些旧工具重新登记成空间进展。

固定对称有界率 \(L\ge0\)，\(H_t=e^{-tL}\)。采用前稿构造的一次分配：
\[
u\ge0,\quad u\in L^1\cap L^2,\quad \Omega=\{u>0\},\quad
\nu_{\rm b}=\nu\mathbf1_\Omega,\quad
\mu_{\rm b}=\nu_{\rm b}-Lu,\quad 0\le\mu_{\rm b}\le\kappa,
\quad \int\nu_{\rm b}=\int\mu_{\rm b}=W_{\rm b}.
\]
域内 \(\mu_{\rm b}=\kappa\)，\(\kappa|\Omega|\le W_{\rm b}\le W\)；good 输入 \(\nu_{\rm g}=\nu\mathbf1_{\Omega^c}\le\kappa\)，\(W_{\rm b}+W_{\rm g}=W\)。写
\(\sigma=\nu_{\rm b}-\mu_{\rm b}=Lu\)。不重定义来源或重新选择 FIRST。

## 2. 必须使用导数/分布 Fourier

非局部情形在域外有
\(\sigma=-\mu_{\rm b}=-\mathcal Ju\)，一般非零。因此
\[
F_{\sigma,x}(v):=H_{e^v}\sigma(x)\longrightarrow-\mu_{\rm b}(x)
\quad(v\to-\infty)
\tag{1}
\]
在相应强 \(L^2\) 意义下成立。它通常不是可积的普通 logtime 函数。有限链的 \(F_{\nu_{\rm b}}\) 在 \(v\to+\infty\) 也可留下 free 零谱常量。不能原样引用热论文“\(F_\sigma\in L^1(dv)\)”或其普通 Mellin 核 (4.1)。

定义
\[
G_\eta(v)=\partial_v H_{e^v}\eta=-e^vLH_{e^v}\eta.
\]
谱定理给
\[
\int_{\mathbb R}\|G_\eta(v)\|_2^2dv
=\tfrac14\|(I-P_{\ker L})\eta\|_2^2\le\tfrac14\|\eta\|_2^2.
\tag{2}
\]
因此使用 \(G_\eta\) 的 \(L^2\) Fourier，或等价地用
\(\widehat G_\eta(\xi)=2\pi i\xi\,\widehat F_\eta(\xi)\) 的分布解释。本文约定 \(\omega=2\pi\xi\)、Fourier 核 \(e^{-2\pi i\xi v}\)。

设 ergodic average \(\widetilde H_t=t^{-1}\int_0^tH_rdr\)，\(R_t\sigma=H_t\sigma-\widetilde H_t\sigma\)。积分微分恒等式给
\[
\widetilde H_t\sigma=(u-H_tu)/t\le0\quad\hbox{于 }\Omega^c.
\tag{3}
\]
这是完整非局部对象，域外负号保留。logtime 中
\[
R_\sigma(v)=\int_0^\infty e^{-a}G_\sigma(v-a)da,\qquad
\widehat R_\sigma(\xi)=\frac{\widehat G_\sigma(\xi)}{1+i\omega}.
\tag{4}
\]
初始常量和 free 零谱常量均在此消去。可从 \(F_\sigma\) 的强 \(L^2\) 有界性与积分分部证明 (4)，无须 ordinary Fourier of \(F_\sigma\)。\(R_\sigma\in L^2(dv;L^2)\)。对几乎处处 \(x\)，其 Fourier 为 \(L^1(d\xi)\)：用 Cauchy–Schwarz 和 \(\int(1+\omega^2)^{-1}d\xi<\infty\)。故 \(R\) 的连续 Fourier 代表及以下 pointwise split 合法。

**几乎处处所有 \(t\) 的代表。** 不能只凭强 \(L^2\) 连续性认同 pointwise 全时间轨道。令正算子 \(M=qI+\mathcal J\)，则 \(|L^k\eta|\le M^k|\eta|\)。对每个整数 \(T>0\) 和 \(\eta\in L^1\)，
\[
\sum_{k\ge0}\frac{T^k}{k!}M^k|\eta|\in L^1,\qquad
\left\|\sum_{k\ge0}\frac{T^k}{k!}M^k|\eta|\right\|_1
\le e^{2QT}\|\eta\|_1.
\]
所以 \(H_t\eta=\sum_{k\ge0}(-t)^kL^k\eta/k!\) 在 \(0\le t\le T\) 几乎处处绝对一致收敛，给出连续代表；相应导数级数也由 \(M e^{TM}|\eta|\in L^1\) 控制。对本文有限个输入和整数 \(T\) 取共同零集即可。另一方面，\(\widehat R=\widehat G/(1+i\omega)\) 属于 \(L^1(d\xi;L^2_x)\)，其反演为强 \(L^2\) 连续函数，并且对几乎处处 \(x\) 为连续标量函数。它与原 \(R\) 在每个有理 logtime 的强 \(L^2\) 代表相同；删去可数零集，再由两侧 pointwise 连续性，得到所有 \(v\)（即所有 \(t>0\)）相同。以下 maximal、域外 sign、cap 与 split 使用这些共同代表；不依赖随 \(t\) 变化的零集。

## 3. fixed Abel 给出的低频费用

Abel Gamma 平均对 logtime 导数同样成立。令
\(\Phi_N(r)=N^Nr^{N-1}e^{-Nr}/\Gamma(N)\)，
\(\mathcal G_N(v)=A_{Ne^{-v}}^N\nu_{\rm b}\)，则
\[
\partial_v\mathcal G_N(v)=\int_0^\infty\Phi_N(r)G_{\nu_{\rm b}}(v+\log r)dr,
\quad
M_N(\xi)=N^{-i\omega}\Gamma(N+i\omega)/\Gamma(N).
\tag{5}
\]
这是 \(L^2\) 导数恒等式；有限链 free 零模不贡献导数。在 \(\Omega^c\)，killed 延拓及其导数为零，所以已证完整 Abel 能量控制 (5) 的外部导数：
\(\int_{\Omega^c}\int|\partial_v\mathcal G_N|^2\le276\sqrt N\kappa W_{\rm b}\)。

无限积给 \(|M_N|^2\ge\exp[-2\omega^2/N]\)。对 \(B\ge1\) 可取
\(N=\lceil256B^2\rceil\)，用 \(\pi<4\) 得 \(|M_N|^2\ge e^{-1/2}>1/2\) 于 \(|\xi|\le B\)，\(\sqrt N<17B\)。于是
\[
E_{\nu_{\rm b}}(B):=\int_{\Omega^c}\int_{|\xi|\le B}|\widehat G_{\nu_{\rm b}}|^2
\le9384B\,\kappa W_{\rm b}.
\]
由 (2)、\(\|\mu_{\rm b}\|_2^2\le\kappa W_{\rm b}\)，以及 \(G_\sigma=G_{\nu_{\rm b}}-G_{\mu_{\rm b}}\)，
\[
E_\sigma(B):=\int_{\Omega^c}\int_{|\xi|\le B}|\widehat G_\sigma|^2
\le18769B\,\kappa W_{\rm b},\qquad B\ge1.
\tag{6}
\]
这些保守常数来自解析界，不由新数值拟合。

定义
\[
T_\Lambda(x)=\sup_v\left|\int_{|\xi|\le\Lambda}
e^{2\pi i\xi v}\frac{\widehat G_{\sigma,x}(\xi)}{1+i\omega}\,d\xi\right|,
\quad\Lambda>1.
\tag{7}
\]
分 \(|\xi|\le1\) 与 \(2^{j-1}<|\xi|\le2^j\)。第一块的 Cauchy–Schwarz 界至多 \(2E_\sigma(1)\)；每个其它块用
\(\int_{\rm block}(1+\omega^2)^{-1}d\xi\le2^{-j}\) 与 (6)，平方范数至多 \(18769\kappa W_{\rm b}\)。约 \(\log(e\Lambda)\) 块的 Minkowski 求和得到
\[
\int_{\Omega^c}T_\Lambda^2
\le C\log^2(e\Lambda)\,\kappa W_{\rm b}.
\tag{8}
\]
这一步不依赖 heat 空间核，不依赖 \(L,K\) 交换，也不假设 \(F_\sigma\) 可积。

## 4. 真的可回用 maximal 合同：尚缺的一条高尾

令
\[
Q_\Lambda(x)=\int_{|\xi|>\Lambda}
\frac{|\widehat G_{\sigma,x}(\xi)|}{\sqrt{1+\omega^2}}\,d\xi.
\tag{9}
\]
它对几乎处处 \(x\) 有限，且 Fourier inversion 给
\(\sup_t|R_t\sigma|\le T_\Lambda+Q_\Lambda\)。但此处已知的 \(L^2\) finiteness 并不支付其空间 \(L^1\)。

由 (3)，在 \(\Omega^c\)，
\[
H_t\nu=H_t\nu_{\rm g}+H_t\mu_{\rm b}+\widetilde H_t\sigma+R_t\sigma
\le\kappa+R_t\sigma.
\tag{10}
\]
这里 \(H_t\nu_{\rm g}+H_t\mu_{\rm b}=H_t\mu_{\rm tot}\le\kappa\)，使用前稿完整目标的同一个 cap；原较宽的 \(2\kappa\) 界也正确。下面仍保留原保守 cap 选择和费用常数，无须重算 guard。
若对同一 \((\nu,\kappa,u,\Omega)\) 证明
\[
\boxed{\ \int_{\Omega^c}Q_\Lambda(x)dx\le C_{\rm hi}W_{\rm b}\ },
\tag{11}
\]
其中 \(C_{\rm hi}\) 一致，取 \(\kappa=\alpha/[4\log(e\Lambda)]\)，则域内 \(|\Omega|\)、低频 Chebyshev (8) 及高频 Markov (11) 给
\[
\alpha|\{\sup_{t>0}H_t\nu>\alpha\}|
\le C(1+C_{\rm hi})\log(e\Lambda)\,W.
\tag{12}
\]
这里来源只拆一次；(12) 是**条件结论**，不把 (11) 当作已证。若 \(\Lambda\) 为 dimension 的固定多项式且 \(C_{\rm hi}=O(1)\)，可得 fixed-family 的 polylog 弱型。它仍不自动给原时变 cube kernel 的正比较、共同尺度 source payment 或实际门回代。

## 5. compound Poisson 的精确高频表示

现在限于常率 \(L=q_0(I-P)\)、\(q_0>0\)，\(P\) 为正保守 \(L^1,L^2\) 收缩且可逆。固定轴向 \(P=n^{-1}\sum_iP_i\) 是特例。

谱变量给 \(\widehat G_\sigma(\xi)=-\Gamma(1-i\omega)L^{i\omega}\sigma\)；
\(L^{i\omega}\) 在零谱定义为零，而 \(\sigma=Lu\) 本已与零谱正交。以 \((I-P)^{1+i\omega}\) 的绝对收敛 binomial 展开（实部为 1），在 \(\Omega^c\) 去掉 \(P^0u=u=0\)，得到
\[
\boxed{\ \widehat R_{\sigma,x}(\xi)
=-i\omega\,q_0^{1+i\omega}
\sum_{k\ge1}\frac{\Gamma(k-1-i\omega)}{k!}(P^ku)(x)\ },\qquad\omega\ne0.
\tag{13}
\]
这是导数/分布 Mellin 的精确式，不是未经许可的普通 \(F_\sigma\) 积分。固定 \(\omega\) 下级数于 \(L^1,L^2\) 绝对收敛，因为系数为 \(O_\omega(k^{-2})\)。\(q_0^{i\omega}\) 只有相位；不能因此消除后面的 \(q_0\|u\|_1\) 费用。

### 5.1 逐 \(k\) 高频三角积分确实发散

对整数 \(k\ge3\)，Euler 无限积给
\[
\frac{|\Gamma(k-1-i\omega)|}{\Gamma(k-1)}
=\prod_{m\ge0}\left(1+\frac{\omega^2}{(k-1+m)^2}\right)^{-1/2}.
\]
利用 \(\log(1+x)\le x\) 和
\(\sum_{m\ge0}(k-1+m)^{-2}\le1/(k-2)\)，在
\(\sqrt k\le\omega\le2\sqrt k\) 上得下界
\(\exp[-2k/(k-2)]\ge e^{-6}\)。因此只取正频率带，也有
\[
\int_{\sqrt k}^{2\sqrt k}
\frac{\omega|\Gamma(k-1-i\omega)|}{k!}\frac{d\omega}{2\pi}
\ge\frac{3e^{-6}}{4\pi(k-1)}
>\frac1{3888(k-1)}.
\tag{14}
\]
最后用 \(e<3,\pi<4\)。当 \(k>(2\pi\Lambda)^2\)，该带完全落在高频。

冻结轴向 \(P\) 的 Fourier multiplier
\(\widehat P(\zeta)=n^{-1}\sum_i\widehat G_{c,\ell}(\zeta_i)\) 在非零 \(\zeta\) 严格小于 1，且模不超过 1。故 \(P^ku\to0\) 于 \(L^2\)。\(\Omega\) 有限体积，于是
\[
\int_\Omega P^ku\le|\Omega|^{1/2}\|P^ku\|_2\to0,\qquad
\int_{\Omega^c}P^ku\to\|u\|_1.
\tag{15}
\]
若 \(u\ne0\)，将 (14) 乘上 (15) 的域外质量，再对 \(k\) 求和，得到**逐层取绝对值的 proposed upper** 为 \(+\infty\)。不能用 Tonelli 声称高频收费有限。该证明没有声称真实复数和 (13) 的高尾发散；跨 \(k\) 抵消恰是必须保留的内容。即使不用 triangle，\(q_0\|u\|_1\) 也不能免费替成 \(W_{\rm b}\)。

## 6. 哪些参数可免费消去，哪个空间接口未付

物理尺度 \(\ell>0\) 是空间 dilation：同时 rescale 输入、cap 和障碍，质量 \(W\) 保持，(11) 的空间积分随密度/体积精确抵消。因此若 scale 1 的相应 family 定理已证，\(\ell\) 不新增费用。乘生成元常数 \(q_0\) 是 logtime translation，亦不改变真实高尾幅度。

但改变 \(c=1-s\) 不仅改变 \(q_0=n/c\)，还改变一维 jump shape \(G_c\)。不能以 (13) 的 \(q_0^{i\omega}\) 相位证明对 \(c\) 一致。若限制 \(c\ge c_0>0\)，仍需证明该整个 shape-family 的 (11)；任意 \(c\downarrow0\) 则还需记录 shape 损失。

一个足够但未证的空间桥是：在真实 obstacle 的域外，
\[
\|\,\mathbf1_{\Omega^c}L_{\rm fr}^{i\omega}\sigma\,\|_1
\le C_0^{\,n}(1+|\omega|)^{C_1n}e^{\theta|\omega|}W_{\rm b},
\quad \theta<\pi/2,
\tag{16}
\]
常数对所需 \(c\)-范围一致。Gamma 因子除以 \(\sqrt{1+\omega^2}\) 有
\(O(|\omega|^{-1/2}e^{-\pi|\omega|/2})\) 衰减；若 (16) 成立，
\(\Lambda=A n\log(n+2)\) 的充分大固定 \(A\) 可支付 (11)，仍使 \(\log\Lambda=O(\log(n+2))\)。
但 (16) 不是一般 \(L^1\) imaginary-power theorem：必须用 capped obstacle 的域外结构和完整原来源。对所有 \(L^1\) 输入的 unrestricted imaginary powers 本就不提供这里的 endpoint；也不能把 (16) 写成仅剩 scalar Gamma 的已证界。

热论文通过其特定空间 potential 收费完成对应桥；jump 中单跳 \(J\) 的外部通量虽已付，多跳 fractional/imaginary kernel 的长尾不能直接由它支配。本轮没有证明 (16) 或其它替代 (11)。

## 7. 三轮新的有限/半解析方法守卫

先写 frozen_jump_mellin_registration_20261007.json，再运行 guard。三轮是 \(m=4,16,64\) 个状态的 uniform-reset 链，\(q_0=1,3/2,7/3\)，\(\kappa=1\)、\(\nu=(m/2)\delta_0\)。明确是有限原密度与 counting measure，**不是空间 dimension 或 actual cube**。

它有显式 \(u\) 仅在 state 0 正，域外吸收密度
\((m-2)/[2(m-1)]\)，内帽为 1，完整质量 \(m/2\)。\(\sigma\) 只在 eigenvalue \(q_0\) 上，故真实
\[
\widehat R_{\sigma,x}(\xi)
=-\sigma(x)\,q_0^{i\omega}\frac{\Gamma(1-i\omega)}{1+i\omega}.
\tag{17}
\]
真实尾指数小；用 \(|\Gamma(1-i\omega)|^2=\pi|\omega|/\sinh(\pi|\omega|)\)，对 cutoff \(\Lambda=m\) 可保守界
\(\int_{\Omega^c}Q_\Lambda<4\cdot2^{-\Lambda}W_{\rm b}\)。
与此同时 \(P^ku\) 在域外对每个 \(k\ge1\) 都是同一个正均值，(14) 的 triangle 为 harmonic 发散。这个例子严格区分“triangle 失败”和“真实高尾失败”。

新 guard 检查：

- 精确 Fraction 的内帽、完整外部吸收、源质量、体积，以及真实尾的有理上界。
- \(\omega=1/2,1,2,4\) 时直接 logtime 积分 \(R\)，与独立 integrand 的 \(\Gamma(1-i\omega)\) 积分交叉核 (17) 及 \(q_0\) 相位。
- \(k=(8\Lambda)^2\) 的 1/2/4 倍、\(\omega/\sqrt k=1,3/2,2\) 的 Gamma 比值 strip 下界；浮点有限乘积仅作诊断，(14) 由解析无限积证明。
- 1/4/16 个 doubling harmonic blocks 的精确费用累加；无穷发散来自解析论证，没有用三点拟合。

结果 PASS_EXACT_AND_NUMERIC，135 个谓词。32 个 panel、logtime \([-40,40]\)，三轮 GL16/32、32/64、64/128；最大 identity error \(8.44\times10^{-15}\)，最大加密差 \(7.65\times10^{-14}\)。这是 floating quadrature，**不是 interval enclosure**。解析 logtime 两尾另界：normalized \(R\) 的遗漏 \(\le1.5e^{-40}\)，Gamma 的遗漏 \(\le e^{-40}+e^{-e^{40}}\)。

首次默认 Python 缺 numpy，在 import 阶段停止，未产生数学结果；随后使用现有 bundled Python 运行同一脚本，未安装包、未改参数/容差，三轮一次通过。保存 script/registration SHA；未声称另一算法独立重构所有结果。

**没有新的 cube 费用。已付/可复用的是固定算子的导数 Mellin 低频接口与精确条件 maximal 合同；真正未付的是 (11) 的 obstacle-restricted 高频空间 \(L^1\) bridge，以及此后原时变 cube 的生成元/参数/实际门回代。**
