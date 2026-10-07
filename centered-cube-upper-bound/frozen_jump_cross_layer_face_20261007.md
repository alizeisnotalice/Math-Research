# 高频跨层接口：proper-face 源预算与未付的 all-visited 核

2026-10-07。仅新增本稿及同前缀 registration/guard/results；不修改主稿、总账或旧数据。本轮的新结论是：冻结轴向过程的 holding 原子及 proper-face 部分在 cutoff \(\Lambda=16n\) 的高 Mellin 尾可用 \(W_{\rm b}/8\) 支付，**对全部 \(c\in(0,1]\)、\(\ell>0\) 一致**。all-visited 完整密度部分仍未支付。本稿不把部分工具费计入实际 cube 主账。

## 1. 查重与 Fourier 对象

已读 [F02 SKILL](/Users/zhengzhihao/.codex/skills/math-f02-tensor-fourier-positivity-soft-square-supremum/SKILL.md)、provenance、method、cube-interface。技能只提供张量/正性审查，未提供 imaginary-power 或 endpoint 定理。总 tex 的空间 Fourier 估计、已有 Mellin 低频接口、full \(L^1\) imaginary-power 限制及逐 \(k\) 方法失败均不重做；本稿不再测试有限 reset 链。

冻结实际正算子及其乘子为
\[
L_{\rm fr}=c^{-1}\sum_{i=1}^n(I-P_i)=q_0(I-P),\quad
q_0=n/c,\quad P=n^{-1}\sum_iP_i,\quad
\widehat P_i(\zeta)=\widehat G_c(\ell\zeta_i)
=\frac1{1+cB(\ell^2\zeta_i^2)},
\]
\[
B(v)=v/\log(1+v)-1,\qquad
\widehat L_{\rm fr}=\sum_i\frac{B(\ell^2\zeta_i^2)}{1+cB(\ell^2\zeta_i^2)}.
\tag{1}
\]
它是求和生成元的乘积 semigroup，不是把一个共同 Brownian clock 免费用于全部坐标。这里只使用已核 \(G_c\) 正概率密度的表示。

保留前稿 obstacle：
\(\sigma=\nu_{\rm b}-\mu_{\rm b}=Lu\)，
\(\nu_{\rm b},\mu_{\rm b}\ge0\)，
\(\int\nu_{\rm b}=\int\mu_{\rm b}=W_{\rm b}\)，
\(0\le\mu_{\rm b}\le\kappa\)，
\(\Omega=\{u>0\}\)，\(u\in L^1\cap L^2\)。
特别
\[
\|\sigma\|_1\le2W_{\rm b}.
\tag{2}
\]
不把 \(q_0\|u\|_1\) 换成来源质量。原域外 ergodic sign 和 a.e 全时间代表沿用 frozen_jump_mellin_interface_20261007.md；不重新 FIRST 或改变来源。

## 2. 明确分离全部轴面，而非丢掉它们

条件于 \(k\) 次轴向跳跃，标签为独立均匀坐标。把正概率卷积 measure \(P^k\) 分成
\[
P^k=M_k+C_k,\qquad k\ge1.
\tag{3}
\]
\(M_k\) 收集至少一个坐标从未被选的全部标签序列，\(C_k\) 收集每个坐标至少被选一次的序列。令 \(M_0=I,C_0=0\)。有
\[
m_k:=\|M_k\|_{\rm measure}
=\mathbb P(\text{some coordinate unvisited at }k)
\le\min\{1,n(1-1/n)^k\}.
\tag{4}
\]
各 \(M_k,C_k\) 都是原过程的未归一化正 measure，不作条件归一化。\(M_k\) 包括所有 proper coordinate faces；\(k=0\) 是 holding 原子。\(C_k\) 在所有 \(G_c\) 为密度时是全维密度。无序列被删除或当作零。

这里的 coupon 标签只来自冻结生成元 \(q_0(I-P)\) 的 Poisson 表示，是内部跳跃所选的坐标。它们不是输入源标签、下界模型的标签数或辅助 killed 算子 \(K\) 的标签；没有假设源测度的坐标独立。任意非产品、XOR、Cantor 或其他联合来源都可以作用这个线性算子。all-visited 只描述该跳跃历史，不能解释为来源复杂度受限。

Fourier 上此分解是逐标签乘积：
\(\prod_i\widehat G_c(\ell\zeta_i)^{N_i}\)，未访问坐标的因子为 1。总 jump kernel 的正空间性与 \(\widehat G_c\ge0\) 是不同事实；这里的收费靠正 measure 的总质量，不靠错误的 cube sinc 正性。

## 3. derivative Mellin 的准确 operator 分解

沿用 \(\omega=2\pi\xi\) 及
\[
\widehat R_\sigma(\xi)
=-\frac{\Gamma(1-i\omega)}{1+i\omega}L^{i\omega}\sigma.
\tag{5}
\]
这是 logtime 导数 Fourier 定义，原 \(H_{e^v}\sigma\) 仍可有非零初始常量。为了严谨处理实部为零的 binomial 级数，先以
\((I-rP)^{i\omega}\) 正则化，\(0<r<1\)，最后 \(r\uparrow1\) 作强 \(L^2\) 极限；\(\sigma=Lu\) 不含 free 零谱。不能假设整个未经正则化的 \(k^{-1}\) 级数在 \(L^1\) 绝对收敛。

由
\((-1)^k\binom{i\omega}{k}=\Gamma(k-i\omega)/[\Gamma(-i\omega)k!]\)，定义 holding/proper-face 部分
\[
\widehat R_{\rm miss}(\xi)
=q_0^{i\omega}\left[
-\frac{\Gamma(1-i\omega)}{1+i\omega}\sigma
+\frac{i\omega}{1+i\omega}
\sum_{k\ge1}\frac{\Gamma(k-i\omega)}{k!}M_k\sigma
\right].
\tag{6}
\]
由 (4)，该级数在每个固定 \(\omega\) 于 \(L^1,L^2\) 绝对收敛，可以解除 \(r\)。令
\(\widehat R_{\rm full}=\widehat R_\sigma-\widehat R_{\rm miss}\)；
剩余项是原 \(C_k\) 的 Abel-limit 跨层和，没有将它标成已付或可逐层绝对收费。

## 4. 新的一致高尾源预算

对所有整数 \(k\ge1\)、\(\omega>0\)，Euler 积给
\[
r_k(\omega):=\frac{|\Gamma(k-i\omega)|}{\Gamma(k)}
=\prod_{j\ge0}\left(1+\frac{\omega^2}{(k+j)^2}\right)^{-1/2}
\le\exp\left[-\frac1{10}\min\{\omega^2/k,\omega\}\right].
\tag{7}
\]
证明可直接核：当 \(\omega\le k\)，每个比例平方 \(\le1\)，\(\log(1+x)\ge x/2\) 和 \(\sum_{j\ge0}(k+j)^{-2}\ge1/k\) 给指数 \(\omega^2/(4k)\)。当 \(\omega>k\)，把递减 log 的和下界为积分，取 \(x\in[\omega,2\omega]\)，使用 \(\log(5/4)\ge1/5\)，给指数 \(\omega/10\)。故 (7) 不需要 Stirling 渐近或数值拟合。

令 \(\rho=1-1/n\)。以 \(k\le\omega\sqrt n\) 和其余分开：
\[
\sum_{k\ge1}\frac{m_k}{k}r_k(\omega)
\le n\log n\,e^{-\omega/(10\sqrt n)}
+n^2e^{-\omega/\sqrt n}
\le2n^2e^{-\omega/(10\sqrt n)}.
\tag{8}
\]
第一段用 \(\sum\rho^k/k=\log n\)，第二段用几何尾及 \(1/k\le1\)。\(n=1\) 时 \(m_k=0\) 对所有 \(k\ge1\)，直接成立。

因为 \(\|M_k\sigma\|_1\le m_k\|\sigma\|_1\)、\(|i\omega/(1+i\omega)|\le1\)，而 \(q_0^{i\omega}\) 只是相位，(8) 在这里允许 Tonelli：
\[
\int_{|\xi|>\Lambda}
\|\text{(6) 的 }k\ge1\text{ 部分}\|_1d\xi
\le\frac{40}{\pi}n^{5/2}
e^{-\pi\Lambda/(5\sqrt n)}W_{\rm b}.
\tag{9}
\]
holding 项另用
\(|\Gamma(1-i\omega)|^2=\pi\omega/\sinh(\pi\omega)\)：
对 \(\Lambda\ge1\)，其高尾 \(L^1\) 至多 \(8e^{-\Lambda}W_{\rm b}\)。

所以完整 proper-face 费用为
\[
\boxed{
\int_X\int_{|\xi|>\Lambda}|\widehat R_{\rm miss}(\xi,x)|d\xi\,dm(x)
\le
\left[\frac{40}{\pi}n^{5/2}e^{-\pi\Lambda/(5\sqrt n)}
+8e^{-\Lambda}\right]W_{\rm b}.}
\tag{10}
\]
它在**全空间**成立，因而也可限制于 \(\Omega^c\)；对全部 \(c,\ell\) 一致。没有使用初始源的特殊形状、层间强 \(L^1\) imaginary powers 或目标空间样本。

### 4.1 预设 cutoff 与常数

固定 \(\Lambda=16n\)。用 \(\pi>3\)、\(40/\pi<e^3\)、\(\log n\le\sqrt n\)，(10) 第一项严格小于
\[
\exp[3+(5/2)\sqrt n-(48/5)\sqrt n]
\le e^{-41/10}<1/16.
\]
holding 项 \(<8\cdot2^{-16n}\le1/8192\)。因此
\[
\boxed{\ \int_X Q^{\rm miss}_{16n}\,dm<W_{\rm b}/8.\ }
\tag{11}
\]
选择与输入/数据无关，未事后挑选有利 cutoff。任意更大 cutoff 也保留此界。

这是已证的 generator-family 费用。其来源是 (2) 的一份 bad 质量；实际 cube 的时变过程并不因此获得共同参数分解或已付新交通。

## 5. all-visited 的精确剩余与弱 endpoint 路线

当前需要控制的是
\[
Q^{\rm full}_\Lambda(x)=
\int_{|\xi|>\Lambda}|\widehat R_{\rm full}(\xi,x)|d\xi
\tag{12}
\]
在原 \(\Omega^c\) 的源预算。它保持原 \(\sigma=Lu\)、完整障碍及全部 \(C_k\)；不能因为每个 \(C_k\) 有密度就免费得到跨 \(k\) 的 \(L^1\) 界。

强界 \(\int_{\Omega^c}Q^{\rm full}_\Lambda\le CW_{\rm b}\) 是充分条件而非必要。也可直接证明
\[
\tau|\{x\in\Omega^c:Q^{\rm full}_\Lambda(x)>\tau\}|
\le CW_{\rm b},
\tag{13}
\]
或证明 full 高频 inverse 的时间最大值具有同样弱界；与旧低频 Chebyshev 和新 (11) 即可结合。本文未证明 (13)，没有对连续参数直接套 \(L^{1,\infty}\) Minkowski。

有一个具体可回代的 CZ 合同，说明尚缺哪一条空间估计。设
\(\mathcal T_\Lambda f=(\widehat R_{\rm full}(\xi,\cdot))_{|\xi|>\Lambda}\)，目标 Banach 空间为 \(L^1(d\xi)\)。由 fixed spectral \(L^2\) 与 (8)，Minkowski 在**强 \(L^2\)** 中给
\[
\|\mathcal T_\Lambda f\|_{L^2_x(L^1_\xi)}
\le B_\Lambda\|f\|_2,\quad
B_\Lambda\le8e^{-\Lambda}
+\frac{20}{\pi}n^{5/2}e^{-\pi\Lambda/(5\sqrt n)}.
\tag{14}
\]
此处 whole \(L^{i\omega}\) 在 \(L^2\) 为 contraction；full 是 whole 减 miss。这个 bound 对 \(c,\ell\) 一致，但 \(\|f\|_2\) 不是源 \(L^1\) 费用。

还缺 full 核的 Banach-valued Hörmander 条件及其参数一致常数：
\[
\sup_{h\ne0}\int_{\|x\|_\infty>2\|h\|_\infty}
\|\mathcal K_\Lambda^{\rm full}(x-h)-\mathcal K_\Lambda^{\rm full}(x)\|_{L^1_\xi}dx
\le H_\Lambda,
\tag{15}
\]
并需合法的离对角 kernel representation。它必须在原 \(\widehat G_c\) family 上证明；低频 multiplier 的正性、有界率或每层 smooth density 都不替代 (15)。

若 (15) 存在，标准 dyadic 分解可以直接用显式阈值处理连续频率，避免弱 Minkowski：对输入 \(f\)，在 height
\(h_0=\tau(3/2)^{n/2}/(2B_\Lambda)\) 分解。三倍 bad cubes 的测度至多 \(3^n\|f\|_1/h_0\)；good 项的 Chebyshev 费用至多 \(4B_\Lambda^22^nh_0\|f\|_1/\tau^2\)；bad 项在外部的 Banach norm 积分至多 \(2H_\Lambda\|f\|_1\)。因此条件性地
\[
\tau|\{\|\mathcal T_\Lambda f\|_{L^1_\xi}>\tau\}|
\le4(6^{n/2}B_\Lambda+H_\Lambda)\|f\|_1.
\tag{16}
\]
对原 \(\sigma\) 可用 (2)。真正需要的是把右侧常数在一个预设多项式 cutoff 下压成 \(O(1)\)；例如 (14) 中的 \(6^{n/2}\) 不能在 \(\Lambda=16n\) 免费略去。允许更大的 poly(n) cutoff，但必须同时核 (15) 的 \(c\)-一致性。式 (16) 是**明确条件接口**，不是新已付弱界。

根端提供的 BF subordination 时间解析 \(4n\) bound 可以辅助研究 (15)，但时间导数 \(L^1\) 控制自身没有空间 Hörmander 差。本文未用该外部工具完成或数值认证 (15)。原先普适强 \(L^1\) 界的限制不在此重复证明。

## 6. 三轮新注册压力：8/32/128

先注册 frozen_jump_cross_layer_face_registration_20261007.json，再用已有 bundled Python 运行同前缀 guard。cutoff 固定 \(16n\)。所有 coupon 质量用正的整数递推
\[
d_{k+1,j}=j\,d_{k,j}+(n-j+1)d_{k,j-1},\quad d_{0,0}=1,
\quad m_k=1-d_{k,n}/n^k.
\]
三轮 \(n=8,32,128\)，逐步到 \(8n\)，检查总字符串数 \(n^k\)、非负性、missing 单调、精确 union bound，保存七个预设步的 mass/hash。

Gamma diagnostics 用精确整数参数 modulus 的有限乘积公式计算 log：
\[
r_k(\omega)^2
=\frac{\pi\omega}{\sinh(\pi\omega)}
\prod_{j=1}^{k-1}(1+\omega^2/j^2).
\]
它们检查 (7) 的上界及四个预设频率的 finite partial sums 对 (8) 的方向。无穷 \(k\) 尾另保留固定频率几何上界；没有把有限 partial sum 当全尾。cutoff 常数的有理余量精确核查，(7)–(11) 的无穷积分由解析证明负责。

结果 PASS_EXACT_AND_NUMERIC，6828 个谓词，三轮一次通过。整数/Fraction 部分精确；log-Gamma 数值为浮点诊断，**非区间认证**，没有拟合 dimension 阶数、核 c 的有限扫选或声称 actual FIRST 样本。对 \(c,\ell\) 的一致性来自标签总质量与相位证明，不来自采样。root 另用 inclusion-exclusion 只读重构全部21个保存 coupon mass/hash，共66检查通过；未重跑本脚本、未独立复核浮点 Gamma。

**终态：proper-face/holding 高频费用 (11) 已证且参数一致；all-visited 高频空间费用仍未付，最小的新空间缺口为 (12)/(13) 或可用的 (15)。不增加 actual cube 主账费用。**
