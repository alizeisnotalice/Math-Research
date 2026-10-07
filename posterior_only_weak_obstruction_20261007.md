# Posterior-only 弱停止引理的精确否定证书

2026-10-07。本文件核验 bernoulli_endpoint_stopping_attempt_20261007.md §2 的新弱策略障碍。**这不是原 G_c 的反例，也不是 reset 模型。** 这里仅反驳“正 Bernstein、层质量、驻点后验和输入高度帽足够给 polylog 弱收费”的放松引理。真实 mask 的共同卷积来源相容性仍必须保留。

## 1. 两个显式有限测度构造

固定 n=8,32,128，W=1，lambda=3/2。令
\[
 p_k=k/n,\quad b_k=p_k^k(1-p_k)^{n-k},\quad
 v_k={n\choose k}b_k/\lambda,\quad A_k=\lambda/b_k.
\]
采用0^0=1。取 E_k=[2k,2k+v_k]×[0,1]^(n−1)，各盒有正间隙、有限体积。令 a_k=A_k1_(E_k)，F_r=sum a_k r^k(1−r)^(n−k)。

每层 ∫a_k=binom(n,k)，因此 **对所有 r**，不是仅测试节点，
\[
 \int F_r=\sum_k{n\choose k}r^k(1-r)^{n-k}=1.
\]
每个 E_k 的连续 r 最大值恰为 lambda。内部 basis 的导数为
 r^(k−1)(1−r)^(n−k−1)(k−nr)，故最大点为 p_k；端点 k=0,n 原样保留。峰值后的代数 posterior 是单点 k，均值 np_k、方差0。

设 D=sum binom(n,k)b_k。严格弱水平集在 t<lambda 时包含所有 E_k，在 t≥lambda 时为空，故
\[
 \|F_*\|_{1,\infty}=D.
\]
这是 supremum，不是在严格 threshold lambda 处取得。

Height-cap 加强版固定
\[
 H=\lambda\max\{1,\max_{k\ge1}[{n\choose k}b_k]^{-1}\}.
\]
只把 E_0 的体积改为1/H、a_0的高度改为H，其余盒和系数不动。此时 f=a_0、∫f=1，全部 a_k≤binom(n,k)H；F_0=f、F_r≤H、每 r 质量仍1。该加强版 weak norm **精确等于**
\[
 D_{\rm cap}=D-1+\lambda/H>1.                  \tag{1}
\]
证明分两段：t<lambda 的 ratio 为 t[(D−1)/lambda+1/H]，趋(1)；lambda≤t<H 仅E_0，ratio=t/H的 supremum 为1；t≥H为空。两种严格端点都已核。

对 m=1,2,4,8,16,32,64，预设 t_m=lambda(1−2^(−m))。两种 ratio 分别恰为 (1−2^(−m))D、(1−2^(−m))D_cap；gap 精确为各自 supremum 的2^(−m)倍。

## 2. 注册与三轮 Fraction 收据

先保存 posterior_only_weak_obstruction_registration_20261007.json，再运行新脚本一次。所有决定性值为 Fraction/整数；显示小数没有进入判断。核每层质量、盒间隙、连续 peak 的导数身份/左右符号、posterior normalization/均值/方差、多个 r 的总质量、逐层 height cap，以及严格阈值逼近和 lambda/H 端点。

|n|通过谓词数|baseline weak supremum 显示|height-cap supremum 显示|
|---|---:|---:|---:|
|8|213|4.2450180054|3.5184555054|
|32|645|7.7740454767|6.9139954108|
|128|2373|14.8553006829|13.9256867751|

终态 PASS_EXACT，共3231个已执行 Boolean 谓词。结果保存完整有理 layer rows、各注册 strict thresholds、test names、分子/分母 binary SHA256，以及注册和脚本 SHA256。有限 r 检查是实现守卫，任意 r 保质量由上述 binomial identity 证明；连续 peak 由导数符号证明。没有拟合维数阶数。

## 3. 适用边界

原 family 的真实系数必须满足 nu_A=G_Af、a_k=sum_(|A|=k)nu_A 和 nu_(A∪i)=G_i nu_A。当前构造不满足它们；特别非零 f 的原 G_[n]f 在所有 receiver 严格正，而这里 a_n 仅在一个有限盒上非零。因此本收据不能用作原 A_ord 的下界，不能称 actual FIRST/CPGP/hardwinner/history 样本。

该证书排除的是新 posterior-only **弱**停止策略，未重跑旧 strong D_n 或原核实验。主账不变；原完整 Bernoulli polylog 弱端点仍未证明或否定。
