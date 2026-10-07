# 原 Bernoulli 弱端点：后验驻点停止策略的相容性检验

2026-10-07。本轮固定 c∈[1/2,1]、物理尺度 1 和一份任意非负 L1 密度 f，W=||f||1。研究完整 T_r=[(1−r)I+rG_c]^(tensor n)，r∈[0,1]，不重做旧强包络、L2 或投影压缩，也未新增温和 shell 实验。本稿没有证明 polylog 弱端点，没有向 cube 主账收费。

## 1. 原核系数必须满足什么

写 G_A=prod_(i∈A)G_c^(i)，nu_A=G_A f。则
\[
 T_rf(x)=\sum_{A\subseteq[n]}r^{|A|}(1-r)^{n-|A|}\nu_A(x)
        =\sum_{k=0}^n r^k(1-r)^{n-k}a_k(x),\qquad
 a_k=\sum_{|A|=k}\nu_A.                         \tag{1}
\]
真实相容条件不止 a_k≥0 和 ||a_k||1=binom(n,k)W：
\[
 \nu_\varnothing=f,\quad
 \nu_{A\cup\{i\}}=G_c^{(i)}\nu_A\quad(i\notin A),\quad
 \widehat\nu_A=\widehat f\prod_{i\in A}\frac1{1+cB(\xi_i^2)}. \tag{2}
\]
等价地，(I+cB_i)nu_(A∪i)=nu_A，按分布/测试函数解释，不假定任意 f 在 B_i 的强算子域内。每条边是同一源的正卷积；不同 i 的边方形 commute。若 0≤f≤H，还必须 0≤nu_A≤H 和 a_k≤binom(n,k)H。

谱上 m_r 的完全单调不使 T_rf(x) 逐点单调；导数含 G_i−I 的 signed 项。因此一个对 r 的普通 Markov stopping/optional-sampling 公式不能仅凭谱单调写出。

## 2. 后验均值及每参数保质量仍不足以支付弱水平集

以下是对一个候选证明步骤的严格否定，不是原 G_c 的反例。放松仅保留 a_k≥0、∫a_k=binom(n,k)W，以及内部最大点的代数后验 E k=np。归一 W=1，给定 lambda>0，令
\[
 p_k=k/n,\quad b_k=p_k^k(1-p_k)^{n-k},\quad
 v_k={n\choose k}b_k/\lambda.
\]
端点用 0^0=1。在 Rn 取互不相交、有限体积盒 E_k，|E_k|=v_k，置 a_k=(lambda/b_k)1_(E_k)。于是每个 r 的响应 F_r 总质量恰一。在 E_k 上 F_r 是单个 Bernstein basis，其最大高度恰 lambda，内部最大 p_k 的后验是单点 k，满足 E k=np_k，甚至 variance=0。严格弱水平集取 threshold t↑lambda：
\[
 \|\sup_rF_r\|_{1,\infty}
 =\lambda\sum_k|E_k|
 =\sum_k{n\choose k}b_k=D_n.                  \tag{3}
\]
所有 E_k 测度有限；严格 > 的端点没有被忽略。旧 D_n 的 sqrt(n) 量级意味着上述信息不足以给 polylog。这里检验的是 **weak** 而不是再次讨论强范数。

这个步骤也不能靠补一个 input-height cap 救回。将 E_0 改成体积 1/H、a_0=H1_(E_0)，其余 k≥1 保持原构造，并选
\[
 H=\lambda\max\left\{1,\max_{1\le k\le n}
                [\,{n\choose k}b_k\,]^{-1}\right\}.
\]
全部 a_k≤binom(n,k)H、a_0=f 及一份初始质量均成立；在 threshold 趋 lambda 下仍有 weak ratio 至少 D_n−1+lambda/H。

**为什么它不是原核输入？** (2) 不成立。特别 f=a_0 非零而原 G_[n]f 在每个 receiver 都严格正（原 G_c 是严格正的从属热密度），上述 a_n 却只支撑 E_n。各 a_k 不能随意独立搬到不交盒；就算每行驻点及每层质量都成立，也不能声称它们来自同一 f。此构造不使用 reset spectrum，不认证 actual FIRST、CPGP 或任何 geom 历史；它只排除忽略真实空间相容性的 posterior-only stopping 引理。

## 3. 保留真实边相容性后的 source-once 公式

令 F_*(x)=max_r T_rf(x)、E={F_*>lambda}；在共同零集之外有限系数均有限，可以取最小 global maximizing r_*(x)。写
 w_A(x)=1_E(x)r_*(x)^|A|(1−r_*(x))^(n−|A|)。
对每个 x，sum_A w_A=1_E；由正 Tonelli 和原 G_A 的对称性，
\[
 \lambda|E|\le\sum_A\int w_A\,G_Af
       =\int f(y)\left[\sum_A G_Aw_A(y)\right]dy. \tag{4}
\]
这是同一来源 f 的一次收款表达，没有把输出所选 source 重领 W。若 E 无限，可先截有限测度 E 并用单调极限。

驻点身份在 receiver 侧是
\[
 \sum_A(|A|-nr_*(x))\,w_A(x)\nu_A(x)=0
\]
（只在内部 winner）。它不把源侧 sum_A G_Aw_A 自动界成 polylog：r_* 是 receiver 选择，卷积内的 r_*(x) 不能换成 r_*(y)。按 (2) 移动一条边后，receiver weight 变为 G_i w_(A∪i)，而不是原 w_A；没有已证的 telescope 或有偿 nonlocal commutator 控制。一般 f 不提供 r_* 的空间正则性。

本轮没有提出足以控制这个 edge-weight 项的新原核引理。对任意 receiver selector 无条件给源侧 pointwise polylog bound 会是一个更强的 strong domination 策略；不能把这一强合同当作弱端点所必需。真正弱证明需要同时利用 f、E={F_*>lambda} 与 (2)，而不是放宽成任意 selector 后付列。

## 4. 终态

已得到明确的方法边界：正 Bernstein、驻点均值、方差和每参数保质量乃至全 input-height cap 都不足以完成弱收费；真实卷积边相容性不能删。上述否定由显式解析构造负责，没有为它再运行代理模型数值，也没有重跑先前原核实验。

当前仍未证明或否定原 A_ord(n)=polylog。下一次若做有限 Fourier/cell 优化，必须从一份 f 生成全部 nu_A，保留原 G_c，并对空间水平集及 Rn lifting 误差负责；本轮 root 建议暂不新增大规模模拟。本文以未闭合状态收口，原主账不变。

后续按 root 明确要求，为§2这个新弱策略障碍另作三轮纯 Fraction 证书：见 [posterior_only_weak_obstruction_20261007.md](posterior_only_weak_obstruction_20261007.md)。n8/32/128共3231项谓词 PASS_EXACT，保留严格阈值、height-cap 版本和连续 peak 的解析身份。它不是原核运行，不是旧 strong D_n 重跑，也没有改变原弱端点未决状态。
