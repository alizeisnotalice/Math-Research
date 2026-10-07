# 同来源最优接受耦合：确定性势、Palm 截断与未付几何

2026-10-07。只新增本稿，不改其它文件或主账。按父任务 6.1-sol high 要求工作；当前工具不公开实际运行模型 ID。本轮自含证明有限核耦合及可测性，不新增数值或阶估计，也不重跑旧守卫。

**准确结论：给定原来源 y，所有接受概率 p_j(y) 的总和 S(y) 可以被一个 balanced 圆周耦合实现。** 此时 N∈{⌊S⌋,⌈S⌉}，来源一次使用的最大 union 为 min(1,S)，最小重复为 (S−1)+。但没有证明原 actual 的 S 高层质量小；只是将最有利的接受次数合同写成了明确的确定性来源势合同。

沿用已读 [E01](/Users/zhengzhihao/.codex/skills/math-e01-occupation-measure-flow-lp/SKILL.md) 的来源守恒约定及 [M04](/Users/zhengzhihao/.codex/skills/math-m04-measurable-supremum-selection/SKILL.md) 的可测版本边界。只用有限 measurable kernels、正 Fubini 和标量凸性。没有 argmax 选择、无限 LP 强对偶或条件 weak 平均。

## 1. 原 finite index joint kernel 的假设

固定完整原输入 μ、0<W=μ(Rⁿ)<∞，先冻结实际 receiver selector、来源切分及所有门。取有限 j=1,…,J。令 Q_j(y,dξ) 是 measurable probability kernel，ξ 属于一个 standard Borel 标记空间，包括 receiver x 及需保留的原 index 内标签/history。原 weak 权和 RN 接受给 measurable α_j(y,ξ)∈[0,1]，定义

\[
 A_j(y,dξ)=α_j(y,ξ)Q_j(y,dξ),\quad
 R_j=Q_j-A_j,\quad
 p_j(y)=A_j(y,\Xi_j),\quad S(y)=\sum_jp_j(y).
 \tag{1}
\]

所以 0≤p_j≤1、0≤S≤J。核积分的 measurability 可先对矩形指示验证，再用 simple functions 和正单调极限，故 p_j、S 都 measurable。若只给 accepted subprobability A_j、未给 rejected prior，可补一个 cemetery state 完成 Q_j；此时不宣称 cemetery 是原 rejected history。以下保持原 rejected joint 的结论以已给 Q_j 为前提。

actual 中 A_j 必须是原 joint 接受核，而非任意自由 a≤1 的代理。保留 FIRST/fullfuture、原 source cell、首跳 t/c_t/落点/退出、hard z、CP/GP、共享 seed、唯一 LCA、失败权、early/allhistory、当前全部删支及弱权；标记必须按原 source y,z 判断。原完整 μ/W 是本稿的来源概率，不借用旧障碍 ν_b 的 cap 或 signed 势。

## 2. 一个圆周变量的可测 balanced 耦合

置 s_0(y)=0、s_j(y)=Σ_{i≤j}p_i(y)，取唯一一个 U∼Uniform([0,1))，独立于原 y。定义

\[
 B_j(y,U)=\sum_{k=0}^{J-1}
   \mathbf1_{\{s_{j-1}(y)\le U+k<s_j(y)\}},\qquad N=\sum_jB_j.
 \tag{2}
\]

每个区间长 p_j≤1，半开约定保证最多包含一个 U+k，所以 B_j∈{0,1}；即使长度为 1 也不会重复计边界。公式是有限个 measurable 不等式的指示之和，故 (y,U)↦B_j measurable。其 marginal 精确为

\[
 \int_0^1B_j(y,u)du
 =\sum_{k=0}^{J-1}|[s_{j-1},s_j)\cap[k,k+1)|=p_j(y).
 \tag{3}
\]

所有区间连续铺满 [0,S)，故

\[
 N(y,U)=\sum_{k=0}^{J-1}\mathbf1_{0\le U+k<S(y)}
 =m(y)+\mathbf1_{U<θ(y)},\quad
 m=\lfloor S\rfloor,\ θ=S-m.
 \tag{4}
\]

S=0、S=J、θ=0 都包含在这个半开公式中。N 的 conditional mean 为 S；有 θ 的概率取 m+1，其余取 m。

圆周上的每个 arc 是原连续区间的 mod 1 像。长度不超过 1 是关键；不对长度大于 1 的单个 index 宣称 Bernoulli。

## 3. 随后抽 endpoints/history，保持原每个 joint

在 p_j>0 处定义 Q_j^+=A_j/p_j，在 p_j<1 处定义 Q_j^−=R_j/(1−p_j)。零分母处任选 Q_j 作为 unused fallback。对每个 measurable D，这些概率核的值是 measurable quotient，加上 measurable 分支；核仍是 measurable。

给定 (y,U)，按 B_j 取 Q_j^+ 或 Q_j^−，再条件独立抽 ξ_j。有限乘积 kernel 可由矩形上的概率乘积及逐次核积分构造；没有需要选择一个不可测 endpoint。于是

\[
 \Pr(B_j=1,\ ξ_j\in D\mid y)=p_jQ_j^+(D)=A_j(D),\quad
 \Pr(B_j=0,\ ξ_j\in D\mid y)=(1-p_j)Q_j^-(D)=R_j(D).
 \tag{5}
\]

每个 index 的 accepted/rejected/source/receiver/history 联合边际都准确保留。p_j=0 或 1 时所用分支也准确，unused fallback 不生成额外 accepted mass。

代价是跨 index copula 被改写。没有保持原 Poisson field、同一 U_sensor 或原路径 histories 之间的联合 law，除非另按第 7 节条件化；也没有把原 FIRST 变成新辅助观察次序的停时。原 shared seed 可以作为每个 index 的标记保留，但其原跨 index 共同随机变量不能凭 (5) 宣称仍有相同 copula。

## 4. 来源一次 union 和重复的全耦合最优性

任意保持 (1) 的同来源 coupling，条件于 y 都有 E N=S。因为 N 为非负整数，

\[
 \Pr(N>0\mid y)\le\min(1,S),\quad
 \mathbb E[(N-1)_+\mid y]=S-\Pr(N>0\mid y)\ge(S-1)_+.
 \tag{6}
\]

后式也由凸函数 (n−1)+ 的 Jensen 不等式得到。(4) 同时达到两界：S<1 时 N∈{0,1}；S≥1 时 N≥1。于是最优 once-used source 与最小重复为

\[
 d\nu_{\rm union}^{opt}(y)=\min(1,S(y))μ(dy)\le μ(dy),\qquad
 D_{opt}=\int(S(y)-1)_+μ(dy).
 \tag{7}
\]

更强地，对每个整数 K≥1，(4) 同时达到

\[
 \mathbb E[(N-K)_+\mid y]=(S-K)_+,\qquad
 \mathbb E[\min(N,K)\mid y]=\min(S,K).
 \tag{8}
\]

这两函数在相邻整数之间是 affine；(4) 的两点分布取其准确线性插值。任意其它 coupling 的 overflow ≥(S−K)+，truncation ≤min(S,K)。balanced coupling 因而最小化所有这类 convex overflow；它**不**自动最小化每个 indicator tail Pr(N>K)，后者不是凸函数。

一个更一般的 convex-order 说明：对定义在整数上的凸 φ，其 linear interpolation φ_lin 凸，φ_lin(N)=φ(N)，所以 Eφ(N)≥φ_lin(S)；(4) 取等号。不得把这个次序推广为所有非凸尾的最优性。

任何 coupling 都保持同一个 S 和 M=∫S dμ，因此优化 copula 不会缩小实际交通或原弱水平集的 volume。最优重复仍有 D_opt≥M−W；只有独立取得 S 的几何控制，才可能把这份最有利分解变成预算。

## 5. 两种 actual 校准和准确 Palm 来源

令

\[
 M=\int S(y)μ(dy),\qquad A=M/W=\mathbb E_μN.
 \tag{9}
\]

M=0 时没有 accepted Palm，下面仅在 M>0 使用。仍由 (5) 保持原 receiver 边际。沿 logistic actual 稿的 normalization C_h=16/3，有两种不同校准：

* 强交通权 g=1_E：C_hM=∫_E T(x)dx；没有逐行 lower threshold 时，不能写为 λ|E|。
* 真正弱归一化 g=τ/T·1_{T>τ}：C_hM=τ|E_τ|，Palm receiver 为 uniform Lebesgue|E_τ|。原 finite index partition不被重判。

这些是 (1) 中原 α_j 的不同输入；不得相互替代。定义 accepted-occurrence Palm

\[
 P^\#(dω,j)=\frac{B_j(ω)}A P_μ(dω),\qquad
 P^\#(dy)=\frac{S(y)}M μ(dy).
 \tag{10}
\]

Palm source 是 S 加权的原来源，通常不是 μ/W。对于 balanced coupling，其 conditional N law 为 m(1−θ)/S 在 m、(m+1)θ/S 在 m+1；m=0 的 0 项为零，S=0 的来源没有 Palm mass。

因此准确

\[
 \mathbb E^\#\frac1N=\frac{\int\min(1,S)dμ}{M},\qquad
 \mathbb E^\#\frac{(N-K)_+}N=\frac{\int(S-K)_+dμ}{M},\qquad
 \mathbb E^\#\min(1,K/N)=\frac{\int\min(S,K)dμ}{M}.
 \tag{11}
\]

K=1 时给最优 stop/deficit 的准确分解 M=∫min(1,S)+∫(S−1)+。第一项≤W；第二项没有来源守恒自动上界。

## 6. Palm 合同等于哪个确定性截断合同

整数 K≥1 的最小 overflow 合同为

\[
 \int(S-K)_+dμ\le δM,\quad0\leδ<1.
 \tag{12}
\]

因为 M=∫min(S,K)+∫(S−K)+≤KW+δM，所以

\[
 M\le\frac K{1-δ}W,\qquad
 τ|E_τ|\le\frac{C_hK}{1-δ}W
 \quad\text{（只在弱校准下）}.
 \tag{13}
\]

(12) 是确定性 source potential 的截断质量合同；它在 balanced coupling 下等价于 (11) 的 Palm truncation≥1−δ。这不是已知预算，只是一条明确充分条件。希望 K=√n n^{o(1)} 仍待证。对强交通校准，(13) 的同一推导结论为 ∫_E T≤C_hK W/(1−δ)，而非原 λ|E| 直接付款。

若坚持原 low-count Palm 形式，则其精确 deterministic tail 为

\[
 H_K(s)=\begin{cases}
 0,&s\le K,\\
 (K+1)(s-K),&K<s<K+1,\\
 s,&s\ge K+1,
 \end{cases}\qquad
 \Pr^\#(N>K)=\frac{\int H_K(S)dμ}{M}.
 \tag{14}
\]

证明是对 (4) 的两点 law 计算 E[N1_{N>K}|y]，再作 Palm 换测度。粗而可验证的 sufficient source condition ∫_{S>K}S dμ≤δM 推出 (14)≤δ。也可直接用 (12)，无需强 pair 或 low-count indicator tail。H_K不是凸函数，不声称 balanced coupling 是任意 copula 下 (14) 的最小值。

## 7. 若要保留原 Poisson fields，最优性须条件化

设 F 是固定的原共同 fields、U_sensor 或其它需保留 copula 的随机环境，且每个 index 有 measurable conditional kernel Q_j(y,F,dξ) 和 conditional accept probability p_j(y,F)。把 (2)–(5) 应用于 (y,F)，得 N∈{⌊S_F⌋,⌈S_F⌉}，S_F=Σ_jp_j(y,F)。先按原 law 抽 (y,F)，再作 conditional circle coupling，保持 **F 与每个 index 内历史的 joint marginal**；fields 的原正增量 law 保留。

这时最小 conditional repeated source为

\[
 D_{F,opt}=\int \mathbb E[(S_F-1)_+\mid y]μ(dy)
 \ge\int(S(y)-1)_+μ(dy),\quad S(y)=\mathbb E[S_F\mid y].
 \tag{15}
\]

一般 K 的 overflow 也有同向 Jensen，conditional union 则 E[min(1,S_F)|y]≤min(1,S)。它是“保持环境且保持每个 index conditional marginal”的耦合类中最优，未必达到 source-only 的全耦合最优。

若保留环境，(12) 应改为 ∫b_K^F(y)μ(dy)≤δM，其中 measurable b_K^F(y)=E[(S_F−K)+|y]；对应 Palm truncation 为 ∫E[min(S_F,K)|y]dμ/M。不能把 b_K^F 换成通常较小的 (S−K)+。

source-only (7) 较有利，但可能破坏原共同 field/cube sensor 的跨 index 几何；conditional-field 版本保留这些几何，却留下随机 source potential S_F 的更大 overflow。不能把 source-only 的低重复费与 conditional-field 的 Markov/copula 性质同时免费使用。原 FIRST 和字段中的辅助 first-accept 仍不同。

## 8. 原 fullfuture 对来源势具体给了什么

在 finite actual 参数 θ_j=(σ_j,L_j)、原不交 receiver cells E_j 下，source-only 的准确势可写

\[
 S(y)=\sum_j\int_{E_j}
       g(x)a_j(y,x)p_{σ_j,L_j}(x-y)dx,
 \tag{16}
\]

其中 a_j 是整份原 accepted history 条件密度，含硬来源 z/seed/LCA 的原条件积分；不能取 arbitrary mask 的 supremum。原 fullfuture 是 receiver 行帽

\[
 F_{s,L}μ(x)\le q(x),\qquad s\geσ(x),\ a\le L\le b.
 \tag{17}
\]

它对 selected 参数给 ∫p_selected(x−y)μ(dy)=q(x)。因 a≤1、g≥0，正 Fubini 仅给

\[
 M\le\sum_j\int_{E_j}g(x)F_{σ_j,L_j}μ(x)dx
      =\int_Eg(x)q(x)dx.
 \tag{18}
\]

有限 rows 在 (17) 的实际参数上时此式准确；finite approximations 若只有冻结原标记，exact cap须另迁移。沿原 C_x logistic 曲线的未来帽仍是一份 receiver 行约束。它没有给 (16) 固定 y 的被不同 receiver/index 接受次数，也没有给 high-source-tail ∫(S−K)+dμ；(18) 右边仍含未知 receiver volume/weight。

一个清楚的正检验：若全部 p_{σ_j,L_j} 恰为同一 probability convolution p，则 disjoint E_j、Σ_jg1_{E_j}≤1、a_j≤1 直接给 S(y)≤∫p(x−y)dx=1。此时 union确可一次支付。moving soft/hard 使 kernel 随 j 变，原分区只在 receiver x 上不交，不能移到同一个 source kernel 下求和。logistic 场的偏序增量或 (Y,U) 几何可能限制这一移动，但该限制尚未证明。

本稿未构造符合全部实际门的 fullfuture 反例，所以不声称 (17) 与 actual geometry 无法推出 (12)。当前精确的最小几何缺口是：在**原 frozen actual accepted kernel** 上，控制 (16) 的 source-weighted overflow，或者在保留 original fields 的版本上控制 E[(S_F−K)+]；需 K=√n n^{o(1)}、δ<1，并对真实 receiver 参数及 grid approximation一致。它不是自由 field 的 path Doob 界，也不是每 C 独立 weak 界的平均。

对 arbitrary RN mask 若用完整自由 kernel 的 source potential作为上界，可以推出一个更强 sufficient lemma；但必须标明丢掉哪些实际门，不能把其未知 strong envelope费登记成 actual预算。原硬 u/q 帽、FIRST、首跳 early、fullfuture、CP/GP/LCA/allhistory 应保留到具体列占用几何证明中。

## 9. 完成状态

已严格证明有限核的 measurable balanced coupling、每 index accepted/rejected/history 边际保持、source-once union与全部整数 overflow 的最优性、Palm deterministic truncation身份，以及 conditional-field 的 Jensen 代价。

未证明 (12) 的实际平方根参数，也未取得原 R_angle 新 paid 费用。此轮没有新阶数或数值候选，按父任务要求不新增 toy 和重复守卫。旧一般 O(n log n) 基线与主账保持原状态。
