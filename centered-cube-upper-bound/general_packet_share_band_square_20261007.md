# 一般完整来源的真实份额 band 平方函数

2026-10-07。本稿独审父任务提出的弱化：保留 packet fullness，去掉固定 dominance 门，按原真实捕获份额分 band。结论是所有 band 的 **平方函数** 常数预算；不是全部 band 相加后的 square 预算，更不是一般 cube 弱上界。沿用已核 [无背景 packet 引理](general_source_packet_square_20261007.md) 的 A01 齐次/source 口径；全部证明在下文给出，不调用新外部定理。

查重：旧 [single capture §4.2](single_capture_source_rearrangement_20261007.md) 给 background 下固定 θ 的 1/(ηθ) packet square；旧 [gated occupancy §7](gated_occupancy_bridge.md) 给有限硬原子 N 的交通费 O(log N) 与 small-share 吸收。二者都不是以下任意完整 μ、fractional packets 的份额平方函数汇总。此前 [覆盖反证](packet_coverage_feasibility_20261007.md) 说明不能以 θ=n^(-1/2) 的固定门覆盖所有来源；本稿保留其 small-share 补集。

## 1. 完整原对象与 band 端点

设 μ 是有限正 Borel 测度，质量 W>0。取合法原 receiver E、原可测边长 R(x)>0，m(x)=μ(Q(x,R(x)))，并仅假设

\[
 m(x)>\tau R(x)^n,\qquad I_1=\tau|E|<\infty. \tag{1}
\]

这包含原真 finite winner E={Mμ>τ}；证明无需改变赢家。预固定 dμ_i=a_i dμ，a_i≥0，Σ_i a_i≤1，M_i=μ_i(R^n)>0。记 m_i(x)=μ_i(Q(x,R(x)))，p_i(x)=m_i(x)/m(x)，于是 Σ_i p_i(x)≤1。

固定 η∈(0,1]。对 j=0,1,…，取 θ_j=2^(-j-1)，并定义互不重叠的真实资格

\[
 E_{i,j}=\{x\in E:m_i(x)\ge\eta M_i,
                        \ \theta_j<p_i(x)\le2\theta_j\}. \tag{2}
\]

右端闭、左端开：p=2^(-k) 属于 j=k，p=1 属于 j=0。所有正份额恰属一个 band；fullness 资格保证 m_i>0，故不存在 p=0 的丢失。对来源 y 取原 incidence profile

\[
 S_{i,j}(y)=\tau\int_{E_{i,j}}
       \frac{\mathbf1_{Q(x,R(x))}(y)}{m(x)}\,dx,
 \quad S_j(y)=\sum_i a_i(y)S_{i,j}(y),
 \quad S_{\rm full}(y)=\sum_{j\ge0}S_j(y). \tag{3}
\]

这些对象用同一 μ、同一原 R、同一 τ。非负 Tonelli 合法处理 countable packets、bands；没有 background β1_D，也未重新领取来源质量。

## 2. 单 packet 的 band 上端给更大分母

对 x∈E_i,j，由 fullness 与份额上端，

\[
 m(x)\ge\frac{\eta M_i}{2\theta_j}=:A_{i,j}. \tag{4}
\]

若 y∈Q(x,R(x))，则 R(x)≥2||x−y||∞。与 (1) 合用，得到

\[
 \frac\tau{m(x)}\le
 \frac{2\tau}{A_{i,j}+\tau(2\|x-y\|_\infty)^n}. \tag{5}
\]

令 V(x,y)=(2||x−y||∞)^n；Lebesgue 体积 |{x:V(x,y)≤v}|=v。因此对任意可测 F，递减径向重排/层蛋糕给

\[
 \int_F\frac{2\tau\,dx}{A+\tau V(x,y)}
 \le2\log\left(1+\frac{\tau|F|}{A}\right). \tag{6}
\]

这里仅把原 incidence 的正上包重排，未重新选择 R 或 query。用 log²(1+t)≤t（由 log(1+t)≤sqrt(t)）并代入 (4)，得到

\[
 \int S_{i,j}(y)^2\,d\mu_i(y)
 \le4M_i\log^2\left(1+
                   \frac{2\theta_j\tau|E_{i,j}|}{\eta M_i}\right)
 \le\frac{8\theta_j\tau}{\eta}|E_{i,j}|. \tag{7}
\]

因此固定 θ 下的 1/θ 损失不再出现；这依赖真实 band 的上端。只有 p_i≥θ 而没有 p_i≤2θ 时，(4) 不成立。

## 3. 所有 band 的平方函数可以一起付

令

\[
 A_j=\tau\theta_j\sum_i|E_{i,j}|.
\]

在每个原 x，band 互斥且 θ_j<p_i 给

\[
 \sum_{j,i}\theta_j\mathbf1_{E_{i,j}}(x)
 \le\sum_i p_i(x)\le1,
 \qquad\sum_{j\ge0}A_j\le I_1. \tag{8}
\]

对每个 y，用 Σ_i a_i(y)≤1 的 fractional Jensen（剩余权重乘零），再用 (7)：

\[
 \int S_j^2\,d\mu
 \le\sum_i\int S_{i,j}^2\,d\mu_i
 \le\frac8\eta A_j.
\]

非负 Tonelli 与 (8) 于是给真正的统一平方函数预算

\[
 \boxed{\int\sum_{j\ge0}S_j(y)^2\,d\mu(y)
 \le\sum_{i,j}\int S_{i,j}^2\,d\mu_i
                 \le\frac8\eta I_1.} \tag{9}
\]

这比逐 band 粗写 ≤8I₁/η 后直接乘 band 数更准确。它对任意维数、任意原来源质量及 countable packets 成立；不需要 input 坐标独立。若一组 labels 的源全部相同地点或同 source 被所有 bands 反复捕获，(9) 仍成立；但这一事实不使 cross terms 消失。

## 4. 总 profile 的同源交叉仍在

S_full=Σ_j S_j，而

\[
 \int S_{\rm full}^2\,d\mu
 =\int\sum_jS_j^2\,d\mu
       +2\sum_{j<k}\int S_jS_k\,d\mu, \tag{10}
\]

等式按扩展非负积分理解。band 在 receiver/label incidence 上互斥，不能据此宣布它们在来源 y 上正交。

若至多 L 个 bands 有效，或对 μ-a.e. y 至多 L 个 S_j(y)>0，则 Cauchy–Schwarz 给

\[
 \int S_{\rm full}^2\,d\mu\le\frac{8L}{\eta}I_1. \tag{11}
\]

是 L 而不是 L²。令 q_full(x)=Σ_i1_{m_i≥ηM_i}p_i(x)，则原 first moment 准确为

\[
 \int S_{\rm full}\,d\mu
       =\tau\int_Eq_{\rm full}(x)\,dx. \tag{12}
\]

若该平均 fullness≥κ>0，(11) 才进一步给

\[
 I_1\le\frac{8L}{\eta\kappa^2}W. \tag{13}
\]

对有限个正质量 packets，可用 m≤W 和 fullness 得 p_i≥ηM_i/W，故若 M_min=min_i M_i>0，安全地取 L≤1+floor(log₂(W/(ηM_min)))。这依赖质量深度，不是维数无关常数；label 数有限也不意味着 L≤label 数。countable packets 的 M_i 可趋零，没有维数控制的 L。微细来源网格能保证平均 fullness 接近1，却也可能把这项质量深度放大。

另一个精确的条件合同是：若 B=Σ_j sqrt(A_j)<∞，用权重 w_j=sqrt(A_j)/B（A_j=0 的 band 由 (9) μ-a.e. 为零）给

\[
 \int S_{\rm full}^2\,d\mu\le\frac8\eta B^2. \tag{14}
\]

仅 ΣA_j≤I₁ 不能控制 B²/I₁；(14) 不是新 paid fee。两种合同分别明确把缺口留在 effective band 数或 band 分布的半阶量，而没有用“平方函数”重命名所缺总能量。

**原幅度带的 packet-relative baseline。** 若进一步保留原 τ<Mμ(x)≤2τ 的幅度带以及 R(x)∈[a,2a]，则每个固定 i 的 fullness incidence 都满足

\[
 \frac{\eta M_i}{2\tau(2a)^n}\le p_i(x)
 \le\frac{M_i}{\tau a^n}.
\]

上下端比≤2^(n+1)/η，故这个 i 自己所用的 band 数至多

\[
 L_\eta=n+3+\lceil\log_2(1/\eta)\rceil. \tag{15}
\]

不同 i 的 band 起点可随 M_i 任意平移；无需所有 i 共用一个有限 band 区间。先逐 i 对 Σ_j S_i,j 用 (15) 的 Cauchy，再对 Σ_i a_i 用 fractional Jensen，最后用 (9) 的 **labelled** 预算，准确得到

\[
 \int S_{\rm full}^2\,d\mu
 \le\sum_i\int\left(\sum_jS_{i,j}\right)^2d\mu_i
 \le L_\eta\sum_{i,j}\int S_{i,j}^2d\mu_i
 \le\frac{8L_\eta}{\eta}I_1. \tag{16}
\]

所以在这个原幅度带内，即使 M_i→0，也有 O((n+log(1/η))/η)I₁ 的 total square baseline；若平均 fullness≥κ，则 I₁≤8L_ηW/(ηκ²)。单凭 (1) 没有上幅度带时不能用此结论。本文此前质量深度合同 (11) 仍合法，但不是原幅度带内唯一可用的基线。要推进平方根，须从真实 geometry/gates 省掉这项同-source 跨 band 的 n 收费；(9) 自身尚未做到。

## 5. 对 actual geom 的边界

本文给了一般完整 μ 的合法新弱化 (9)，并保留真实 small-share 来源。它使用 capture、fullness、原 m>τR^n；未使用 FIRST/fullfuture/CPGP/LCA/all-history。若另有真证明把 actual 余项的来源-incidences 映到 (2) 且保留 coefficient，再有优于 (16) 的可付跨 band 预算或 (11) 的 L=O(sqrt(n)n^{o(1)})，才有主目标级费用。当前都未由 (9) 自动获得。

本稿没有新的 endpoint 阶数猜测，未执行 toy 或重跑旧数值。前一覆盖稿的三轮49有理守卫只核其 packing 反证，不充当本稿空间/gate 或跨 band 预算的数值证书。原主账与一般目标状态不变。
