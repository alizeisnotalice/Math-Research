# 保留最低尺度的 capture-block 重排：固定维数界及径向互斥损失

2026-10-07，独立新稿。只读此前 single_capture_source_rearrangement 和 true_winner_log_gram_counter，不改旧稿或主账，不重跑旧数值。使用已读 L03 的预登记/Fraction 组件证书方法；解析自含。

结论分清两个归一化。保留 a 的 block proxy 对固定 n、固定 b/a 有统一界，不能再像旧 log proxy 那样仅令来源标签数增加而发散。但是这个 proxy 相对 I1=β|E| 仍可在真实连续 centered-cube winner 上指数增长，因而不能单凭 capture-overlap Gram 和逐 block 独立重排得到 √n polylog(n)·I1。这不是 proxy/W 或 exact-source-square 的反例，也不认证 original actual FIRST/history 余项。

## 1. 精确 scale-aware 重排系数

沿此前完整背景合同，μ=ρ+β1_Ddx，β>0；当前 receiver 的 selected Q(x,R(x))⊂D，R∈[a,b]，M=β+m_A/R^n。在 capture-block E_A 上完整捕获固定来源集合 A，m_A=ρ(A)>0，v_A=|E_A|。精确系数为
\[
 s_A=\int_{E_A}\frac{\beta\,dx}{\beta R(x)^n+m_A}.
\tag{1}
\]
对任一 y∈A，捕获同时保证 R^n≥a^n 及 R^n≥V_y(x)=(2||x-y||∞)^n；V_y 的 sublevel 体积为 u。和前稿相同的层饼重排，自含给
\[
 s_A\le F_a(v_A,m_A):=
 \int_0^{v_A}\frac{\beta\,du}{\beta\max(a^n,u)+m_A}.
\tag{2}
\]
记 A0=a^n、z=m/β，则准确闭式
\[
 F_a(v,m)=
 \begin{cases}
 v/(A0+z),&0\le v\le A0,\\
 A0/(A0+z)+\log[(v+z)/(A0+z)],&v>A0.
 \end{cases}
\tag{3}
\]
两端在 v=A0 连续；v=0 为0。F 对 v 单调递增、凹，对 m 递减；且
\[
 F_a(v,m)\le\frac{\beta v}{\beta A0+m},
 \qquad F_a(v,m)\le\log(1+\beta v/m).
\tag{4}
\]
因此此前 log 反例的小体积 block 已降到体积线性，并不会复制旧反例。

令
\[
 Q_a=\sum_{A,B}m_{A\cap B}F_a(v_A,m_A)F_a(v_B,m_B)
 =\int\left[\sum_{A:y\in A}F_a(v_A,m_A)\right]^2\rho(dy).
\tag{5}
\]
仅完整捕获 packet/原子可写这一个共同 s_A；partial source 仍按原 kernel 保留，不自动替换为全来源质量。

## 2. 固定维数无 block-count 发散：完整来源一次求和

对每个固定原来源 y，所有包含它的 E_A 不交，且都在 Q_b(y)，因此
\[
 \sum_{A:y\in A}v_A\le b^n.
\tag{6}
\]
用 (4) 和 weighted Cauchy，再对原ρ正 Tonelli：
\[
\begin{aligned}
 Q_a
 &\le b^n\sum_A v_A
          \frac{\beta^2m_A}{(\beta a^n+m_A)^2}\\
 &\le\frac14(b/a)^n\,\beta\sum_Av_A
 \le\frac14(b/a)^n I1,\qquad I1=\beta|E|.
\end{aligned}
\tag{7}
\]
最后一步用 t/(1+t)^2≤1/4。没有来源重启、atom-count、树深度或 block-count 费用；无限 countable 分解可由有限部分与单调收敛。式 (7) 只依 capture/尺度/不交 receiver，甚至不需要 winner 极大性。因此固定 n 与固定 b/a 的旧 J→∞ 反例被排除，不能把“未找到”冒充证明。

每个 diagonal 仍由 mF_a²≤m log²(1+βv/m)≤βv 支付。式 (7) 是交叉的一个一般上界，但 b/a≤2 时为2^n/4，不是 √n 或 n。

## 3. 真实连续赢家的径向分块反例机制

取 n≥2、a=1、b=2、β=τ=1，
\[
 J=2^n,\quad \delta=\frac1{100nJ^2},\quad
 t_j=1+j/J\ (0\le j\le J),\quad
 \rho=\delta_0+\frac{\delta}{J+1}\sum_{j=0}^J\delta_{t_je_1}.
\tag{8}
\]
背景取完整固定盒 D=[-2,4]×[-2,2]^{n-1}，μ=ρ+1_Ddx，全部来源质量
\[
 W=1+\delta+6\cdot4^{n-1}.
\tag{9}
\]
这里的 central mass 为1，所有辅助来源真实质量总和为δ，不能删去或把背景费用免费化。

研究 receiver cone
\[
 x_1>0,\quad |x_i|<x_1\ (i>1),\quad
 R0=2x_1\in(1,2).
\tag{10}
\]
所有尺度 query 在 D 内。中央来源第一次被捕获的尺度恰 R0。在其前，spike 捕获质量≤δ，故 signal≤δ<2^{-n}；中央到达的 signal≥R0^{-n}≥2^{-n}，严格更大。

对 j=0,…,J-1，只保留同一实际输入的 receiver 子区间
\[
 t_j<R0< u_j:=t_{j+1}/(1+\delta).
\tag{11}
\]
这些区间非空。R0 时恰完整捕获中央来源及辅助点 t_0e1,…,t_je1，记此真实集合 A_j，质量
\[
 m_j=1+(j+1)\delta/(J+1)\le1+\delta.
\tag{12}
\]
若后续还捕获一个新辅助点，其位置 t_k≥t_{j+1}>R0(1+δ)，它第一次到达需尺度
\[
 R'=2t_k-R0>R0(1+2\delta).
\]
即使把全部 spike mass1+δ 都计入，后续 signal≤(1+δ)/(R')^n<R0^{-n}。其中 (1+2δ)^n>1+δ。两次到达间固定捕获质量的响应严格下降。于是整个连续 [1,2] 的严格唯一 winner 就是 R0；其 query 完整捕获 A_j。不是自由 selector，也不是仅在一个离散候选点发现局部峰值。闭面约定使中央点在 R0 被捕获，receiver 间隔端点均可忽略。

在 cone 中，以 R0=s 作变量，截面体积为 s^{n-1}、dx1=ds/2。因此 (11) 的真实 receiver 体积准确为
\[
 v_j'=\frac{u_j^n-t_j^n}{2n}>0.
\tag{13}
\]
其对应全局 capture-block E_{A_j} 的体积至少 v_j'。由均值公式及 t_{j+1}-t_j=1/J，
\[
 v_j'\le\frac{t_{j+1}^n-t_j^n}{2n}
 \le\frac{2^{n-2}}J=\frac14.
\tag{14}
\]
用 1-(1+δ)^{-n}≤nδ，
\[
\begin{aligned}
 \sum_jv_j'
 &\ge\frac{J-1}{2n}
       -\frac{\delta}{2}\sum_jt_{j+1}^n\\
 &\ge\frac{J-1}{2n}-\frac1{200n}
 \ge\frac{J-1}{4n}.
\end{aligned}
\tag{15}
\]
由 (3)，对这些小体积及 m_j≤1+δ<2，
\[
 F_1(v_j',m_j)=v_j'/(1+m_j)\ge v_j'/3.
\]
全局 block 的体积可能更大，但 F 递增，所以不会减少这一个已核下界。所有 A_j 共有同一个中央原来源，质量1；故
\[
 Q_a\ge\left[\sum_jF_1(v_j',m_j)\right]^2
 \ge\frac{(J-1)^2}{144n^2}.
\tag{16}
\]

完整严格水平集 E={max_[1,2] h_L*μ>1} 不是人为限制为 cone。背景响应处处≤1，所以 E 中必须有 spike 被捕获。所有 spike 位于线段[0,2]e1，必有
\[
 E\subset[-1,3]\times[-1,1]^{n-1},\quad
 I1=|E|\le2^{n+1}=2J.
\tag{17}
\]
此空间盒中所有 query 也在 D，背景实际为1；因此完整 E 实际等于来源的 side2 cubes 的并集。来源第一坐标0、1、1+1/J,…,2 的间距均≤1，side2 区间覆盖整个 [-1,3]，其余坐标都是 [-1,1]。所以原子闭面约定下精确 E=[-1,3]×[-1,1]^{n-1}、I1=2J；本稿只需 (17) 的方向。完整 E 还保留所有 cone 外输出，不删原源。

结合 (16)/(17)：
\[
 \boxed{\frac{Q_a}{I1}\ge
 \frac{(J-1)^2}{288n^2J}
 \ge\frac{2^n}{1152n^2}.}
\tag{18}
\]
因而真实 centered-cube 连续赢家并不能让这个 scale-aware block proxy 相对 I1 满足 O(√n polylog n)，甚至不存在多项式维数统一界。本反例不是固定 n 后 J→∞；它与 (7) 完全相容。

## 4. L1 范围及不被否定的预算

式 (8) 为合法完整有限测度输入，原核/连续赢家/阈值完全真实。若需要严格 L1，可把所有点换为半径 ε 的均匀小方盒，并仍保留整个背景。固定 n,J,δ 后，先在 (11) 的每个 receiver cone 内取一个闭紧子集，避开 |x_i|=x1、径向两端及任意新辅助捕获的边界。这些子集体积可从下逼近 v_j'。

在此紧子集，中央小方盒被部分捕获时，其他 n−1 方向已完整包住；中央捕获量对 L 的导数为1/(4ε)。其 signal 导数为 [L·mass'(L)-n·mass(L)]/L^{n+1}，而总质量≤1+δ，L≥1。只要 ε<1/[4n(1+δ)]，中央 partial-band signal 严格上升，完成后直到下一 helper 到达严格下降。中央完整端为 Rε=R0+2ε；全部 helper 捕获状态在紧子集稳定，后续候选仍有 (11) 给出的严格统一 signal gap，选择 ε 充分小即保持中央完整端为严格全 winner。

于是同一固定背景及真正 L1 来源，可恢复几乎全部 v_j'，并保持共同中央 packet mass1。所有 F 对 v,m 连续，故选择一个充分小但有限的 ε，可保留 (16)/(18) 的任意预定稍小常数（例如下界的一半）。E 的空间盒增加 ε，体积上界趋于2J；没有只以原子点态极限假称赢家稳定。此 L1 迁移是上述解析稳定性，不在新有限守卫中模拟原 packet 流。

需要清楚区分：

- 背景 W≈6·4^{n-1} 真实存在。(16) 不否定 Q_a≤√n polylog(n)·W；本构造下 Q_a/W 的已核下界很小。尚无这种 W 归一化反例或证明。
- 精确 source profile 仍≤1+n log2；完整 ∫S dμ=I1，所以 ∫S²dμ≤(1+n log2)I1。独立重排后的指数 proxy 大不说明精确平方指数大。
- 本构造的中央来源在 (12) 占主要份额，可能完全落入旧已付 dominance 子支；本稿不认证 actual R_angle 的非集中/lowS/lowcoin/FIRST/history 门。
- 不宣称共同有限尺度族版本；(11) 内的赢家尺度 R0 随 receiver 变动，本轮连续合同足够排除 proxy 通用维数猜测。

## 5. 剩余一般接口

最低尺度恢复了固定维数的标签稳定性，但每个 block 可以各自把 v_j' 放到来源周围最内层 a-cube；这些相互冲突的径向重排并不对应原真实 E_A 的共同空间位置。反例中真实 block 位于不同 R0 壳，proxy 在 (3) 中却把每一块都按 R=1 计算。辅助来源总质量可以任意小，仍能产生真实不同 capture 标签。

所以若目标是 source-tail 的 √n·I1，一般下一步必须保留同一 y 的联合径向占用、或者直接保留原 R(x) 的 exact c_A；capture-overlap 与独立体积系数不足。逐源只删质量分母并联合积分原尺度，得到上面的已知 n 级包络，仍不自动改善到 √n。本稿没有把这个必要修复写成已证明 √n 新合同。

## 6. 新证书预登记与范围

新三轮 n=4,8,12，J=2^n，无随机种子；注册后执行 Fraction。只检查新的原径向分块构造：连续赢家的前后 signal margin、真实 v_j'、小体积界、质量/公共来源、总 retained volume 及 (16)–(18) 下界。所有参数 t_j、δ、u_j、v_j' 为有理数；没有 nth-root 误差、receiver 网格、旧 log 反例或旧 fixed-K 数据。

这是原 hard 连续赢家的解析组件守卫，不是抽象 incidence toy；L1 迁移及 all-n 证明由上文解析给出，数值不冒称 actual FIRST。

## 7. 终态与计算范围

三轮新守卫 52,453 项 Fraction 检查 PASS（n4/8/12 各204/3,084/49,164，另注册1项），0.111秒。无 live handle、无随机数，无旧 log 反例/旧 fixed-K 结果重跑。

| n | J | 已证 Q_a 下界除以完整 I1=2J |
|---|---:|---:|
| 4 | 16 | 25/8192 |
| 8 | 256 | 7225/524288 |
| 12 | 4096 | 207025/2097152 |

这些小维数下界本身不足展示渐近 √n 失效；排除统一多项式维数合同的是 (18) 的所有 n 解析指数下界，而不是数据外推。结果还保存实际 retained cone 总体积的有理摘要、每个真实前后 signal gap 和公共来源平方的下界摘要。不会为逐 block 不同分母建立无必要巨大 LCM。

第一次新执行尝试在构造输出摘要时触发 Python 的4300位整数转字符串限制，尚未写终态 JSON；不是数学残差失败。修正为保留逐 coefficient≥v/3 的已核残差，再对同一原来源用 [Σv/3]² 作有理 lower certificate，避免计算全部 prefix 分母的巨大 LCM。最终执行和脚本 SHA 如下；没有覆盖旧任务结果。

- 脚本：[scale_aware_capture_rearrangement_20261007_exact_guard.py](scale_aware_capture_rearrangement_20261007_exact_guard.py)，SHA-256：7dd20ff095ef80a23abd6bc5b177b74f89881666bb4c5f499b63e4185b210f38。
- 预登记：[scale_aware_capture_rearrangement_20261007_registration.json](scale_aware_capture_rearrangement_20261007_registration.json)。
- 终态：[scale_aware_capture_rearrangement_20261007_results.json](scale_aware_capture_rearrangement_20261007_results.json)。

新结果只关闭“最低尺度足以使 independent block proxy 有 √n·I1 预算”这个过强路线；完整 actual 空间付款及 W 归一化 proxy 仍未闭合。
