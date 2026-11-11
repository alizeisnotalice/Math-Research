# r088 reflection-center columns独立审计

结论：`core88.md` 的 affine patch 上下界与一般 center-compatible residual 公式均通过审计。校准结论只针对给定 affine source 的 patch `X`；它不提供一般输入的 full-soft budget 上界，原稿已正确保留该限制。

## 1. 算子与 lower bound

`A_{c,r}` 的另一端为 `(1+r)c-rx`，故对任意以 `x` 为角点的矩形都有唯一 `c=(a+rx)/(1+r)`；切分一般包含 `x` 的矩形给 corner pieces，因而 `sup_c A_{c,r}=M_S`。固定 `c` 时只是 r028 fixed-ratio operator 的平移，是真实 common-source column，不是 r075 那种随 receiver 重锚定的 source 变换。

在 affine example 中，`t=log M_SF=1/2-sum z_i/N`、
`C(z)=(v+rx)/(1+r)` 正确。`dx/dv<0` 的绝对值小于 `2/N`，所以
\[
|C'(z)|=\frac{g(v)|1+r,dx/dv|}{1+r}>\frac13
\]
对所有 `r in [0,1]` 一致；因此从 center ball 拉回 `z` 的 Jacobian 损失至多 `3^N`。

深度 slab 的平移和 central-unimodal density 计算正确：`P_b` 上
\[
\frac32N^{-b}\le t\le\frac52N^{-b},\qquad |P_b|\ge4^{-N}N^{-b}.
\]
`b=2` 确实达到 flat saturation；`b=3` 的 `chi` 下界为
`q_N>=e^{-2N/3}/2`。固定中心成功时，`a_i<1`；clip 后的平方亏损界给
\[
\sum_i(a_i^+-v_i)^2<6N^{1-b}.
\]
这小于任何一个 clip 到零所需的 `v_i^2>=1/16`，故实际上没有 clip；再得实际边长至少 `1/8`、`J(G)>8^{-N}`。单列成功截面是一个 center ball 的 inverse image，体积至多 `(C N^{-b/2})^N`。Tonelli 对 countable weights 与实际单列费用给出饱和的 `N log N+O(N)` 和 flat-soft 的 `3N log N/2+O(N)` lower bounds；没有混入原 r084 的 fixed-endpoint 论证。

## 2. whole-patch upper construction

mesh `sqrt(delta)/64` 的 center grid 对任意 `t>=delta/2` 给 `|a_i-v_i|<=sqrt(delta)/32`。精确 loss 满足 `sum l_i<delta/300`，故
\[
\log A_{c,r}F(x)\ge t-\delta/150>0.
\]
这有真正严格余量，覆盖的不是边界等号。取 `delta=log(1+eta)` 覆盖整个 saturation target。

dyadic `delta_j=2^{-j}/4` 覆盖 `X` 的全部正深度，而非仅一个 slab；各层 target 由 cutoff 单调性支配。大深度层数只为 `O(log N)`，可吸收进 `C^N`。小深度时
\[
q_j\le e^2e^{-\eta/(2\delta_j)},
\]
并以 `v_j=eta/delta_j` 分离 `v^{N/2}e^{-v/4}` 的最大值，得尾和常数对 N 一致。故 eta=`N^{-2}` 时的 upper exponents与 lower bounds匹配。每个 grid 有限、dyadic family 可数且总费有限，符合 budget 定义。

## 3. 一般残差公式

公式 (5) 正确：center shift 使原 witness rectangle 在新平均盒内，体积扩大至多 `prod(1+u_i)<=1+delta`，严格 margin 保留成功。以 `z_i=log(1+u_i)` 计算得到
\[
|C_j(x)|=\frac{|R_j(x)|}{(1+r)^N}
\int_0^{\log(1+\delta_j(x))}\frac{e^ss^{N-1}}{(N-1)!}ds.
\]
它只使用 outward containment，未把这个 certified center opportunity 误认为全部 success opportunity。

在 Borel witness data 下 incidence set 为 Borel；即使 `pi_j` 没有 density，section probability `p_j(x)` 仍可测。取有限 packets、每 packet 独立抽取 `k_j` 个 centers 后，Fubini 给 expected residual
\[
\int\rho(x)\prod_j(1-p_j(x))^{k_j}dx.
\]
target integrability使该式合法，故有一个确定 realization 达到 (6)，预算正是实际 `sum k_jJ(G_j)`。这不是 atomic receiver pricing，也不要求 centers 的绝对连续分布。

## 4. 限定范围

本轮的 patch upper 不可推广成一般 `F` 的 uniform upper：后者仍需为每个输入共同构造 packets、proposals 和 overlap payment，在 `N^{o(1)}J(F)` budget 与物理 residual尺度下同时工作。反过来，patch 可按物理体积以 `|X|<=N^{-N}` residual 删除，故其 superexponential pointwise budget不反驳允许物理残差的一般分解，也不更新 `C_N`。

## 主代理与Luna代码验收补记

Luna只读审查指出精确Log与近似数混用可能在Exp中生成巨型整数幂；最终实现先数值化每个对数项，再作有误差界的log-sum-exp。静态推断未冒充实际求值路径分析。主代理已用单项内存限制和系统余量监控运行全部独立核验，结果和精确区间互相包含。最新用户要求是约预留1GiB内存并尽量使用安全可分配算力；本批次采用12核且测得余量充足。
