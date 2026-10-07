# 冻结 all-visited Hörmander 稿的独立解析审计

审计对象：`frozen_jump_allvisited_hormander_20261007.md` 最新正文，重点 §4、§5、§7。结论：在该稿及已审 fixed-obstacle / proper-face 接口的范围内，未发现阻塞性错误。低频估计已经明确限定在原障碍域外 \(\Omega^c\)，这一限定是 §7 回代必需的。没有运行或重启任何数值，本收据不认证原时变 cube 的 actual FIRST 转移。

## 1. 矩形传播与移位

保留完整 all-visited 核后，复时间 strip 的全空间范数 \(M=2\cdot3^n\) 与实轴空间差估计共同作用；不能先对各跳数取绝对值再求和。宽 \(2\theta\)、高 \(\theta\) 的矩形中，底边 harmonic measure 在中心及 \(0\le y\le\theta/2\) 由
\[
\cos\frac{\pi x}{2\theta}\,
\frac{\sinh[\pi(\theta-y)/(2\theta)]}{\sinh(\pi/2)}
\]
从下控制，其中心最小值为 \([2\cosh(\pi/4)]^{-1}>1/4\)。这个函数在侧边及顶边为零，底边不超过一，故比较方向正确。

对 \(L^1\) 对偶函数逐个用 two-constants，再取上确界，确实给出 \(M^{3/4}m^{1/4}\)。若实轴上界 \(m>M\)，直接用 \(M\le M^{3/4}m^{1/4}\)；不需要错误地坚持 \(m\le M\)。相邻底边时间的 \(e^\theta\) 因子可被正文的 \(32n\) 吸收。两个 logtime 尾积分是 \(4+8\)，因而正文 \(J_n=48\,3^{3n/4}n^{1/4}\) 相符。

中间高度的两个尾都一致趋零，足以使竖边在 Banach 范数中消失。向 \(-\operatorname{sign}(\omega)i/256\) 移位的符号正确：\(e^{-i\omega(v+iy)}\) 的模是 \(e^{\omega y}\)。于是平移差取得 \(e^{-|\omega|/256}\)，连续频率 Tonelli 给正文 Hörmander 常数。未发现这里把实轴衰减直接冒充复时间衰减的交换。

## 2. Mellin 与局部核的匹配

小时间 full mass \(\le(t/c)^n\)，大时间 \(\|P_t^{\rm full}\|_2\le10^{n/2}t^{-n/4}\)。故对 \(f\in L^1\cap L^2\)，正文 (12) 在 logtime 中是实际 \(L^2\) Bochner 积分。中间时间区间允许依赖固定 \(c>0\)；定义核的局部常数不被当作统一来源费用。

正文“有限 inclusion-exclusion”应按下列正则解释，不能将各项本来发散的 \(dt/t\) 积分逐项作普通积分。对非零空间频率，记
\[
\lambda_S=\sum_{i\notin S}\psi_c(\zeta_i^2)+|S|/c>0.
\]
先在 \(\Re z>0\) 使用
\[
\int_0^\infty t^{z-1}\prod_i(e^{-t\psi_i}-e^{-t/c})\,dt
=\Gamma(z)\sum_S(-1)^{|S|}\lambda_S^{-z}.
\]
完整乘积的小时间 \(O(t^n)\) 消去使左边延拓到 \(\Re z>-n\)，右边在 \(z=0\) 的极点由 \(\sum_S(-1)^{|S|}=0\) 消去。随后取 \(z=-i\omega\)，配合
\(i\omega\Gamma(-i\omega)=-\Gamma(1-i\omega)\)，即得 (13)。这给出了正文表述的合法解释，不新增假设。\(\omega=0\) full 为零，空间 \(\zeta=0\) 仅是欧氏 Fourier 的零测集合；有限链常数模不能由此处理。

也可用正文引用的 Abel 正则作相同匹配：对 \(0<r<1\)，full 正则核满足
\[
H_t^{{\rm full},r}=e^{-(1-r)nt/c}H_{rt}^{\rm full}\le H_t^{\rm full}
\]
（不等式为核的正项逐项比较），因此上述小/大时间界支配 \(L^2\) Mellin 极限。没有从这里推出全跳数的强 \(L^1\) 绝对可和性。

对有界空间集 \(E\)，正文局部质量界的小/大时间指数均为正可积指数。应用 §4 的局部 Banach 传播得到 logtime 的一致尾衰减，再移位得到
\(\int_{|\xi|>\Lambda}\|1_E\mathcal K_{2\pi\xi}^{\rm full}\|_1d\xi<\infty\)。因此 \(L^1_\xi\)-值 kernel 是空间局部可积的。对紧支撑测试输入，接收区与来源区的差集有界，局部 Fubini 使卷积与 multiplier 匹配；继而以 \(L^2\) 连续性延拓。用于坏块的离对角 kernel representation 已有来源，不能只凭乘子估计省略它。

## 3. CZ 与固定 maximal 的回代

原稿使用强 Banach 值 \(L^2_x(L^1_\xi)\) 界与 Hörmander 差界，未调用 weak-\(L^1\) Minkowski。三倍坏块外的几何满足 \(\|x-c_Q\|_\infty>2\|y-c_Q\|_\infty\)，因此取消项确实落入已证空间差域。正文 height 的两项优化给 \(4\,6^{n/2}B_\Lambda\)，坏块给 \(4H_\Lambda\)。full 弱尾与 miss 强尾按各一半阈值合并为 (20)，其常数 \(5W_b/16\) 正确。

§7 的低频界必须、且最新正文已经写为
\[
\|T_\Lambda\sigma\|_{L^2(\Omega^c)}^2
\le C_0\log^2(e\Lambda)\kappa W_b.
\]
在 \(\Omega^c\) 才使用 \(u=0\)、\(H_tu\ge0\)，得到 (21)。域内另按 \(\kappa|\Omega|\le W_b\) 支付体积，不能把低频界扩大到全空间。设 \(\kappa=\alpha/[2\log(e\Lambda_n)]\)，域外以 \(\alpha/4\) 分别分配低/高尾，费用分别为 \(8C_0\log(e\Lambda_n)W_b\)、\(5W_b/4\)；域内为 \(2\log(e\Lambda_n)W_b\)。这正是 (22)。

非负 \(L^1\) 密度的递增有界紧支撑截断，先对有理时间单调极限，随后利用固定有界生成元的幂级数取得 a.e. 紧时间连续代表，是合法的全时间延拓。该步骤没有把原子输入的 holding 或轴面输出免费转成密度。

## 4. 审计范围与当前 cube 账

结论是每一个固定 \((c,\ell)\) 生成元自己的 \(O(\log(n+2))\) 弱 maximal 定理，常数对其参数一致；没有证明对 \((c,\ell)\) 同时取上确界，也没有构造原时变 cube actual kernel 的正支配。各冻结算子的障碍域和来源清算不能合并为共同一份 \(W\)。因此原 cube 当前仍保留 \(R_\dagger\) 未付空间余项，不增加冻结 family 费用、不宣称一般目标闭合。

数值范围：本次只读审查解析正文和已审接口；未重跑任何旧或新守卫，未另认证保存的浮点组件、原 FIRST 样本或连续参数数值范围。
