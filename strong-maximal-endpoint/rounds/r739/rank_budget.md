# 保留逐秩 KL 亏损：精确中间式与逆检验障碍

本稿只处理 r710 负责混合熵证明中被舍去的量。得到一个严格的平均 KL／反馈熵亏损预算，并说明将其投影到 NP 事件时损失了什么。没有证明一般近中央逆检验比较，也没有声称存在满足真实盒几何与全部近中央假设的反例。结论不需要新增小维目录或低率模型。

## 1. 目标与原谱的精确秩参数化

令 \(P=P_\pi\) 是有限真实矩形目录上的负责 LF 混合，\(Q\) 是完整并上的均匀律，\(w=dP/dQ\)，\(M=D(P\Vert Q)\)。所有秩用 \(Q\)-质量归一化为 \((0,1)\)。选择常数 \(C_0\) 足够大，使 r710 对所有反馈策略的统一熵界为
\[
ES_\tau-h(S_\tau)\le B_n:=\mu n-\tfrac12\log n+C_0.
\]
该常数也吸收原负责混合熵界的常数。记非负、相差绝对常数无碍的亏损
\[
A=B_n+1-M\ge0,
\qquad \mathcal T=V-M.
\]
题设不含常数的 \(D=\mu n-\frac12\log n-M\) 可能为负；以下使用 \(A=D+C_0+1\)，所以出现 \(\log(1+A)\) 时没有符号问题。所需联合补偿等价于 \(V\le\mu n+o(\log n)\)。例如更强的
\[
\mathcal T\le\tfrac12\log n+O(\log(1+A))+o(\log n)
\]
会闭合，因为 \(O(\log(1+A))-A\) 有统一上界；本稿没有得到这项估计。

负责权保证正权行的并覆盖 \(Q\)-几乎所有点：对正权标签全集应用 Hall 条件，左端为 1。由于成员原子有限，\(w\) 在 \(Q\)-几乎处处有正下界和有限上界。因此可以直接使用 \(\log w\)，无须正则化后交换极限。

设 \(w^*(r)\) 是递减重排，并定义概率密度
\[
k(t)=e^{-t}w^*(e^{-t}),\quad t\ge0;
\qquad H(u)=\int_u^\infty k(t)\,dt.
\]
以 \(T\) 表示密度为 \(k\) 的**原混合秩变量**。这是原空间密度的精确参数化，不是把 \(T\) 宣称为反馈和 \(S_\tau\)。直接换元给
\[
\boxed{M=ET-h(T).} \tag{1}
\]
NP 最优预算为 \(Qf=e^{-V}\)、精确命中为 \(Pf=q\)，故
\[
\boxed{H(V)=q.} \tag{2}
\]
在密度平台上均匀打破秩并列即可；秩变量本身连续，平台上的源随机化没有造成额外原子。这里必须使用最优 top-rank 预算，不能把任意可行源的体积代入 (2)。

## 2. 逐秩证明中保留的全部 KL 项

r710 对同一个 \(g=\log w+\text{常数}\) 的负责／重排比较给
\[
M\le\int_0^\infty e^{-u}\sup_\tau
E\log w^*(e^{-(u+S_\tau)})\,du. \tag{3}
\]
每个 \(S_\tau\ge0\) 是原证明允许的可预测动作和，具有有限均值、有限微分熵及密度。最后一步的独立指数卷积保证这些性质。每个秩 \(u\) 允许选择不同策略，这仍只是上比较，不是原盒核的等式表示。

定义
\[
J_\tau=B_n-[ES_\tau-h(S_\tau)]\ge0,
\]
以及
\[
K_\tau(u)=D\!\left(\mathcal L(S_\tau)\;\middle\Vert\;
\frac{k(u+\cdot)}{H(u)}\right)\ge0.
\]
右侧参考律是原秩变量条件在 \(T\ge u\) 后减去 \(u\) 的律。保留归一化项，有**逐策略精确恒等式**
\[
\boxed{
E\log w^*(e^{-(u+S_\tau)})
=u+B_n-J_\tau-K_\tau(u)+\log H(u).
} \tag{4}
\]
这是因为 \(\log w^*(e^{-t})=t+\log k(t)\)，而
\(K_\tau(u)=-h(S_\tau)-E\log k(u+S_\tau)+\log H(u)\)。

代入 (3)，利用 \(\int_0^\infty ue^{-u}du=1\)，得到本稿最主要的严格中间式：
\[
\boxed{
\int_0^\infty e^{-u}
\left[-\log H(u)+\inf_\tau\{J_\tau+K_\tau(u)\}\right]du
\le A.
} \tag{5}
\]
三项费用的方向都不可反转：缺额 \(A\) 给的是可以选到的策略家族的**上预算**。它不直接给 \(\mathcal T\) 的上界，也没有推出原谱正态逼近。

(5) 不需要预先选一个全秩最优策略。若需要以下联合实验，可选择可测的近最优策略家族：有限成员原子使 \(\log w^*(e^{-t})\) 为有限阶梯函数；\(S_\tau\) 有密度，故各策略的 (3) 被积期望是 \(u\) 的连续函数。这些函数组成紧集上一致收敛拓扑下的可分函数族，取可数稠密子族，再取第一个达到给定误差的策略，即得到可测近选择。于是对任意 \(\epsilon>0\)，可使
\[
E_U[J_{\tau(U)}+K_{\tau(U)}(U)-\log H(U)]\le A+\epsilon,
\quad U\sim\mathrm{Exp}(1). \tag{6}
\]
这里的 \(U\) 是均匀接收秩的对数，不是原空间后验标签。

## 3. 可以抽出的两个有意义控制

### 3.1 近比较反馈和的平均动作损失

设 \(a_\tau=\mu n-ES_\tau\ge0\)。r710 的分位耦合和分箱熵估计给
\[
J_\tau\ge a_\tau-\log(1+a_\tau+C\sqrt{a_\tau})-C'.
\]
从 \(x-\log(1+x+C\sqrt x)\ge x/2-C''\) 得
\[
\boxed{E_Ua_{\tau(U)}\le C(A+1+\epsilon).} \tag{7}
\]
再调用原证明实际构造的耦合
\[
W_1(\mathcal L(S_\tau),\mathcal L(S_n^*))
\le a_\tau+C\sqrt{a_\tau},
\]
其中 \(S_n^*\) 是固定最优动作的独立增量和。令
\[
R=U+S_{\tau(U)},\qquad R^*=U+S_n^*,
\]
\(R^*\) 中两项独立；对 \(U\) 条件耦合后得到
\[
\boxed{W_1(\mathcal L(R),\mathcal L(R^*))\le C(A+1+\epsilon).} \tag{8}
\]
(7)—(8) 控制的是近比较反馈律 \(R\)，没有把原秩 \(T\) 偷换成它。

### 3.2 原秩与反馈秩之间的正确 KL 方向

对 (6) 的联合分布 \((U,S_{\tau(U)})\)，参考**次概率密度**是
\[
e^{-u}k(u+s)\,1_{u,s\ge0}.
\]
其 \(t=u+s\) 边际为
\[
b(t)=(1-e^{-t})k(t),\qquad
Z:=\int_0^\infty b(t)dt<1.
\]
补上一个质量 \(1-Z\) 的墓地状态，KL 数据处理给
\[
\boxed{
D\!\left(\mathcal L(R)\middle\Vert b/Z\right)-\log Z
\le E_U[K_{\tau(U)}(U)-\log H(U)]
\le A+\epsilon.
} \tag{9}
\]
这是普通概率 KL 加上正确的次概率归一化费用；不能把 \(b\) 误写为概率密度，也不能遗漏 \(-\log Z\)。

令 \(a=P(R\ge V)\)、\(b_V=\int_V^\infty b(t)dt\)。由 (2)，
\[
(1-e^{-V})q\le b_V\le q.
\]
把墓地状态与阈值下侧合并，记二元相对熵为 \(d(\cdot\Vert\cdot)\)，则
\[
d(a\Vert b_V)\le A+\epsilon,
\]
进一步有
\[
\boxed{
d(a\Vert q)\le A+\epsilon+
\log\left(1+\frac{q e^{-V}}{1-q}\right).
} \tag{10}
\]
最后一步由
\[
d(a\Vert q)-d(a\Vert b_V)
=a\log(b_V/q)+(1-a)\log[(1-b_V)/(1-q)]
\]
直接估计。近中央时最后一项指数小，但 \(A\) 本身未被近中央得分假设控制为小量。

## 4. 为什么这些已证控制尚不给所需逆检验式

**第一处是 KL 方向及其尺度。** (9) 是反馈律到原秩律的 KL，而非两律的位移控制。对任何概率律 \(R\) 和任意远离其主要质量的律 \(B_L\)，
\[
D\bigl(R\Vert(1-q)R+qB_L\bigr)\le-\log(1-q),
\]
不随 \(B_L\) 移动距离增长。这个简单事实不是实际矩形反例，也没有宣称满足原混合的均值／熵／Hall 条件；它只说明不能从这一方向的 KL 上预算单独恢复远尾位置。

投影式 (10) 在 \(q=1/2\) 尤其明显：所有 \(a\in[0,1]\) 都满足 \(d(a\Vert1/2)\le\log2\)。所以若可用缺额预算已超过 \(\log2\)，该二元不等式对 \(P(R\ge V)\) 完全没有约束。即使预算是更小的固定常数，也仅迫使尾概率接近到常数精度；宽度 \(\sqrt n\) 的中央分布中，这仍对应 \(\sqrt n\) 量级的分位移动，而不是 \(O(\log n)\)。

**第二处是耦合对象。** (8) 的小平均位移发生在 \(R\) 与 \(R^*\) 之间，不发生在 \(T\) 与 \(R^*\) 之间。\(R^*\) 的密度峰高为 \(O(n^{-1/2})\)，故仅由 (8) 与标准平滑换算得到的 Kolmogorov 误差至多是
\[
O\bigl((A+1)^{1/2}n^{-1/4}\bigr),
\]
尚且不是所需的 \(o(\log n/\sqrt n)\) 精度；更不能跨过 (9) 把它赋给 \(T\)。该估计不证明更强耦合不可能，而是说明当前这组估计本身不够。

**第三处是可预测反馈放松。** 原证明允许对每个接收秩分别取策略上确界。不同秩的 \(\tau(u)\) 不必由同一个真实目录、同一个 NP 源或同一套负责运输实现。(5) 正确保留了这个上界的亏损，却没有恢复被放松的全局兼容性。若直接把 \(S_{\tau(u)}\) 当成原谱 \(T-u\) 的生成机制，就恰好丢掉了 \(K_\tau(u)\) 所测量的差异。

这些论点排除的是“把 (5) 的 KL 预算配上 Pinsker／二元事件比较便自动推出对数检验价”的推理；它们**不证明** (5) 连同所有真实几何信息永远无法得到更强结果，也不构成一般近中央反例。

## 5. 尚可尝试但本稿没有补上的几何接口

若能在 r718 选定的近中央目录中证明 \(\mathcal T\) 对预算 \(A\) 的更强约束，仍可闭合。现在具体需要的是以下一类额外输入：

- 把 (9) 的反向 KL 转化为对原谱 \(T\) 的尾位置控制，同时利用其均值熵恒等式 (1) 和来自同一盒核的限制；单独替换成 Pinsker 不够。
- 证明同一个实际负责运输使各秩的比较策略具有额外兼容性，从而严格改进 (3) 中逐秩独立的上确界。
- 直接对原密度重排证明临界尾补偿。r737 的正态逼近条件只是足够强的一种实现，不是必须证明的唯一形式。

原 r710 的统一熵界仍完整有效；若最终得到 \(\mathcal T\le\tfrac12\log n+o(\log n)\)，半对数项已被该熵上界抵消，无需额外容量增益。本稿没有再次引入这个已纠正的假缺口。

**核验记录。** 已重读 r710 的逐秩／分位耦合／熵证明，以及 r731 的检验—容量接口。本机既有 Python/SymPy 实际对 (4) 的归一化代数与 (10) 的二元 KL 差恒等式化简，两个残差均精确为零。另用原生 Astra 独立审查了 (1)—(10) 的方向、次概率墓地补全、NP 平台处理与适用前提；没有进行新小维数值试验。符号计算只核验代数，(3) 的几何比较与控制过程估计仍依赖 r710 已有证明。

## 6. 精确 Gamma 混合：对数尺度上的抽象谱障碍

以下构造只针对“熵上界＋合法递减秩＋矩控制＋小反向 KL 已足够推出所需谱位置”这组**抽象**推理。它不是矩形目录，没有认证 Hall／LF 条件，也不声称满足 (3) 的全部真实反馈比较。尤其此处的比较密度是明确给定的 Gamma 密度，不冒充真实最优动作和 \(S_n^*\)。

取固定 \(L>0\)，令
\[
m_n=\mu n-1,\qquad \theta_n=1+n^{-1/4},\qquad
\alpha_n=m_n/\theta_n,\qquad
\epsilon_n=L\frac{\log n}{\sqrt n}.
\]
以下取充分大 \(n\)，使 \(0<\epsilon_n<1/2\)、\(\alpha_n>8\)。用 \(\operatorname{Gam}(\alpha,\theta)\) 表示 shape 为 \(\alpha\)、scale 为 \(\theta\) 的概率密度。定义
\[
f_0=\operatorname{Gam}(\alpha_n,\theta_n),\quad
f_- =\operatorname{Gam}\!\left(\alpha_n,\theta_n(1-2/\sqrt{\alpha_n})\right),
\quad
f_+ =\operatorname{Gam}\!\left(\alpha_n,\theta_n(1+1/\sqrt{\alpha_n})\right),
\]
\[
g_n=\tfrac13 f_-+\tfrac23 f_+,
\qquad
\boxed{k_n=(1-\epsilon_n)f_0+\epsilon_n g_n.} \tag{11}
\]

**均值精确保留，秩条件严格成立。** 因
\[
\tfrac13(1-2/\sqrt\alpha)+\tfrac23(1+1/\sqrt\alpha)=1,
\]
每个 \(k_n\) 的均值恰为 \(m_n\)。三个 scale 都严格大于 1：最小者相对 1 的正增量是 \(n^{-1/4}-O(n^{-1/2})\)。对 shape \(\alpha>1\)、scale \(\theta>1\) 的 Gamma 密度，
\[
\frac{d}{dt}[e^t f(t)]
=e^tf(t)\left(\frac{\alpha-1}{t}+1-\frac1\theta\right)>0.
\]
故 \(e^t k_n(t)\) 严格递增，确实定义了合法的递减秩密度
\(w_n^*(r)=k_n(-\log r)/r\)，且 \(\int_0^1w_n^*=1\)。这是抽象可测空间上的合法密度；它不是有限盒混合的构造。具体地，因 scale 大于 1，所得 \(w_n^*(r)\) 在 \(r\downarrow0\) 时无界，且在 \(r\uparrow1\) 时趋于零，不能将它直接纳入第 1 节的有限正密度类。混合方差为 \(\mu n(1+o(1))\)，三阶绝对中心矩为 \(O(n^{3/2})\)。

**反向 KL 甚至趋于零。** 对 \(f_d=\operatorname{Gam}(\alpha,\theta(1+d/\sqrt\alpha))\)、\(d=-2,1\)，直接积分得
\[
\int\frac{f_d^2}{f_0}
=(1-d^2/\alpha)^{-\alpha}.
\]
因此
\[
\chi^2(g_n\Vert f_0)
\le\tfrac13(1-4/\alpha_n)^{-\alpha_n}
+\tfrac23(1-1/\alpha_n)^{-\alpha_n}-1=O(1).
\]
记该上界为 \(C_n\)，则
\[
D(k_n\Vert f_0)\le\epsilon_n^2C_n,
\qquad
D(f_0\Vert k_n)\le\frac{\epsilon_n^2C_n}{2(1-\epsilon_n)^2}
=O(\log^2n/n)=o(1). \tag{12}
\]
第二式对 \(-\log(1+\epsilon_n(g_n/f_0-1))\) 作带下界 \(1-\epsilon_n\) 的 Taylor 估计；线性项因积分为零而消失。

**熵变化一般是一阶 \(O(\epsilon_n)\)，不是 \(O(\epsilon_n^2)\)。** 这个正混合保留均值，却没有保留 \(E\log T\)。精确恒等式为
\[
h(k_n)-h(f_0)=\epsilon_n A_n-D(k_n\Vert f_0),
\]
\[
A_n=-(\alpha_n-1)
\left[\tfrac13\log(1-2/\sqrt{\alpha_n})
+\tfrac23\log(1+1/\sqrt{\alpha_n})\right]
=1+O(\alpha_n^{-1/2}).
\]
所以 \(h(k_n)-h(f_0)=\epsilon_n+O(\epsilon_n/\sqrt{\alpha_n}+\epsilon_n^2)=o(1)\)，这已经足够。Gamma 熵公式与 Stirling 给
\[
h(k_n)=\tfrac12\log n+\tfrac12\log(2\pi e\mu)+o(1).
\]
于是对应的原谱熵量满足
\[
\boxed{M_n=E_{k_n}T-h(k_n)
=\mu n-\tfrac12\log n-1-\tfrac12\log(2\pi e\mu)+o(1).} \tag{13}
\]
它满足 r710 形式的统一熵上界，并且 \(A_n^{\rm budget}=B_n+1-M_n=O(1)\)。此预算记号与上段的熵线性系数 \(A_n\) 是不同的量。

**中位数却有正的 \(\log n\) 位移。** 令 \(V_n\) 为 \(k_n\) 的中位数，即 \(\int_{V_n}^\infty k_n=1/2\)。中心 \(m_n\) 处，Gamma 标准化给
\[
F_{g_n}(m_n)\longrightarrow
\tfrac13\Phi(2)+\tfrac23\Phi(-1)<\tfrac12.
\]
基 Gamma 的中心 CDF 为 \(1/2+O(n^{-1/2})\)，中心附近密度为
\(\phi(0)/\sqrt{\mu n}\,(1+o(1))\)。在 \(O(\log n)\) 窗口内一致成立，因而
\[
\boxed{V_n=m_n+c_0L\log n+o(\log n),} \tag{14}
\]
\[
c_0=\left[\tfrac12-\tfrac13\Phi(2)-\tfrac23\Phi(-1)\right]
\sqrt{2\pi\mu}>0.
\]
正性不依赖小数核验：方括号等于
\[
\tfrac13\left(\int_0^1\phi(x)dx-\int_1^2\phi(x)dx\right)>0,
\]
因为正半轴上的正态密度严格递减。CDF 展开的 \(O(n^{-1/2})\) 基误差小于 \(\epsilon_n\)，而 \(\epsilon_n\sqrt n=L\log n\)，正好产生 (14) 的对数位移。

合并 (13)—(14)，得到
\[
V_n-\mu n=c_0L\log n+o(\log n),
\qquad
V_n-M_n=\left(\tfrac12+c_0L\right)\log n+O(1)+o(\log n).
\]
因此即使保均值、\(O(n)\) 方差、三阶矩自然尺度、合法递减秩、锐到常数的熵上界，以及对一个同均值、同渐近方差、同锐熵尺度比较律的 \(o(1)\) 反向 KL 都具备，仍不能仅据这些标量性质推出目标半对数逆检验式。

这里没有宣称 (12) 的 \(f_0\) 就是原 r710 允许的 \(S_\tau\)，也没有验证构造满足 (5) 的真实反馈全族约束。若要加强到真实 \(S_n^*\) 的比较，需要另证明相应 KL 收敛；若要加强到真实盒反例，则更须构造同一目录、负责运输和共同 LF 源。当前构造只精确定位了抽象谱推理的缺口，不承担这两项更强结论。

本节额外用本机 SymPy 实际核验了均值系数恒等式，以及熵线性系数展开
\(A_n=1+2/(3\sqrt{\alpha_n})+1/(2\alpha_n)+O(\alpha_n^{-3/2})\)。没有为此进行小维目录或昂贵数值搜索。

第 6 节又经原生 Astra 独立审查：均值、秩单调、双向 KL 上界、熵一阶项和对数中位数系数均通过；审查特别指出上述无界端点，已明确记入适用范围。父任务另实际以本机 Wolfram 15、12 个内核完成 12 组高维参数的计算，输出为 [gamma_obstruction_wolfram.json](gamma_obstruction_wolfram.json)，状态为 `complete`；这些数值只核验解析谱机制，不认证真实盒可实现性。

## 附录 A：逐维化为有界正的有限阶梯秩密度

第 6 节的端点无界不是上述抽象障碍所必需的。下面证明：对每个固定 \(n\)，可以先截断，再取有限递减阶梯近似，使中位数、熵、矩和指定方向的 KL 都任意接近；最后沿 \(n\) 作对角选择。这里不要求显式截断尺度，也不声称截断后的同一函数继续**精确**保均值。

### A.1 截断并归一化

本节固定 \(n\)，暂省下标。令第 6 节的秩密度为
\[
w(r)=k(-\log r)/r,\qquad0<r<1.
\]
对 \(0<\ell\le1\le H<\infty\)，定义
\[
c_{\ell,H}(r)=\min\{H,\max\{\ell,w(r)\}\},\quad
Z_{\ell,H}=\int_0^1c_{\ell,H}(r)dr,\quad
w_{\ell,H}=c_{\ell,H}/Z_{\ell,H}.
\]
它仍递减、正、有界且积分为 1。其秩密度为
\(k_{\ell,H}(t)=e^{-t}w_{\ell,H}(e^{-t})\)。

令 \(\ell\downarrow0\)、\(H\uparrow\infty\)。由于
\(c_{\ell,H}\to w\) 且 \(c_{\ell,H}\le w+1\)，支配收敛给 \(Z_{\ell,H}\to1\) 和 \(L^1(dr)\) 收敛。对每个固定整数 \(j\le3\)，还有
\[
\int_0^1w_{\ell,H}(r)(-\log r)^jdr
\longrightarrow E_kT^j. \tag{A1}
\]
支配函数是 \((w+1)(-\log r)^j\)；第一部分是 Gamma 混合的有限矩，第二部分积分恰为 \(j!\)。

同样，所有至三阶的有符号及绝对对数矩收敛：
\[
\int_0^1w_{\ell,H}(\log w_{\ell,H})^jdr
\longrightarrow\int_0^1w(\log w)^jdr. \tag{A2}
\]
理由不是笼统地由 \(L^1\) 收敛推出熵收敛：当值小于 1 时，\(x|\log x|^j\) 在 \([0,1]\) 有界；当值大于 1 时，截断后的 \(x(\log x)^j\) 被原 \(w(\log_+w)^j\) 控制。后者可积，因为在 \(t\) 坐标中
\[
|\log w(e^{-t})|\le C_n(1+t+|\log t|), \tag{A3}
\]
且 \(k\) 是固定参数的有限 Gamma 混合。归一化只增加 \(\log Z_{\ell,H}\to0\) 的常数项。由 (A1)—(A2)，
\[
h(k_{\ell,H})=E_{k_{\ell,H}}T-
\int_0^1w_{\ell,H}\log w_{\ell,H}\,dr
\longrightarrow h(k).
\]
原信息密度 \(L=\log w(e^{-T})\) 的均值、方差及三阶绝对中心矩也连续收敛。

### A.2 指定反向 KL 的尾部确实可控

固定比较密度仍是第 6 节的 \(f_0\)，不随 \((\ell,H)\) 改变。关键精确式为
\[
\log k_{\ell,H}(t)
=-t+\operatorname{clip}\bigl(\log w(e^{-t}),\log\ell,\log H\bigr)
-\log Z_{\ell,H}.
\]
因截断区间包含 0，
\[
\left|\operatorname{clip}(y,\log\ell,\log H)\right|\le|y|.
\]
(A3) 在 \(f_0\) 下可积：Gamma 律有有限一阶矩与 \(E|\log T|<\infty\)。因此可以在 \(f_0\) 下使用支配收敛，得到
\[
\boxed{D(f_0\Vert k_{\ell,H})\longrightarrow D(f_0\Vert k).} \tag{A4}
\]
这不是由总变差或熵收敛自动推出的。上端截断会把很远尾部改为速率 1 的指数尾，而 \(f_0\) 的尾部较慢；但对每个固定 \(n\)，其对数费用被可积的 \(C_n(1+T+|\log T|)\) 控制，所以尾部不会留下固定 KL 缺口。此论证没有给出关于 \(n\) 统一的截断效率，也不需要它。

### A.3 有限递减阶梯化

固定一组 \((\ell,H)\) 后，\(w_{\ell,H}\) 连续延伸至 \([0,1]\)，并在两个端点附近分别为常数。于是可以在足够细的有限均匀分割上以单侧取值构造递减正阶梯函数，再归一化为 \(w^{\rm step}\)，使其对 \(w_{\ell,H}\) 一致收敛。

由于固定的正下界存在，进一步有
\[
\sup_{0<r<1}|\log w^{\rm step}(r)-\log w_{\ell,H}(r)|\longrightarrow0.
\]
令 \(k^{\rm step}(t)=e^{-t}w^{\rm step}(e^{-t})\)。于是
\[
|D(f_0\Vert k^{\rm step})-D(f_0\Vert k_{\ell,H})|
\le\sup_r|\log w^{\rm step}(r)-\log w_{\ell,H}(r)|
\longrightarrow0. \tag{A5}
\]
同时，\(t^j\) 矩的一致误差由 \(j!\|w^{\rm step}-w_{\ell,H}\|_\infty\) 控制；对数矩在固定有界正值区间上连续。因此 (A1)—(A2) 及熵的收敛在阶梯化后仍成立。若需要有理有限数据，还可在保持递减性和正下界的前提下对各阶梯高度作足够细的有理舍入，再归一化；这不改变上述论证。

### A.4 中位数稳定及逐维同时选择

\(L^1(dr)\) 收敛等价于相应 \(t\)-密度的 \(L^1(dt)\) 收敛，故尾分布函数一致收敛。原 \(k\) 在中位数 \(V\) 邻域连续且严格为正；因此对任意 \(\eta>0\)，其在 \(V-\eta\) 与 \(V+\eta\) 的尾概率分别严格大于、小于 \(1/2\)。近似尾分布保留这两个严格不等式，就使其唯一中位数落在 \((V-\eta,V+\eta)\)。这直接证明中位数连续性，不要求阶梯密度在中位数处可微。

综上，对每个固定 \(n\)，先取充分小 \(\ell_n\)、充分大 \(H_n\)，再作充分细的有限阶梯化，可以同时使下列各量同原构造的误差小于 \(1/n\)：

- \(ET\)、\(h(T)\)、\(M=E\log w(e^{-T})\)；
- \(\operatorname{Var}T\)、原信息密度的方差及三阶绝对中心矩；
- 指定方向 \(D(f_{0,n}\Vert k_n^{\rm step})\)；
- NP 中位数 \(V\)。

这里方差等组合量按其连续性单独控制；不能只把一阶和二阶原始矩各误差取成 \(1/n\)，便未经检查声称方差误差也是 \(1/n\)。采用更小的原始矩容差即可满足有限多个同时要求。

于是有限阶梯版本仍满足
\[
M_n^{\rm step}=\mu n-\tfrac12\log n-1-
\tfrac12\log(2\pi e\mu)+o(1),
\]
\[
D(f_{0,n}\Vert k_n^{\rm step})=o(1),\qquad
V_n^{\rm step}=\mu n+c_0L\log n+o(\log n).
\]
其均值为 \(\mu n-1+O(1/n)\)，**不是**未经修正的精确 \(\mu n-1\)。所有阶梯数、上下界和精度都允许依赖 \(n\)，本构造没有高效或统一复杂度断言。

### A.5 信息密度的自然方差尺度

除了秩变量 \(T\) 的方差，还可核对 r737 使用的原信息密度方差。第 6 节的三个 Gamma 分量均有 \(\operatorname{Var}T=\Theta(n)\)，其混合满足 \(\operatorname{Var}_kT=\mu n(1+o(1))\)。在一个 scale 比为 \(r_d=1+d/\sqrt\alpha\) 的分量下，Gamma 的精确矩公式给
\[
\operatorname{Var}[\log f_0(T)]
=(\alpha-1)^2\psi_1(\alpha)+\alpha r_d^2-2(\alpha-1)r_d=O(1),
\]
且分量之间 \(E\log f_0(T)\) 相对基分量的差为
\[
(\alpha-1)\log r_d-\alpha(r_d-1)=O(1).
\]
故混合下 \(\operatorname{Var}(\log f_0(T))=O(1)\)。另外，写 \(X=(T-m_n)/(\theta_n\sqrt{\alpha_n})\)，三个密度比的显式指数形式给
\[
|\log(k/f_0)|\le C(1+|X|),\qquad E_kX^2=O(1),
\]
从而 \(\operatorname{Var}_k\log k(T)=O(1)\)。原信息密度为
\(L=T+\log k(T)\)，所以
\[
\boxed{\operatorname{Var}_kL=\mu n(1+o(1)).} \tag{A6}
\]
上述阶梯逼近保留该方差至 \(1/n\)，故端点修复没有破坏自然的 \(\sqrt n\) 标准差尺度。

本机 SymPy 实际核验了上述 Gamma 分量方差与均值差展开，分别为
\(d^2+1/2+2d/\sqrt\alpha+O(1/\alpha)\) 与
\(-d^2/2+O(1/\sqrt\alpha)\)。这只核验矩代数，逼近及尾部控制由 (A1)—(A5) 的支配收敛和一致对数控制证明。

**适用范围没有扩大成盒反例。** 本附录确实消除了“密度无界／非有限阶梯”这一抽象层面的逃口；但它仍没有构造满足真实负责 Hall 约束及共同 LF 最优性的矩形目录，也没有验证整个反馈比较族。有限阶梯秩表示本身不补上这些几何条件。

### A.6 重要交叉核对：真实几何的矩约束已经排除此 Gamma 候选

本轮的[真实几何矩约束证明](lp_spectrum.md)比本稿保留的均值／熵／KL 标量条件更强，且**严格排除**第 6 节的 Gamma 候选。此处据其已核验结果记录结算，并交叉引用 [r741 的有理区间证书](certify_lp.json) 和 [Wolfram 展开核验](lp_expansion_wolfram.json)。这项排除不能被有界阶梯化修复。

对于任意实际负责混合，取 \(p=(1+t)/t\)、\(t>0\)。负责平均与强极大函数的 \(L^p\) 界，经对偶给
\[
\|w\|_{L^{1+t}(Q)}\le c_p^n,
\qquad
E_Pe^{tL}=\int_0^1w(r)^{1+t}dr\le e^{nK(t)},
\]
其中 \(K(t)=(1+t)\log c_p\)，\(c_p>1\) 是一维锐范数的相应根。父任务核验的局部展开是
\[
K(t)=\mu t+\tfrac\mu2t^2-\frac{\mu(\mu+3)}{24}t^3+O(t^4).
\]
第 6 节基 Gamma 谱的同一矩率则为
\[
F_1(t)=\mu(1+t)\log(1+t),
\qquad
F_1(t)-K(t)=\frac{\mu(\mu-1)}{24}t^3+O(t^4)>0
\]
对充分小的固定 \(t>0\) 成立。\(\theta_n\to1\) 不改变领先率；又 \(k_n\ge(1-\epsilon_n)f_{0,n}\)，所以正混合保留此矩下界。特别地，上述有理证书在 \(t=1/10\) 给出严格率差
\[
F_1(1/10)-K(1/10)
>0.00001411891677682855533678003335204521785516464876468147205.
\]
这是固定参数上的严格正率差，因而对充分大 \(n\) 造成 \(O(n)\) 的对数矩违反，不是端点微小误差。

附录的逼近也能同时保留这项违反。固定 \(n\) 和 \(t=1/10\)，原 Gamma 混合满足 \(\int_0^1w^{1+t}<\infty\)（充分大 \(n\) 时所有 scale 均趋于 1，远尾指数仍可积）。截断后
\[
c_{\ell,H}^{1+t}\le w^{1+t}+1,
\]
故归一化后的 \((1+t)\)-矩由支配收敛趋于原值；固定截断后的一致阶梯逼近也保留该矩。将这一个额外连续量纳入 A.4 的有限个要求，可以选择
\[
\left|\log\int(w_n^{\rm step})^{1+t}
-\log\int w_n^{1+t}\right|<1/n.
\]
因此存在本附录所述同时保持谱障碍、又保持 \(O(n)\) 几何矩违反的有界正有限阶梯版本。结论是：端点修复不会把此候选变成合法负责盒谱；真正排除它的是已保留真实几何的矩约束。这个排除是本轮的新约束证据，而不是一般近中央结论已闭合。

## 附录 B：完整逐秩预算给出的密度层数下界与三点联合排除

本附录得到一个普适必要条件：若实际负责混合的密度只有 \(K\) 个不同正值，则其熵缺额预算 \(A=B_n+1-M\) 必须满足
\[
\boxed{K\ge c\,\frac{\sqrt n}{A+1}.} \tag{B1}
\]
这里 \(c>0\) 是绝对常数。计数的是**混合密度 \(w=dP/dQ\) 的不同密度层数**，不是输入源 \(f\) 的层数、矩形行数或成员原子总数。

这严格排除了 [几何矩约束稿](lp_spectrum.md) 附录中的三点候选：它虽满足全部正参数矩界，但只有三个密度层且 \(A=O(1)\)，与 (B1) 不相容。因此该特定候选不仅尚未认证负责实现，而且不可能满足已经证明的完整逐秩比较。本结论仍不是一般近中央结算。

### B.1 有限阶梯秩密度的指数距离矩

将递减秩密度 \(w^*\) 的相邻同值层合并，得到 \(K\) 段。在 \(t=-\log r\) 坐标，每一段是区间
\([s_j,s_j+\ell_j)\)，其中 \(\ell_j>0\)，最后的无界区间允许 \(\ell_j=\infty\)。所有起点构成 \(K\) 点集
\[
\mathcal S=\{s_1,\ldots,s_K\},
\]
且其中包含 0。对应秩密度在第 \(j\) 段上形如
\[
k(t)=C_j e^{-(t-s_j)}.
\]
条件于该段，残差 \(Y=T-s_j\) 的律就是截断在 \([0,\ell_j)\) 的单位指数律。因此对任意 \(0<\lambda<1\)，
\[
E[e^{\lambda Y}\mid T\text{ 在第 }j\text{ 段}]
=\frac{1-e^{-(1-\lambda)\ell_j}}
{(1-\lambda)(1-e^{-\ell_j})}
\le\frac1{1-\lambda}. \tag{B2}
\]
无界区间按极限解释，恰取等。由于
\(d(T,\mathcal S)\le T-s_j\)，对各段平均得到
\[
\boxed{E_k e^{\lambda d(T,\mathcal S)}\le(1-\lambda)^{-1}.} \tag{B3}
\]
这个估计不依赖各段的长度、位置和质量，也不需要三点构造的特殊参数不等式。

### B.2 完整 KL 预算控制比较律到这些起点的距离

由第 2—3 节，对任意 \(\epsilon>0\)，存在近最优比较律
\[
R=U+S_{\tau(U)},\qquad
\beta(t)=(1-e^{-t})k(t),
\]
满足
\[
D_{\rm sub}(\mathcal L(R)\Vert\beta)\le A+\epsilon,\qquad
W_1(\mathcal L(R),\mathcal L(R^*))\le C(A+1+\epsilon), \tag{B4}
\]
其中 \(R^*=U+S_n^*\)，其 \(U\sim\mathrm{Exp}(1)\) 与最优独立动作和独立，\(C\) 为绝对常数。

因 \(\beta\le k\)，(B3) 同时控制 \(\int e^{\lambda d}\beta\)。对次概率 KL 直接用 Gibbs 变分，得
\[
\lambda E d(R,\mathcal S)
\le A+\epsilon-\log(1-\lambda).
\]
特别取 \(\lambda=1/2\)，
\[
\boxed{E d(R,\mathcal S)\le2(A+\epsilon+\log2).} \tag{B5}
\]
这里没有遗漏次概率归一化，也无须证明 \(h(R)\) 有限。距离函数至多线性增长，\(ER\le1+\mu n\)，指数积分由 (B3) 有限；可先截断测试函数再取极限。等价地，先用 \(\beta\le k\) 得 \(D(R\Vert k)\le A+\epsilon\)，再应用概率版变分式也成立。

### B.3 最优独立比较律的密度峰高强迫足够多起点

最优独立动作具有精确表示
\[
S_n^*\overset d=\Gamma(n,1)+\mu\operatorname{Bin}(n,p_*),
\qquad p_*=(\mu-1)/\mu,
\]
两项独立。因此
\[
R^*\overset d=\Gamma(n+1,1)+\mu\operatorname{Bin}(n,p_*).
\]
Stirling 下界给 \(\Gamma(n+1,1)\) 的密度峰高
\[
\frac{n^ne^{-n}}{n!}\le\frac1{\sqrt{2\pi n}}\le\frac1{\sqrt n}.
\]
平移混合不增加峰高，所以 \(\|f_{R^*}\|_\infty\le n^{-1/2}\)。任意 \(K\) 个中心的半径 \(r\) 邻域总长度至多 \(2Kr\)，故
\[
P\{d(R^*,\mathcal S)\le r\}\le2Kr/\sqrt n.
\]
取 \(r=\sqrt n/(4K)\)，得到
\[
\boxed{E d(R^*,\mathcal S)\ge\sqrt n/(8K).} \tag{B6}
\]
距离函数为 1-Lipschitz，故由 (B4)—(B6)，
\[
\frac{\sqrt n}{8K}
\le2(A+\epsilon+\log2)+C(A+1+\epsilon).
\]
令 \(\epsilon\downarrow0\)，得到精确的常数形式
\[
\boxed{
K\ge\frac{\sqrt n}
{8[(C+2)A+C+2\log2]}.
} \tag{B7}
\]
这证明 (B1)。尤其 \(A=O(1)\) 强制 \(K=\Omega(\sqrt n)\)；更一般地 \(A=O(\log n)\) 强制 \(K=\Omega(\sqrt n/\log n)\)。这一结论不以 LF 最优性为前提，适用于满足负责条件的有限混合。

### B.4 三点候选的严格联合排除

三点构造的参数为
\[
m=\mu n,\quad\delta=\tfrac12\log n+4,\quad
h=m-\delta,\quad p=e^{-m},
\]
其信息密度层为 \(h+a,h-a,b\)，相应 P 质量为 \(1/2,1/2-p,p\)。其精确均值已经计算为
\[
M=h-p(h-a-b)
=\mu n-\tfrac12\log n-4+O(ne^{-\mu n}).
\]
因此
\[
A=B_n+1-M=C_0+5+O(ne^{-\mu n})=O(1),
\qquad K\le3.
\]
这与 (B7) 矛盾。故所有上述三点参数变体均不可能满足完整逐秩比较，尽管它们满足全部正参数 \(L^p\) 矩族。

这个排除展示了两个已有真实几何约束的互补性：正参数矩族排除了先前的 Gamma 候选；完整逐秩 KL／耦合预算又排除了通过全部矩界的三点候选。两者都没有穷尽一般谱。特别是附录 A 的有限阶梯逼近允许阶梯数随 \(n\) 任意增大，(B7) 不因“有限”二字就自动排除它；该 Gamma 逼近仍由 A.6 的矩约束排除。

本附录的 KL 变分、条件截断指数矩、密度峰高与 Wasserstein 步骤均为解析证明；没有以新增数值试验替代论证。原生 Astra 已独立审查并确认本附录的截断指数矩强化、常数方向和三点应用；这一证明避免了通过最大熵估计产生的不必要 \(\log n\) 费用。


### B.5 对实际负责混合的直接熵改进

由 \(C+2\log2\le C+2\)，(B7) 还给
\[
K\ge\frac{\sqrt n}{8(C+2)(A+1)},
\qquad
A\ge\frac{\sqrt n}{8(C+2)K}-1.
\]
代入 \(A=B_n+1-M\)，得到某些绝对常数 \(c,C_1>0\) 下的严格一般子类上界
\[
\boxed{D(P_\pi\Vert Q)
\le\mu n-\tfrac12\log n-c\frac{\sqrt n}{K}+C_1.} \tag{B8}
\]
适用条件是：有限实际矩形目录、可负责权 \(\pi\)，以及该混合密度恰有至多 \(K\) 个不同正值。不要求 LF、单相位或固定几何形状。\(K\) 较大时额外节省可能没有渐近作用；本式不提供一般 \(V\) 的上界，因而不宣称已闭合精确端点系数。

(B8) 是完整逐秩预算在有限混合密度层数上的新定量后果。特别地，\(K=O(1)\) 的实际负责混合不能达到只差 \(O(\log n)\) 的锐熵水平；这比单独使用全参数矩界更强，并精确解释三点抽象谱为何被联合排除。
