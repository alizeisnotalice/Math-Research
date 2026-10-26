# A095 独立证明核查：全空间小幅平滑修正

审查对象固定为 C095 manifest
`e35b524a3424ab3ed3f688380850194f5e6ef966e943c7e06bdaf08cfdb24f6a`，
theorem SHA `d8ff73cbc86ae25ce827308185c86b8423f55390bb6e196f77eda09213e84c16`。
结论为 **PASS_PROOF_WITH_SCOPE**：该精确对象给出的正幅度修正存在性与全空间源
支配成立。下面给出独立的论证链及易错处说明；不凭三个小 fixture 代替全称证明。

## 1. 真实模型、输入和所需假设

固定 V={(1,±1,…,±1)}，32 个方向各权 1，未单位化，|v|²=6。
生成元 A=−Σ_v|v·D| 的二差分系数是 1/(2π)；在函数局部为零或常数时，
合并 s 与 −s 得到单差分系数 1/π。这与标准 Cauchy 密度 t/[π(t²+z²)] 一致。
没有移入 augmented frame、其他 α、其他 jump rate 或非 Lebesgue 背景。

已逐式读取并锁定实际输入：六维 theory/001 的全空间正部源与六维障碍证明，
paper_finite_source.md §A–C 的一般有限源核证明，以及旧 U10 bump 接口。
001 对 n=6 使用 Ct^-6 的真实热核界、H¹ 生成元域、频率尾与大 cutoff 回收质量；
其证明不依赖“只在三维已接受”的标签。所用准确基的 |det B|=160 是已有输入，
本轮没有重跑该基或旧扫描。

因此可用：G lower semicontinuous、S2 上为无穷、S2 外有限连续；K0=(G−c)+∈L¹；
AK0=ν0−δ0 于全空间，ν0 正概率、ν0({0})=0；以及 G 的 −5 次齐次性。
本轮的结论明确依赖这些已建立的同模型输入，不能推广到任意换过的 generator。

## 2. 核非负及支撑 buffer

支撑两盒位于 ±p、halfwidth=1/100，第一坐标绝对值至少 119/100，故与 Q 不交。
六个绝对坐标之间任意相邻差至少 9/50。任意三个 body 方向的线性组合，其坐标
至多有四种绝对值，因为八种 ± 三符号组合按整体反号配成四对。因此每个 S3
向量至少有两坐标绝对值相等。反三角不等式与二维 Cauchy–Schwarz 给

\[
\operatorname{dist}(z,S_3)\ge(9/50)/\sqrt2>1/10=\delta
\quad(z\in\operatorname{supp}b).
\]

对 x=z+s v、任意实际 body pair plane W，有

\[
\operatorname{dist}(x,W)\ge\operatorname{dist}(z,W+\mathbb Rv)
\ge\operatorname{dist}(z,S_3)\ge\delta.
\]

所以有限个完整圆柱 C_b=∪_v(supp b+Rv) 全部避开 S2，而非只支撑点避开。
每个圆柱是紧集与闭子空间之和，因而闭；Ab 在 C_b 外为零。所有 jump-in 远尾
仍在这个未截断的圆柱集合中。只要求支撑 off-S2 的弱化版本是错误的，见新 C2。

对 off-S2 的 cube 边界点 y，应使用实际连续性：从 cube 内且避开 essential
例外零集的点趋近 y，得到 G(y)≥g0。单独 lower semicontinuity 不能推出这一点，
作者没有犯该错误。对支撑点 x，x/||x||∞ 仍 off-S2，故齐次性给
G(x)≥g0/||x||∞⁵>g0/3=2c。因此 K0>c 于 supp b、0≤b≤1，

\[
b\le K_0/c=C_0K_0\quad\text{于全空间}.
\]

支撑紧且包含于开集 Ω={G>c}，所以 b 在闭集 E=Ω^c 的一个邻域恒零。
E 与 S2 不相交，因为 S2⊂Ω。

## 3. uniform G 尾：不能偷用接近奇异方向的角紧性

不同的三个 unoriented body 方向确实独立。若 w=au+bv 且 w,u,v 都为 ±1
向量，u/v 坐标中两种符号均出现，则 |a+b|=|a−b|=1，推出 ab=0，从而 w
等于 ±u 或 ±v，与第三方向不同矛盾。

固定 |x|≥R≥1、dist(x,S2)≥δ，δ≤1。设 d=δ/(2m√6)、ell=R/(m√6)。
任何映到 x 的标量参数都至少有三个绝对值>d；否则其余小参数和距一个 body
pair span 至多 m√6 d=δ/2。最大参数绝对值至少 ell，且 ell≥d。因此可由
N=m C(m−1,2)=14880 个指定三参数事件覆盖。

对每个事件独立选包含该三列的实际六列基。三条 restricted 密度的上界分别为
t/(π ell²)、t/(πd²)、t/(πd²)，另三条峰值为 (πt)^−1。每个非零整数 det
至少 1，未解的 26 个 Cauchy 变量只积分为 1。这是正 pushforward/coarea 的
未归一事件界，不是对 X=x 条件概率的计算，也不混合不同基中的独立变量。
故 a.e. 的上界为 A0/(R²δ⁴)，A0=16N32⁶6³/π⁶。

连续 p_t 将它延拓为所需逐点界。闭边界也可直接取 x_j=(1+1/j)x：两项距离
都严格增大，先用同一 R、δ 的开域界，再令 j→∞，避免 R=1 时参数单侧极限
的形式问题。该补充只说明已用的连续延拓，不增加假设。

结合既有 Ct^-6 界，任意 T>0 给

\[
G(x)\le\frac{\bar A T}{R^2\delta^4}+\frac{\bar C}{5T^5},
\quad\bar A=16N32^6 6^3/3^6,\quad\bar C=1/(160\,3^6).
\]

作者整数 T≥max(1,4 Cbar/(5c)) 保证第二项≤c/4；
R²≥4 Abar T/(cδ⁴) 保证第一项≤c/4。ceil 与 ceil_sqrt 只向安全方向取整。
从而 Ω∩C_b⊂{|x|<R,dist(x,S2)≥δ}。这里允许圆柱趋向奇异方向，所用的是
固定正 S2 距离的 uniform 尾，不能用不适用的角紧集论证代替。

## 4. Ω 中真实源密度及显式正下界

h=min(G,c) 全空间连续有界。S2 点附近，G=∞ 与 lower semicontinuity 使 h
局部恒 c；其余点用 G 的实际连续性。Ω 内 h 局部恒 c，所以测试 −Ah 或直接
交换远 jump 积分给

\[
\nu_0|_\Omega=m_0(x)dx,\quad
m_0(x)=\frac1\pi\sum_v\int_\mathbb R
\frac{(c-G(x+s v))_+}{s^2}\,ds.
\]

每个 Ω 内紧集具有统一 small-s 零段，其外由 c/s² 支配，故该密度局部有限连续。
这不意味着 m0 在所有轴/S2 上严格正；本证明只需 Ω∩C_b。

固定一条 v0，取 L=12R。若 pair plane W 含 v0，沿 v0 的距离保持≥δ。
若不含 v0，正整数 Gram determinant 及 det Gram(w1,w2)≤36 给
dist(v0,W)²≥1/36。对 |x|<R、s∈[L,2L]，有

\[
\operatorname{dist}(x+s v_0,W)\ge s/6-|x|\ge R\ge\delta,
\qquad |x+s v_0|\ge s-|x|\ge11R.
\]

于是整条同一参数段都满足 G≤c/2。仅一个 v0 项已给

\[
m_0(x)\ge c/(4\pi L)\ge c/(16L)=\bar m>0
\quad(x\in\Omega\cap C_b).
\]

没有未识别的紧集最小值或源采样。本段只在写出的解析 jump 积分中积分 s^-2，
本轮没有实际评价 m0、Ab 或其他 source density 积分。

## 5. Ab 是普通 L¹ 连续函数，没有新增界面奇异源

η 在端点 flat，两个 bump 的闭支撑不交。直接微分可用
||η′||∞≤4、||η″||∞≤148；Hessian 对角≤148 hbox^-2、非对角≤16 hbox^-2。
对称 Hessian 的 row-sum 给 H=2280000；不需把互不相交的两支撑再乘 2。
小 |s| 二差分由 s²|v|²H 控制，大 |s| 由 4||b||∞ 控制。按原系数积分给

\[
\|Ab\|_\infty\le(192H+128)/\pi<145920043=M.
\]

同样的 L¹ 二差分界（D²b 与 b 均 L¹）使积分绝对可交换；逐平移项的空间积分
为零，故 ∫Ab=0。支配收敛给连续性，且 Fourier 中 ψ b̂∈L¹ 给
(P_tb−b)/t→Ab 于全空间 L∞。这些是 source 作为普通密度的理由，而非“紧支撑
所以没有非局部外尾”。因此没有新增平面/支撑边界原子；原正源的奇异部分不会
被 Ab dx 的减法抹掉。

Ω 上，(Ab)+ 在 C_b 外为零，在 C_b 上≤M，故

\[
(Ab)_+dx|_\Omega\le(M/\bar m)\nu_0|_\Omega
=(16LM/c)\nu_0|_\Omega.
\]

## 6. 闭补集 E 的测度支配：Portmanteau 方向正确

取

\[
\nu_t=t^{-1}(P_tK_0-K_0+U_t)dx,\quad U_t=\int_0^t p_sds.
\]

这里括号按 L¹ 密度等价类理解；不在 S2 上作点态 ∞−∞。正部半群不等式给
ν_t≥0，质量守恒给 ν_t(R⁶)=1。光滑测试的生成元极限给 AK0+δ0=ν0；由于
总质量一致有界，紧支撑连续测试可由光滑函数一致逼近，所以升级为 vague 收敛。

E 上 K0=0 a.e. 且 b=0，b≤C0 K0 与半群保正给
ν_t(x)≥P_tb(x)/(C0 t) a.e.。对 φ∈C_c、φ≥0，φ1_E 确为上半连续，因 E 闭。
无需引用黑箱方向，可令

\[
F_n(x)=\phi(x)\max(1-n\operatorname{dist}(x,E),0)\downarrow\phi1_E.
\]

各 F_n 都连续、紧支撑，且≥φ1_E。先固定 n 使用 vague 收敛，再令 n→∞，得

\[
\limsup_{t\downarrow0}\int_E\phi\,d\nu_t
\le\int_E\phi\,d\nu_0.
\]

另一方面，b=0 于 E，加上全空间一致生成元极限，给

\[
\liminf_{t\downarrow0}\int_E\phi\,d\nu_t
\ge C_0^{-1}\int_E\phi Ab\,dx.
\]

这里 φ 支撑体积有限，所以一致收敛足以传递积分。b 在 E 邻域为零且非负，
故 Ab≥0 于 E。合并两式得到 C0 Ab dx|E≤ν0|E，保留了 E 上一切界面奇异源。
没有使用开放集的错误方向，也没有假定 E 的 Lebesgue 体积为零。新 C1 的严格
反向不等式反例说明此方向不可交换。

取 C1=16LM/c≥C0，Ω 与闭 E 分解给全空间

\[
\boxed{(Ab)_+dx\le C_1\nu_0.}
\]

## 7. 正幅度、质量和一般输入/全尺度

ε=c/[32L(M+1)] 是严格正有理数，εC0<1/2、εC1=M/[2(M+1)]<1/2。
于是 Kε=K0−εb≥K0/2、νε=ν0−εAb dx≥ν0/2，且
AKε=νε−δ0、νε(R⁶)=1。Kε 是非零偶 L¹ 核；b 在 Q 上为零，cube 支配完全
保留。负 δ0 系数不变，正奇异部分也不变。没有依据“核仍正”省略源条件。

paper_finite_source.md §C 的假设只是一般 K≥0、K∈L¹、AK=a(μ−δ0)、μ 正概率，
并没有要求 K 必须是未修正的正部。六维 H¹ 障碍输入已由 001 建立，故应用到
K=Kε/N0、N0=64(g0−c)=320c、a=1/N0 合法。Fourier/L² 严格给

\[
K_r*Au=\frac a r(\mu_r*u-u),\qquad K_r(x)=r^{-6}K(x/r).
\]

在同一个障碍集 D={u>0} 外，u=0 且 μ_r*u≥0 a.e.，后者即使有奇异测度也由
dx dμ_r 的 Fubini 得到。K_r*f≤κ∫K。先合并正有理尺度的零集，再用
||K_r*f−K_s*f||∞≤||f||∞||K_r−K_s||1，把有界紧支撑输入扩至所有实 r>0。
最后 f_j=min(|f|,j)1_{[-j,j]^6} 单调增至 |f|，正卷积/MCT 和 sup_j sup_r
交换给一般 L¹，并保证可测性。没有按 r 重新选择障碍集，没有只证明有限尺度。

Q 上 Kε/N0≥1_Q/64，故上述 r 正是中心、轴平行、六坐标共同半径的 cube 比较。
本轮没有把不同生成元的 source、mass 或 lower floor 混用。

## 8. 非空参数及收益的准确范围

gstar 的既有 det160 基给每个逆元≤5!/160=3/4。限制其余 26 参数于 [-1,1]，
对 x∈Q 已解坐标≤6·(3/4)·27=243/2<122。1≤t≤2 时各已解 Cauchy 密度
至少 1/[4(4+122²)]，各未解参数的区间概率至少 1/10。因此固定正 coarea
积分给 p_t≥[160(4(4+122²))⁶10²⁶]^-1；长度一的时间带给 G≥gstar。
这是实际模型的解析正下界，不是虚构一个尚未支付的 g0，也没有执行 body 评价。
从此 c,T,R,L,C0,C1,ε 都由有限有理运算/整数上取整指定，幅度集合非空。

inner-half 两盒总容积为 2 hbox⁶，逐因子 η≥e^-2>1/9，故 ∫b≥Bmin>0。
质量下降严格等于 ε∫b/N0，至少 Bmin/[10240L(M+1)]。
原 R_c 的 radial/face 质量式仍是 I/c^(6/5)；只是在同一归一化下减去这个正量。
对已认证 I_U 的 slack 改善公式正确，但没有由此推出 slack>0、成本<64 或新 C6。
非空可行幅度和竞争性幅度大小是不同结论。

## 9. 新旧机制的边界

旧三维 U09/U10 已属于空间 bump subtraction，U10 也已有 L¹/连续生成元、Ω 内
密度、Ω 外自动结算和远尾归约。不能把本轮称作全新的修正类别。
本轮真实新增的六维内容是：选在 cube 外而仍位于 Ω 的 off-S3 盒，控制其所有
方向圆柱避开 S2；给 fixed-distance uniform Green 尾及同一线段上的显式源下界；
从而对全部 Ω∩C_b 支付一个可计算 C1，并得到实际正 ε，而非留下未认证紧区域。
闭 E 的测度论证也不需要旧 U10 的 level-set 零体积结算。
这只是本地已有源的比较，不支持文献范围的首创声明。
