# 复倾斜与 anchor 倾斜特征函数：一个可直接使用的比较命题

本稿在 `complex_probe.py` 实际运行并保存 `RESULTS.json` 后提出。核心结论是式 (A)：在真实 anchor 倾斜后，原 GF source 与共享相位的人工解耦模型的 n 坐标特征函数接近。误差在频率零点精确消失；不需要对 n 坐标矩母函数取复对数。它不推出事件概率下界。

## 1. 模型和唯一使用的有限域输入

令 q=2^d，K=q−1，g 为 GF(q) 本原元，b≥1 个 F₂ 独立基元素定义 B=2^b 个 trace 标签。端点平滑函数为

\[
P_\delta((j+s)/K)=\operatorname{clip}((s-U_j+\epsilon/2)/\epsilon,0,1)-s,
\quad U_j=(y_j+1/2)/B,\quad 0<\epsilon=K\delta<1/B.
\]

固定共同相位 θ，令 X=(j+θ)/K，尺度 L_ℓ>0，

\[
H_\ell(j,\theta)=P_\delta(X+L_\ell/2)-P_\delta(X-L_\ell/2),\quad
B_\ell(j,\theta)=\sum_{m=1}^M c_m\sin(\pi mL_\ell)
\cos(2\pi m(j+\theta)/K),
\]

其中 c_m≥0、V=Σc_m，U=2(b+1)、W=U+V、C=2+V。将 H 视为域上函数时，采用同一 trace 公式给零态赋值。归一化 Walsh 范数与循环 DFT 系数绝对和范数满足

\[
\|H_\ell\|_W\le U,\quad |H_\ell|\le2,\qquad
\|B_\ell\|_A\le V,\quad |B_\ell|\le V. \tag{1}
\]

理由简述：clip 关于自然整数标签单调；每个二进制前缀/后缀指标可拆为至多 b 个二进制柱集合，每个柱集合 Walsh 范数1。单调 [0,1] 函数是这些指标的正组合，总权≤1，所以 clip 范数≤b；减去常数 s 增至≤b+1，两个端点相减给 U。端点的粗格偏移使 trace 参数乘固定 g 次幂，保持拉回字符的系数绝对和上界。每个余弦的 DFT 范数≤1，频率合并不增范数。这些范数均有乘积次乘法性。

对任意域函数 F 和循环函数 Ψ，j 在 0,…,K−1 均匀时，

\[
|E F(g^j)\Psi(j)-E F(g^j)E\Psi(j)|
\le\lambda\|F\|_W\|\Psi\|_A,\qquad \lambda=\sqrt q/K. \tag{2}
\]

这是复函数的双线性协方差估计，不取复共轭。展开两套有限 Fourier 级数，常数频率在协方差抵消；非平凡加性/乘性 Gauss 和的模为 √q。其平方换元 x=ty，t=1 内和为 K，其余为−1，得 K+1=q。加性字符在非零态均值为−1/K，已包括在常数频率的精确抵消中。这里不用 Sidon 或任何不同尺度独立性。

置 a_ℓ=p/(KL_ℓ)>0，a_*=max a_ℓ，并假定

\[
r=Wa_*<1. \tag{3}
\]

定义两个正实随机向量：

\[
Y_\ell=1+a_\ell(H_\ell(j,\theta)+B_\ell(j,\theta)),
\qquad
\widetilde Y_\ell=1+a_\ell(H_\ell(j_1,\theta)+B_\ell(j_2,\theta)).
\]

原模型共享 j；比较模型在给定同一 θ 后 j₁、j₂独立均匀，全部 H 尺度仍共用 j₁、全部 B 尺度仍共用 j₂。θ 可服从任意共同分布；真实 Haar 情形取 Uniform[0,1)。由 (3)，Y、Ỹ≥1−Ca_*>0。

一切复幂只定义为 `Y^τ=exp(τ log Y)`，其中 log Y 是正实数上的实对数。定义 μ(τ)=E ΠY_ℓ^{τ_ℓ}，μ̃ 同理。

## 2. 一般复倾斜及安全的局部复对数

对任意有限复向量 τ，记

\[
S_1(\tau)=\sum_\ell|\tau_\ell|a_\ell,\quad
S_2(\tau)=\sum_\ell|\tau_\ell|a_\ell^2,
\]

\[
D(\tau)=\lambda UV\frac{S_1(\tau)^2+S_2(\tau)}{(1-r)^2}
\exp\!\left(\frac{WS_1(\tau)}{1-r}\right).
\]

则

\[
|\mu(\tau)-\widetilde\mu(\tau)|\le D(\tau). \tag{4}
\]

**证明。** 复二项式系数仍满足
`|binom(τ,k)|≤(|τ|)_k/k!`，因为 `|τ−j|≤|τ|+j`。绝对系数因此被

\[
\prod_\ell[1-a_\ell(Uu+Vv)]^{-|\tau_\ell|}
\]

支配。纯 H、纯 B、常数项在两模型中相同；对混合项用 (2)。正系数级数的混合系数和等于单位方块上的混合二阶导数积分；逐项求导后用分母≥1−r，得到 (4)。r<1 确保范数绝对收敛且对 θ 一致，故二项式展开、逐项期望、θ 积分和求和交换均合法。复数性没有改变这一论证。

复数 μ 不一定非零；不能直接沿用实对数的下界。一个充分条件是

\[
Z(\tau)=\frac{CS_1(\tau)}{1-Ca_*},\qquad
\rho=e^{Z(\tau)}-1<1.
\]

因为每个 |Στ_ℓ log Y_ℓ|≤Z，得 |μ−1|、|μ̃−1|≤ρ。二者及连线都在不含零的圆盘 D(1,ρ)，主值 Log 在该盘解析，积分导数 1/z 给

\[
|n\operatorname{Log}\mu-n\operatorname{Log}\widetilde\mu|
\le\frac{nD(\tau)}{1-\rho}. \tag{5}
\]

这里 n 为正整数，式 (5) 比较的是指定的单坐标局部 Log 乘 n。虽然 `exp(n Log μ)=μ^n`，一般不能把 `n Log μ` 写成 `Log(μ^n)`；后者可能绕分支。条件 ρ<1 若失效，(4) 仍成立，(5) 不作主张。

## 3. 主要命题：e₀ anchor 下的特征函数

选择家族中一个 anchor 尺度0，记 e₀=(1,0,…,0)。线性给出

\[
A:=E Y_0=E\widetilde Y_0>0. \tag{6}
\]

这对任意共同 θ 分布精确成立，因为两模型保留每个纯 H、纯 B 均值。真实 Haar 情形更有 A=1：周期函数的两个平移差 H₀ 均值0，每个整数余弦 B₀ 的 Haar 均值0。固定 θ 或有限 θ 混合时 A 未必1，必须保留分母。

在两个模型上分别定义概率倾斜 `dP₀=(Y₀/A)dP`、`dP̃₀=(Ỹ₀/A)dP̃`。对 η∈R^J，令

\[
\varphi(\eta)=E_{P_0}\exp\!\left(i\sum_\ell\eta_\ell\log Y_\ell\right)
=\frac{\mu(e_0+i\eta)}A,
\]

φ̃ 同理。取 n 个 iid 倾斜坐标后，联合 log-score 和的特征函数分别为 φ(η)^n 与 φ̃(η)^n。

把 S₁=Σ|η_ℓ|a_ℓ、S₂=Σ|η_ℓ|a_ℓ²，并置

\[
D_0(\eta)=\lambda UV\exp\!\left(\frac{WS_1}{1-r}\right)
\left[
\frac{(1+a_0W)(S_1^2+S_2)}{(1-r)^2}
+\frac{2a_0S_1}{1-r}
\right].
\]

**结论：**

\[
\boxed{
|\varphi(\eta)^n-\widetilde\varphi(\eta)^n|
\le \min\{2,\ nD_0(\eta)/A\}.}
\tag{A}
\]

尤其 D₀(0)=0，而直接把 τ=e₀+iη 代入 (4) 的较粗界在 η=0 不消失。

**证明。** 使用精确分解

\[
\prod_\ell Y_\ell^{(e_0+i\eta)_\ell}-Y_0
=Y_0\left(\prod_\ell Y_\ell^{i\eta_\ell}-1\right).
\]

令 w=Uu+Vv、Q(w)=Π(1−a_ℓw)^{−|η_ℓ|}。上式的绝对系数被

\[
R(u,v)=(1+a_0w)(Q(w)-1)
\]

支配。仍只有全局同时含 H 与 B 的项贡献两模型差；减去的 Y₀ 因 (6) 对比较差无贡献。令

\[
T_1(w)=\sum_\ell\frac{|\eta_\ell|a_\ell}{1-a_\ell w},\quad
T_2(w)=\sum_\ell\frac{|\eta_\ell|a_\ell^2}{(1-a_\ell w)^2}.
\]

Q′=QT₁，Q″=Q(T₁²+T₂)，因此

\[
\partial_u\partial_vR
=UVQ\big[(1+a_0w)(T_1^2+T_2)+2a_0T_1\big].
\]

在 0≤u,v≤1 上，用 w≤W、T₁≤S₁/(1−r)、T₂≤S₂/(1−r)²、Q≤exp(WS₁/(1−r))，由 (2) 得

\[
|\mu(e_0+i\eta)-\widetilde\mu(e_0+i\eta)|\le D_0(\eta).
\]

由于两个 anchor 分母是同一个 A，故 |φ−φ̃|≤D₀/A。φ 和 φ̃ 是真正概率分布的特征函数，模均≤1；代数恒等式
`z^n−w^n=(z−w)Σ_{k=0}^{n−1}z^k w^{n−1−k}`
直接给 n 倍界，另有上界2，得到 (A)。该步骤没有复对数、分支或取指数导致的额外损失。

共同的确定性中心化仅使两个特征函数同乘单位模相位，不改变界。若各自按不同的倾斜均值中心化，还需要另外控制两套中心之差；本命题不自动提供此控制。

## 4. 归一化 η=t/√v 与尺度数

若 L_ℓ≥L_min>0、Σ|t_ℓ|≤T，置 η=t/√v。式 (A) 中

\[
S_1\le a_*T/\sqrt v,\qquad S_2\le a_*^2T/\sqrt v.
\]

真实 Haar A=1，固定 p,L_min，U=O(log B)、V 有统一上界，且 r→0 时

\[
nD_0(t/\sqrt v)
=O\!\left(
\frac{n\log B}{K^2\sqrt K}
\left[\frac{T^2}{v}+\frac{T}{\sqrt v}\right]
\right).
\]

在接受基线 `n/K²=Θ_p(v)`、`B=Θ(v)`、`K/v^r→∞`（每个固定 r>0）下，这成为

\[
O_{p,L_{\min}}\!\left(
\frac{\log B}{\sqrt K}(T^2+T\sqrt v)
\right)=o(1) \quad(T\text{ 固定}). \tag{7}
\]

因此对所有有限共同尺度包、所有 Σ|t|≤T，该倾斜特征函数比较统一成立，不另付尺度数 J；不能把总 ℓ¹条件换成每个频率坐标单独有界。增长 T 只有在右侧及指数因子仍受控时才可用，不能据此直接放宽为任意增长频率域。

对于一般实 anchor α，可用 (4) 和正 anchor 下界处理两个不同分母；本稿优先保留 e₀ 的精确相同分母及零频消失界，不用假设两个非线性 anchor 归一化常数相同。

## 5. 实际数值证据、复现与尚未解决的环节

`complex_probe.py` 以 GF(2⁷)、模0x83、本原元2、b3、真实 M4 Fourier 系数及三尺度 10/9、7/6、11/9 构造端点 clip。两相位1/7与4/9各自保留127个原态及127²个比较态，再作相同的有限混合。这不是 Haar θ 积分近似的认证；其共同 anchor A 是按有限模型直接计算的，未强设为1。

实际完成15个 anchor 倾斜特征函数案例（5组 t、n=1/100/10000）及2个一般复倾斜局部 Log 案例；全部满足相应界。最大单坐标矩差约1.25×10⁻⁷，最大 n 坐标特征函数差约.00123114。一般复 Log 例中 ρ 约.0143–.0156，符合近1圆盘条件。

大参数部分只用 Decimal 代入解析界：n=10³²、K=2⁴⁷−1、p=.99、b11、V≤40/9，并按接受基线 `v=9p²n/(29K²)` 得 v≈1535.660835。t=(1,0,0)、(0,1,−1)、(3,−2,1) 的 (A) 界约.00278235、.00509919、.01667354。这些是80位 Decimal 求值，不是向外舍入端点、n维模拟或概率认证。

复现：

```sh
/opt/homebrew/opt/python@3.15/bin/python3.15 complex_probe.py
```

依赖只有已有 NumPy 与 Python 标准库；脚本不读取或修改46轮产物。

式 (A) 是频域分布比较工具。有限频率上的小误差不能直接给罕见事件下界，也不提供最大事件的总变差界。若要做反演，需要合适的频率截断、平滑误差和比较模型的反集中/局部下界；若 J 增长还需支付反演核及尺度重叠成本。这里没有推出 C_n 的新 log 增益。
