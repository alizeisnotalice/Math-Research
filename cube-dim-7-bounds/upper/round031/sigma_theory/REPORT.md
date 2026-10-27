# R031 · 固定 t 的 σ 可行域、唯一边界与合法一维接口

这里只分析 R030 已证明的 general-σ two-region 公式，不改变覆盖算法，不声称全 t 的全局最优，也不启动下界研究。所有结论适用于整数 n≥1 和实数 t>1，包括 n=1。该报告的证明为解析证明；附带的精确有理小检查只诊断代数实现。

## 1. 输入与完整实数域

固定 n,t，记

\[
\beta=1+1/t,\quad A=\beta^n,\quad
\ell_t(\rho)=1+\frac1{2t}-\frac1{2\rho},\quad
h=\ell_t(1)=\frac{1+1/t}{2},\quad K=1-h^n.
\]

K 是旧 cap 包络在单位核心内的质量，0<K<1。记

\[
B=K+n\int_1^t
\left[1-\tfrac12\ell_t(\rho)^{n-1}\right]\frac{d\rho}{\rho}.
\]

对 σ∈[1/t,1] 定义

\[
D(\sigma)=\sigma^n-
\left(\frac{\sigma+1/t}{2}\right)^n,
\]

\[
P(\sigma)=n\int_{t\sigma}^{t}
\left[1-\tfrac12\ell_t(\rho)^{n-1}\right]\frac{d\rho}{\rho},
\quad E(\sigma)=B-D(\sigma),
\]

\[
g(\sigma)=P(\sigma)-D(\sigma)-\sqrt{A E(\sigma)}.
\]

R030 的严格原算子接口要求 1/t<σ<1、E>0、g>0，并给出

\[
C_n\le F_{\rm cut}(t,\sigma)
=\left(\sqrt A+\sqrt{E(\sigma)}\right)^2.
\]

### 全域正性：不存在 E≤0 造成的缺域

端点与导数为

\[
D(1/t)=0,\quad D(1)=K,
\]

\[
D'(\sigma)=n\sigma^{n-1}
-\frac n2\left(\frac{\sigma+1/t}{2}\right)^{n-1}>0. \tag{1}
\]

严格正性在整个闭区间成立：当 n=1 时 D'=1/2；当 n>1 时 (σ+1/t)/2≤σ，故 D'≥(n/2)σ^(n−1)>0。于是 0≤D≤K。

记

\[
P_0=P(1/t)=B-K.
\]

由于 0<ℓ≤1，有

\[
P_0\ge\frac n2\log t>0.
\]

因此对全部 σ∈[1/t,1]，

\[
E(\sigma)=B-D(\sigma)\ge B-K=P_0>0. \tag{2}
\]

所以 g 是完整闭区间上的连续实函数，并在开区间 C1；以下不存在跨越复平方根区域后重新出现的可行岛。

## 2. 唯一边界定理

由可变下限积分直接得到

\[
P'(\sigma)=-\frac n\sigma
\left[1-\tfrac12\ell_t(t\sigma)^{n-1}\right]<0. \tag{3}
\]

所以 0≤P(σ)≤P0<B。又 E'=−D'，故

\[
g'(\sigma)=P'(\sigma)-D'(\sigma)
\left[1-\frac{\sqrt A}{2\sqrt{E(\sigma)}}\right]. \tag{4}
\]

**定理 1（全 t,n 的可行域分类）。** 令

\[
g_0=g(1/t)=B-K-\sqrt{AB}.
\]

若 g0≤0，则严格可行域为空。若 g0>0，则存在且仅存在一个 σbar∈(1/t,1) 满足 g(σbar)=0，且严格可行域恰为

\[
\{\sigma\in(1/t,1):g(\sigma)>0\}=(1/t,\bar\sigma).
\]

每个零点都严格向下穿越；事实上所有 g≥0 的点均满足 E>A，并有 g'<0。

**证明。** 因 P<B，总有 P−D<E。若 g≥0，则

\[
\sqrt{AE}\le P-D<E.
\]

用 (2) 除以 √E，得 E>A（严格）。更全局地，E 随 σ 严格下降，故

\[
J=\{\sigma\in[1/t,1]:E(\sigma)\ge A\}
\]

是一个初始区间，可能为空或等于全区间。在 J 上，(4) 中方括号≥1/2，因此 g'≤P'−D'/2<0。在 J 外 E<A，由 P−D<E<√AE 直接得 g<0；在 E=A 的边界同样由 P−D<E=√AE 得 g<0。于是全部非负值都只能位于一个 g 严格下降的初始区间。这一步排除了所有 later positive islands，而非只检查左端符号。

若 g0≤0，在初始区间 J 上 g 不会增至正值，在其外也严格负，故没有严格可行点。若 g0>0，则

\[
g(1)=-K-\sqrt{A(B-K)}<0,
\]

连续性给出至少一个内部零点。零点必在 E>A 的初始区间内，那里 g 严格下降，故零点唯一。其前严格正，其后在 J 内严格负、J 外也严格负。式 (4) 还给零点处 g'<0。证毕。

该定理不声称 g 在所有负值区域都单调；结论所需的是在 E≥A 的初始区间单调，以及其余区域完全不可行。

## 3. 固定 t 的最佳边界与 closure 的原算子意义

对整个闭区间，由 (1)–(2)，

\[
\frac{\partial F_{\rm cut}}{\partial\sigma}
=-D'(\sigma)
\left[1+\frac{\sqrt A}{\sqrt{E(\sigma)}}\right]<0. \tag{5}
\]

**推论 2（固定 t 的最优扣减）。** 若 g0>0，定理 1 的唯一根 σbar 满足

\[
\inf_{1/t<\sigma<\bar\sigma}F_{\rm cut}(t,\sigma)
=F_{\rm cut}(t,\bar\sigma).
\]

该边界值仍是原算子常数的严格有效上界：

\[
C_n\le\Psi_n(t):=
\left(\sqrt{A_n(t)}+
\sqrt{B_n(t)-D_n(t,\bar\sigma_n(t))}\right)^2. \tag{6}
\]

**证明。** (5) 与连续性证明标量下确界结论。选择任意 σj↑σbar 且 1/t<σj<σbar。每个 σj 均有严格 g>0，故 R030 原算子定理分别给 Cn≤Fcut(t,σj)。取这一组已经有效的实数上界的极限，得到 (6)。这里并未在 g=0 的参数上直接实施一个缺少 guard 余量的有限 a 算法；每个 σj 都先独立满足严格可行条件，然后才对数值上界取极限。选择历史、有限族和任何 selector 均不需要收敛。证毕。

## 4. 一维 Ψ 接口与边界认证

令

\[
\mathcal T_n=\{t>1:B_n(t)-[1-h(t)^n]-\sqrt{A_n(t)B_n(t)}>0\}.
\]

对 t∈Tn，Ψn(t) 由唯一根和 (6) 明确定义，并有

\[
C_n\le\inf_{t\in\mathcal T_n}\Psi_n(t)
\]

（Tn 非空时）。也可在 Tn 外把 Ψ 定义为 +∞，作为纯粹的公式族下确界约定；若想为单一不可行 t 保留原旧上界，则应另行显式使用 Fn(t)，不能把不存在的 σbar 填入 Ψ。

数值接口可以先准确判定 g0 的符号，再在 [1/t,1] 内夹逼唯一根。虽然 g 未被证明全域单调，唯一根及其正/负两侧符号足够保证有符号二分夹逼的合法性。若候选使用实际严格 σ，则只需认证该 σ 的 E>0 和 g>0；这通常比准确使用 closure 更直接。若使用 σbar 的区间，目标函数随 σ 下降，所以对根的下端 σlower≤σbar 计算 Fcut 可给向外上界；但仍需 B、D 和根式的区间误差同向传递。根的上端代入目标一般给下包围，不能冒充算子上界。

本报告没有证明 Tn 的全局区间结构、t 上唯一极小点、全 t 根完整性或参数全局最优；没有把有限区间 [2,32] 的搜索升级为上述任何断言。后续若研究外壳后序区或多个共同前序区域，应单独证明原实例映射，不能沿用单位核心区的体积公式替代新几何。

## 5. 实际运行与证据边界

`verify_algebra.py` 在本机使用标准库 `Fraction` 对 n=1,2,7、t=2,7 的六个有理实例检查 D 的端点、正导数和 P' 方括号。检查点只诊断公式实现；定理 1 和推论 2 的全域结论由上面的解析证明给出，不由有限采样推出。没有进行新的 t、σ 参数搜索，没有共享 MATLAB 会话，也没有修改任何 round030 文件。

