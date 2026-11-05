# 第 5 轮：单调双曲见证机制的首系数障碍与有界核心归约

日期：2026-11-06（Asia/Shanghai）。状态：受限机制的严格定理；不作为一般 C_n 上界。

## 1. 机制的精确定义与结论

每个维数允许任意比例向量 r_i∈[0,1]，T=ΣZ_(r_i)，及任意非负、非减可测轮廓 h:[0,∞)→[0,∞)。
记 J(h)=∫_0∞e^-t Phi_n(h(t))dt，定义
\[
 Q_n=\sup_{\mathbf r,h:\ E h(T)\ge1,\ 0<J(h)<\infty}J(h)^{-1}.
\]
Q_n 是第2轮“单个双曲接收集合+一组固定矩形比例”机制的最优极限比值，不是全部实际 halo 的最优值。
对等号期望可先乘1+eta、再截断，最后 eta↓0，故第2轮归约给 Q_n≤C_n。

**本轮定理**：
\[
 \boxed{\frac{n\mu-\log Q_n}{\sqrt{n\log n}}\longrightarrow\sigma.}\tag{1}
\]
它覆盖任意层数、连续软轮廓和任意非等参数 r_i。第3轮的 log log n 收益与本定理兼容：多层可改善后续项，但在这个机制内不能降低首个负项系数。
它不排除加权对数坐标、非单调轮廓、多中心、相关结构、有限截断深度 L 或全空间实际 halo 给出新机制。

## 2. 轮廓的点态控制与远尾

固定 epsilon∈(0,sigma)，令 d=(sigma−epsilon)sqrt(n log n)、s=nmu−d。
假设 J≤e^-s。非减性给 J≥e^-t Phi_n(h(t))，因此
\[
 \Phi_n(h(t))\le e^{t-s}.\tag{2}
\]
于 t≤s−2log n，h(t)≤n^-2；于 t≤s，h≤1。
若 s≤t≤s+n/4 且 h>1，令 x=t−s、y=log h，则
\[
 (n-1)y/(1+y)\le\log\Phi_n(h)\le x,
\]
所以 y≤x/(n−1−x)≤2x/n（n≥4），进而 h−1≤4x/n。

记 M=E T≤nmu、V=Var T≤C n。中间部分的期望超额至多
\(4E(T-s)_+/n\le4(d+\sqrt V)/n=O(\sqrt{\log n/n})\)。
远尾 D_+={T>s+n/4} 对充分大的n包含于 {T>nmu+n/8}。
把其上 h 分成 h<exp(4T/n) 与补集。

第一部分一致地指数小：由精确 mgf
\[
 E e^{\theta Z_r}=(1-\theta)^{-1}(1+r^{1-\theta})/(1+r),\qquad\theta<1,
\]
在固定小正参数区间内，log mgf 的二阶导数对 r∈[0,1] 一致有界。
故 log E exp(eta Z_r)≤eta m(r)+C eta²≤eta mu+C eta²。
取固定充分小的 theta>0、eta=theta+4/n，Chernoff 给
E[exp(4T/n)1_(D_+)]≤exp(−c1 n)。

第二部分利用所有 r_i 的统一密度比上界
\[
 R(t)=e^tf_T(t)\le2^n t^{n-1}/(n-1)!.
\]
对 r_i>0 条件于 Bernoulli 集合 A，位移 a_A=Σ_(i∈A)(−log r_i)，密度比是
Σ_A P(A)e^(a_A)(t−a_A)_+^(n−1)/(n−1)!；每项系数≤1，共2^n项。r_i=0 直接删去退化位移项或取极限。
当 h≥exp(4t/n) 时，Phi_n(h)≥h(4t/n)^(n−1)，故
\[
 E[h(T)1_{D_+\cap\{h\ge e^{4T/n}\}}]
 \le J\frac{2^n(n/4)^{n-1}}{(n-1)!}=e^{-c_2n}.
\]
括号 n 次根趋于 e/2<B，J≤exp(−nmu+o(n))，所以 c2>0 可一致选择。
结合中间部分，令 A_n={T≤s−2log n}，得到
\[
 E h(T)\le1-(1-n^{-2})P(A_n)+O(\sqrt{\log n/n})+e^{-c n}.\tag{3}
\]

## 3. 统一低尾下界导致矛盾

写 D=nmu−M=Σ(mu−m(r_i))≥0。
如果 D>2d，则阈值 s−2log n 比 M 至少高 d−2log n。Chebyshev 和 V≤Cn 给 P(A_n)=1−O(1/log n)，与 (3) 及 E h≥1 矛盾。

如果 D≤2d，由 m 在 r=B^-1 唯一、非退化最大，以及区间[0,1]的紧性，存在绝对C使
\[
 |v(r)-\mu|\le C\sqrt{\mu-m(r)}.
\]
局部证明是 Taylor 展开：mu−m(r) 与 (r−B^-1)² 可比，v 连续可微；远处均值缺口严格为正，端点 m=v=1。
Cauchy–Schwarz 给 |V−nmu|≤C sqrt(nD)=o(n)。

三阶中心矩也一致有界：用 Z=Y+cJ、c=−log r，可取
E|Z−m|³≤4(E|Y−1|³+c³E|J−p|³)≤4(16+sup_(0<r≤1)r|log r|³)<∞。
且 v(r)≥1。因此独立非同分布的 Berry–Esseen 不等式给分布函数绝对误差≤C/sqrt(n)。
阈值低于均值的标准化距离至多
\[
 z_n\le(1-\varepsilon/(2\sigma))\sqrt{\log n}
\]
（阈值高于均值时该低尾下界更容易）。于是以 a=1−epsilon/(2sigma)∈(0,1)，由正态 Mills 下界
\[
 P(A_n)\ge \Phi_{N}(-a\sqrt{\log n})-C/\sqrt n
 \gg\sqrt{\log n/n}.
\]
普通绝对误差在此足够，因为高斯尾为 n^(-a²/2)/sqrt(log n)，a<1，严格大于 n^-1/2 log因子；没有在临界a=1处误用 Berry–Esseen。
这再次与 (3) 矛盾。因此充分大n时 J>exp(−nmu+(sigma−epsilon)sqrt(nlogn))，一致于所有 r_i 和 h。
得到 (1) 的上界方向；第3轮指数轮廓给匹配的下界方向。

## 4. 一般有界核心到示性 halo 的严格归约

更一般地，若0≤f≤H、0<tau<H、delta∈(0,1)，令 E={f>tau(1−delta)}。
逐点 f≤tau(1−delta)+(H−tau+tau delta)1_E，故
\[
 |\{M_Sf>\tau\}|\le
 \frac{T_n(\tau\delta/(H-\tau+\tau\delta))}{\tau(1-\delta)}\int\Phi_n(f).\tag{4}
\]
这里只用了 |E|≤∫Phi_n(f)/(tau(1−delta))。T_n 为示性输入 Tauberian 常数。
这不是新的维数一致估计，但把核心几何缺口明确地归到 T_n(alpha) 的维数依赖，保留了阈值和代价。
若 tau=1、H=1+epsilon，选delta=sqrt(epsilon)时 alpha=1/(1+sqrt(epsilon))（epsilon<1）。

## 5. 常幅单盒障碍与本机校准

对 Q=[0,1]^n，强矩形的独立性给 M_S1_Q=∏M_1 1_[0,1]。
每坐标的 −log M_1 1_[0,1] 在0有质量1、正半轴有密度2e^y。
设 f=e^A 1_Q，A>0，则实际 halo 体积精确为
\[
 V_n(A)=1+\sum_{k=1}^n{n\choose k}\frac{2^k}{(k-1)!}\int_0^A e^s s^{k-1}ds.
\]
用积分≤e^A A^k/k、二项系数≤n^k/k!，得到
\[
 V_n(A)\le\sum_{k\ge0}\frac{(2e^AnA)^k}{(k!)^2}\le e^{2\sqrt{2e^AnA}}.
\]
最后一步由 I_0(x)≤e^x，也可从 (2k)!≤4^k(k!)² 和 cosh(x)≤e^x 逐项证明。
因此 A与sqrt(logn/n)可比时，log V_n=o(nA)，
\[
 \log L_{\rm box}=\log V_n-A-(n-1)\log(1+A)=-nA+o(nA),
\]
无法产生 B^n 的主指数底数。此障碍限于常幅单盒及其真实 halo。

|n|log真实halo体积（数值）|log模数（数值）|log归一化比（数值）|
|--:|--:|--:|--:|
|8|3.95932315|3.39382913|0.56549402|
|32|7.17048638|9.14856795|-1.97808157|
|128|11.88754642|22.78693681|-10.89939039|
|512|18.82214517|53.61426468|-34.79211951|

MATLAB R2026a 在 A=sqrt(logn/n) 以定向无关的确定性积分算出此表；它是公式的数值核对，不是新数值证书。
理论障碍由上面的解析界证明。另对2001个r网格的方差稳定性比值做实际抽查，最大约0.84634；常数存在性由紧性证明，不靠采样。

## 6. 当前范围与下一轮

一般 C_n 保持第3轮严格软轮廓下界与 sqrt(B)B^n 上界。Q_n 的首系数现在闭合，但 Q_n≤C_n 不能反向用于 C_n。
本轮没有解决原附件全空间双曲示性族的有限L上界。

下一轮最有区分力的两项动作：① 转向加权对数和与分块聚合，改变接收几何或允许更多独立见证；② 对一般有界核心的 Tauberian 归约跟踪维数常数，以第3轮软轮廓作反例检验。
复现单盒：MATLAB `rows=sm_box;`。Berry–Esseen 的可直接核验来源为 Vershynin arXiv:2602.06234v1 Theorem 0.1，详见 sources.md；该标准定理是本轮唯一概率逼近输入。
本轮无已跟踪冗余文件需要清理。
