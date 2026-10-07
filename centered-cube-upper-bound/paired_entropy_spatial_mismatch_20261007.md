# 原有序传播的时间配对：两个空间组件的错配障碍

2026-10-07，gated_radial_energy_audit。固定原 \(c=1\)、物理尺度 1、\(Z=1\)、\(\lambda=1/2\)。本稿只研究
\[
\mathfrak P(f)=\int_0^1\sqrt{d_f(z)a_f(z)}\,dz
\]
这一**先分别积分全空间、再按时间配对**的上界。结论：有真正的非负 \(L^1\cap L^\infty(\mathbb R^n)\) 来源使 \(\mathfrak P(f)\ge \sqrt n\,\|f\|_1/2048\)。因此一般 \(\operatorname{polylog}(n)W\) 合同也不能由这项时间配对得到。这不否定更细空间／边配对、signed commutator 或原 weak endpoint，更不认证原 actual geom 的全部门。

已读并使用 L03 SKILL、provenance、method、cube-interface：只使用有理残差和认证区间工作流，不引用技能作为本定理。fast 单 mode 使用 root 已审的 ordered_capped_test_energy_20261007.md §3–5，未重跑其数值。新的 slow 组件、同一欧氏来源拼接及配对下界在这里推导。

## 1. 泛函与已有 fast 组件

对一个完整输入 \(f\ge0\)，定义
\[
v_z=T_{1-e^{-z}}^1f,\quad J_z=\sum_i(G_{e^{-z},i}-I),\quad
M(x)=\max_{0\le z\le1}v_z(x),\quad E=\{M>\lambda\}.
\]
取最早 winner \(\tau\)，令
\[
t_z=\mathbf1_E\frac{\lambda+v_z}{M}\mathbf1_{\{z<\tau\}},\qquad 0\le t_z<2.
\]
准确的两个非负速率是
\[
\begin{split}
d_f(z)&=\frac{\lambda}{2}\sum_i\int dx\,G_{e^{-z},i}(dh)
\frac{(\Delta v_z)^2}{(\lambda+v_z(x+h))(\lambda+v_z(x))},\\
a_f(z)&=\frac{\lambda}{2}\sum_i\int dx\,G_{e^{-z},i}(dh)(\Delta t_z)^2.
\end{split}\tag{1}
\]
一般 signed commutator 有 \(|C|\le\mathfrak P(f)\)，但上界本身未必由一次 \(W\) 的 polylog 控制。

引用已审 fast 组件：整数 \(K\ge64\)，\(0<\varepsilon\le10^{-3}\)，周期来源
\[
F_{\rm fast}(x)=1+\varepsilon\cos(K\sum_i x_i).
\]
以单位周期体积归一，原 response 是
\[
v_z^{\rm fast}=1+\varepsilon m_n(z)\cos(K\sum_i x_i),\quad
m_n(z)=[e^{-z}+(1-e^{-z})g_1(K^2)]^n.
\]
原 winner 为正 cosine 的 0、负 cosine 的 1，且完整 \(E\) 是全周期。已有每轴 Dirichlet 下界给
\[
a_{\rm fast}^{\rm per}(z)\ge n/4,\qquad 0<z<1.
\tag{2}
\]
这是真实 winner test，不是外加 arbitrary mask。其熵在晚期极小并不影响 (2)；本稿不重证这个单 mode。

## 2. 另一处空间的慢熵来源

取同一个 \(c,\lambda\)，但只在第一坐标使用低频：
\[
F_{\rm slow}(x)=1+\eta\cos x_1,\qquad\eta=1/4.
\]
原符号 \(B(v)=v/\log(1+v)-1\) 给
\[
g=g_1(1)=\log2,\quad
b(z)=e^{-z}+(1-e^{-z})g,\quad
v_z^{\rm slow}=1+\eta b(z)\cos x_1.
\tag{3}
\]
这里其它 \(n-1\) 轴保持常数，不引入维数衰减。记周期相位平均为 \(\langle\cdot\rangle\)，\(A_0=\lambda+1=3/2\)。由精确熵导数（不是从浮点估算反推）：
\[
\begin{split}
d_{\rm slow}^{\rm per}(z)
&=-\eta b'(z)\left\langle\cos x_1\,
\Phi'_\lambda(1+\eta b(z)\cos x_1)\right\rangle\\
&=\frac{\lambda\eta^2 b(z)e^{-z}(1-g)}{A_0}
\left\langle\frac{\cos^2x_1}{A_0+\eta b(z)\cos x_1}\right\rangle.
\end{split}\tag{4}
\]
第二行使用 \(\langle\cos\rangle=0\) 和恒等式
\[
-\frac{q}{A_0+aq}
=-\frac q{A_0}+\frac{a q^2}{A_0(A_0+aq)}.
\]
因此 \(b(z)\ge g\)、\(\langle\cos^2\rangle=1/2\)、\(A_0+\eta b\cos\le A_0+\eta\) 给
\[
d_{\rm slow}^{\rm per}(z)
\ge\frac{\lambda\eta^2 g e^{-z}(1-g)}
{2A_0(A_0+\eta)}
>\frac1{3024},\qquad 0\le z\le1.
\tag{5}
\]
常数链只需 \(2/3<\log2<3/4\) 与 \(e<3\)：
\[
\frac{(1/2)(1/16)(2/3)(1/3)(1/4)}
{2(3/2)(7/4)}=\frac1{3024}.
\]
取共同晚期时间段 \(I=[1/2,3/4]\)，长度 \(1/4\)。在这里 fast 的 \(a\) 与 slow 的 \(d\) 同步存在，却位于两处不同空间。这正是待查的错配。

## 3. 必须是一个 \(\mathbb R^n\) 输入：双盒与互相污染

固定 \(n\ge8\) 后才取极限。令 \(L\to\infty\) 为 \(2\pi\) 的整数倍，\(Q_L=[-L/2,L/2]^n\)、\(D_L=L^2\)，定义同一个来源
\[
f_L(x)=\mathbf1_{Q_L}(x)F_{\rm fast}(x)
+\mathbf1_{D_Le_1+Q_L}(x)F_{\rm slow}(x-D_Le_1).
\tag{6}
\]
两盒不交；\(f_L\ge0\)、\(f_L\le5/4\)，且每个 cosine 的盒积分为零，故
\[
W_L=2L^n.
\tag{7}
\]
对 (6) 重新计算一个共同 \(v^L,M_L,E_L,\tau_L,t^L,d_L,a_L\)；不把抽象直和中的两个 winner 当成 (6) 的精确 winner。

在入口 \(c=1\)，\(G_1=w\) 是尺度 \(\sqrt t\le1\) 的 Laplace 混合：
\[
g_1(\xi^2)=\int_0^1(1+t\xi^2)^{-1}dt,\qquad
\Pr_w(|Y|>R)\le e^{-R}.
\]
原 \(T_{1-e^{-z}}\) 的各轴是 holding 与这个固定 \(w\) 的混合。取 bulk 厚度 \(R_L=\sqrt L\)。在第一盒 bulk，对所有 \(0\le z\le1\)：
\[
|v_z^L(x)-v_z^{\rm fast,per}(x)|
\le (5/4)n e^{-R_L}+(5/4)e^{-(D_L-L)}.
\tag{8}
\]
第一项控制自身截断；第二项控制第二盒的输入：进入第二盒必须在第一坐标跳过至少 \(D_L-L\)。第二盒 bulk 对 slow 周期响应有同样的界。于是**两盒互相污染在整个紧时间区间统一消失**；不需要假设真实跳路不能跨盒。

bulk 体积比 \((1-2R_L/L)^n\to1\)，因为先固定 \(n\)。所有逐时响应均由 \(L^\infty\) 压缩性界于 \(5/4\)。fast 和 slow 周期轨道在 cosine 非零处的最大点唯一。由 (8) 的一致轨道收敛及紧区间 argmax 稳定性，原 \(\tau_L\) 分别趋向各自周期 winner；同时 \(M_L\) 收敛，bulk 最终属于 \(E_L\)。对 \(0<z<1\)，\(t_z^L\) 在避开零相位的 bulk 中趋于相应周期 test。并未要求有限 \(L\) 的 winner 恰为 0 或 1。

极小 fast winner gap 的控制允许依赖固定 \(n\)：先丢弃 \(|\cos|<\kappa\) 的相位，再令 \(L\to\infty\)，最后 \(\kappa\downarrow0\)。不得把 \(L\) 的选择声明为随 \(n\) 一致的尺寸。

## 4. 正边积分提升与配对下界

对固定 \(z\in I\)，先将每轴 jump 限于 \(|h|\le H\)，并只在相应来源盒更小的 bulk 内积分。这是 (1) 中的正子积分；同时两端仍在同一 bulk。利用 (8)、winner 稳定和 \(t^L\le2\)，相位周期平均与支配收敛给：

- 在 slow 盒中，正 entropy edge 积分／\(L^n\) 的下极限至少为相应周期的截断积分；
- 在 fast 盒中，正 test edge 积分／\(L^n\) 的下极限至少为相应周期的截断积分。

按 \(L\to\infty\)、相位排除宽度 \(\kappa\downarrow0\)、\(H\to\infty\) 的顺序，用正 Tonelli 放开 jump 截断，得几乎处处 \(z\in I\)：
\[
\liminf_L\frac{d_L(z)}{L^n}\ge d_{\rm slow}^{\rm per}(z)\ge1/3024,\qquad
\liminf_L\frac{a_L(z)}{L^n}\ge a_{\rm fast}^{\rm per}(z)\ge n/4.
\tag{9}
\]
轴向核可支撑在低维面；这里只用其概率核和正积分，不假设全维密度。两个下界是**同一个 (6) 输入的两个全空间速率**，所以可取非负乘积的下极限并对时间用 Fatou：
\[
\begin{split}
\liminf_L\frac{\mathfrak P(f_L)}{L^n}
&\ge\int_I\sqrt{d_{\rm slow}^{\rm per}(z)a_{\rm fast}^{\rm per}(z)}\,dz\\
&\ge\frac{\sqrt n}{8\sqrt{3024}}.
\end{split}\tag{10}
\]
因此结合 (7)：
\[
\liminf_L\frac{\mathfrak P(f_L)}{W_L}
\ge\frac{\sqrt n}{16\sqrt{3024}}.
\tag{11}
\]
对每个固定 \(n\) 可以选择一个有限 \(L=L_n\)，使比值至少为 (11) 的一半。由于 \(32^2\cdot3024<2048^2\)，得到具体真实有限 \(L^1\) 输入
\[
\boxed{\mathfrak P(f_{L_n})\ge \frac{\sqrt n}{2048}W_{L_n}.}
\tag{12}
\]
这排除任意固定幂 polylog \(W\) 的全局时间配对合同。不是用改变输入的极限偷去完整来源：每个 \(n\) 的反证见证都是单个有限 \(f_{L_n}\)，质量就是 (7)。

## 5. 初始已付点也不能修好这种全空间配对

若先删初始集合并只对 \(S=E_{\rm new}=\{f\le\lambda/2,M>\lambda\}\) 做 winner test，仍可使用同一双盒机制。仅将 fast 来源换成已审的 \(1+\cos(K\sum_i x_i)\)。其 gap test 在 \(I\) 上已有
\[
a_{{\rm fast},S}^{\rm per}(z)\ge n/12
\]
（每轴 test Dirichlet \(>1/3\)，scaled factor \(\lambda/2=1/4\)）。slow 盒初始高度至少 \(3/4>\lambda/2\)，不属于 bulk 的 \(S\)，但其真实熵仍计入全空间 \(d(z)\)。重复正积分提升给
\[
\liminf_L \mathfrak P_S(f_L)/W_L
\ge \frac{\sqrt n}{8\sqrt{36288}}.
\tag{13}
\]
取半极限并用 \(16^2\cdot36288<4096^2\)，有限 \(L^1\) 见证满足 \(\mathfrak P_S(f_L)\ge\sqrt n\,W_L/4096\)。此 gap 版本上界 \(f_L\le2\)，(8) 的常数改成 2，不改变提升论证。它没有重新计算 fast 单 mode。

## 6. 真正失败的交换与仍可用的局部接口

对于空间不交的理想极限，真实 signed edge commutator是两个本地 commutator之和；跨组件边在远移极限中消失。但
\[
\sqrt{(d_1+d_2)(a_1+a_2)}
\]
额外加入了 \(d_{\rm slow}a_{\rm fast}\) 的人工交叉。式 (12) 仅证明这种已做全空间 Cauchy 的时间配对太强。它不证明 \(|C|\) 本身为 \(\Omega(\sqrt nW)\)，也不否定
\[
\sum_P\sqrt{D_PA_P}
\]
在保留空间／坐标／jump 局部划分后可能有更好的共同来源付款。尤其不能把本例的两个周期相位谱当成原 actual \(K\) 标签或完整 FIRST/history。

本例在原有序核、同一输入与真实 winner 上成立；它没有 cube hard winner、full future cap、LCA／CP／GP、失败币等全 actual 资格。因此不作为 \(R_{\rm angle}\) 或原 geom 反例，也不新增主账费用。

## 7. 新守卫登记与范围

paired_entropy_spatial_mismatch_registration_20261007.json 在执行前登记三轮 \(n=8,32,128\)。仅运行同前缀新 guard：用正 artanh 级数给 \(\log2\) 的有理包围、正指数级数给 \(e^z\) 的有理包围，核 (4) 的有理差身份、(5) 下界、(10)–(13) 平方常数链。来源和时间参数预先固定，不拟合，不重跑 fast mode。

这些 Fraction 证书认证新解析常数；它们不模拟双盒原 jump/history，也不单独认证 finite-\(L\) 几何提升。该提升由 §3–4 的解析论证承担。

新 guard 已运行一次并终态 exit0。paired_entropy_spatial_mismatch_results_20261007.json 保存三轮各 300 项、共 **900 项 exact_PASS**：\(\log2\) 和 \(e^z\) 正级数有理包围，四个预定时间节点的 slow 熵下界、逐有理相位值的代数差身份，以及 full/gap 配对平方常数链。相位取值仅验代数身份，不把有限点当周期积分；周期 \(\langle\cos^2\rangle=1/2\) 由解析计算。没有 Monte Carlo、浮点拟合、原 fast 核重跑、live handle 或待完成采样。
