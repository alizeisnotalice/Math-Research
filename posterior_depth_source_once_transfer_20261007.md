# 实际捕获后验的径向深度与一次来源删除预算

2026-10-07。根节点自含推导。使用 H05 的尾界、矩母函数及依赖结构审计流程；已读 SKILL、method、provenance、cube-interface，不引用任何未核外部定理。本稿的概率是原完整捕获后验，不假定来源坐标独立。主目标仍未证明。

## 1. 原连续赢家强制的随机变量序

令完整有限正 Borel 来源为 μ，W=μ(R^n)>0；Q(x,s) 为全边长 s 的闭立方体，0<a≤s≤b。设 M(x)=sup_{a≤s≤b} μ(Q(x,s))/s^n。给定可测候选 L(x)∈[a,b] 且 m_L(x)=μ(Q(x,L))>0，写

\[
g(x)=\log\frac{M(x)L(x)^n}{m_L(x)}\ge0,\qquad
\pi_{L,x}(dy)=\frac{\mathbf1_{Q(x,L)}(y)\mu(dy)}{m_L(x)}.
\]

对该后验中的 Y 定义

\[
A(x,Y)=n\log\frac{L(x)}{\max(a,2\|x-Y\|_\infty)},\qquad
B(x)=n\log(L(x)/a).
\]

于是 0≤A≤B。对 0≤t≤B，尺度 s=L exp(−t/n) 合法，包含边界原子在内有准确身份及上界

\[
\pi_{L,x}(A\ge t)=\frac{\mu(Q(x,s))}{m_L(x)}
\le\min\{1,e^{g(x)-t}\}.                 \tag{1}
\]

t>B 时左端为零。(1) 在 t=0、t=B 也成立；当 B=0 时 A=0，不能声称 A 无原子或有连续指数分布。因而 A 被 min(B,g+Z) 一阶随机支配，其中 Z 是参数 1 的辅助指数变量；只作分布比较，不把 Z 当实际独立的来源坐标。

由非负 layer cake，所有 p>0、0<σ<1 满足

\[
\mathbb E_{\pi_L}A\le g+1,\qquad
\mathbb E_{\pi_L}e^{\sigma A}\le\frac{e^{\sigma g}}{1-\sigma},
\quad
\mathbb E_{\pi_L}(A-g)_+^p\le\Gamma(p+1).       \tag{2}
\]

最后一个 Γ 记号只表示显式积分 p∫_0^∞ t^{p−1}e^(−t)dt；整数 p 时即 p!。这些式不要求多个 receiver 或多个样本独立。

精确原赢家 g=0 时，保留窗口截断还能得到

\[
\mathbb E A\le1-e^{-B},\qquad
\mathbb E e^{\sigma A}\le
\frac{1-\sigma e^{-(1-\sigma)B}}{1-\sigma}.       \tag{3}
\]

证明分别把 (1) 积分于 [0,B]，使用 E e^(σA)=1+σ∫_0^B e^(σt)P(A>t)dt。闭端点的质量不改变这个积分。普通 finite-scale winner 不满足 (1)，除非另证其 maximum 等于整个连续窗口的 maximum。有限原子来源可通过全部真实 arrival 加 a,b 的精确最大值获得这种资格。

## 2. 原弱型 joint 下按来源只计一次

现在 L=R 为原精确赢家。设 E={M>τ}，I=τ|E|，并定义非归一化真实 joint

\[
J(dx,dy)=\tau\mathbf1_E(x)dx\,\pi_{R,x}(dy).
\]

其总质量 I，来源边缘为 S(y)μ(dy)，其中

\[
S(y)=\tau\int_E\frac{\mathbf1_{Q(x,R)}(y)}{\mu(Q(x,R))}\,dx.
\]

对 t>0 定义 S_{\rm deep,t}(y) 为同一积分增加 1_{A(x,y)>t}。Tonelli 和逐 receiver 的 (1) 严格给

\[
0\le S_{\rm deep,t}\le S,\qquad
\int S_{\rm deep,t}\,d\mu\le e^{-t}I.           \tag{4}
\]

深度 A≤t 的保留部分至少有 (1−e^(−t))I 交通；被删部分没有按尺度、原子数、方向或历史再领取 W。该深度条件等价于来源在原赢家的厚度约 t/n 外壳内，另含 R接近a时的整个floor区；不能删掉 max(a,D) 或把floor区当空。

若一个非负实际贡献在同一 joint 上有密度 q(x,y) 且已独立证明 0≤q≤C，则其深部部分 ≤C e^(−t) I。更一般，Hölder 在同一个有限 joint 上给

\[
\int q\mathbf1_{A>t}\,dJ
\le \left(\int q^p dJ\right)^{1/p}
       (e^{-t}I)^{1-1/p},\quad p>1.             \tag{5}
\]

这只是一份条件迁移接口。尚未核验原 R_angle 的 q 有这样的统一 C 或 Lp预算，原软核可在硬捕获之外带质量，原历史还可能偏置后验；因此 (4) 本身不能给 R_angle 标记 paid。

## 3. 不是一般上界的原因

(1)–(4) 控制每个原 receiver 的已捕获来源在该 cube 内有多深。它们未限制有多少 receiver 选择相似外壳，未控制 source column S 的高度或外壳的全空间重叠。单原子来源的赢家恰在捕获面上，A=0，(1)–(4) 完全合法，却不给面 trace 的密度上界。这个例子只说明接口范围，不是目标弱型上界的反例。

既有包络 joint 的 KL 恒等式中 d=A，故可加入 E_P d≤1；这不把宽径向参考 law 改成集中 law，也不自动产生角 score 对 I/W 的下界。独立审计另见 posterior_radial_age_entropy_20261007.md。

数值压力由独立注册的 posterior_radial_age_probe_20261007 系列负责。本稿不以有限样本推导次数，也不把 receiver 离散检测当作 Lebesgue 积分证明。
