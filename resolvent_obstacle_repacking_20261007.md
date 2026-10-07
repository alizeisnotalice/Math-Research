# 同一障碍的正预解重装：保留域外符号和封顶，一次来源费用至多翻倍

2026-10-07。一般原核代数工具；未证明 ordered 弱端点或 geom 余项。使用用户 J03 非局部退出工具的域外数据要求，所有域外条件均在整个 Ωc 上，不把退出落点压到几何边界。

## 1. 假设与正预解

沿用 `ordered_frozen_spectral_difference_20261007.md` 的原固定入口：Gi 为原坐标保质量对称正算子，S=Σ(I−Gi)，Hq=exp(−qS)，Kr=∏(I−r(I−Gi))。原障碍 u≥0 属于 L1∩L2，Ω={u>0}，σ=Su=νb−μb，0≤μb≤κ，μb=κ 于 Ω；νb≥0 支撑 Ω，∫νb=∫μb=Wb，κ|Ω|≤Wb。这些假设来自旧固定障碍构造。

任意 t>0 定义

\[
R_t=(I+tS)^{-1}=\int_0^\infty e^{-a}H_{ta}\,da,\qquad
h_t=R_t\sigma=\frac{u-R_tu}{t}.
\]

Bochner 积分在 L1 和 L2 上合法，Rt 正、保质量且是收缩。该恒等式来自 (I+tS)Rt=I，S 有界；没有用未证的联合乘子 L1 定理。

由 σ≥−κ（全空间），得 h_t≥−κ；可把保常数正算子作用于 σ+κ，或在 L∞ 上处理下界，不能说常数属于 Rn 的 L1。在 Ωc，u=0，所以 h_t=−Rtu/t≤0。另有

\[
\int h_t=0,\quad \|h_t\|_1\le2W_b,
\quad \int h_t^-\le W_b.
\]

这些是同一来源的完整时段界，不按 r 或 t 逐份重新收费。

## 2. 新封顶障碍，仍是同一个 u 和 Ω

定义有界 Markov 生成元

\[
B_t=(I-R_t)/t,
\quad \nu'_t=\mathbf1_\Omega(\kappa+h_t),
\quad \mu'_t=\kappa\mathbf1_\Omega+\mathbf1_{\Omega^c}R_tu/t.
\]

逐点检查：ν'≥0，支撑 Ω；μ'≥0，且 Ω 上 μ'=κ，Ωc 上 μ'=−h_t≤κ。于是

\[
\nu'_t-B_tu=\mu'_t,\qquad 0\le\mu'_t\le\kappa,
\qquad \mu'_t=\kappa\quad\hbox{on }\{u>0\}.
\]

质量由 ∫h_t=0 给

\[
W'_t=\int\nu'_t=\int\mu'_t
=\kappa|\Omega|+\int_{\Omega^c}h_t^-
\le2W_b.
\]

故这是一个合法的新障碍数据组，保持原 u 和整个域外条件。它是从原障碍导出的重装，不宣称 ν'=νb，也不宣称 W'≤Wb。无需再次假设新障碍的存在性；上述显式式子已经给出一组解。Rt 的核一般包括 holding 与所有坐标面，因此不能称它为处处绝对连续的 full-dimensional 核。

## 3. 差算子的简单有界 L1 定义与 L2 工具

固定 t=1，记 h=h1。准确恒等式

\[
(K_r-H_r)\sigma=\mathcal N_rh,
\qquad \mathcal N_r=(K_r-H_r)(I+S).
\]

它在 L1 和 L2 均有定义，且单个 r 的粗界为 ||N_r||1→1≤2(1+2n)。这个粗界不提供极大弱型所需的 polylog n；不可遗漏外面的 sup。

已有谱差积分对 d_r=k_r−e^{-rs} 给

\[
\int_0^1(|d_r|^2+|r d'_r|^2)\,dr/r
\le40\min(1,s^{-2})\quad(s>0).
\]

乘以 (1+s)^2，并用 (1+s)^2 min(1,s^-2)≤4，可得与旧对数参数 Sobolev 论证相同的

\[
\left\|\sup_{0\le r\le1}|\mathcal N_rh|\right\|_2
\le\sqrt{160}\|h\|_2.
\]

s=0 时差乘子为0。轨道共同几乎处处连续，沿用有限 Kr 多项式与有界 S 指数级数的支配，不只诉诸强 L2 连续。这里的正预解来源避免了 S/(1−e^-S) 的额外 L1 函数演算问题，但 h 的正部仍未封顶，不能用 h≥−κ 和 ||h||1≤2Wb 推断 ||h||2²≤CκWb。

## 4. 数值守卫与尚缺接口

先登记 `resolvent_obstacle_repacking_registration_20261007.json`，随后一次执行同名前缀 guard。三轮 n=2,3,4，有限离散立方体 Gi=(I+坐标翻转)/2；单尖峰、双尖峰、全支撑三个 u；t=1/4,1,4，共27组，Fraction 精确消元，234项恒等式/符号/质量检查通过。结果保存在同前缀 results JSON，含注册及脚本 SHA256。

这些有限 Markov 模型只检查上述一般代数，明确不是原 Rn 卷积核实验、不是原 geom 模型、不是阶数拟合或端点证据。一般性结论由 §§1–3 的正算子证明承担。旧原核数值不重跑。

仍需估计 Ωc 上 sup_r[N_r h]+ 的水平集，利用 h=(I−R1)u 的障碍结构支付 polylog(n)Wb；或直接支付原 geom 联合余项。新封顶 μ' 和来源 W'≤2Wb 不自动给出联合 Kr 与新 killed 生成元的交换、Abel 泄漏或相同 FIRST 标签。未闭合，主费用没有减少。
