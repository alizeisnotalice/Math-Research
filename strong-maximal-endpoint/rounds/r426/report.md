# 第426轮　最小捕获集的逐盒门槛传递

- A 已覆盖全部等权深度轮廓。对任意实际接收分区，平均命中约束可传递为一族严格 halo 约束；稳定门槛处的体积损失至多为平方根乘多对数因子。
- 若最小捕获集产生的行平均具有非消失方差，则不降低门槛即可保住正比例接收体积。
- 仍缺近乎常值但在门槛两侧极不平衡时的无损传递。一般系数差仍为 $(\mu-1)/8$。

2026-11-13。沿用第425轮 $A=\sum |U_i|$、$D$ 的体积 $v_q$。按概率 $|U_i|/A$ 抽行，令 $X=|D\cap R_i|/|R_i|$。由排序构造 $\mathbb E X=q$。对任意 $0<t<1$，$A\Pr(X>t)\le v_q H_n(t)$，因为每个成功行的原负责区域都被 $D$ 的严格 halo 覆盖。因此
\[
 q\le\int_0^1\min(1,v_q H_n(t)/A)dt,\qquad
 v_q\ge\frac{A\delta}{(1-q+\delta)H_n(q-\delta)}\quad(0<\delta<q).
\]
第二式用 $q\le q-\delta+(1-q+\delta)\Pr(X>q-\delta)$。特别在第349轮未分带的稳定近最优目录，取其门槛为 $q$，$A\ge(1-1/n)mH_n(q)$，$H_n(q-\delta)/H_n(q)=1+O_I(1/\log n)$，$\delta=1/(\sqrt n\log^2 n)$，得到 $v_q\gtrsim_I m/(\sqrt n\log^2 n)$。此处不能直接把未分带的稳定性移至第421轮移动后的带门槛。

在同一门槛，$\operatorname{Var}(X)\le\mathbb E|X-q|=2\mathbb E(X-q)_+\le2(1-q)\Pr(X>q)$。故 $v_q\ge A\operatorname{Var}(X)/[2(1-q)H_n(q)]$。更强地，对于任意 $\epsilon>0$，若 $v_q H_n(q)/A\le\epsilon$，则 $\mathbb E|X-q|\le2(1-q)\epsilon$，即平均约束到逐盒约束的缺陷必表现为绝大部分接收体积上的行密度集中。有限目录下捕获集可由单元分割实现；即使使用一般可测捕获集，以上比较亦适用于指标 halo 定义。下一步应利用门槛等值单元内可自由选取部分质量的性质，检验能否同时改善多行严格门槛，而不增大 $v_q$。
