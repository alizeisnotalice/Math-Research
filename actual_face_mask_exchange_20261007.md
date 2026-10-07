# 实际 hard 面与 future 诊断 mask：后验、合法行支与未付项

2026-10-07。本轮只新增本稿及同前缀脚本/JSON，不改主账、不运行旧数据。结论：不能把先验 K/n 直接用于原捕获后的 mask；盒外坐标给原核上的确定性相关障碍。但对**不是原 soft 盒外坐标的 hard 面**，原核的 inside φ 权重可比较，给一个保留原标签的后验界 ≤8K/(3n+5K)。它能支付 unforced-hit 的固定低 K 行支，仍不覆盖 forced-hit、nonhit 或一般高 K 的 R_dagger。没有取得完整平方根 W 预算。

## 1. 查重与原诊断分配

已读 `actual_future_capacity_bridge_20261007.md` 的 P^A/P 诊断、κ_C 与已付高 K；`double_cutoff_actual_geometry_interface_20261007.md` 的 k_ε、O_s 及原盒支持；`short_shell_jump_face_audit_20261007.md` 的真实 first-jump 坐标、hard 面标签和后验匹配。以下 A 仍为**额外诊断**，不是原 early continuation 的活动坐标；不重证旧首跳均匀性障碍或全 mask peak envelope。

保留原完整 μ、共同 q∈[3λ/4,λ]、FIRST/future cap、原 L_s,R_h、hardband、原 soft y/hard z、共享森林和唯一 LCA、失败币、strict CP/GP、early 全历史、lowS、K 双 cutoff 和 adaptive-core 补集。对每一原实际 (x,y,z,θ,history) 标签，再以

\[
 \pi_x(A\mid y)=P^A_{\sigma,L_s}(x-y)/P_{\sigma,L_s}(x-y)
 \tag{1}
\]

作正分配。权重和一，原交通与完整 W 不改变。原 history 条件在分配前固定，不能把这个诊断 A 改为实际 history mask；P=0 的原交通为零。

对 0<σ<1，令 v_i=(x_i−y_i)/L_s、h_i=1_{|v_i|≤1/2}、φ_i=φ(v_i)。则真实来源条件后验是坐标 Bernoulli 的乘积，活动概率

\[
 p_i(x,y)=\frac{\sigma\phi_i}{(1-\sigma)h_i+\sigma\phi_i}.
 \tag{2}
\]

这是在固定原 y 下的诊断核分布；把 y 按完整 soft posterior 混合后一般不再独立。把来源对、K 截门或 history 再条件化也不能把 (2) 改成先验 σ。

当前已登记的新门对诊断 A 仅通过 K=|A| 使用；原其它门在原 history 上，来源核心门在原 z 上。因此固定全部原标签、再条件 K=k 与当前截门接受时，A 的角分布仍为 (1) 在 |A|=k 内的归一化。若未来新增依 A 角位置的门，这句话须重新审计；不能因总接受≤1而给接受后的条件分布继承本界。

## 2. fixed K 的准确带强迫坐标分布

设 O_s(x,y)={i:|x_i−y_i|>L_s/2}，o=|O_s|，I=[n]\\O_s，d=n−o。原 h 使用闭边界，|v_i|=1/2 属 I。φ 在全空间正，所以每个有正权的 mask 必含 O_s，K≥o。当 k≥o，写 ℓ=k−o；固定 K=k 后

\[
 A=O_s\cup B,\quad |B|=\ell,\qquad
 \Pr(B\mid K=k,x,y)=\frac{\prod_{i\in B}\phi_i}{e_\ell(\phi_I)}.
 \tag{3}
\]

e_ℓ 是 elementary symmetric polynomial。σ 系数、O_s 上共同的 φ 因子在条件化后恰消去，没有将来源重新归一化作为新的输入。σ=0 仅有 K=0；σ=1 仅有 A=[n]，确定端点另处理。

hard cone 面仍取原 z：J(x,z)=min argmax_i|x_i−z_i|。即使先验条件 K=k 的 mask 在坐标上均匀，若 J(x,z)∈O_s(x,y)，式 (3) 给

\[
 \Pr\{J(x,z)\in A\mid K=k,x,y,z,history\}=1.
 \tag{4}
\]

因此无条件 K/n 不是此实际后验的上界，且 K 小不消除该 forced-hit。它完全保留 y,z 的不同角色。

## 3. 原核/正 L¹ 来源的相关支持证据及范围

固定 L_s=R_h=1、σ=2/n。取 ε=1/(100n)、c_y=(1/(2n),…,1/(2n))，完整 μ 是以 c_y 和 c_z=c_y+(3/5)e_1 为中心、半边长 ε 的两份均匀正盒，各质量 1/2，W=1。接收盒中心 c_x=c_y+e_1，同半边长 ε。对这些盒中每个 x,y,z，原 soft 来源 y 在第一盒、hard 来源 z 在第二盒时：

- |x_1−y_1|≥1−2ε>1/2，其它坐标位移≤2ε<1/2，故 O_s={1}；
- 2/5−2ε≤|x_1−z_1|≤2/5+2ε<1/2，其余位移≤2ε，故原 hard 面唯一 J=1，且 z 严格在 hard cube 内；
- 原 P^A 核正性给 K=1 时唯一有正权的 A={1}。相交后验为1，而均匀 prior 给1/n。

这是同一完整正 L¹ 输入和原核上的正体积支持证据，不是任意抽象后验数组。第一 soft 盒位于原 origin-0 细格内。它**没有认证**该 μ 的完整 FIRST/future cap、实际 maximal 两赢家/硬带、far、strict CP/GP、coins、GOOD、birth/owntrace、全部历史或最新核心补集；不能称 R_dagger 反例。它只排除仅凭 prior 交换所声称的点态通用 K/n。若这些完整门能排除 forced-hit，仍须给出相应真实不等式，不能从诊断 K 小自行推断。

## 4. unforced hard 面的可证后验界

对于 j=J(x,z)∉O_s，原 inside 权重满足

\[
 m=3/8<\phi_i\le1=M\quad(i\in I).
 \tag{5}
\]

下界由原 φ(1/2)>3/8 和偶递减性给出（`short_shell_jump_face_audit_20261007.md` §7 已核原 tex 5800–5811）；上界因为 φ=h*w、w 是概率、h≤1。保留闭面不损失下界。

对任意正 inside 权重，记 E_r=e_r(φ_{I\\{j}})。ℓ≥1 时

\[
 \Pr(j\in A\mid K=k,x,y)
 =\frac{\phi_jE_{\ell-1}}{E_\ell+\phi_jE_{\ell-1}},
 \quad
 E_\ell\ge\frac{m(d-\ell)}\ell E_{\ell-1}.
 \tag{6}
\]

后一式来自每个 ℓ−1 子集补入一个未选坐标，再除每个 ℓ 子集被数 ℓ 次；不是独立性猜测。因此

\[
 \boxed{\Pr(j\in A\mid K=k,x,y)
 \le\frac{8(k-o)}{3(n-o)+5(k-o)}
 \le\frac{8k}{3n+5k}\le\frac{8k}{3n}.}
 \tag{7}
\]

ℓ=0 时概率0；ℓ=d 时为1；以上包含端点。第二个不等式用 (k−o)/(n−o)≤k/n。式 (7) 对固定原来源对和 history 成立，故可在原完整输入上积分，无需把 hard z 替成 soft y，也无需新的 W。

## 5. 能支付的条件行支，不能支付的互补部分

取输出无关、预先固定整数 M_0≤n，只在当前 R_dagger 内切

\[
 \mathcal H_{\rm unforced,low}=
 \{J(x,z)\notin O_s(x,y),\ J(x,z)\in A,\ K\le M_0\}.
 \tag{8}
\]

保留所有原标签及 K 双 cutoff。固定 K 后用 (7)，再对原尚未加入本次 angular hit 的低 K 基础交通积分；其完整归一化行质量≤1，hard 幅度≤4λ。因此

\[
 R_{\rm unforced\ hit,\ K\le M_0}
 \le\frac{8M_0}{3n+5M_0}\,R_{\rm current,\ K\le M_0}
 \le\boxed{\frac{32M_0}{3n+5M_0}\lambda|E|}.
 \tag{9}
\]

这是一个合法可吸收的**行子支**，费用没有 J 或 atom count；在依次互补账中才能登记一次。若给定吸收余量 ε_abs，只在明确 32M_0/(3n+5M_0)≤ε_abs 时使用。它不是整个 R_dagger 的 √n W 预算，也不在本稿擅改主账或增用已有吸收额度。

完整剩余至少有 forced-hit（式4）、所有 nonhit（J∉A）和未被其它已付 K 门删除的 K>M_0。κ_C=nσ+√(2nσ log(1/p_C))+2log(1/p_C)/3 可超过 n；k_ε 也未在所有 current 行被证明为 O(√n)。因此当前双 cutoff 不供应整份 (9) 所需的 √n 计数上限。就算额外证明低 K，nonhit prior 往往接近1；不能拿 hit 的小概率乘到整份 traffic。

nonhit 仅保留原 soft kernel 在该坐标的 h_{L_s}(x_j−y_j) 支持。它不证明 hard z 与 y 近，也不给重新冻结该坐标的 FIRST/future cap；旧真实 continuation 冻结审计已明确该转接失败。若面敏感度公式只依 hit 起作用，仍须从原未截门核推导该公式，并支付被门截去的差分负部；不能由 (9) 免费消除 holding/nonhit collar。

## 6. 三轮新守卫与最终范围

`actual_face_mask_exchange_20261007.py` 只运行本轮新守卫，n=512、4096、32768：原两盒共同正 L¹ 来源/接收盒的精确支持余量；forced K=1 支持及1对1/n相关；inside φ∈[3/8,1] 的 elementary symmetric marginal 恒等式和 (7) 的有理系数界。权重测试数组只核通用 interval lemma，不声称任意选定数组就是实际 φ 空间值。全部几何和系数决定用 Fraction，无 Monte Carlo、kernel积分误差或阶数拟合。

运行已结束，结果为 `actual_face_mask_exchange_20261007_results.json`：三轮各27份 marginal 记录，共81份；每份核边缘界和独立 partition 恒等式，另每维5项来源盒支持条件，共177项精确检查全部通过。实际运行2.75秒；三个 positive-source 支持证据均为 forced-hit=1、prior=1/n，未填写全 actual 资格为通过。没有 running handle 或待完成采样。

本轮没有新的全覆盖阶数/下界猜测，所以不另启动下界族搜索，不重跑旧 mask、FIRST、shell 或 packet 实验。新的有效接口止于 (7) 的真实后验比较与 (9) 的条件行支；原 forced-hit 和 nonhit 的 source-once 空间费、或其完整门排除/有偿面差分转接仍欠。一般 √n n^{o(1)} 目标与主账 R_dagger 保持未决。
