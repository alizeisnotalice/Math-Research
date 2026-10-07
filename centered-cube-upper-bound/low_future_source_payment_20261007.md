# 低平均未来响应之后的逐原来源支付

2026-10-07，root。证明状态：下述正支配已证；剩余空间预算未证。使用 A02 外放松与误差预算 skill 的原共同选择/正支配/收费检查，不引用外部优化定理。

## 1. 保留原来源的实际子后验

保持 `future_softness_moment_budget_20261007.md` 的所有原合同，特别是原 FIRST 的 σ(x),L_s(x),q，不重新选择。写

\[
 p_x(y)=P_{\sigma(x),L_s(x)}(x-y),\quad
 k_x(y)=\overline P_{\sigma(x),L_s(x)}^{(\alpha)}(x-y),\quad
 \rho_x(dy)=p_x(y)d\mu(y)/q,\quad Z_x(y)=k_x(y)/p_x(y).
\]

剩余 σ>1/n，所以原 φ 的严格正性给 p_x(y)>0。ρ_x 是概率，EρZ=Fbar/q。y 是原软源，不是硬源 z、首落点 w 或终点 V。

实际历史子后验 Π_x 的原 y 边际由 ρ_x 支配。这比只知道 Π_x 的总质量≤1强：对任意原来源可测集合 B，将原输入限制为 μ|B，真实 first-exit/early 正分解给其密度≤Pσ,L*(μ|B)/q；由任意 B 得测度支配。原 sourcepair 门依赖 y,z,历史不影响此结论，因为门介于0和1。必须在扩大到固定 b 的共同参考之前使用；不能从参考总质量≤2W反推此逐行原后验支配。

令 r_x 为硬赢家的归一化来源概率；实际硬子源由 r_x 支配。当前剩余交通的任何非负子支可写为或被下式支配：
\[
 \int_E u(x)\int \!\int g_x(z,\zeta)\,r_x(dz)\Pi_x(d\zeta)\,dx,\quad 0\le g\le1.
\]
仅需正支配，不假定软硬源在带门后独立。

## 2. 高 Z 来源支的空间支付

对任意预先固定 τ>0，在上式中添加 `Z_x(y(ζ))≥τ`。逐行删除 g、积分硬子概率，再用原 y 边际支配：
\[
 R_{Z\ge\tau}\le\int_E\frac{u(x)}q\int_{Z_x(y)\ge\tau}p_x(y)d\mu(y)dx
 \le\frac{C_h}{\tau}\int_E\int k_x(y)d\mu(y)dx
 \le\frac{C_h K_{n,\alpha}}{\tau}W. \tag{1}
\]
最后一步是上一轮已经证明的全输入平均核列费，原完整 μ 只计一次。Tonelli 适用于全部非负量；比值及原选择器可测，分支也可测。等号 Z=τ 计入已付支，未付支为严格 `<τ`。

取 α=ceil√n，τ=ε_n=log^-4(n+2)，得到 `O(√n log^5(n+2))W`，不消耗任何 λ|E| 吸收。即使输出本来满足 EρZ<ε，也可能有少量来源 Z≥ε；所以输出筛选和逐源筛选不是同一事件。二者依次取互补，费用至多两次 C_hK/ε，仍可并入既有 log^6 费用。

本引理本身不需要低平均输出条件；此处在上一轮 low-average 补支上应用，为下一步同时保留两项约束。没有用 Markov 舍弃而增加吸收项；此前考虑的此类吸收预算不进入最终账。

## 3. 与已付工具连接

对每一份未付原软来源，已有
\[
 Z_x(y)<\epsilon_n,\qquad \overline F_\alpha(x)/q<\epsilon_n.
\]
记 a_i=φ((x_i−y_i)/L_s)/ψσ((x_i−y_i)/L_s)>0，D=−∑log a_i。加权 AM–GM 与 Jensen 给
\[
 Z=\alpha\int_0^1v^{\alpha-1}\prod_i(1-v+va_i)dv
 \ge\alpha\int_0^1v^{\alpha-1}e^{-vD}dv
 \ge e^{-\alpha D/(\alpha+1)}.
\]
因此 `D>(1+1/α)log(1/ε_n)`。此式与同伴的 replacement 谱推导交叉一致，不证明来源之间独立。

若需要再支付 D≤B_n 的来源，则 p_x(y)≤e^{B_n}P_{1,L_s}(x-y)，原 y 边际的同一支配证明费用≤C_h e^{B_n}K_{n,soft}W，其中 K_{n,soft}≤1+(log(b/a)/2)√(10n)。例如 B_n=8loglog(n+2) 可付 √n log^8。此项是原全软列工具的有偿应用，非本轮新的空间定理；可选择使用，但本稿主账尚未将其叠加。

## 4. 更新后的精确缺口

原全部 FIRST/future cap/真实winner/出生/CPGP/LCA/first exit/full continuation，以及 nonconcentration、固定短壳、远离已选顶点等门不变；残余再逐软来源满足 Z<ε。此前吸收系数49/8192保持不变。

尚需证明这份实际残余具有 √n n^{o(1)} W 空间费用，或进一步可支付分解。上面的 D 下界及坐标捕获计数只是必要条件，不能当作空间费用。一般主定理仍未完成。本稿数值守卫仅检查有限正测度和相关门的支配代数，不认证实际 FIRST，也不拟合阶数。
