# r089 一般输入 cap 相变独立审计

审计结论：在 r073、r021、r022、r011、r028 与 r013/r014 所引用的已验收定理范围内，`core89.md` 的三段 cap 结论闭合；未发现把标量构造误作一般 physical upper、遗漏实际 Orlicz 费用或临界误差大于主项的漏洞。它仍不确定 unrestricted `C_N` 的未知 `log N` 系数。

## 1. cap-two clipping 与 (2\sqrt\mu) 上阈值

`f=F wedge a`、`g=F-f` 有 pointwise exact identity
\[
J(F)=J(f)+\Delta J,
\]
即使两者原函数支撑重合，也没有使用错误的 Orlicz additivity。`g` 的 `L^p` 费用来自 `Phi(a+u)-Phi(a)`，而非丢弃原高度；`(log Phi)'>=cN` 给
\[
\|g\|_p^p\le\frac{(Cp/N)^p}{\Phi_N(a)}\Delta J.
\]
与 strong `L^p` 估计结合后，低/high 两个 halo 的费用是 `A J(f)+H Delta J`，因此可由 `max(A,H)J(F)` 支付。这是使用实际 modular 的合法 max ledger。

取 `p=b sqrt(N/log N)` 后，式 (3) 中的根号项在
\[
b=\frac{\sigma}{\sqrt{2P+1}},\qquad
K=\sigma(1+\sqrt{2P+1})
\]
恰为 `-sigma`。与 cap-two lower benchmark 相比，`z<sqrt(log N)` 产生的负 `log log N` 修正及 `(p/N)^p` 次级项使 high term 严格更小，故 (1) 没有循环。

在 `K=2sigma`，`s=N^{-1}(log N)^{-1/4}` 使 low dilation 是 `1+O((log N)^{-1/4})`；同一 high-term 计算留下 `-(sigma/4)sqrt(N/log N)log log N`，仍低于 benchmark。临界常数处理正确。`delta/rho -> infinity` 的 (3a) 也正确：`p=b/delta`、`s=e^{-cNdelta^2}` 在 `Ndelta^2>>log N` 时满足所有 `L^p` 条件，且 bracket 的严格负性来自 `c<1/(2mu)`；fixed-small-delta 版本须保留原稿列出的“小端点先选定”限制。

## 2. (\sigma\le K<2\sigma) 的一般 physical upper/lower

r021 的全预测 CDF 是单侧 stochastic inequality。对有界非减 test，Stieltjes/layer-cake 确实把
`F_pi>=F_G'-2N^-4` 转为 `E_pi V<=E V(G')+2N^-4 osc(V)`；不需要 policy 独立。这里必须使用整个连续且有界的 capped envelope
\[
1_{x\ge s}+N^{-1}\min(N\epsilon,(x-s)_+)+e^{x-s}1_{x<s},
\]
原稿正是如此处理，未孤立低侧指数项。

由单调 modular budget 得 `Phi(h(x))<=e^{x-s}`；低侧、inverse-Phi ramp 与真实 cap 能拼成上述 envelope。`G'` 的低侧指数期望为 `O(L/N)`，normal mean/variance shift 的 cdf 误差为 `O(L^2/N)`，ramp coupling 同量级，故 (5) 对整个 policy/profile 类一致。该 remainder 小于所有固定 `K<2sigma` 的 root signal `f_0D_K`：`K=sigma` 时为 order `N^{-1/2}`，`sigma<K<2sigma` 时仍多项式地大于 `L^2/N`。

Gaussian mean identity使用
\(mathbb E(a-G)_+=\sigma\sqrt N[\phi(y)+y\Phi(y)]\)。在上端 standardized point `w`，这等于 `sigma sqrt(N)[w+Psi(w)]`，正好抵消 `epsilon-F_G(q_0)` 中的 `sigma w/sqrt N`；因此原稿的
\[
U(q_0)=1+\mathbb E(q_0-G)_+/N-\sigma\Psi(w)/\sqrt N
\]
没有漏项。`U'`、局部 density comparability 与 root bracketing均正确，给
\[
\log C_N(1+K\rho)=q_0-(1+o(1))D_K.
\]
真正 lower 使用 `H_s` 本身、fixed iid `r_*` increments、绝对 Edgeworth remainder和
`j(H_s)<= (1+N epsilon)e^{-s}`；`O(L)` modular overhead 被 `D_K>>L` 吸收。因此这不是只对 Gaussian comparison 的数值结论。

## 3. 更小 caps

对 `epsilon=N^{-1}`、`(N log N)^{-1/2}` 和固定 `K<sigma`，upper 的两层 budget 使低层 CDF 必须小于 `epsilon`，再用 r021 predictive CDF 给
\[
\log C_N(1+\epsilon)\le q_\epsilon+O(L).
\]
lower two-level profile 的真实 modular 至多 `2e^{-s}`，严格可行条件恰为
`F_*(s)<epsilon[1-F_*(s+ell)]`。`ell<=N epsilon`，r013/r014 fixed-iid relative moderate deviations 的适用范围覆盖这里的 `z=O(sqrt L)`；所选 `C_1L` decrement 压过 cap-ramp 项。因此 (10)--(12) 的量词和实际 `Phi_N` 费用正确。特别是 (12) 仅固定 `K<sigma`，原稿没有错误地宣称它在 `K up to sigma` 一致。

## 4. 必须保留的依赖与范围

本稿的 general upper 依赖 r021 的已验收 all-predictable-policy CDF inequality和 r011 的 physical-to-feedback reduction；若这些历史定理被削弱为固定 vector 或仅 scalar model，§2 将不成立。r022 在这里用于真实 lower profile/known fixed-iid Edgeworth，不能替代 r021 的 general direction。当前引用的作用域匹配。

所有相变陈述均为 fixed `K`/`P` 的渐近；不包含随 N 收缩的 K transition window、critical polynomial exact rate或 unrestricted `log N` coefficient。r089 正确没有作这些额外断言。
