# 实际面／诊断 mask 交换及主账吸收：独立审计

2026-10-07，gated_radial_energy_audit。只读审计 actual_face_mask_exchange_20261007.md，并核对最新原提示词的完整行表示及 current_joint_budget_20261007.md 当前账。不修改这些文件，不重跑已结束的 177 项守卫。结论：在下述单一诊断分配和角向门合同下，式 (3)、后验界及新的互补行支可合法登记；主账的 \(147/143\) 因子正确。一般 geom 目标仍有明确剩余。

## 1. 固定原标签后进行一次诊断分配

完整原行为

\[
p_x(dy)=q^{-1}P_{\sigma,L_s}(x-y)\mu(dy),\quad
r_x(dz)=u(x)^{-1}h_{R_h}(x-z)\mu(dz),\quad u(x)\le4\lambda.
\tag{1}
\]

两者是完整输入上的概率。实际子交通的原 FIRST、futurecap、原 soft \(y\)／hard \(z\)、出生、种子、唯一 LCA、失败币、strict CP/GP、实际 early/full continuation 及此前所有删除门，先保留在原接受权中。原提示词明确原实际交通或其正支配可由 \(0\le G\le1\) 表示；这里没有调用参考放大的 \(C_h\) 行。

在固定 \(x,y,z,\vartheta,\mathrm{history}\) 上，以

\[
\pi_x(A\mid y)=\frac{P^A_{\sigma,L_s}(x-y)}
{P_{\sigma,L_s}(x-y)}
\tag{2}
\]

装饰一个额外诊断标签。该抽样是**条件于原 \(x,y\)、独立于原 history 的正分配**；它不是自然路径的 continuation mask。原 \(P=0\) 的来源没有原 soft 交通，任意补全分配均无影响。

当前已经登记的 \(K\) 双 cutoff 与 \(w_K\) 门使用同一个诊断 \(A\)，不可重新抽样后悄悄把旧 \(K\) 当新 \(K\)。独立性只指分配前与原 history 的条件独立，不指 \(A\) 与 \(y\)、\(z\)、\(K\) 或原接收行无条件独立。

原 history 及其全部门在 (2) 分配前冻结。当前新增的诊断接受权仅依赖 \(K=|A|\) 及固定原标签；强 \(S\)、来源核心、父组、币等门并不读取 \(A\) 的角位置。这一条件足以在固定 \(K=k\) 后保持 (2) 的角分布。若后来增加依 \(A\) 坐标位置的门，则接受后的角分布可能偏置，不能自动继承本审计的后验概率界；仍可在新增角门之前使用未归一化正支配。

## 2. 式 (3) 与 unforced 面后验界

当 \(0<\sigma<1\)，设

\[
O_s=\{i:|x_i-y_i|>L_s/2\},\quad o=|O_s|,\quad
I=[n]\setminus O_s,\quad d=n-o.
\]

原硬因子的闭边界 \(|x_i-y_i|=L_s/2\) 属 \(I\)。逐坐标原 mixture 分解给权重

\[
\sigma^{|A|}(1-\sigma)^{n-|A|}
\prod_{i\in A}\phi_i\prod_{i\notin A}h_i.
\]

故正权 mask 必含 \(O_s\)。给定 \(K=k\ge o\)，\(\ell=k-o\) 时 \(A=O_s\cup B\)，且

\[
\Pr(B\mid K=k,x,y)
=\frac{\prod_{i\in B}\phi_i}{e_\ell(\phi_I)},\qquad
B\subset I,\quad |B|=\ell.
\tag{3}
\]

原 \(\sigma\) 因子以及 \(O_s\) 上公共的 \(\phi\) 因子完全消去。因此原稿 (3) 正确，且不要求不同来源间独立。把完整 soft posterior 对 \(y\) 混合后通常失去坐标独立，不能由此退回先验 \(\sigma\) 或 \(k/n\)。

令 \(j=J(x,z)\) 为按固定 tie rule 定义的原 hard 面。若 \(j\in O_s\)，则后验 hit 概率恰为 1；本支不支付它。若 \(j\notin O_s\)，inside 原核满足 \(m=3/8<\phi_i\le1=M\)。对 \(1\le\ell\le d\)，删去 \(j\) 的 elementary symmetric 系数满足

\[
\ell E_\ell\ge m(d-\ell)E_{\ell-1}.
\]

这是逐 \(\ell-1\) 子集补坐标的有限正计数。于是

\[
\Pr(j\in A\mid K=k,x,y,z,\mathrm{history})
\le\frac{8(k-o)}{3(n-o)+5(k-o)}
\le\frac{8k}{3n+5k}.
\tag{4}
\]

第二个比较等价于 \(n(k-o)\le k(n-o)\)，即 \(o(k-n)\le0\)。当 \(\ell=0\) 概率为 0；当 \(\ell=d\) 为 1。若 \(d=0\)，不存在 unforced \(j\)，所以不取 \(0/0\)。在 \(\sigma=0\) 端点仅 \(A=\varnothing\) 有正权、hit 空；\(\sigma=1\) 仅 \(A=[n]\)，本次 \(M_0<n\) 的低 \(K\) 支空。原实际 \(\sigma>1/n\) 的范围无需修改。

## 3. 当前余项中的互补行支

仅在当前 \(R_\dagger\) 的实际接受交通内定义

\[
\mathcal H=\{j\notin O_s,\ j\in A,\ K\le M_0\}.
\tag{5}
\]

在未加本次角向 hit 门的基础交通中，先固定全部原标签并给定 \(K=k\)。此前诊断门只依赖 \(k\)，因而 (4) 可用于未归一化交通。再按原完整概率 \(p_xr_x\) 积分，有

\[
R_{\mathcal H}
\le\frac{8M_0}{3n+5M_0}R_{\rm base,low}
\le\frac{32M_0}{3n+5M_0}\lambda|E|.
\tag{6}
\]

最后一界仅使用完整行质量 \(\le1\) 与 \(u\le4\lambda\)。没有额外 \(C_h\)、森林深度或同顶组二倍。若将实际 history 先积分，早期子核比值已由 \(Q/P\le1\) 纳入原权，不能再乘一次该比值。

令 \(R_{\rm new}\) 为 \(R_\dagger\) 在 \(\mathcal H^c\) 上的原正子交通，便有精确分割

\[
R_\dagger=R_{\mathcal H}+R_{\rm new}.
\tag{7}
\]

仅 unforced-hit 且 \(K\le M_0\) 被新支删除。forced-hit、所有 nonhit，以及原双 cutoff 仍允许的 \(K>M_0\) 全部保留。不能将 (6) 的系数乘到整份 \(R_\dagger\)；更不能删除此前已付门再重新清算这些历史交通。

## 4. 主账常数及端点

取预先固定的

\[
M_0=\lfloor n/65536\rfloor,\qquad
\alpha_n=\frac{32M_0}{3n+5M_0},\qquad C_0=\frac{8192}{49}.
\tag{8}
\]

若 \(M_0>0\)，由 \(M_0/n\le1/65536\) 及单调性，

\[
\alpha_n\le\frac{32}{196613}<\frac1{6144},
\qquad C_0\alpha_n\le\frac4{147}<1.
\tag{9}
\]

\(n<65536\) 时 \(M_0=0\)，(5) 因 \(K=0\) 不可能 hit 而为空；可仍写统一的宽上界，但没有实质新删除。含等号 \(K=M_0\) 属新 paid 子支，补集不是无条件 \(K>M_0\)，因为 forced/nonhit 仍允许低 \(K\)。

将当前既有费用统一记为 \(F_nW\)，根给出的账是

\[
\lambda|E|\le4\mathcal B_{\rm paid}
+C_0F_nW+C_0R_\dagger.
\]

代入 (6)–(7)，得到较精确的版本

\[
\lambda|E|
\le\frac1{1-C_0\alpha_n}
\left[4\mathcal B_{\rm paid}+C_0F_nW+C_0R_{\rm new}\right],
\]

以及统一可登记的

\[
\boxed{\lambda|E|\le
\frac{147}{143}
\left[4\mathcal B_{\rm paid}
+\frac{8192}{49}F_nW
+\frac{8192}{49}R_{\rm new}\right].}
\tag{10}
\]

该因子放大整个原右端一次，没有新增一份来源 \(W\) 收费。它不改写上游三尺度证明依赖，不将当前系数改回旧 HR/MM 的 \(128/63\)，也不重用旧吸收额度。因此可以作为**真正的当前主账互补更新**，但仅更新余项资格和有限常数，未证明 \(R_{\rm new}\) 的来源一次空间预算。

## 5. 收据范围

未运行或重跑数值。原稿 177 项已终态的 Fraction 守卫验证其有限正分解、inside 权区间 lemma 和来源盒支持；其中正体积 forced-hit 模型未认证完整实际 FIRST/history。式 (9)–(10) 的主账常数在本审计中直接精确推导，不依赖那些有限模型成为 actual 样本。

可执行的新补集仍是全部原门、原两赢家与来源、原 actual history、lowS、lowcoin、诊断双 cutoff、adaptive-core 补集，加上 \(\neg\mathcal H\)。其 forced-hit/nonhit 的空间清算仍欠。该更新有效，但没有取得一般 \(\sqrt n\,n^{o(1)}\) 完整结论。
