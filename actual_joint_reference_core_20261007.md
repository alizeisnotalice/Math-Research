# Actual 同门 reference 的规范审计：高标签 reference 完全抵消

2026-10-07。本稿只审查当前实际 \(R_{\rm angle}\) 的联合 reference/signed 核心；不改主账、不重跑任何旧守卫。使用 A02 的正支配/误差预算工作流，已读其 SKILL 与 provenance；只用其中要求的真实可行集映射、支配方向和误差清算，不将方法说明当新的数学定理。实际原件及记号沿 actual_signed_bridge_reassessment_20261007.md、原 tex 4875–4935/5522–5585、actual_gate_contract 和 current_joint_budget 最新节。

本轮结论是一个严格的路线排除：**同门全 reference 与高标签 signed 项联合后，高标签上的 reference 无论如何选都消失。** 优化其 Poisson 时钟、物理尺度或正混合本身不能给出高标签实际交通的新空间费用。没有用易尾支代替这一核心，也没有把未知当反例。

## 1. 在原 joint measure 上写身份

保留原完整输入 \(\mu=f\,dx\)、\(\sum_C\nu_C\le\mu\)、硬来源 \(\nu_h\le\mu_{\rm hi}\le\mu\)、全部 FIRST/fullfuture、出生、strict CP/GP、重捕获、唯一 distinct-child LCA、共享 seed/\(\theta\)、失败币，以及最新 lowS/lowcoin/双诊断 cutoff/source-core/角门补集。

先冻结原接收 selector，再在原 \((x,\ell,z,\theta)\) joint measure 上作接受密度 RN 分配。首跳标签
\[
 \ell=(C,y,t,i,\rho,w),\quad w=y+\rho e_i,\quad
 c_t=1-t,\quad p_x=(\sigma(x)-t)/(1-t),\quad L=L_s(x)
\]
仍在原 \(y,z,C\) 上判门；\(t>\sigma(x)\) 时所有接受函数置零。这里的 \(A\) 仅为原正 endpoint 分解标签，既非原物理路径访问集合，也非已经筛过的辅助诊断 \(K\)。特别是最新角门切分后，不重新独立抽样旧角后验。

记 \(h_L^{(n)}=L^{-n}\mathbf1_{[-L/2,L/2]^n}\)，
\[
 P_A(x;\ell)=\beta_A(p_x)
       (h_L^{(n)}*G_{c_t,L}^{A})(x-w),\qquad
 \beta_A(p)=p^{|A|}(1-p)^{n-|A|}.
\]
accepted 标签密度准确写为
\[
 F_A(x;\ell,z,\theta)=a_A(x;\ell,z,\theta)P_A(x;\ell),
 \qquad 0\le a_A\le1.                                     \tag{1}
\]
这里 \(a_A\) 是完整历史门相对于 joint endpoint measure 的条件密度，不是一个独立先验 gate。零 prior 分量的 accepted 密度为零，取 \(a_A=0\)；先固定 selector 后定义它，避开逐参数任意 RN 代表沿不可数 graph 取值的问题。

令 \(\mathfrak L_x\) 表示原首跳测度、原 hard factor \(h_{R_h(x)}(x-z)/q\)、\(\nu_h(dz)\) 与共享 \(\theta\) 的同一正积分。它在本稿所有项中保持不变。canonical actual 表示满足
\[
 T(x)=\mathfrak L_x\sum_A F_A.
\]
若只有保持门的正支配接口，则以上表示定义 \(T_{\rm rep}\ge T_{\rm act}\)。后续全部代数对该表示成立；扩大后的表示不自动有 \(4\lambda\) cap。只有原 \(T_{\rm act}\) 继承该 cap，不能将两者混写。

固定整数 \(m\ge1\)，分空标签 \(A=\varnothing\)、低标签 \(1\le |A|\le m\)、高标签 \(|A|>m\)。\(m\) 与原诊断 \(K\)、face \(M_0\) 完全区分。

## 2. 任意正 reference 的准确 gauge 公式

任取非负 endpoint reference 密度 \(R_A(x;\ell)\)，可依原 selector、首标签和其它冻结标签；只要求涉及的正积分有限。它不必是 Markov kernel。冻结同一个 \(a_A\)，定义
\[
 F_A^R=a_AR_A,\qquad D_A^R=F_A-F_A^R,\qquad
 B_R=\mathfrak L_x\sum_A F_A^R .
\]
则
\[
 B_R+\mathcal D_{\rm high}^R
 =T-\mathcal D_{\rm low}^R-\mathcal D_{\varnothing}^R
 =T_{\rm high}+B_{R,\rm low}+B_{R,\varnothing}.
                                                               \tag{2}
\]
其中每项都在相同原 joint 标签上积分；例如 \(T_{\rm high}=\mathfrak L_x\sum_{|A|>m}F_A\)。最后三个项均非负。

这是逐 \((x,\ell,z,\theta)\) 的有限 mask 代数身份，先成立再作正积分；没有交换实际 future cap、停止时间或原来源。对于 canonical actual 表示，右侧的 \(T_{\rm high}\) 就是原实际 accepted 高 endpoint 标签正交通。对于正支配表示，它是该表示的高标签项，不能据此反称原 actual 高标签已经有同样大小。

因此：

* 改动任意高标签 \(R_A\) 完全不改变 (2)，哪怕其 high-mask 差有很大的负 repeated 核。
* 改动低/空标签只改变 \(B_{R,\rm low}+B_{R,\varnothing}\)，同时改变其被单列的 signed 缺陷。误差费必须重新核，不可保留旧费而免费最小化新 reference。
* (2) 自身非负。原高标签 reference 的负校正全部与 \(B_{R,\rm high}\) 逐标签抵消，没有留下一个可额外利用的高标签负项。

更准确地，对两个 reference，
\[
 (B_R+\mathcal D_{\rm high}^R)
 -(B_{\widetilde R}+\mathcal D_{\rm high}^{\widetilde R})
 =\mathfrak L_x\sum_{|A|\le m}a_A(R_A-\widetilde R_A).
                                                               \tag{3}
\]
它的自由度只在未合回的低/空块。称全部 reference 为完全无关的 gauge 也不准确；**高块**是完全 gauge，低/空块存在有偿替换自由度。

如果要求空标签仍可免费按非正差丢弃，则必须满足 \(a_AR_\varnothing\ge F_\varnothing\)，而不能任意取 \(R_\varnothing=0\)。在这样的 reference 类中
\[
 B_R+\mathcal D_{\rm high}^R
 \ge T_{\rm high}+T_\varnothing.                              \tag{4}
\]
原 \(H_p\) 的空块满足这一条件，因为 \(e^{-np}\ge(1-p)^n\)。任意 reference 未必满足，空块正缺陷不能沿用 H 的旧符号结论。

## 3. 原 H 的联合核心究竟是什么

同一 \(c_t,L_s(x)\) 下取
\[
 R_A=h_{L_s}^{(n)}*H_{p_x,A},\quad
 H_{p,A}=e^{-np}\sum_{\operatorname{supp}\mathbf j=A}
                  \frac{p^{|\mathbf j|}}{\prod_i j_i!}G_c^{\mathbf j}.
\]
每个 support 下 repeated 系数非正的差已在上一稿准确保留。可是全 \(B_H\) 与高 signed 项 **联合**后，
\[
 V_H:=B_H+\mathcal D_{>m}^H
 =T_{>m}+B_{H,\le m}\ge T_{>m},                              \tag{5}
\]
这里 \(B_{H,\le m}\) 包含空 mask。

上一稿实际 source-once 费用仍有效：
\[
 \int_E[\mathcal D_{\rm low}^H]_+
 \le2C_hN_h\,\frac{m(m+1)(m+2)}{3(n+1)}W,
 \quad C_h=16/3,\quad N_h=1+n\log(b/a).
                                                               \tag{6}
\]
取 \(m=\lfloor n^{1/6}\rfloor\)，该费用不超过 \(4C_h\sqrt n\,W\)。本稿不重复其 6025 项新守卫。

因为空差非正，准确的后续接口为
\[
 T_{\rm act}\le V_H+e,\qquad
 e=[\mathcal D_{\rm low}^H]_+,\qquad
 T_{\rm act}\le\min(4\lambda,V_H)+e.                          \tag{7}
\]
这比先分开要求 \(B_H\) 与 \([\mathcal D_{>m}]_+\) 强 L1 费用更贴近原对象；但 (5) 揭示其本质：**它至少保留全部原 accepted 高标签正交通**。给 \(V_H\) 新弱预算确有用，却不能仅通过改变高 reference 获得。

原 capped_reference_absorption 工具仍可用于 (7)，只在已经证明 \(V_H\) 整体 weak 费用后应用。现没有这样的证明。不得将原固定 free \(H\) 的 weak 费用直接写给 \(V_H\)，因为后者还含原 positive 高标签本身。

## 4. 多 Poisson 时钟、实测时间与正混合均不突破 (2)

让 \(R_A^{(j)}\) 为不同 Poisson 时钟、不同实测半群时间或其它正 reference。允许有限 \(j\)，以及任意非负权
\[
 \gamma_{j,A}(x,\ell,z,\theta),\qquad \sum_j\gamma_{j,A}=1.
\]
定义 effective reference \(\overline R_A=\sum_j\gamma_{j,A}R_A^{(j)}\)。无论权是否依来源、接收点或历史，
\[
 \sum_j\gamma_{j,A}
   \{a_AR_A^{(j)}+[F_A-a_AR_A^{(j)}]\}=F_A
\]
逐标签成立；高块合回后仍只有 \(T_{\rm high}\)。输出后验选择时钟、source-only 分配时钟与正混合对此没有区别。

source-only 权、\(\sum_j\nu_j\le\nu\) 对 **独立成立的各来源费用** 可防止重复收 W；它不会改变上述 high-block gauge。若每个 reference 自有 weak 费用，则求和还须有限弱混合的有偿合同，不能用 \(\sum_j O(\log n)\|\nu_j\|\) 免费证明一个 sum 的 weak 范数。若时钟权依接收点，用作门可以，但不能当预先来源质量分配。

原 \(c_t\) 窄窗的 fixed-clock 网亦只能减少 parameter 困难：对 \(\rho=c_{\rm grid}/c_t\)，
\[
 \rho G_{c_{\rm grid}}-G_{c_t}
   =(\rho-1)G_{c_{\rm grid}}G_{c_t}\ge0,\qquad
 H_p^{c_t}\le e^{np(\rho-1)}H_{\rho p}^{c_{\rm grid}}.
\]
当 \(\rho-1\le2/n\)、\(p\le1\) 时系数不超过 \(e^2\)。它不能削掉 (5) 中的实际 \(T_{>m}\)，也不能将 moving \(h_L^{(n)}\) 当成同一个 source。

另外，不同正半群时间的核不存在可直接用的点态单调序。固定原 \(c,L\)、\(0\le\tau_1<\tau_2\)，两份 \(h_L^{(n)}*H_{\tau_j}^{c,L}\) 都有质量1。它们的 Fourier transforms 在原 \(h_L\) transform 非零的零频邻域满足
\[
 \widehat h_L(\xi)\,
 e^{-\tau_j\sum_i(1-\widehat G_c(L\xi_i))}.
\]
原非退化 \(G_c\) 在某非零小 \(\xi\) 上有 \(\widehat G_c<1\)，故两核不同。若一份点态支配另一份，由相等质量只能相等，矛盾。因此时间差的正负部分都非零；换到更晚时间本身没有全空间负符号。这里使用原核和完整 L1 密度，没有抽象有限链或 toy 反例；也没有否定任何弱预算。

## 5. 真正可缩小核心的合同必须针对原高标签，而非 reference 名字

(2) 直接给出最小缺口：要付 \(V_R\)，至少要对原保持全部门的 \(T_{>m}\) 证明空间/来源预算。reference 的高块不存在可优化掉的误差。

把 canonical 原高标签交通按原软来源 \(y\) disintegrate，可写
\[
 T_{>m}(x)=\int k_{\rm hi}(x,y;\mu)\,\mu(dy),\qquad k_{\rm hi}\ge0.
                                                               \tag{8}
\]
核 \(k_{\rm hi}\) 内仍包括原 \(\eta_C(y)\)、原首跳标签/退出、原硬 \(z\)、LCA/\(\theta\)、所有门、真实 continuation 参数和实际两赢家；不是新输入 \(h_L\mu\)，不重启 FIRST，也不是只剩一份 free kernel。

为了明确 fullfuture 对偶缺的量，下面是一种**充分证书格式，尚未获得**。若在同一原输入上能构造 \(\Phi(x,y)\ge0\) 和一个不依 \(y,C,\ell,z\) 的非负物理 future multiplier measure \(b_x(ds,dL)\)，使
\[
 \begin{split}
 k_{\rm hi}(x,y;\mu)&\le
       \Phi(x,y)+\int_{\substack{s\ge\sigma(x)\\L\in[a,b]}}
                      p_{s,L}(x-y)b_x(ds,dL),\\
 \int\mu(dy)\int_E\Phi(x,y)dx&\le B_nW,\\
 q\|b_x\|&\le\delta\lambda ,
 \end{split}                                                  \tag{9}
\]
则原完整 future cap、同一 \(\mu\) 与正 Tonelli 准确给
\[
 \int_E T_{>m}\le B_nW+\delta\lambda|E|.                         \tag{10}
\]
这是一份 source-once 空间缺陷，加可吸收项；不是一个纯计数/行 MGF。主账若用它，\(\delta\) 须小于 \((\gamma_n\,8192/49)^{-1}\)，不能任意按已吸收的旧角概率再次使用。

目前没有证明 (9) 的 \(\Phi\) 空间费用。把 \(b_x\) 选为每个首标签或每个捕获子源自己的倍率，不满足其共同性：fullfuture 只保证
\(\int p_{s,L}(x-y)\mu(dy)\le q\)，不保证带 \(y\)-依赖倍率后的同一个 q cap。把原 \(z\)/LCA 门移到 landing \(w\) 或 reference history 也不满足 (8) 的合同。现有 FIRST-MGF 支已支付的是原辅助 Bernstein 诊断尾，它没有给 (8) 的 after-first endpoint 高标签这个 \(\Phi\)。

这不是重新提交 fractional_actual_pair LP 作为改善；(9) 仅准确显示从 fullfuture 行约束通往本轮真实剩余空间费所缺的共同源证书。没有构造到这样的证书，就不能把任意 clock 混合称为 paid branch。

若为了简化而把 original diagnostic \(K\) 换成实际 visited-axis 或 after-first \(A\)，会改变新 joint law。此前 fullfuture/来源高权切分或可在另一耦合下重新证明，但当前最新角门的条件后验仅在原独立诊断分配下单次使用；不能保持原 \(R_{\rm angle}\) 名字与吸收常数，免费重标记所有历史。因此 \(K\) 的上限没有被本文当作 \(|A|\) 的上限。

## 6. 本轮边界与执行状态

本轮完成的都是一般解析审计：

1. (2)/(3) 在保持全部 actual joint 标签后成立，高 reference 的规范自由度被准确定位。
2. 原 H 联合高 signed 核心是正的 \(T_{>m}+B_{H,\le m}\)，并非一个仍有自由高负项的谱差。此前 signed 低校正费保留，不能把它当 low positive traffic 的整份付款。
3. 多时钟/多 reference 的正混合不改变实际高块；共同 source 拆分能管理已有费用，不能生成高块新预算。
4. 不同半群时间没有可免费使用的正核单调性；已核 \(c_t\) 网不解决原高块或 moving hard seed。

没有得到 (9) 或别的真正新高标签 source-once 空间估计。没有构造完整实际门下的反例，不能将未知高标签预算写成失败定理。late-p/低 mask 的 Binomial 下尾虽然可直接正包络支付，仅是已有 coefficient 工具的易支推论；按本轮范围不扩大、不作为 headline、不为它跑新试验。

由于没有形成新的可付高标签策略，本轮**不注册、不运行新数值**，不以有限 toy gauge 检查充当 actual 样本；既有 6025、523、316 等证书全部保持原终态。只新增本 md。当前一般 cube 目标及 \(R_{\rm angle}\) 仍未闭合，主账没有新费用或新裁剪。
