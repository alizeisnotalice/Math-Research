# 来源死亡过程的适应停止、质量倾斜与永生格子标签

2026-10-07。父任务指定 6.1-sol high。仅新增本文件，不改主笔记、原稀疏化稿或他人文件。本轮解析，没有运行数值。

**严格所得。** 已核定时 \(16uW\) 预算可推广到任意有界适应停时 \(T\)，费用精确等于 \(W\mathbb E_QT\)。完整质量倾斜 \(Q\) 可直接构造成“按原格子质量抽一个永生标签，其余格子独立指数死亡”。有限 \(N\) 格停到 \(k\ge1\) 格的 \(Q\) 平均时间准确为 \(H_{N-1}-H_{k-1}\)，与权重无关；输入复杂度仍未消除。停止到全灭给出不能无条件推广终值质量守恒的真实边界。

使用已读 [E04](</Users/zhengzhihao/.codex/skills/math-e04-snell-envelope-optimal-stopping/SKILL.md>)、[J01](</Users/zhengzhihao/.codex/skills/math-j01-inhomogeneous-jump-generator/SKILL.md>) 及其 method、cube-interface、provenance 的滤过、域和停止核验流程。本稿自证有限网格停止、密度合同和指数寿命表示，不援引未验收的连续 Snell／无界可选抽样定理。[来源稀疏化稿](source_cell_thinning_budget_20261007.md) 的碰撞、加权碰撞、真实删除及定时 \(16uW\) 已独审，不再次记为新预算。目录检索未找到此前已写出的本稿质量倾斜／永生标签合同。

## 1. 同一完整来源、真实滤过与定时构件

沿用原全连续窗口 \([a,b]\)、\(0<a<b\)、阈值 \(\tau>0\)、完整有限正来源 \(\mu\)、\(W>0\)。固定原半开网格 \(C_i\)，
\[
 0<d\le\min\{a/(8n),(b-a)/4\},\qquad
 \mu_i=\mu|_{C_i},\quad w_i=\mu_i(\mathbb R^n)>0,\quad\sum_iw_i=W.
\]
忽略零质量格；先设正质量格子数 \(N<\infty\)。原完整功能为
\[
 M_\nu(x)=\max_{a\le r\le b}r^{-n}\nu(Q(x,r)),\qquad
 I(\nu)=\tau|\{M_\nu>\tau\}|,\qquad
 \Phi(\nu)=\tau\int\log_+(M_\nu/\tau)dx,\qquad Z=1+n\log(b/a).
\]
已核构件对每份当前完整来源成立：
\[
 0\le\Phi(\nu)\le(Z/e)W(\nu),\qquad
 \sum_iD_i(\nu)\le I(\nu)+16W(\nu),\qquad
 D_i(\nu)=\Phi(\nu)-\Phi(\nu-\nu|_{C_i})\ge0.
\tag{1}
\]
后者来自原赢家后验碰撞 \(\tau\int\chi\le6W(\nu)\)、加权碰撞 \(\tau\int\chi\log(M_\nu/\tau)<2.5W(\nu)\) 及真实删除帽。每次死亡必须对当前完整来源重新使用它。

取独立 \(P\)-寿命 \(L_i\sim\mathrm{Exp}(1)\)，
\[
 A_t=\{i:L_i>t\},\qquad
 \mu_t=e^t\sum_{i\in A_t}\mu_i,\qquad
 W_t=e^tw_{A_t},\qquad w_A=\sum_{i\in A}w_i.
\tag{2}
\]
死亡时取删除后的右连续状态，空状态保持空。滤过 \(\mathcal F_t\) 只观察到时刻 \(t\) 的死亡及标签，不能提前知道未来寿命。可加入独立外部随机种子。

为构造 \(Q\)，采用允许 \(L_i=\infty\) 的规范路径空间和原始观察滤过。**不能把整个无限时域所有 \(P\)-零测事件先加入 \(\mathcal F_0\)**：“有一个格子永不死亡”在 \(P\) 下为零、在 \(Q\) 下为一，预先全局补齐会破坏有限时域密度。各有限时域的共同零测补齐可另行进行。

确定 \(s\le t\) 时，指数剩余寿命及条件 Tonelli 给
\[
 \mathbb E_P[W_t\mid\mathcal F_s]=W_s.
\tag{3}
\]
条件于 \(\mathcal F_s\)，已核定时 \(16u\) 预算应用于当前来源，给
\[
 \mathbb E_P[\Phi(\mu_t)\mid\mathcal F_s]
 \ge\Phi(\mu_s)-16(t-s)W_s.
\tag{4}
\]
当前为空时两边均零。这使用当前完整赢家，不冻结初始赢家。

## 2. 有界适应停时的自含有限网格证明

置 \(B_t=\int_0^tW_sds\)、\(X_t=\Phi(\mu_t)+16B_t\)。由 (3)、(4) 及条件 Fubini，
\[
 \mathbb E_P[X_t\mid\mathcal F_s]\ge X_s,\qquad
 \mathbb E_P[W_t\mid\mathcal F_s]=W_s.
\tag{5}
\]
固定 \(h<\infty\)，逐路径 \(W_t\le e^hW\)、\(\Phi(\mu_t)\le(Z/e)e^hW\)、\(B_t\le he^hW\) 对 \(t\le h\) 成立。有限格路径有有限次死亡，\(\mu_t\) 在 TV 下右连续；\(\Phi\) 的 TV-Lipschitz 界使 \(X_t,W_t\) 右连续。

阈值平台处右增长导数是 \(\tau|\{M_\nu\ge\tau\}|\)，严格 \(I(\nu)\) 只是沿无跳路径的 \(dt\)-几乎处处导数。这里用的是已正确处理平台的积分合同 (4)，不假设全状态点态生成元为 \(I-\sum D_i\)。

令 \(T\le h\) 为适应停止时。向上有限网格取整得 \(T_m\downarrow T\)。对相邻网格时刻 \(t_j,t_{j+1}\)，事件 \(\{T_m>t_j\}\in\mathcal F_{t_j}\)，逐路径
\[
 X_{T_m}=X_0+\sum_j\mathbf1_{\{T_m>t_j\}}(X_{t_{j+1}}-X_{t_j}).
\]
条件期望逐项非负，故 \(\mathbb E_PX_{T_m}\ge X_0\)。同样展开 \(W\) 的增量得 \(\mathbb E_PW_{T_m}=W\)。右连续与确定有界包络允许直接支配收敛，得到
\[
 \boxed{\mathbb E_P\Phi(\mu_T)\ge\Phi(\mu)
 -16\mathbb E_P\int_0^TW_sds,\qquad \mathbb E_PW_T=W.}
\tag{6}
\]
看见全部未来寿命后再挑时刻不属于该适应停止类。

## 3. 质量倾斜的全路径构造与准确 h-transform

对有限 \(h\)，定义 \(Q_h(F)=\mathbb E_P[\mathbf1_FW_h]/W\)，\(F\in\mathcal F_h\)。由 (3) 这些概率一致。直接构造共同 \(Q\)，无需仅凭一致性引用扩张定理：

1. 初始按 \(Q(J=i)=w_i/W\) 抽一个格子标签。
2. 置 \(L_J=\infty\)，其余寿命独立为 \(\mathrm{Exp}(1)\)。
3. 保持 (2) 的来源形状及补偿因子 \(e^t\)。

固定标签 \(J=i\) 后，有限时域观察分布等于 \(P\) 条件于 \(L_i>h\) 的分布。因此
\[
 Q(F)=\sum_i\frac{w_i}W e^hP(F\cap\{L_i>h\})
 =\mathbb E_P\left[\mathbf1_F
 \frac{e^h\sum_{i\in A_h}w_i}W\right]=Q_h(F).
\tag{7}
\]
这是同一来源的共同全路径构造，标签保护整个格子的原测度，不以一个抽样点替换该格。也可先按 \(\mu/W\) 抽来源点再取所属格，但不产生来源坐标独立性。

隐藏标签的真实后验为
\[
 Q(J=i\mid\mathcal F_t)
 =\frac{w_i\mathbf1_{\{i\in A_t\}}}{w_{A_t}}.
\tag{8}
\]
在 \(Q\) 下存活集永不为空。当前非空集合 \(A\) 中格子 \(i\) 的死亡率为
\[
 \boxed{q_i(A)=1-\frac{w_i}{w_A}
 =1-\frac{\mu_t(C_i)}{W_t}\quad(i\in A).}
\tag{9}
\]
若标签是 \(i\)，该格不死；否则其剩余寿命率为 \(1\)。有限集合 \(A\) 的总率为
\[
 \sum_{i\in A}q_i(A)=|A|-1.
\tag{10}
\]
小步密度也直接核验该公式：\(P\) 无死亡概率 \(1-|A|\Delta+o(\Delta)\) 乘无死亡的质量比 \(e^\Delta\)，给 \(Q\) 无死亡概率 \(1-(|A|-1)\Delta+o(\Delta)\)。格 \(i\) 死亡概率 \(\Delta+o(\Delta)\) 乘 \(e^\Delta(w_A-w_i)/w_A\)，给 (9)。确定增长仍为 \(d\mu_t/dt=\mu_t\)，没有额外来源漂移。

## 4. 费用与停止时密度合同

对任意停止时（费用允许无穷），非负 Tonelli 与 \(\{T>s\}\in\mathcal F_s\) 给
\[
 \boxed{\mathbb E_P\int_0^TW_sds
 =\int_0^\infty\mathbb E_P[W_s\mathbf1_{\{T>s\}}]ds
 =W\int_0^\infty Q(T>s)ds=W\mathbb E_QT.}
\tag{11}
\]
中间使用有限时域密度 (7)，不预设终值质量守恒。

若 \(T\le h\)、\(F\in\mathcal F_T\)，简单停止时由 (3) 逐格得到
\(\mathbb E_P[W_h\mathbf1_F]=\mathbb E_P[W_T\mathbf1_F]\)。向上网格逼近仍有 \(F\in\mathcal F_{T_m}\)，有界支配给
\[
 Q(F)=\mathbb E_P[W_T\mathbf1_F]/W.
\tag{12}
\]
所以有界停止的完整合同为
\[
 \mathbb E_P\Phi(\mu_T)\ge\Phi(\mu)-16W\mathbb E_QT,
 \qquad\mathbb E_PW_T=W.
\tag{13}
\]
在 \(Q\) 下 \(W_T>0\)。其真实终值比值 \(G_T=\Phi(\mu_T)/W_T\) 满足
\[
 \mathbb E_QG_T=\mathbb E_P\Phi(\mu_T)/W
 \ge\Phi(\mu)/W-16\mathbb E_QT.
\tag{14}
\]
若原输入 \(\Phi(\mu)/W\ge K-\varepsilon\)，故存在正质量停止实现保留比值至少 \(K-\varepsilon-16\mathbb E_QT\)。合法来源类须对删除和放大封闭，完整正 Borel、完整非负 \(L^1\) 类都满足。

更直接的 \(Q\)-条件合同为
\[
 G_t=\Phi(\mu_t)/W_t,\quad0\le G_t\le Z/e,\qquad
 \mathbb E_Q[G_t\mid\mathcal F_s]
 =\frac{\mathbb E_P[\Phi(\mu_t)\mid\mathcal F_s]}{W_s}
 \ge G_s-16(t-s).
\tag{15}
\]
它在 \(Q\)-几乎处处的正质量状态解释。故 \(G_t+16t\) 有条件次鞅合同，有界停止仍由有限网格自证。不能把这个来源过程称成原几何 Snell 预算已经付款。

## 5. 可数格的 TV 路径与有界停止

若正质量格子可数无穷，仍取独立 \(P\)-寿命并按 \(w_i/W\) 抽永生标签。非负 Tonelli 使 (3)、(7)、(8)、(11) 原样成立。固定 \([0,h]\)，同寿命有限截断有确定误差
\[
 \sup_{t\le h}\|\mu_t-\mu_t^{(N)}\|_{\rm TV}
 \le e^h\sum_{i>N}w_i\longrightarrow0.
\tag{16}
\]
此界在 \(P,Q\) 下均成立。因有限截断为 TV 右连续、有左极限的路径，完整 \(\mu_t\) 也具有这些性质。TV-Lipschitz 界使 \(\Phi_t\) 右连续，仍有共同确定包络 \(W_t\le e^hW\) 等。

当前完整可数来源已有定时 (4)，故 (5) 及有界网格停止直接推广，无需让截断赢家一致。也可在全寿命滤过中对有限截断证明：遗漏寿命是与截断独立的外部随机性，条件合同仍成立，再用 (16) 取极限。不能将全滤过停时未经说明当成只观察有限截断的停时。

单格边际率 (9) 仍有效；无限存活格子的总率为无穷，不能用有限总率 (10) 或声称紧时域只有有限次死亡。质量加权 TV 路径、定时合同和有界网格停止不需要这项错误假设。

## 6. 无界停止、class-D 与未停止尾项

任意停止时的截断 \(T\wedge h\) 都有 (13)。由 (12) 对 \(\{T\le h\}\) 应用有界停止，
\[
 \mathbb E_P[W_T\mathbf1_{\{T\le h\}}]=WQ(T\le h),
\]
从而在有限停止事件上准确有
\[
 \boxed{\mathbb E_P[W_T\mathbf1_{\{T<\infty\}}]
 =WQ(T<\infty).}
\tag{17}
\]
故 \(P(T<\infty)=1\) 单独不保证终值质量为 \(W\)，还须核 \(Q(T<\infty)=1\)。

若 \(Q(T<\infty)=1\)、\(\mathbb E_QT<\infty\)，则
\[
 0\le\mathbb E_P[\Phi(\mu_h)\mathbf1_{\{T>h\}}]
 =W\mathbb E_Q[G_h\mathbf1_{\{T>h\}}]
 \le(Z/e)WQ(T>h)\longrightarrow0.
\tag{18}
\]
把
\[
 \mathbb E_P\Phi(\mu_{T\wedge h})
 =\mathbb E_P[\Phi(\mu_T)\mathbf1_{\{T\le h\}}]
 +\mathbb E_P[\Phi(\mu_h)\mathbf1_{\{T>h\}}]
\]
代入截断合同，第一项与 (11) 的费用单调收敛，给出
\[
 \mathbb E_P[\Phi(\mu_T)\mathbf1_{\{T<\infty\}}]
 \ge\Phi(\mu)-16W\mathbb E_QT,\qquad
 \mathbb E_P[W_T\mathbf1_{\{T<\infty\}}]=W.
\tag{19}
\]
若同时 \(P(T<\infty)=1\)，就是通常的完整终值陈述；否则本稿仅对有限停止贡献陈述 (19)，不擅造未停止实现的终值来源。正质量实现选择仍在有限停止事件上成立。

若 \(Q(T=\infty)>0\)，(18) 尾项不能删除，(17) 显示质量流失。若 \(\mathbb E_QT=\infty\)，费用是无穷，所得下界通常平凡。无界停止没有上述条件或终值一致可积性时，只能保留截断合同及尾项。

有限格全灭时间 \(T_0=\max_iL_i\) 是具体反例：\(P(T_0<\infty)=1\)，但 \(\mu_{T_0}=0\)、\(\mathbb E_PW_{T_0}=0\)；在 \(Q\) 下永生标签使 \(T_0=\infty\)。原 \(P\)-过程 \(W_h\to0\) 几乎处处，却每个 \(h\) 有均值 \(W\)，故无界时域质量族不一致可积，不能免费宣称全局 class-D。它与有界停止证明不矛盾。

## 7. 停到剩 \(k\) 格：准确调和时间与复杂度

有限初始正质量格子数 \(N\)，\(1\le k\le N\)，定义
\[
 T_k=\inf\{t:|A_t|\le k\},\qquad H_0=0,\quad H_m=\sum_{j=1}^m1/j.
\]
在 \(Q\) 下固定永生标签，其余 \(N-1\) 个寿命独立为 \(\mathrm{Exp}(1)\)。从 \(j\) 个存活格子降到 \(j-1\) 个时，有 \(j-1\) 个非永生指数剩余寿命，最小值的存活函数为 \(e^{-(j-1)t}\)，均值 \(1/(j-1)\)。死亡后重复同一记忆无关计算，得到
\[
 \boxed{\mathbb E_QT_k
 =\sum_{j=k+1}^N\frac1{j-1}
 =H_{N-1}-H_{k-1}.}
\tag{20}
\]
不用外部纯死亡链定理。时间分布也与权重无关；死亡标签分布仍由 (9) 依赖权重，不是均匀标签删除。

这些停时在 \(P,Q\) 下都几乎处处有限，且 \(Q\) 平均有限，故 (19) 给
\[
 \mathbb E_P\Phi(\mu_{T_k})
 \ge\Phi(\mu)-16W(H_{N-1}-H_{k-1}),\qquad
 \mathbb E_PW_{T_k}=W.
\tag{21}
\]
近优输入有正质量停止实现，比值至少
\[
 K-\varepsilon-16(H_{N-1}-H_{k-1}).
\tag{22}
\]
\(k=N\) 停时为零，\(k=1\) 费用为 \(16H_{N-1}W\)；\(k=0\) 不属于这些公式。

若初始有可数无穷个正质量格子，任意有限时刻仍有无穷个非永生格子存活：独立正保留率使每个无限尾列仅有限存活的概率为零，再对整数时刻取可数交。于是 \(P,Q\) 下 \(T_k=\infty\) 对每个有限 \(k\) 成立。有限质量不等于有限格子数量，不能将截断调和公式变成统一复杂度预算。按质量定义新停止规则可能有用，但须另控它的实际 \(\mathbb E_QT\)。

## 8. 对一般主接口的准确位置

本稿消除了定时稀疏化不能依已见信息停止的缺口，给出完整来源一次的精确费用和永生标签法则；没有消除输入复杂度。想压到常数格，有限来源仍要付可能任意大的调和费用，可数来源按数量停止甚至永不到达。

仍须从真实中心立方体赢家、实际高值资格及完整来源几何构造停止规则，控制 \(\mathbb E_QT\) 的目标维数费用，并证明停止后的真实输入属于已有预算分支。该过程不是原 FIRST/fullfuture/CP--GP/LCA/allhistory 的跨输入耦合；永生标签也不是 receiver 的免费 Palm 产品化。一般 \(\sqrt n\,n^{o(1)}\) 闭合仍未完成。
