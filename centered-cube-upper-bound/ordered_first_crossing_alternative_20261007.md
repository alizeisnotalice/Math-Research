# 固定 receiver 阈值的首次触达：hinge 付款与回落余项

2026-10-07，gated_radial_energy_audit。本稿固定原 \(c>0\)、物理尺度、有限 \(Z\)，任意非负 \(f\in L^1(\mathbb R^n)\)，\(W=\int f\)。它给出完整一般的 receiver 弱接口：共同确定阈值 \(\eta>0\) 的 first-touch gate 消除双方 alive 的正 commutator，并使**当前跨过 \(\eta\) 的全部定向边**由一份 \(W\) 支付，包含任意小 contrast。剩余是 inactive 端已在阈值下的回流；尚未闭合一般 polylog ordered weak 预算或原 cube 主账。

查重已读 ordered_propagation_endpoint_20261007.md §4（原 first-event 出口强范数反例，不是 weak 反例）、stopped_budget_audit/counterexample.md（实际 middle 正账的共同尺度取消障碍）、receiver_crossing_budget_next.md、continuous_crossing_next.md，以及当前 ordered_oriented_small_contrast_20261007.md。未重跑其模型。使用 L03 的精确常数／诊断分离流程、J02 的杀死与 Duhamel 合同；均读 SKILL、provenance、method、cube-interface，技能不是本稿定理的外部证明。

## 1. 阈值触达与保留未完成 receiver

原传播及真实生成元为
\[
v_z=T_{1-e^{-z}}^cf,\qquad
J_z=\sum_i(G_{ce^{-z},i}-I),\qquad v'_z=J_zv_z.
\tag{1}
\]
有限 mask 展开给同一满测集上的连续、实际为有限解析系数的逐点轨道。取规范代表；固定零测异常集不改 Lebesgue 或 \(L^1\) 来源。

定义 receiver 的确定时间
\[
\tau_\eta(x)=\inf\{z\in[0,Z]:v_z(x)\ge\eta\},
\]
无触达时取 \(\infty\)。若 \(f(x)\ge\eta\)，置 \(\tau_\eta=0\)。触达指闭水平，允许相切和 \(z=Z\)；不要求单调响应或正导数。它可测，因为连续轨道的紧区间最大值可由可数有理点给出。

定义全 receiver prehit gate
\[
H_z(x)=\mathbf1_{\{f<\eta\}}\mathbf1_{\{z<\tau_\eta(x)\}},\qquad
t_z(x)=H_z(x)(1+v_z(x)/\eta).
\tag{2}
\]
这里**包括截至 \(Z\) 尚未触达的 receiver**。所以不需要未来 hit 集合作为空间门；before-hit \(v_z<\eta\)，故 \(0\le t_z<2\)。空间上 \(H,t\) 可能有无限支撑，不能称它们是有限 receiver 集合的 \(L^2\) test。

记 \(F=\{f<\eta,\tau_\eta\le Z\}\)、\(N=\{f<\eta,\tau_\eta=\infty\}\)。逐 receiver 微积分与 \(\int_0^Z\|J_zv_z\|_1dz\le2nZW\) 给精确财富式
\[
\eta|F|+\int_N v_Z-\int_{\{f<\eta\}}f
=\int_0^Z\!\int H_zJ_zv_z.
\tag{3}
\]
终端未完成项 \(\int_Nv_Z\) 为正，必须先保留，再为上界舍去。finite mask 上包络使 \(F\) 有限测度。对 \(E=\{\max_{[0,Z]}v_z>\eta\}\)，\(E\subset\{f\ge\eta\}\cup F\)，且两初始来源区域互斥。因此
\[
\eta|E|\le W+\int_0^Z\!\int H_zJ_zv_z-\int_Nv_Z.
\tag{4}
\]
初始高、低质量只收一份 \(W\)，不是分别上界后收两份。

## 2. 无限支撑 test 的合法交换与精确正余项

令 \(L(v)=\log(1+v/\eta)\)，\(\Phi_\eta(v)=v-\eta L(v)\)。已有准确 Bregman 式
\[
\frac{\eta}{\eta+v}J_zv=\eta J_zL(v)+b_z,\qquad b_z\ge0,
\quad D=\int_0^Z\!\int b_z\le W.
\]
由于 \(t\eta/(\eta+v)=H\)，(4) 的 flow 等于
\[
\int t_zb_z+C,\qquad
C=\eta\int_0^Z\!\int t_zJ_zL(v_z)
=-\frac{\eta}{2}\int d\mathfrak m\,\Delta t\,\Delta L(v),
\tag{5}
\]
其中 \(d\mathfrak m=\sum_i dz\,dx\,G_{ce^{-z},i}(dh)\)。

绝对 Fubini 不使用 \(t\in L^1\) 或 \(L^2\)：\(|t|\le2\)、\(L(v)\le v/\eta\)，因此 \(\int|tJL|\le4nW/\eta\) 逐时；对称边积分亦由 \(|\Delta t|\le2\)、\(\int G|\Delta L|\le2W/\eta\) 控制。也可用 \(|Jt|\le4n\) 与 \(L(v)\in L^1\) 解释交换。outside-source 常数部分的差分为零，不能对其裸 \(L^2\) 范数收费。

两端 alive 时 \(t=1+v/\eta\)，\(\Delta t\,\Delta L\ge0\)，所以整支 commutator 非正。两端 inactive 时 \(t=0\)，费用为零。正支只可能由 alive 端 \(x\) 指向 inactive 端 \(x+h\)。令 \(v=v_z(x)\)、\(u=v_z(x+h)\)，其 signed cross 部分精确为
\[
C_{\rm cross}=\int d\mathfrak m\,
H_z(x)[1-H_z(x+h)](\eta+v)
\log\frac{\eta+u}{\eta+v}.
\tag{6}
\]
系数没有 \(1/2\)：反向边给相同 signed cross。正费用需要 \(u>v\)。这是 receiver gate，不是 Markov 轨迹 stopping law。

## 3. 当前 level 跨边由 hinge 熵支付

设 \(\Psi_\eta(s)=(s-\eta)_+\)。它在非负 \(L^1\) 上是 Lipschitz，且 \(0\le\int\Psi_\eta(f)\le W\)。沿 (1) 的轨道，逐点绝对连续链式与 Fubini 给
\[
\int\Psi_\eta(f)-\int\Psi_\eta(v_Z)
=Q_\eta,\qquad
Q_\eta=\int d\mathfrak m\,
\mathbf1_{\{u\ge\eta>v\}}(u-v)\le W.
\tag{7}
\]
其中 \(u\ge\eta\) 的等号处理如下：对几乎处处 \((z,x)\)，\(v_z(x)=\eta\) 时 \(v'_z(x)=0\)（绝对连续标量函数的 level-set 性质）。所以在链式微分中可选 \(\Psi'_\eta=\mathbf1_{s\ge\eta}\) 而非 \(\mathbf1_{s>\eta}\)，不改变积分。对称化后正反方向抵消 \(1/2\)，得到 (7)。有限时间的 edge 积分以 \(u+v\) 支配；或用左移 hinge 的光滑逼近并保留闭上端点。不能丢弃有正空间质量的水平面。

在 (6) 的 alive 端必有 \(v<\eta\)。若另端当前 \(u\ge\eta\)，则它必 inactive，并且
\[
(\eta+v)\log\frac{\eta+u}{\eta+v}\le u-v.
\tag{8}
\]
因此当前高端的整个定向支 \(C_{\rm high}\) 非负，且
\[
\boxed{C_{\rm high}\le Q_\eta\le W.}
\tag{9}
\]
没有 \(1/\delta\) 的 contrast 代价，即使 \(u,v\to\eta\) 也成立。这比“只付 regularized 大 contrast”增加了一个真正通用同边付款：只要求当前边跨越共同固定水平。

## 4. 准确回落余项与一般弱接口

保留 signed 的低 inactive 端交叉量
\[
R_\eta^{\rm signed}
=\int d\mathfrak m\,H_z(x)[1-H_z(x+h)]
\mathbf1_{\{u<\eta\}}(\eta+v)
\log\frac{\eta+u}{\eta+v}.
\tag{10}
\]
它可正可负；其正定向上包为
\[
R_\eta^+
=\int d\mathfrak m\,
\mathbf1_{\{H_z(x)=1,\ H_z(x+h)=0,\ \eta>u>v\}}
(\eta+v)\log\frac{\eta+u}{\eta+v}.
\tag{11}
\]
inactive 高端分类是：

- 初始 \(f(x+h)\ge\eta\)，目前已经降到 \(u<\eta\)；
- 初始低，但在此前 receiver first-touch 之后，当前回落至 \(u<\eta\)。

这不是“第二次 hit”：(11) 本身是一条空间 jump 边，可能没有后续 receiver 第二次上穿。未触达 receiver 保持 alive；没有原 future-event 空间门造成第三类 inactive。

由双方 alive 非正、(9) 和 \(t b\le2b\)，得完整任意输入的接口
\[
\boxed{\eta|E|
\le W+2D+Q_\eta+R_\eta^{\rm signed}
\le4W+R_\eta^{\rm signed}
\le4W+R_\eta^+.}
\tag{12}
\]
其中可继续保留 (4) 的负 terminal nonhit 项。所有阈值、来源、receiver 都是一份原 \(f\)。这个精确余项比原 max-winner 小 contrast 分类具有新结构，但其空间费尚未由 \(W\) 的 polylog 给出。

若对输出越阈高度使用 \(2^j\eta\) 等确定分层，(7) 只逐固定阈值给 \(Q_{2^j\eta}\le W\)。同一来源可能跨多层，不能无限相加这一 \(W\)。本稿没有将首次触达或分层自动变成来源熵的可和预算；只用一个共同水平即可表示 \(\{\max v>\eta\}\)，并不需要先知道 \(M\) 的大小。

## 5. 单独 test-energy 仍不是一般可付合同

保留未hit点修复了未来 membership，却不自动支付 test-energy。作为已审 fast 模型的解析推论（不重跑数据），取同一个原 \(c=1,Z=1\)、\(f=1+\varepsilon\cos(64\sum_i x_i)\)，但本节 first-touch 阈值取 \(\eta=1\)。正 cosine 处初始已高，所以 \(H=0\)；负 cosine 处 \(v_z<1\) 在整个有限区间内成立，故 \(\tau_\eta=\infty\)、\(H=1\)。于 \(0<z<1\)：
\[
t_z=\mathbf1_{\{\cos<0\}}\,[2+\varepsilon m_n(z)\cos].
\]
已有原核高频谱界给基准 \(2\mathbf1_{\{\cos<0\}}\) 的每轴 Dirichlet 至少 \(2(99/100)\)；扰动的 Dirichlet 半范数至多 \(\varepsilon\)。因此每轴真实 test Dirichlet \(>1\)，scaled \(a(z)>n/2\)。

用已审周期大盒提升可得有限非负 \(L^1\) 输入同样有 \(\int a=\Omega(n)W\)。这里 bulk 的负 cosine 从未触达仍正确：自身截断响应不超过原周期正输入的响应；正 cosine 初始高域也精确保留。该例 actual weak 输出已由初始高域支付，不能用作 weak 反例。它只否定本 first-touch 版本的 standalone polylog test-energy 合同，说明仍须保留 (9)／(10) 的同边符号而不是回到粗 Cauchy。

## 6. 人为构造合法 killed process：一次吸收不等于一次回流

可作一个**新的辅助路径构造**，而不是把固定 receiver 的 \(\tau_\eta(x)\)直接称为 Markov 停时。令 \(A_z=\{H_z=1\}\)，随时间递减。原跳过程在有限 \(Z\) 内总率 \(n\)，几乎必然只有有限次跳。对每个 piecewise-constant 路径位置，显式设置该位置的 deterministic deadline \(\tau_\eta(x)\)；landing 不在 \(A_z\) 时立即杀死；在保持段达到 deadline 也杀死；初始 \(f\ge\eta\) 的路径在 0 杀死。这才是以已知 deterministic Borel 场构造的路径 stopping rule。deadline、landing 与 terminal equality 均保留。

令 \(\ell_z\) 是 surviving 密度，\(\kappa(ds,dy)\) 是一次首次吸收的正时空来源，包含时间 0 的初始高源。原 bounded-rate path conditioning 给
\[
v_z=\ell_z+\int_{[0,z]}U(z,s)\kappa(ds),\quad
\ell_z=0\ \hbox{于 }A_z^c,\quad
\|\ell_z\|_1+\kappa([0,z]\times\mathbb R^n)=W.
\tag{13}
\]
这里 \(U\) 是 (1) 的真实非齐次正 cocycle。可直接由有限次 Poisson 跳与 holding 的条件分解建立，不调用未经核对的无限杀势极限。首次吸收确有一次质量 \(\kappa\le W\)，但它是辅助 barrier 的来源，不继承原 actual FIRST/LCA/history 语义。

记 \(\beta_z=v_z-\ell_z\)。在 inactive 位置 \(y\)，\(v_z(y)=\beta_z(y)\)。因此 (8) 的 \(\log\le\) 线性支配也给
\[
R_\eta^+\le
\int_0^Z dz\sum_i\int dx\,H_z(x)
\int G_{ce^{-z},i}(dh)\,
\mathbf1_{\{H_z(x+h)=0,\ v_z(x+h)<\eta\}}\beta_z(x+h).
\tag{14}
\]
对 (13) 做正 Tonelli，把右端写为每个首次吸收来源的 free continuation **从 inactive 返回 alive 的期望 jump 计数**。同一个已吸收 packet 可以返回、再次离开，再返回；不可将这个计数当成一份首次吸收概率。当前严格粗界仅为
\[
R_\eta^+\le n\int_0^Z\|\beta_z\|_1dz
=n\int_{[0,Z]}(Z-s)\,\kappa(ds,dy)\le nZW.
\tag{15}
\]
这是一份来源的真实有偿表示，但阶数仍不足。进一步的 polylog 付款需要利用这些真实回流的空间结构或 signed 取消；仅“inactive 不可逆增长”与 \(\kappa\le W\) 尚不能关闭 (14)。原 first-event 强范数反例亦禁止把这种 free continuation 的强最大积分免费压成 polylog。

## 7. 历史 deficit 的准确恒等式：仍不能循环付款

root 提示的历史 deficit 可在原 receiver 轨道上严格核对。令 \(B_z=1-H_z\)，
\[
p_z(x)=B_z(x)(\eta-v_z(x))_+.
\]
初始 \(B_0=1\) 处 \(f\ge\eta\)，所以 \(p_0=0\)；新 inactive 出生恰在 \(v=\eta\)，其 deficit 亦为零。因此 \(p_z\) 沿每条 receiver 轨道绝对连续，且几乎处处
\[
\partial_zp_z=-B_z\mathbf1_{\{v_z<\eta\}}J_zv_z.
\]
定义原未归一化 signed 线性回流及正 refill：
\[
\begin{split}
R_{\rm hist}&=\int d\mathfrak m\,
\mathbf1_{\{B_z(x)=1,\ v_z(x)<\eta,\ H_z(x+h)=1\}}
[v_z(x)-v_z(x+h)],\\
Q_{\rm refill}&=\int d\mathfrak m\,
\mathbf1_{\{B_z(x)=1,\ v_z(x)<\eta,\ v_z(x+h)\ge\eta\}}
[v_z(x+h)-v_z(x)].
\end{split}
\]
两个 inactive 低端内部的线性边精确取消；当前高端必 inactive。由正/绝对 Fubini：
\[
\boxed{\int p_Z=R_{\rm hist}-Q_{\rm refill},\qquad
0\le Q_{\rm refill}\le Q_\eta\le W.}
\tag{18}
\]
全时间导数由 \(2nZW\) 支配；closed level 或出生端点不额外产生财富。

这还给 (10) 与 deficit 的精确关系。将空间边定向为 active 低端值 \(v\)、inactive 且低于阈值的另一端值 \(u\)，令
\[
B_{\rm rec}=\int d\mathfrak m\,H_z(x)[1-H_z(x+h)]
\mathbf1_{\{u<\eta\}}(\eta+v)
\left[\frac{\eta+u}{\eta+v}-1-\log\frac{\eta+u}{\eta+v}\right].
\]
其非负且 \(B_{\rm rec}\le2D\le2W\)，并且
\[
R_\eta^{\rm signed}=R_{\rm hist}-B_{\rm rec}
=\int p_Z+Q_{\rm refill}-B_{\rm rec}.
\tag{19}
\]
所以仅改名历史 deficit 并没有付款。事实上 \(B_Z\) 是闭 receiver 曾触达集合（含初始高端），且
\[
\int p_Z\le\eta|B_Z|
\le\int p_Z+\int v_Z\le\int p_Z+W.
\tag{20}
\]
因此支付终端 deficit 已经与支付这一 receiver 水平集体积等价，最多差一份 \(W\)。式 (18)–(20) 是反循环合同，不是新的空间费用。

## 8. 新原核探路与证据界限

新数据只使用一个未运行过的双频正输入：
\[
f=\frac45+\frac7{20}\cos\theta-\frac25\cos\varphi,\quad
\theta=4x_1,\quad \varphi=64\sum_i x_i,\quad \eta=1,\ c=1,\ Z=1.
\tag{16}
\]
完整 torus 平均质量 \(W=4/5\)，且 \(f\ge1/20\)。对 \(n\ge2\)，两相位的 Haar 推前为二维均匀相位；这是同一相关来源的两变量表示，不是两个抽象互不相关 response。原响应准确为
\[
v_z=\frac45+\frac7{20}b_4(z)\cos\theta
-\frac25b_{64}(z)^n\cos\varphi,\quad
b_k(z)=e^{-z}+(1-e^{-z})g_1(k^2).
\]
在相位 Fourier 标签 \((j,l)\) 上，真实 \(J_z\) 符号是
\[
g_{e^{-z}}((4j+64l)^2)-1
+(n-1)[g_{e^{-z}}((64l)^2)-1].
\tag{17}
\]
尤其第一轴不能独立 reset 两相位。receiver first-touch 用此精确双频响应求根，所有未hit点保留；记录 hinge-drop/current-level 通量、(3) stopped wealth、完整 commutator、高当前端和低 inactive 端 signed 分支。它允许先上穿后回落，避免只重跑旧单 mode。

first-touch 的峰值搜索没有使用时间网格漏峰。令 \(t=e^{-z}\)，(16) 的响应是常数加一个 \(t\) 的线性项再减一个正仿射函数的 \(n\) 次幂（两个 cosine 系数各有自己的符号）。导数最多一个内驻点；只有两个 cosine 系数均正时该驻点为最大点，其余符号情形只需端点或一个最小点。脚本解析计算全部合法内最大点，与两个端点比较，再在包含唯一首次上穿的 bracket 中做 64 次 bisection；两个 cosine 系数均负、响应先降后升时，谷前始终低于阈值，所以谓词 \(v_z\ge\eta\) 在这个 bracket 中仍只由假变真一次。closed 触达包括最大点恰等于1的相切和 \(z=Z\) 触达；浮点未对这些等号作 interval 认证，但算法未把它们定义为拒绝。已执行脚本注释把该 bracket 称为 “increasing interval” 用语偏强，但算法不变；脚本未改，保留执行哈希 \(824b70654a450b41523380a2b2fc8267e7cbb1c895cd49f98dac41dee6bef715\)，不重跑。

按 L03，已在执行前写同前缀 registration。只执行新 ordered_first_crossing_alternative_guard_20261007.py 一次，终态 exit0，约1.21秒。results 同前缀保存三轮 \(n=8,32,128\)，相位点每轴32/64/128、时间 Simpson 128/256/512 panels。**2766项 Fraction exact_PASS** 核 (8) 的精确 \(\log r\le r-1\) 后常数、含等号 hinge scalar 正性和正来源最小值。另 **54项双精度 screen 全通过**；固定预注册绝对容差0.03未事后修改。

|n|新触达比例|hinge 通量积分|当前高端 commutator|低 inactive signed 回流|终端 deficit|
|---:|---:|---:|---:|---:|---:|
|8|0.0751953|0.0749725|0.0640521|0.0355550|0.0345285|
|32|0.1115723|0.0750304|0.0634800|0.0391655|0.0375417|
|128|0.1243286|0.0750099|0.0631459|0.0408972|0.0389940|

三轮初始质量均 \(W=0.8\)。财富式时间积分与精确相位网格 terminal 式的误差绝对值分别 \(1.19,2.63,2.27\times10^{-5}\)；hinge 积分／drop 的误差分别 \(2.16\times10^{-5},2.16\times10^{-6},1.22\times10^{-5}\)；历史 deficit 身份误差分别 \(9.70\times10^{-6},2.84\times10^{-5},1.06\times10^{-5}\)。这些是非区间诊断，不保证 sharp 相位边界已被误差认证；不同维数与不同网格也不能作为阶数拟合。

原符号 FFT、first-touch 求根、Simpson 与 sharp gate 相位积分属于双精度诊断；周期输入不是 \(\mathbb R^n\) actual geom 样本。一般 (3)–(20) 由上文解析证明承担。未运行旧 mode、旧 first-event 或其他终态数据；没有 live handle。当前可接受的新工具是 (9)/(12) 的 hinge 同边付款与准确回落补集，仍是通用弱端点的等价余项接口，没有逼近闭合或新 cube 主账费用。
