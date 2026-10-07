# 全中间尺度真赢家：无 atom-count 的径向输运及全局方波测试

2026-10-07，nearflat_cube_geometry，父任务继续指定 6.1-sol high。仅新增本稿，不改主 TXT、已核 nearflat 稿或他人文件。读 D01 SKILL、method、provenance、cube-interface，采用正概率边缘、成本方向及迁移量词审计；下面仅用直接耦合下界、CDF/Fubini，不调用未读 Kantorovich 强对偶或最优计划达值定理。

**严格进展：全部原中间尺度约束给出 arbitrary μ 的径向传输下界；固定、与 receiver 无关的随机周期方波将它合法代入 nearmax 对称切换预算。不存在 atom 数量费用。该具体全局化有显式 n 级坐标损失，加上窄尺度的 1/n 几何距离，目前信号可仅为 n^-2；不能由此声称 sqrt(n) 付款。**

已检索本目录及 cube_general_20261003 的 radial transport、posterior CDF、periodic threshold/global test 相关内容，读旧 sparse_pair_cover/l1_neighbor_transport 和 conditional_cube_transport/band_matrix/cylinder_pair/theorem。旧随机 grid 处理来源点对的近 l1 graph retention 与 Schur；它不等于这里的两个不同捕获后验之 CDF 差及同一个全局来源扰动。本稿不重新记旧 graph 支付，不运行新数值。

## 1. 完整 μ、连续原窗口与全部中间尺度

μ 是完整有限正 Borel 来源，W>0，h_s=s^-n 1_[−s/2,s/2]^n，原全部连续尺度 a≤s≤b≤2a。M=sup_s h_s*μ，E={M>τ}，I=τ|E|，Z=1+n log(b/a)，receiver 都用 Lebesgue。

固定 E 上一个可测原真赢家 R(x)，另取原候选 L(x)≥R(x)，U_L=h_L*μ>0。令

\[
g(x)=\log(M/U_L),\quad \theta(x)=n\log(L/R),\quad
\pi_R=\mu|Q(x,R)/\mu(Q(x,R)),\quad
\pi_L=\mu|Q(x,L)/\mu(Q(x,L)).
\tag{1}
\]

嵌套与原 R 的 truewinner 性准确给

\[
\beta:=\frac{\mu(Q_R)}{\mu(Q_L)}=e^{g-\theta}\le1,
\qquad0\le g\le\theta.
\tag{2}
\]

对每一个原合法中间尺度 s∈[a,L]，完整源响应 U_s≤M，故

\[
\boxed{\pi_L\{y:2\|y-x\|_\infty\le s\}
=\frac{\mu(Q(x,s))}{\mu(Q(x,L))}
\le e^g(s/L)^n.}
\tag{3}
\]

这是由原 μ、原 winner 与全部中间尺度一起强制的后验径向帽，没有给来源坐标独立，也没有自由构造 posterior。closed cube 和边界原子均保留。界可能大于1，此时仅取它的有效部分。

## 2. 局部输运证书及精确积分

π_R 支撑在 Q_R。任意正概率耦合 κ∈Π(π_R,π_L) 都有

\[
\int\|y-z\|_\infty d\kappa
\ge\int\operatorname{dist}_\infty(z,Q_R)\,\pi_L(dz)
=\frac12\int_R^L\pi_L\{2\|z-x\|_\infty>s\}\,ds.
\tag{4}
\]

证明是逐耦合点 y∈Q_R 的距离下界及正 layer cake；不需要最优耦合存在。来源可奇异、坐标可相关；两个后验的有界 cube 支撑保证一阶矩有限。设 u0=L e^-g/n≥R，(3) 给

\[
\boxed{\Lambda(R,L,g)
:=\frac12\int_R^{u_0}[1-e^g(s/L)^n]ds
=\frac12\left[\frac n{n+1}u_0-R+
\frac R{n+1}e^{g-\theta}\right]\ge0.}
\tag{5}
\]

因此每一个 κ 的 max-norm 传输成本≥Λ。g=θ 时 u0=R、Λ=0；此时两个捕获后验相同。g=0 且 θ固定、n→∞、L固定时

\[
\Lambda=\frac L{2n}(\theta-1+e^{-\theta})+O_{\theta}(L/n^2).
\tag{6}
\]

这只是明确的距离量，不是 source budget，也不是弱型常数。

一个便于使用的粗形式：若 θ≥g+c，c>0，取
s=L exp[−(g+c/2)/n]，则 s−R≥a c/(2n)，π_L(Q_s^c)≥1−e^-c/2，因此

\[
\Lambda\ge\frac{a c}{4n}(1-e^{-c/2}).
\tag{7}
\]

局部有界测试 f_x(z)=min(dist∞(z,Q_R),(s−R)/2) 满足 π_R f_x=0、π_L f_x≥(s−R)(1−e^-c/2)/2。**f_x 随 receiver 变化，不能直接拿它作为全局来源扰动 v。** 下一节消除这个量词缺口，保留明确维数损失。

## 3. CDF 代价与原几何的联系

令 F_R,i、F_L,i 是两后验第 i 个坐标的 CDF，定义有限可测量

\[
D_i(x)=\int_\mathbb R|F_{R,i}(u)-F_{L,i}(u)|du,
\qquad D(x)=\sum_{i=1}^nD_i(x).
\tag{8}
\]

不必使用一维 W1 等距公式。直接按两边 tail 的 Fubini：π_R 的第 i 坐标支撑在 [x_i−R/2,x_i+R/2]，所以

\[
D_i\ge\int (|z_i-x_i|-R/2)_+\,\pi_L(dz).
\]

对 i求和，Σ_i(|z_i−x_i|−R/2)_+≥dist∞(z,Q_R)，由 (4)–(5) 得

\[
\boxed{D(x)\ge\Lambda(R(x),L(x),g(x)).}
\tag{9}
\]

这里只用原真 winner 径向质量帽。没有将一个任意 source pair 的距离当作两个不同 receiver 的共同捕获概率。

## 4. 与 x 无关的固定随机周期方波

固定 T=b。先选 i uniform于1,…,n，c uniform于[0,2b)，定义一个作用于整个空间的来源测试

\[
v_{i,c}(y)=
\begin{cases}+1,&(y_i-c)\bmod 2b\in[0,b),\\
-1,&(y_i-c)\bmod 2b\in[b,2b).
\end{cases}
\tag{10}
\]

这个 auxiliary law 不依赖 μ、x、R、L、g、source coordinate law 或 winner。每次抽取后的同一个 v用于所有 receiver；|v|=1，故是已核对称变分允许的完整空间 source test。

两个后验都支撑于同一 Q(x,L)，每个坐标支撑跨度≤L≤b；方波 jump之间距离 b。因此除 phase 零测集外，任一 coordinate support 内至多有一个 jump。对 σ_i=π_R,i−π_L,i，总质量零：若该 jump位于 u，便有
|∫v_{i,c}dσ_i|=2|F_R,i(u)−F_L,i(u)|；若无 jump则该积分为零。一个 phase period 内有两个 jump，uniform c下 jump 位置密度是1/b。正 Fubini给**准确身份**

\[
\boxed{\mathbb E_{i,c}|\pi_Rv_{i,c}-\pi_Lv_{i,c}|
=\frac2{nb}\sum_i D_i(x).}
\tag{11}
\]

source atoms落在某 phase jump 时可能改变选取的左右 CDF，但仅影响 c零测集，du积分不变。L=b时两个 jump同时在两端的 phase也为零测。不需要来源支撑整体有界，因为这里只对原捕获后验使用有界支撑；v本身在整个 R^n上定义。

这是 exact source-independent/random-global-test 转换。与 (7) 的局部 f_x不同，(11) 可以合法平均同一个完整来源变分上界。

## 5. 对称 nearmax 预算的无 atom-count 几何左端

沿用 nearflat_cube_geometry_20261007.md 自证的对称式：若 Φτ(μ)≥(K−ε)W，K是同一个完整窗口、全部完整正来源与阈值的 Φ/W 上确界，则对任意 |v|≤1、0<t<1，

\[
\tau\int_E(t|\pi_Rv-\pi_Lv|-g)_+dx
\le2\epsilon W+\frac{t^2I}{(1-t)^2}.
\tag{12}
\]

两新质量 W±=W±t∫vμ 之和是2W，所以 (12) 不要求这个随机 v保持质量；不能把它无条件套到仅在固定质量子类中的 nearmax。对 (12) 平均 (10)，用 hinge凸性、(11)、正 Tonelli和 (9)，得到

\[
\boxed{\tau\int_E\left(\frac{2tD(x)}{nb}-g(x)\right)_+dx
\le2\epsilon W+\frac{t^2I}{(1-t)^2},}
\tag{13}
\]
\[
\boxed{\tau\int_E\left(\frac{2t\Lambda(R,L,g)}{nb}-g\right)_+dx
\le2\epsilon W+\frac{t^2I}{(1-t)^2}.}
\tag{14}

全部式对 arbitrary finite μ有效；没有 atomcount、source-square、pairGaussian 或 source坐标独立假设。

两真赢家 g=0 时，令 m=I/W、Bε(t)=2ε/t+t m/(1−t)^2，则

\[
\tau\int_ED(x)dx\le\frac{nb}2 B_\epsilon(t)W.
\tag{15}
\]

特别，A={x∈E:两个已选真赢家满足 θ≥c} 上 (7) 给

\[
\boxed{\tau|A|\le
\frac{2n^2(b/a)}{c(1-e^{-c/2})} B_\epsilon(t)W.}
\tag{16}
\]

与先前 physical-atom collision 的 √N不同，(16) 的全部费用仅含显式 n、窗口比例与 ε。它给定量 a.e.唯一性的 arbitrary μ版本；仍只是 nearmax必要条件。若 0<ε≤Z/4取 t=√(ε/(4Z))，Bε≤6√(εZ)，则固定 c>0、b/a≤2的 tie 体积≤C_c n²√(εZ)W/τ。要跨维数使此量小，须显式选择 ε_n n^4 Z_n→0，不能把 ε固定。

近赢家也有量化。若 A上 θ≥g+c、g≤η，设
d_c=a c(1−e^-c/2)/(2n²b)。由 (7)、(13)得

\[
\tau|A|(td_c-\eta)_+
\le2\epsilon W+\frac{t^2I}{(1-t)^2}.
\tag{17}
\]

例如 η≤td_c/2 时，τ|A|≤2 Bε(t)W/d_c。其 gap 资格相当窄，不能免费认定原 nearwinner都落入该域；可直接保留无阈值的 (13)/(14)以避免隐藏这一限制。

## 6. 连续中间尺度不可被有限端点免费替代

本稿 (3) 的关键是每个 s∈[R,L] 都是原合法尺度。仅含有限目录 J 时，定义 s^+(u)=min{s∈J:s≥u}，对 u∈[R,L]只有

\[
\pi_L\{2\|y-x\|_\infty\le u\}
\le e^g(s^+(u)/L)^n.
\tag{18}
\]

于是替代 Λ的有效 finite-family下界为
Λ_J=½∫_R^L[1−e^g(s^+(u)/L)^n]_+du；它可能为0，不能冒充 (5)。全局方波身份 (11)及 (13)仍适用，因为仅需捕获后验支撑跨度≤b。

此前三轮 saved-cell probe仅验证 finite原目录的对称切换，不验证本稿 continuous中间尺度帽或 (5)/(16)。本轮没有重新执行那些数据，没有数值迁移声称。要新数值检验必须先登记新的完整 continuous arrival oracle或仅检验 finite替代式，交由父任务协调。

## 7. 实际完整来源显示本转换的 n^-2 信号损失

以下不是自由 posterior 或 weak常数反例；它是一个完整 finite正 μ，验证 (11) 这条具体全局化不能被无依据升级为 √n信号。

固定原完整窗口 [a,2a]，b=2a。给定 c>0和 n≥c/log2，令 L0=a exp(c/n)≤b。来源 Y按以下规则产生：I uniform于1,…,n，符号±uniform，r∈[0,L0]密度 n r^(n−1)/L0^n，Y=±(r/2)e_I。μ为这个完整概率law，支持2n条轴向射线，无来源坐标独立。

在 receiver x=0，准确有 μ(Q(0,s))=(s/L0)^n（s≤L0），其余 s≥L0 捕获全部质量。因此

\[
U_s(0)=L_0^{-n}\ (a\le s\le L_0),\qquad
U_s(0)=s^{-n}\ (L_0\le s\le2a).
\tag{19}
\]

这是对原全部中间尺度的真赢家陈述。R=a、L=L0均是真赢家，g=0、θ=c。取 τ<L0^-n使0∈E；响应在0附近连续（该 ray来源在相关截面无边界原子），故 E含正体积邻域，未用一个不合法 receiver离散模型。

π_R的 ray半径为 n r^(n−1)/R^n，π_L的为 n r^(n−1)/L0^n。每个 coordinate有共同原点质量1−1/n，余下1/n是正负对称半径。由一维CDF直接积分，

\[
D_i(0)=\frac{L_0-a}{2(n+1)},\qquad
D(0)=\frac{n(L_0-a)}{2(n+1)},
\]
\[
\boxed{\mathbb E_{i,c'}|\pi_Rv_{i,c'}-\pi_Lv_{i,c'}|
=\frac{L_0-a}{b(n+1)}\sim\frac c{2n^2}.}
\tag{20}
\]

这里 phase c'与固定 log-scale gap c不同。式(20)是严格实际来源的点态检测值；没有声称这个输入是 nearmax，也没有声称 positive-volume区域全部 exact ties。它只否定对这套固定周期坐标测试的 dimension-improved 点态检测下界。不能据此排除其它全局多坐标测试、更强近极大约束或最终平方根目标。

## 8. 年龄后验：root 结论的独立端点核验及 nearflat 能量删除

本节已完整读 root 的 posterior_depth_source_once_transfer_20261007.md，仅作为独立核验；年龄尾/MGF及 source-once一阶删除不重复记为本稿新成果。

记 B=n log(L/a)，D=2||x−Y||∞，A=n log[L/max(a,D)]。对 0≤t≤B，准确 {A≥t}={D≤L exp(−t/n)}，所以 (3) 给 PπL(A≥t)≤min(1,e^(g−t))。**t=B仍含D≤a的来源，可能有正质量；只能在t>B置零。t<0时尾为1。** B=0时A=0，不能称连续指数law。

在本稿 nested L≥R设置中 0≤g≤θ≤B，因此直接积分还给

\[
\mathbb E_{\pi_L}A\le g+1-e^{g-B},
\quad
\mathbb E_{\pi_L}e^{\sigma A}\le
\frac{e^{\sigma g}-\sigma e^{g-(1-\sigma)B}}{1-\sigma}
\le\frac{e^{\sigma g}}{1-\sigma},\quad0\le\sigma<1.
\tag{21}
\]

σ<0时最后这个式一般不成立；σ≥1不享有上述维数无关上界。root稿对一般正响应L也有效，但其g可能大于B；那时不能未经分支使用本稿(21)的truncated闭式，只能使用其一般尾积分或g+1/MGF粗界。

取L=原truewinnerR即g=0，在原selected joint
J(dx,dy)=τ 1_E(x)dx π_R,x(dy)，质量I。逐receiver条件尾与Tonelli给 J{A>t}≤e^-t I、E_{J/I}A≤1、E_{J/I}exp(σA)≤1/(1−σ)。这无需nearmax。这里A就是winner_envelope_entropy_bridge里的depth d；可把其P-mean上界改进为1，但原宽径向参考Q未变，source/angle independence也未产生。

以下相对二阶推论使用新的 nearflat约束，而非重新提出 all-input source-square合同。令
S(y)=J_Y/μ、Sd,t(y)=τ∫E 1_{Q_R}(y)1_{A>t}/μ(Q_R)dx，Δ=∫|S−K|dμ，I2=∫S²dμ。root的一阶尾给∫Sd,t dμ≤e^-t I。原包络 k* 与准确 h_R=e^-A k*(x−y)（在incidence上）同时给

\[
0\le S_{d,t}\le e^{-t}Z,\qquad
|I_2-KI|\le Z\Delta.
\]

因此从同一个完整来源正积分，

\[
\boxed{\int[S^2-(S-S_{d,t})^2]d\mu
\le2e^{-t}(KI+Z\Delta)
\le2e^{-t}I_2+4e^{-t}Z\Delta.}
\tag{22}
\]

精确flat时Δ=0，固定常数t即可删掉固定比例的深度source交叉能量；近flat只多保留显式ZΔ误差。root近优归约Δ≤3√(εZ)W给ZΔ≤3√ε Z^(3/2)W，需要按n显式选ε。此处没有付浅层source能量，没有凭未知K得到它的小维数阶，也不把此相对能量删除登记为完整weak预算。

floor不能删除：A≤t含R接近a时的整个floor区。若需要将保留incidence解释为真正外壳，可先保留原底端winner输出 {R≤a exp(t/n)}，其既有窄窗包络付τ体积≤(1+t)W（超出b时取b，费用更小）；其余receiver的A≤t确为D≥R exp(−t/n)的浅外壳，边界按≤/>一致选择。这个narrow-window envelope是旧工具，不重复记新paid。

## 9. 完成状态与剩余问题

局部径向运输→同一个 global source test 的量词缺口已由 (11)–(14) 严格消除；atomcount也已去除。代价明确为坐标平均1/n与窄尺度距离1/n，(20)说明该具体转换的损失确有实际来源资格。

要由这一路径控制未知 K，还需至少一个新的几何桥：大 K 的全局 nearmax完整输入必须有足够的、原合法候选尺度 nearwinner，并在允许的 gap下使 D的总积分产生不受 (20) 模型机制压低的信号；或者构造一个不同的 source-independent多坐标测试族，结合全部中间尺度帽和 nearmax自一致性改善检测预算。本稿未证明这种桥，亦未推断来源产品化、免费所有尺度或原 actual FIRST/history资格。

已证明工具是 arbitrary μ 的 (3)/(5)/(11)/(13)/(14)/(16)，以及 nearflat 下的相对深度能量删除(22)。年龄source-once结论已独立核验root稿。没有新增 sqrt(n)n^{o(1)} weak(1,1)结论，没有新的主账维数付款。
