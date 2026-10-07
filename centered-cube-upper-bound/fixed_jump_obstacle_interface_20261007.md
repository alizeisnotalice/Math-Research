# 任意密度输入的固定 jump 障碍接口

2026-10-07。只新增本稿及同前缀 registration/guard/results；不修改总稿、主账、旧 Abel 收据或他人代码。

**结论：一般“对称有界率”假设不够。** 有限可逆链的准确条件是每个封闭连通分量的质量不超过其 cap 容量；\(\mathbb R^n\) 还需有限占用的稳定化条件。下面证明：有正 \(L^1\cap L^2\) supersolution 时存在最小饱和障碍；固定的实际轴向卷积生成元可用有限均值 coupon-collector 停时构造该 supersolution。域外 good 来源必须先拆出。随后可合法去掉旧 Abel 工具的有界 \(V=\kappa/u\) 假设。以上均是单个固定算子的结果，未闭合时变 cube。

## 1. 查重与证据范围

已读 [J03 SKILL](/Users/zhengzhihao/.codex/skills/math-j03-nonlocal-balayage-exit-measure/SKILL.md) 及 provenance/method/cube-interface；它要求保留整个外部落点和退出/杀死/无穷停时质量，未提供这里的任意输入障碍存在性定理。

旧 nonlocal_balayage_exit_route.md 是有序软度、来源小格与共享早退出 measure，不是任意输入的 stationary obstacle。旧总 tex 8353 后已证明：不同 cube 尺度不能共用非平凡饱和帽；本稿不重复该结论。外部 [Spector–Stockdale §2](https://arxiv.org/html/2609.05377v1) 的构造针对 Laplacian 和其势论，不能直接替换成 bounded-rate jump。以下结论由本稿证明。

沿用 nonlocal_abel_absorption_review_20261007.md 的正算子符号
\[
L=qI-\mathcal J,\qquad
\mathcal Jf(x)=\int f(y)j(x,dy),\qquad
m(dx)j(x,dy)=m(dy)j(y,dx),\qquad q(x)=j(x,X)\le Q.
\tag{1}
\]
底测度为有限可逆链的 \(m\) 或 Lebesgue。固定 \(\kappa>0\)、\(\nu\ge0\)、\(\nu\in L^1\cap L^2\)、\(W=\int\nu\)。

## 2. 有限链：准确容量、存在性和临界常数模式

设 \(C\) 是有限封闭连通分量。必要条件为
\[
W_C:=\int_C\nu\,dm\le\kappa m(C).
\tag{2}
\]
因为 \(\int_C Lu=0\)，任何 \(0\le\nu-Lu\le\kappa\) 都满足 (2)。仅检查全链的总容量不够。

**严格容量。** 若 \(W_C<\kappa m(C)\)，在 \(u\ge0\) 上最小化
\[
\mathcal F_C(u)=\tfrac12\langle u,Lu\rangle_C+
\langle\kappa-\nu,u\rangle_C.
\tag{3}
\]
写 \(u=c\mathbf1+v\)，\(v\) 均值零、\(c\ge0\)。连通性给零均值谱隙 \(g_C>0\)，而常数方向的系数是
\(\delta_C=\kappa m(C)-W_C>0\)。故
\(\mathcal F_C\ge (g_C/4)\|v\|_2^2+\delta_Cc-O_C(1)\)，在正锥上 coercive，最小值存在。KKT 为
\[
Lu+\kappa-\nu\ge0,\qquad
u(Lu+\kappa-\nu)=0.
\tag{4}
\]
于是 \(\mu_{\rm tot}=\nu-Lu\le\kappa\)，且在 \(\Omega=\{u>0\}\) 上 \(\mu_{\rm tot}=\kappa\)。在 \(u=0\) 的状态，
\(\mu_{\rm tot}=\nu+\mathcal Ju\ge0\)。因此完整目标非负、cap、守恒。

严格容量下最小点唯一：两个最小点之差在 \(\ker L\)，因而在 \(C\) 为常数；每个最小点至少有一个零状态，否则沿减常数方向降低 (3)。两者非负且都取到最小值零，常数差只能为零。

**临界容量。** 若 \(W_C=\kappa m(C)\)，\(\nu-\kappa\) 均值零，先解
\(Lv=\nu-\kappa\)，再令 \(u=v-\min_Cv\ge0\)。则 \(\mu_{\rm tot}=\kappa\) 于整个 \(C\)，满足目标合同。它是非负 Poisson 解中的最小代表。仍可加非负常数而失去唯一性；这里选择 \(\min u=0\)。所以临界情况不是“不存在”，但不能使用严格容量的 coercivity/唯一性证明。

逐分量应用即可。有限链上 \(\kappa/u\) 在正集自动有界；其 bound 可以依输入而变，不能将它视为统一常数。

### 2.1 可核算法

严格分量从 \(A=\{\nu>\kappa\}\) 开始，解
\(L_{AA}u_A=\nu_A-\kappa\)、\(u_{A^c}=0\)。每轮加入所有域外 \(\nu-Lu>\kappa\) 的状态再解。初始解正；每次新增状态的 residual 为正，Dirichlet \(M\)-matrix 的非负逆使新解逐坐标增加，新状态也正。严格容量保证 active set 不会等于整个分量：否则此前内态饱和、所有余态均超 cap，违反质量恒等式。每轮严格扩张，有限步终止，满足 (4)。这也给出最小非负 supersolution，而非某个任意 stationary 解。

## 3. \(\mathbb R^n\)：仅有界率为何不够

最简单的失败是 \(L=0\) 且 \(\nu>\kappa\)。即使令率处处为 1 也不够：用不相交单位格 \(C_z\)，令 \(P\) 在每格内均匀平均、\(L=I-P\)。取 \(\nu=2\kappa\mathbf1_{C_0}\)；封闭格 \(C_0\) 的质量超过容量，虽全空间体积无限也不可能。

**不可约、紧支撑有界输入也不足保证 \(L^1\) odometer。** 取不交 \(E_k\) 位于一个有界盒内，\(|E_k|=2^{-3k}\)，选严格正、有界、可积 \(a\)，使 \(a=2^{-4k}\) 于 \(E_k\)。令
\[
J(x,y)=a(x)a(y),\qquad A=\int a,\qquad q(x)=Aa(x),
\qquad \nu=2\kappa\sum_{k\ge1}\mathbf1_{E_k}.
\tag{5}
\]
核对称、总率有界，对每个正测度集合均有正跳概率；\(\nu\) 紧支撑且属于 \(L^1\cap L^2\cap L^\infty\)。但任何非负 \(u\) 若满足 \(\nu-Lu\le\kappa\)，必有
\[
q(x)u(x)=\nu(x)-\mu(x)+\mathcal Ju(x)\ge\nu(x)-\kappa,
\quad
\int u\ge(\kappa/A)\sum_{k\ge1}2^k=\infty.
\tag{6}
\]
因此连不可约性也不能补足工具要求的有限占用域。必须排除无限 holding-time 代价，不能把无限容量或“会跳走”当作存在性证明。

## 4. 正 supersolution 的充分定理与完整证明

假设存在 \(w\ge0\)、\(w\in L^1\cap L^2\)，满足
\[
Lw\ge\nu-\kappa\quad\text{几乎处处}.
\tag{7}
\]
这里常数 \(\kappa\) 不必可积，(7) 是点态不等式。定义 \(u_0=0\)，在 \(q>0\) 上递推
\[
u_{k+1}=\left(\frac{\nu+\mathcal Ju_k-\kappa}{q}\right)_+;
\tag{8}
\]
在 \(q=0\) 上令 \(u_k=0\)。此处 (7) 自动要求 \(\nu\le\kappa\)，且该处跳核零。

算子 (8) 保序，且由 (7) 映 \(w\) 到不超过 \(w\)。故
\(0\le u_k\uparrow u\le w\)，并于 \(L^1,L^2\) 收敛。核的单调收敛给固定点 (8)。因此
\[
\mu_{\rm tot}:=\nu-Lu\in[0,\kappa],\qquad
\mu_{\rm tot}=\kappa\ \text{于 }\Omega=\{u>0\}.
\tag{9}
\]
在域外 \(\mu_{\rm tot}=\nu+\mathcal Ju\ge0\)，并由固定点不等式不超过 \(\kappa\)。
\(\nu,Lu\in L^1\)，对称有界率使 \(\int Lu=0\)，于是
\[
\int\mu_{\rm tot}=W,\qquad \kappa|\Omega|\le W.
\tag{10}
\]
任意非负 supersolution 都逐轮支配 \(u_k\)，所以 \(u\) 为最小非负 supersolution。这里不用 \(L^{-1}\)，没有将全空间谱的零点删除。

### 4.1 一个显式有限步 barrier

若 \(L=q_0(I-P)\)，\(q_0>0\) 固定，且某有限 \(r\) 满足 \(P^r\nu\le\kappa\)，则
\[
w=q_0^{-1}\sum_{k=0}^{r-1}P^k\nu,\qquad
Lw=\nu-P^r\nu\ge\nu-\kappa.
\tag{11}
\]
正 \(L^1,L^2\) 收缩使 \(w\) 有限。若 \(P\) 为有界概率密度 \(p\) 的卷积且 \(\|p^{*r}\|_\infty\to0\)，任意当前 \(\nu\) 都可用 \(W\|p^{*r}\|_\infty\le\kappa\) 选择 \(r\)。

该衰减有简短充分证明：\(p\in L^1\cap L^\infty\subset L^2\)，其特征函数在非零频率严格小于 1（密度不能支撑于一个相位恒定的超平面族）。Plancherel、反演及对 \(|\widehat p|^{2r}\) 的支配收敛给 \(\|p^{*2r}\|_\infty\to0\)。这不适用于有原子或始终缺坐标的有限步核。

## 5. 固定轴向卷积：完整停时平滑的 barrier

固定 \(c=1-s>0\)、物理尺度 \(\ell>0\)。沿旧正核表示，令 \(G_{c,\ell}\) 为偶、正的一维概率密度，且 \(G_{c,\ell}\in L^\infty\)。固定算子是
\[
L_{\rm fr}=\frac1c\sum_{i=1}^n(I-P_i)=q_0(I-P),\qquad
q_0=n/c,\quad P=n^{-1}\sum_iP_i,
\tag{12}
\]
其中 \(P_i\) 只在坐标 \(i\) 卷积 \(G_{c,\ell}\)。有限步 \(P^r\) 仍有漏跳坐标的奇异部分，不能直接称它为全维 bounded density。

令跳标签独立均匀于 \(1,\ldots,n\)，各跳位移独立使用对应 \(G\)。令 \(\tau_r\) 为每个坐标至少跳 \(r\) 次的第一离散时刻。逐轮收集全部坐标的 \(r\) 个 coupon-collector 轮次支配它，故
\[
\mathbb E\tau_r\le rnH_n<\infty.
\tag{13}
\]
终核 \(S_r\) 条件于标签计数 \(N_i\ge r\) 后，其密度为
\(\prod_iG^{*N_i}\)。卷积额外概率密度只会降低 \(L^\infty\) norm，因此
\[
\|S_r\|_\infty\le\|G^{*r}\|_\infty^n\longrightarrow0.
\tag{14}
\]
它是总质量 1 的完整终核；没有将漏坐标路径删掉。选择有限 \(r\) 使 \(W\|G^{*r}\|_\infty^n\le\kappa\)，则 \(S_r\nu\le\kappa\)。

令 \(\eta_k\) 为第 \(k\) 跳之后仍未停的未归一化来源密度（保留原初始来源），定义
\[
w=q_0^{-1}\sum_{k\ge0}\eta_k.
\tag{15}
\]
每个 \(\eta_k\) 是总核质量 \(\mathbb P(\tau_r>k)\) 的正平移混合，所以
\[
\|w\|_1=(\mathbb E\tau_r)W/q_0,\qquad
\|w\|_2\le(\mathbb E\tau_r)\|\nu\|_2/q_0.
\tag{16}
\]
下一跳标签独立，故 \(P\eta_k=\eta_{k+1}+\zeta_{k+1}\)，其中 \(\zeta\) 是该步终止密度。对 \(k\) 求和、尾部 \(L^1,L^2\) 消失，得到
\[
L_{\rm fr}w=\nu-\sum_{k\ge0}\zeta_{k+1}
=\nu-S_r\nu\ge\nu-\kappa.
\tag{17}
\]
这提供 §4 所需的真实有限占用 barrier，适用于任意当前非负 \(L^1\cap L^2\) 密度。

旧 frozen 正性表示未在此重证。其一维 Fourier 密度
\(\widehat G_c(\xi)=1/[1+c(\xi^2/\log(1+\xi^2)-1)]\) 在固定 \(c>0\) 时可积（大频率为 \(O_c(\log|\xi|/\xi^2)\)），所以已接受的正概率表示确有 bounded density；由 §4.1 得卷积幂峰值衰减。尺度伸缩不破坏这些固定参数性质。

(16) 的 \(r,H_n,c\) 只用于证明 odometer 有限；没有把占用总量作为新主账费用。它们没有提供不同 \((s,\ell)\) 的共同 \(u,\Omega,\mu\)，更没有共同 FIRST/history 门。

## 6. 不能漏掉域外 good 输入

一般障碍的总 \(\nu\) 未必支撑于 \(\Omega\)。必须一次性拆
\[
\nu_{\rm bad}=\nu\mathbf1_\Omega,\quad
\nu_{\rm good}=\nu\mathbf1_{\Omega^c},\quad
\mu_{\rm bad}=\nu_{\rm bad}-Lu
=\mu_{\rm tot}-\nu_{\rm good}.
\tag{18}
\]
在域外 \(\mu_{\rm bad}=\mathcal Ju\ge0\)，在域内为 \(\kappa\)。于是
\[
0\le\nu_{\rm good}\le\mu_{\rm tot}\le\kappa,\quad
0\le\mu_{\rm bad}\le\mu_{\rm tot}\le\kappa,\quad
\int\mu_{\rm bad}=\int\nu_{\rm bad}=W_{\rm bad},
\quad \kappa|\Omega|\le W_{\rm bad}\le W.
\tag{19}
\]
质量式使用 \(u\in L^1\)。并且 \(W_{\rm bad}+W_{\rm good}=W\)，源只拆一次。旧 Abel 的 \(Ku=\nu\)、\(Tu=\mu\) 应当分别替换为 \(\nu_{\rm bad}\)、\(\mu_{\rm bad}\)。good 来源的任意固定 free Abel/半群响应都 \(\le\kappa\)。

## 7. 合法去掉 \(V=\kappa/u\) 有界假设

在 \(\Omega\) 上设 \(V=\kappa/u>0\)，通常无界。令
\[
K=L_\Omega+V,\qquad D(K)=D(V)=\{f\in L^2(\Omega):Vf\in L^2\}.
\tag{20}
\]
乘法 \(V\) 自伴非负，\(L_\Omega\) 为有界自伴非负扰动；亦可用闭非负形式定义同一个 \(K\)。Trotter 给其保正、sub-Markov 半群。由 (19)，
\[
Vu=\kappa\mathbf1_\Omega\in L^2,\qquad
u\in D(K),\qquad Ku=\nu_{\rm bad}|_\Omega.
\tag{21}
\]
这里 \(K\) 可以无 gap，但 \(\ker K=0\)，因为 \(V>0\) 几乎处处。零集外无须除以 \(u\)。

\(T\) 仍用整个外部项定义。它不再是全局 bounded \(L^2\) 算子，但从 \(D(K)\) 的 graph norm 到 \(L^2\) 有界；且在正序区间 \(0\le z\le u\) 上，
\[
0\le Tz\le Tu=\mu_{\rm bad}\in L^1\cap L^2.
\tag{22}
\]
对 \(g_t=e^{-tK}\nu_{\rm bad}\)，\(t>0\) 有 \(g_t\in D(K)\)、\(g_t\ge0\)。谱定理给
\(\int_0^b g_tdt=u-e^{-bK}u\uparrow u\) 于正序和 \(L^2\)。
乘法及跳核项的 Tonelli 于是给
\[
\int_0^\infty Tg_tdt=Tu=\mu_{\rm bad}.
\tag{23}
\]
不使用“无界 \(T\) 与 Bochner 积分可任意交换”。

设 \(z_j=s^j(s+K)^{-j-1}\nu_{\rm bad}\)。它在 \(D(K)\)，且
\[
z_0=u-B_su\ge0,\qquad z_0\le u,\qquad
z_j=B_s^jz_0\le B_s^ju\le u.
\tag{24}
\]
最后一条迭代使用 \(B_su\le u\)，来自 \(Ku=\nu_{\rm bad}\ge0\)。因此层 \(Tz_j\) 定义良好。(23) 的正积分重新给 Poisson 表示；所有 cap、质量、Beta、导数和 resolvent 泄漏等式均成立。导数可在 graph domain 内核查，或直接用正积分。原 abstract scalar Schur 只作用于固定 free \(L\)，故旧
\[
\int\|s\partial_s(A_s^N\nu_{\rm bad}-EB_s^N\nu_{\rm bad})\|_2^2\,ds/s
\le276\sqrt N\,\kappa W_{\rm bad}
\tag{25}
\]
继续成立，不受 \(V\) 的最大值影响。

若要把总输入写入全空间平方函数，可另加 free good 项。谱变量代换给
\[
\int_0^\infty[N a_s^N(1-a_s)]^2\,ds/s
=N/[2(2N+1)]\le1/4
\]
（free 零谱导数为 0）。因此 good 的能量 \(\le\kappa W_{\rm good}/4\)。两项合并可用 \(552\sqrt N\kappa W\)；外部的原总 \(A_s^N\nu\) 同样被覆盖。该可选强化仍不是时变空间费用。

## 8. 逃逸、不可逆与原子输入的费用

- 全空间 \(L\) 一般没有有界逆；上述证明只用有限占用 barrier 和 killed resolvent，不假设 \(L^{-1}\)。过程 transient 不保证全 Green occupation 属 \(L^1\)；recurrence 也不排斥 §5 的有限停止占用。若只得到非 \(L^1\) 的 \(u\)，不能将形式上的 \(\int Lu=0\) 用作无穷远质量守恒。
- 在保守 bounded-rate 模型中，§5 的停时有限且均值有限，完整终核质量 1，无逃逸残量。对有杀死/爆炸或无限停时质量的替代模型，终核缺质量必须单列，不能沿用 (17)/(19) 的守恒。
- 如“不可逆”指非可逆链而非 \(L\) 不可逆，正 active-set/toppling 可能另有扩展，但本稿未证明；fixed selfadjoint spectral Schur 更不能直接套非自伴算子。
- Lebesgue 底测度下，原子输入不属于本定理。普通函数 \(Lu\) 不能消去输入原子成为 bounded density。可先取 \(\nu_\varepsilon=\beta_\varepsilon*\nu\)，正概率 smooth \(\beta_\varepsilon\)：质量仍为 \(W\)，无额外源质量倍数，\(\|\nu_\varepsilon\|_2\le W\|\beta_\varepsilon\|_2\) 随 \(\varepsilon\) 可发散。障碍、体积及 (25) 的最终常数均不依该 \(L^2\) norm，但极限中的域、响应与 actual FIRST/gates 仍需独立稳定性论证。平滑不是免费正支配原输入，也没有自动给出原门的误差界。

## 9. 新三轮精确压力收据

预注册 fixed_jump_obstacle_interface_registration_20261007.json 后执行同前缀 guard.py，states 7/11/19，所有算术为 Fraction。每轮 10 个预设 profile：6 个严格容量、2 个临界容量、全局超容量拒绝、全局容量足够但封闭分量超容量拒绝。检查 active 解正性/扩张单调、目标 cap、域内饱和、完整质量、域外实际吸收、(18)–(19) 和域体积由 bad 质量支付。

结果：PASS_EXACT，348 个谓词，三轮均通过。实际正外部吸收被强制检查；临界常数模式没有被丢弃。结果保存 script/registration SHA、各 profile 的 \(u,\mu\) SHA 与 active trace。没有重跑上一 Abel guard，也没有将有限链测试称作 \(\mathbb R^n\) 存在性证明；连续充分结论来自 §§4–7。

**终态：固定轴向生成元的任意 \(L^1\cap L^2\) 密度障碍接口已给充分证明；更广任意对称有界率存在性不成立。共同 \((s,\ell)\) 分解、原实际 FIRST/history 回代及一般空间账仍未建立，不向主账加费。**
