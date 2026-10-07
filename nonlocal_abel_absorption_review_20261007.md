# 固定非局部算子的 Abel 吸收层：独立理论审核

日期：2026-10-07。范围：有限可逆跳链和有界总跳率的对称跳核，包括奇异跳测度作用于密度函数。本稿只证明固定算子的抽象引理，不认证原 cube 的时变 FIRST、winner、来源对或历史门。未运行或重跑数值。

**结论：候选成立。** 必须固定正算子的符号、包括完整外部吸收，并把 Abel 差的求和范围写成 \(0\le j<N\)。不要求 \(L,K\) 交换，也不要求 \(K\) 的严格谱隙。以下采用有界 \(V\) 的充分域假设。

## 1. 明确的充分假设和完整外部项

底空间为 \(\mathbb R^d\) 与 Lebesgue 测度，或有限状态空间与可逆测度 \(m_x>0\)。写底测度为 \(m\)。可测非负跳核 \(j(x,dy)\) 满足

\[
M(dx,dy):=m(dx)j(x,dy)=M(dy,dx),\quad
q(x)=j(x,X),\quad {\rm ess\,sup}\,q\le Q<\infty.
\]

令 \((\mathcal Jf)(x)=\int f(y)j(x,dy)\)、\(L=qI-\mathcal J\ge0\)。这里 \(L\) 是负的 Markov 生成元。\(L\) 在 \(L^2(m)\) 上有界自伴非负；\(e^{-tL}\)、\(A_s=s(s+L)^{-1}\) 保正且为 \(L^1,L^\infty\) 收缩。对称性和总率有界保证 \(L^1\) 中积分守恒。普通对称密度 \(J(x,y)dy\) 是特例；允许跳核没有环境空间 Lebesgue 密度。

给定 \(u\ge0\)、\(u\in L^1\cap L^2\)、\(\Omega=\{u>0\}\)，以及 \(\nu\ge0\)、\(\nu\in L^1\cap L^2\)、\(\nu=0\) 于 \(\Omega^c\)。假设

\[
\mu=\nu-Lu,\quad0\le\mu\le\kappa,\quad
\int\mu\,dm=\int\nu\,dm=W<\infty,\quad
V=\mu|_\Omega/u\in L^\infty(\Omega).
\tag{1}
\]

因此 \(\mu\in L^2\)、\(\|\mu\|_2^2\le\kappa W\)。\(E\) 为零延拓、\(P=E^*\)，
\(L_\Omega=PLE\)、\(K=L_\Omega+V\)、\(B_s=s(s+K)^{-1}\)。
\(K\) 有界、自伴、非负，半群保正；所有 \(s>0\) 的 resolvent 都存在，不要求 \(K^{-1}\) 存在。有限链的密度、积分、cap 都相对于 \(m\)，且在 \(u_x>0\) 的有限域上 \(V\) 自动有界。

定义从域内到**整个底空间**的正算子

\[
(Tz)(x)=
\begin{cases}
V(x)z(x),&x\in\Omega,\\
\displaystyle\int_\Omega z(y)j(x,dy),&x\notin\Omega.
\end{cases}
\tag{2}
\]

\(T:L^2(\Omega)\to L^2(X)\) 有界；外部项由双边总率的 Schur 界控制。不能只记局部“边界”或跳离概率。逐点计算给出
\[
LEz=EKz-Tz,\qquad Ku=\nu|_\Omega,\qquad Tu=\mu.
\tag{3}
\]
外部使用 \(u=\nu=0\)，故 \(\mu(x)=-Lu(x)=\int_\Omega u(y)j(x,dy)\)。这部分必须纳入完整源质量账。

## 2. 无谱隙的吸收恒等式

记 \(q_{\rm out}(y)=j(y,\Omega^c)\)。域内二次型为
\[
\langle f,Kf\rangle
=\tfrac12\int_{\Omega\times\Omega}|f(x)-f(y)|^2M(dx,dy)
+\int_\Omega(q_{\rm out}+V)|f|^2dm.
\tag{4}
\]
若 \(f\in\ker K\)，则 \(Vf=0\)，跨外域的 \(|f(y)|^2\) 积分也为零。对外部跳项应用 Cauchy–Schwarz 得 \(T\ker K=0\)。

令 \(g_t=e^{-tK}\nu\ge0\)。由 \(Ku=\nu\) 和谱定理，
\[
\int_0^b g_tdt=u-e^{-bK}u\uparrow u-P_{\ker K}u
\quad\text{于正序且于 }L^2,\qquad
\int_0^\infty Tg_tdt=\mu.
\tag{5}
\]
正序单调来自 \(g_t\ge0\)，不是假设任意半群轨道单调。末式由 \(T\) 有界、杀死 \(\ker K\) 及 Tonelli 得到。无吸收闭合分量可留在 \(\ker K\)，但不贡献 \(\nu\) 或 \(T\)。

## 3. 层、质量、导数及 Beta 恒等式

对 \(j\ge0\)，令 \(\rho_j(s)=T s^j(s+K)^{-j-1}\nu\)、\(p_j(r)=e^{-r}r^j/j!\)、\(F(t)=Tg_t\)。Laplace resolvent 公式与 (5) 给出
\[
\rho_j(s)=\int_0^\infty p_j(st)F(t)dt,\qquad
\sum_{j=0}^\infty\rho_j(s)=\mu.
\tag{6}
\]
求和成立于非负点态、\(L^1\) 和 \(L^2\)。\(p_0\le1\)；对 \(j\ge1\)，Stirling 下界给
\(\sup p_j\le(2\pi j)^{-1/2}\le(j+1)^{-1/2}\)。于是
\[
0\le\rho_j(s)\le(j+1)^{-1/2}\mu,\qquad
\int_0^\infty\|\rho_j(s)\|_1\frac{ds}{s}=\frac Wj\quad(j\ge1).
\tag{7}
\]
质量式用 \(\int_0^\infty p_j(r)dr/r=1/j\)。**不能用于 \(j=0\)**，其积分一般发散。固定层在 \(s>0\) 可微，且
\[
s\partial_s\rho_j=j\rho_j-(j+1)\rho_{j+1}.
\tag{8}
\]
对 \(1\le j<N\)，直接积分指数核得
\[
j\rho_j(s)=N\int_0^1\rho_N(s/\theta)
\frac{(N-1)!}{(j-1)!(N-j-1)!}
\theta^{j-1}(1-\theta)^{N-j-1}d\theta.
\tag{9}
\]
不需要 \(K\) 严格可逆。支配收敛给 \(\rho_0(s)\to\mu\) 当 \(s\downarrow0\)，\(\rho_0(s)\to0\) 当 \(s\to\infty\)，均在 \(L^2\)。

## 4. 全空间 Abel 差

把 \(\nu\) 理解为零延拓，设 \(D_N(s)=A_s^N\nu-E B_s^N\nu\)。由 (3) 的 resolvent 恒等式逐项消去，
\[
D_N(s)=\sum_{j=0}^{N-1}A_s^{N-j}\rho_j(s).
\tag{10}
\]
这是有限和。用 \(sA_s'=A_s(I-A_s)\) 和 (8)，中间项抵消，得
\[
sD_N'(s)=
N\sum_{j=0}^{N-1}A_s^{N-j}(I-A_s)\rho_j(s)-N A_s\rho_N(s).
\tag{11}
\]
在 \(\Omega^c\)，\(D_N=A_s^N\nu\)，外部响应没有被省略。

## 5. 完整的 abstract spectral square-function 核验

只用固定 \(L\) 的自伴非负性、Abel 保正 \(L^\infty\) 收缩及 (6)–(9)，可得
\[
\int_0^\infty\|sD_N'(s)\|_2^2\frac{ds}{s}
\le276\sqrt N\,\kappa W\le280\sqrt N\,\kappa W,\qquad N\ge1.
\tag{12}
\]
证明组织参考 [Spector–Stockdale, arXiv:2609.05377v1, §3.4](https://arxiv.org/html/2609.05377v1)。以下独立列出迁移所需标量计算；原文球极大函数结论不据此转用于 cube。late Schur 是该节最后的标量核步骤，不采用先前通信误写的“式 (3.24)”编号。

由 (7)，对 \(j\ge1\)，
\[
\mathcal B_j:=\int_0^\infty\|\rho_j(s)\|_2^2\frac{ds}{s}
\le\kappa W\,j^{-1}(j+1)^{-1/2}\le\kappa W j^{-3/2}.
\tag{13}
\]

**零层。** 设 \(E_N(s)=\langle\rho_0,A_s^N\rho_0\rangle\)。则
\[
sE_N'=N\langle\rho_0,A_s^N(I-A_s)\rho_0\rangle
-2\langle\rho_1,A_s^N\rho_0\rangle.
\]
\(E_N(\infty)=0\)，但一般 \(E_N(0)=\langle\mu,P_{\ker L}\mu\rangle\ge0\)。积分时该端点以负号出现；由 \(A_s^N\rho_0\le\kappa\) 和 (7)，
\[
N\int\langle\rho_0,A_s^N(I-A_s)\rho_0\rangle\,ds/s\le2\kappa W.
\]
标量 \(0\le N a^N(1-a)\le1\) 的平方不超过自身。谱定理给
\[
\int\|N A_s^N(I-A_s)\rho_0\|_2^2\,ds/s\le2\kappa W.
\tag{14}
\]

**早层。** 令 \(J_N=\lceil N/2\rceil\)。对 \(1\le j<J_N\)，
\(\|N A_s^{N-j}(I-A_s)\|_{2\to2}\le2\)。
Minkowski、(13) 与 \(\sum_{j=1}^M j^{-3/4}\le4M^{1/4}\) 给早层范数至多 \(8N^{1/4}\sqrt{\kappa W}\)，平方至多 \(64\sqrt N\kappa W\)。

**晚层。** 在 \(L\) 的固定谱纤维 \(\lambda\ge0\) 上记 \(a_s=s/(s+\lambda)\)。由 (9)，晚层等于对 \(N\rho_N(t)\) 作用的非负标量核
\[
H_{N,\lambda}(s,t)=N a_s(1-a_s)
\sum_{j=J_N}^{N-1}\binom{N-1}{j}(s/t)^j
[a_s(1-s/t)]^{N-1-j}\mathbf1_{t>s}.
\tag{15}
\]
Beta 积分给行界
\[
\int_s^\infty H(s,t)\,dt/t
=N a_s(1-a_s)\sum_{j=J_N}^{N-1}a_s^{N-1-j}/j\le2.
\tag{16}
\]
列积分固定 \(t\)，用 \(\theta=s/t\)、
\(\theta+a_s(1-\theta)=a_s/a_t\) 和
\(a_s(1-a_s)ds/s=da_s\)，把 binomial 部分和扩大为全和，得
\[
\int_0^t H(s,t)\,ds/s
\le N\int_0^{a_t}(a/a_t)^{N-1}da=a_t\le1.
\tag{17}
\]
\(\lambda=0\) 时直接 \(H=0\)，不作除零运算。Schur 范数至多 \(\sqrt2\)，对谱纤维中的向量同样成立。对纤维积分得晚层平方范数至多 \(2N^2\mathcal B_N\le2\sqrt N\kappa W\)。

**终层。** \(A_s\) 为 \(L^2\) 收缩，故平方范数至多 \(N^2\mathcal B_N\le\sqrt N\kappa W\)。四块相加再平方，用 \(\|\sum_{i=1}^4f_i\|^2\le4\sum_i\|f_i\|^2\)，得 \(4(2+64+2+1)\sqrt N\kappa W=276\sqrt N\kappa W\)。\(N=1\) 的早、晚层为空，也满足该界。全程未把 \(K,L\) 同时对角化。

## 6. 域边界与真正接入缺口

1. 充分域是对称有界率跳核/跳测度、\(u,\nu\in L^1\cap L^2\)、有界 \(V\)。无 gap 已证明。无界 \(V\) 需闭形式与 \(T\) 的域/正序积分核查，不是默认范围；无界总跳率或无限活跃奇异核也未覆盖。
2. 允许奇异跳测度作用于密度，不等于允许 Lebesgue 底测度下的原子 \(\nu\)：普通函数 \(Lu\) 不能把输入原子消去成有界密度。后者需正则化或另一套 measure-valued 域。
3. 此处给定 \(u,\mu\)，没有为任意输入构造 capped 分解。仅 \(\mu\le\kappa\) 不能推出 \(|\Omega|\le W/\kappa\)；若需 bad-set 体积账，必须另有域内饱和/障碍问题性质。外部吸收也可有正质量。
4. 只读旧 tex 4878–4913 确认：其固定 \((s,L)\)、\(s\le S<1\) 的轴向正概率卷积表示若成立，\(\mathsf A_{s,L}=(1-s)^{-1}\sum_i(I-G_{1-s,L}^{(i)})\) 属于这里的有界率对称跳测度类型，率为 \(n/(1-s)\)。本稿未重证旧从属过程构造，也未为 frozen 算子构造 (1)。原过程随 \((s,L)\) 改变：逐时有界/自伴不足以使用固定谱 Schur；各时刻重做 obstacle 不能反复领同一 \(W\)。
5. 还需回代原实际目标：保留 FIRST、winner、CP/GP、来源分配及全部原门，控制高频/最大响应，并核主账系数。平方函数估计自身不是原空间/交通费用结论。

## 7. root 三轮 guard 的只读审查

全文只读了 nonlocal_abel_exact_numeric_guard_20261007.py：states 5/9/17，正来源和完整内外 cap、守恒、\(Ku=\nu\)、\(Tu=\mu\)，36 组有理 \((s,N)\) 的 (10)/(11) 与删 outside 的严格缺陷；还检查九层和的精确余项。有限零谱采用连通图常数态，只将数值极小特征值重置为 0。

数值全平方函数有 15 组 \(N=1,2,4,8,16\)，在 \(s\in[2^{-24},2^{24}]\) 的 16 个对数 panel 上 GL 加密。这是正数值交叉检查，**不是 interval 认证**；解析省略尾界另记。独立核尾界：完全图边权下界 \(c\) 给自由非零谱 gap \(\ge mc\)，\(K\ge(\min V)I\)，\(\|L\|\le2\max q\)、\(\|K\|\le2\max q+\max V\)。小 \(s\) 导数范数至多 \(Ns\|\nu\|_2(g_L^{-1}+g_K^{-1})\)，大 \(s\) 至多 \(N\|\nu\|_2(\Lambda_L+\Lambda_K)/s\)。积分后两条尾界分别为相应平方常数除以 \(2^{49}\)，与代码一致；自由零谱在导数中贡献 0。安全常数 280 来自解析 276 界，没有拟合。

未复跑 guard，也未独立重构全部保存矩阵哈希。结果中的 PASS/438 与加密差是作者运行收据；定理来自 §§1–5 的证明。

**终态：抽象一般引理及 \(\sqrt N\kappa W\) transfer 接受；原时变 cube 实际接入未证明。本稿不添加 cube 已付分支，不重复 coin 账。**
