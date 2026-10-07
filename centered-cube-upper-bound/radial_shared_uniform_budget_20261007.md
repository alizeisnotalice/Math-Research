# 共享 uniform 的径向占用：精确身份与真实中心硬核的小 count 障碍

2026-10-07。只新增本稿及同前缀证书，不改主账。读取 logistic_chain_actual_interface 的共同 U fullsoft场、acceptance_count_projection 与原 source-cone身份。本稿保留同一输入、真实 finite winner及 threshold normalization；不将 σ=0模型写成全部实际 FIRST/fullfuture 门反例。

结论有两层。共享 U 的径向 mean 是原 weak体积的精确换元，单独没有新付款。更强的共同 U Palm小 count合同，即使在真正中心硬核、完整超水平域与 L1输入中，也因罕见同时激活而失败；问题不来自任意 sensor或任意 receiver gate。原 weak体积本身在该输入中有常数预算，故这是该 coupling尾合同过强，不是原 weak目标反例。

## 1. 纯 hard 的径向身份，完整来源只一次

先取 finite尺度 \(L_j=ae^{\beta_j/n}\in[a,b]\)，\(j=0,\ldots,J\)，
\[
 F_j=h_{L_j}^{(n)}*\mu,\quad M=\max_jF_j,\quad
 E=\{M>\tau\},\quad
 g_j(x)=\frac{\tau}{M(x)}\mathbf1_{A_j}(x),
\tag{1}
\]
其中 A_j 为 E上真实最小索引 winner不交分区，μ≥0、W=μ全质量。若指定外域，所有 A_j取同一外域交集。来源 y∼μ/W抽一次、共享 \(U\sim\mathrm{Unif}[-1/2,1/2]^n\)，令 \(Y_j=y+L_jU\)。给定(y,U)，独立 coin以 \(g_j(Y_j)\) 接受，计 N。正 Fubini仍给
\[
 \mathbb EN=W^{-1}\sum_j\int g_jF_j=\tau|E|/W.
\tag{2}
\]
这个 mean与 independent sensors版本相同，其 union/尾部不同。

将 U精确写成 \(U=\tfrac12e^{-d/n}\omega\)，\(\omega\)为原 cube-cone概率（均匀面、符号与其余面坐标），\(d\sim\mathrm{Exp}(1)\)，独立。除U=0及面tie零测外无例外。记 \(\rho_t=\tfrac a2e^{t/n}\)。于是
\[
 \begin{split}
 \mathbb EN
 &=W^{-1}\int d\mu(y)\int d\varsigma_n(\omega)
       \int_0^\infty e^{-d}\sum_jg_j(y+\rho_{\beta_j-d}\omega)\,dd\\
 &=\int_{-\infty}^{B}\Psi(t)\,dt,\quad B=\max_j\beta_j,\\
 \Psi(t)
 &=W^{-1}\int d\mu(y)\int d\varsigma_n(\omega)
 \sum_j e^{t-\beta_j}\mathbf1_{t\le\beta_j}
          g_j(y+\rho_t\omega).
 \end{split}
\tag{3}
\]
每个固定(y,ω,t)的 x=y+ρ_tω只有一个winner，所以求和等于
\[
 e^{t-\beta(x)}\frac{\tau}{M(x)}
 \mathbf1_E(x)\mathbf1_{t\le\beta(x)}.
\tag{4}
\]
不存在网格个数 J费用，但 \(\int\Psi=\tau|E|/W\)恰是目标；宣布其为“原质量加权占用”本身不是估计。粗界 \(\int_{t<0}\Psi\le1,\ \int_{0}^{B}\Psi\le B\)只是旧一源 cone预算。

连续 hard winner（允许全部[a,b]）还能使用原已知 inner捕获比 \(m_{\mathrm{inner}}/m\le e^{-v}\) 与低β支付；finite网不含 \(Re^{-v/n}\)时不能套。此处不把旧短β/长delay再收一次费。fullfuture在C=0沿本曲线只重复 hard帽；其所有非零softness测试不自动控制给定(y,U)的接受次数。特别 (3) 没有把 receiver帽转换成一次命中后的条件财富。

## 2. 完整超水平域的真实 finite hard模型

取任意 n≥1、J≥2，
\[
 q=1+\frac1{2J},\quad L_j=q^j\ (0\le j\le J),\quad
 \mu=\delta_0,\quad
 \tau=\frac12L_J^{-n}.
\tag{5}
\]
所有尺度在[1,2]：\(q^J<e^{1/2}<2\)。真实响应为 \(F_j=L_j^{-n}\mathbf1_{Q_{L_j}}\)。完整 E正是 \(Q_{L_J}\)（边界零测不影响），不是只截outer annuli：
\[
 A_0=Q_1,\qquad A_j=Q_{L_j}\setminus Q_{L_{j-1}},\qquad
 g_j=\tau L_j^n\mathbf1_{A_j},\quad j\ge0.
\tag{6}
\]
整个E严格超过τ，所有g_j≤1/2。取
\[
 H=\{2\|U\|_\infty>1/q\},\quad
 p=\Pr(H)=1-q^{-n}.
\]
对所有j≥1，\(Y_j=L_jU\in A_j\)当且仅当H；\(Y_0=U\in A_0\)始终。因而准确
\[
 N=B_0+\mathbf1_H S,\quad
 B_0\sim\mathrm{Bern}(\tau),\quad
 S=\sum_{j=1}^J B_j,\quad B_j\sim\mathrm{Bern}(\tau L_j^n),
\tag{7}
\]
全部coins彼此及H独立。S为Poisson-binomial，不误写为相同参数的Binomial。
\[
 \mathbb EN=\tau+p\tau\sum_{j=1}^JL_j^n
 =\tau L_J^n=\frac12,\quad
 s:=\Pr(N>0)=\tau+(1-\tau)p\left[1-\prod_{j=1}^J(1-\tau L_j^n)\right].
\tag{8}
\]
此 mean恰为τ|E|，weak量常数有偿；没有 weak反例。

对每固定n先取J→∞，\(\tau\to\tfrac12e^{-n/2}>0,\ p\to0\)，于是
\[
 \mathbb E^\#(1/N)=s/\mathbb EN=2s\longrightarrow e^{-n/2}.
\tag{9}
\]
所以若指望这个共享U coupling在所有网格上有 \(C_n=\sqrt n\,n^{o(1)}\) 的Palm倒数下界，原centered hard模型已排除它。

更直接地，对固定n、任意有限K=K_n≥1，
\[
 \Pr^\#(N\le K)=2\mathbb E[N\mathbf1_{1\le N\le K}]
 \le 2(1-p)\tau+2pK\,\Pr(S\le K).
\tag{10}
\]
因为 \(\Pr(S\le K)\le1\)，(10)已足以给上界 \(2\tau+2pK\)；每固定n,K都有p→0，无需再调用尾工具。因此
\[
 \limsup_{J\to\infty}\Pr^\#(N\le K_n)\le e^{-n/2}.
\tag{11}
\]
对任意固定δ<1，先选固定n令 \(e^{-n/2}<1-\delta\)，再令J大，合同
\(\Pr^\#(N\le K_n)\ge1-\delta\)失败；这里K_n可依n但不能依网格J。

有限J的源law中H稀有，但其S均值随J增长，导致 accepted-occurrence Palm看到非忽略的大count。真正的 receiver Lebesgue winner不交性没有防止同一射线多次跨不同annuli命中。

还有精确的尾不稳定：每固定n，N依分布趋Bern\((e^{-n/2}/2)\)，而所有J的mean都是1/2，极限law的mean仅 \(e^{-n/2}/2\)。因此这个count族不 uniformly integrable；不能用 marginal/finite网逼近就推该coupling的count尾一致性。

## 3. 真正 L1 输入，显式共同coupling误差

现在取 \(\mu_\varepsilon=(2\varepsilon)^{-n}\mathbf1_{[-\varepsilon,\varepsilon]^n}dx\)，W=1，同一尺度和τ，重算全部 \(F_j^\varepsilon,M^\varepsilon,A_j^\varepsilon,g_j^\varepsilon\)。没有冻结atomic winner当L1winner。来源写 \(y=\varepsilon V\)、V∈[−1,1]^n。sharedU及coins沿用，以便比较完整N law。

若
\[
 \frac1{2q}+2\varepsilon<\|U\|_\infty<
 \frac12-2\varepsilon,
\tag{12}
\]
则任意y∈来源支撑、每个j≥1，有
\[
 \frac{L_{j-1}}2+\varepsilon<\|y+L_jU\|_\infty
 <\frac{L_j}2-\varepsilon.
\tag{13}
\]
因此小于j的所有cube与来源不交，j的cube完整捕获来源，大于j的响应≤L_k^{-n}<L_j^{-n}，真实uniquewinner为j且 \(g_j^\varepsilon(Y_j)=\tau L_j^n\)。j0也在Q1内部完整捕获，真实winner0。

若 \(\|U\|_\infty<1/(2q)-2\varepsilon\)，则j≥1的输出在其前一尺度完整捕获区，前一响应较大，故不能以j为winner，\(g_j^\varepsilon(Y_j)=0\)；j0仍准确取τ。与(7)完全相同。

设 \(4\varepsilon\le1-1/q\)。两个good区的补集仅包括 norm的内临界层与outer面层，其概率至多
\[
 [(1/q+4\varepsilon)^n-(1/q-4\varepsilon)^n]
       +[1-(1-4\varepsilon)^n]\le12n\varepsilon.
\tag{14}
\]
同coins保证good区全部接受向量一致，故
\[
 \operatorname{TV}(\operatorname{Law}N^\varepsilon,
                         \operatorname{Law}N^0)\le12n\varepsilon,
 \qquad |\mathbb EN^\varepsilon-\tfrac12|\le12n(J+1)\varepsilon.
\tag{15}
\]
取预固定 \(\varepsilon_J=1/(10^6nJ^2)\)，(14)资格对J≥2成立，且mean/任意[0,J+1]count statistic误差≤18/(10^6J)。N包含j=0..J共J+1项，不能漏掉A0。由 (15) 与阈值实际winner的Fubini身份，整个L1输入的 \(\tau|E^\varepsilon|=\mathbb EN^\varepsilon\to1/2\)。Palm倒数 numerator是stop概率，误差≤12nε；smallcount numerator范围≤J+1，误差≤12n(J+1)ε。因此(9)/(11)对这份真实L1序列同样成立。没有atomic-only漏洞，也不需要未经认证的winner浮点稳定性。

## 4. fullfuture与actual资格的范围

模型用了真实中心cube、同一完整来源、finite全域winner、严格阈与正确maxweight；共享U几何真实，A0未删除。它没有原完整allsoft future帽、FIRST、σ>1/n、CP/GP、两来源/LCA、seed/history等资格。C=0亦被原early资格排除。因此它只否定欲对纯hard/任意来源先证明、再免费移入actual的这个**coupling-specific强尾合同**；不能据此断言actual余项的相同count合同也失败。

有价值的较窄非循环合同必须使用这些遗漏的actual结构，或改变coupling而不要求sharedU的Palm小count。可直接改为加权rareactivation费用
\[
 \mathbb E[N\mathbf1_{N>K_n}]
 \le B_n+\delta\,\mathbb EN,\qquad B_n=O(\sqrt n\,n^{o(1)}),
\tag{16}
\]
结合 \(\mathbb E[N\mathbf1_{N\le K_n}]\le K_n\Pr(N>0)\le K_n\)即给mean预算。它允许(7)中罕见大count贡献由source预算支付，不把它错当必须小Palm概率。**目前未证(16)**；对任意模型取B=EN仍会循环，实际证明需独立column/几何charge。本文不把它登记为新付款，仅明确旧强小count合同何处需要减弱。

## 5. 新三轮守卫登记及终态范围

先注册n2/4/8、J8/16/32，全部Fraction计算(5)–(8)的完整联合law、Palm、threshold/winner margin和(14)/(15)的L1误差界。守卫不抽receiver网格、不模拟actual history，不重跑任何旧count fixture。全J渐近失败由(9)–(15)解析给出，有限三轮只检查新表示/系数，不能由样点替代渐近证明。

第一版执行后独审发现bounded-statistic误差系数漏A0：N最大为J+1，须把12nJε改为12n(J+1)ε。首版script/registration/results按同前缀v1分别保存，未覆盖；首版执行script SHA为14633da9fda815cef640e0ad5c0c1d02f3a57be3f5a7dc98f99f63dbabecae91。修订注册后只重跑本轮新guard；atomic count law与全部Palm值未变化，更正的是L1误差范围。

修订版 **300/300 Fraction PASS**（三轮57/89/153，加注册1），约0.51秒，无live handle。原hard原子模型的count law完整精确；L1迁移由uniform margin与bad-U有理体积上界认证，不把它们写成通用toy或实际history样本。

|n|J|checks|mean|Palm倒数显示值（仅显示）|
|---:|---:|---:|---:|---:|
|2|8|57|1/2|0.558038|
|4|16|89|1/2|0.352549|
|8|32|153|1/2|0.248382|

大Fraction以无符号整数bytes的numerator/denominator SHA、bit长度及显示近似保存，全部判断使用原Fraction。当前 [script](radial_shared_uniform_budget_20261007_exact_guard.py) SHA256为4bcc7683aeb8656bd7a601b339674645fb7fca24d627afd57753bbb0dbf0e2fe；[registration](radial_shared_uniform_budget_20261007_registration.json)、[results](radial_shared_uniform_budget_20261007_results.json)保留每轮完整law记录。无随机种子、无actual FIRST/fullfuture样本、无新增主账费用。
