# Round 042：固定 histogram selector 的独立数学审计

本文件的数学结论自含。先实际运行 `audit.py`，再整理本证明；`AUDIT.json` 是该运行输出。使用本机 Python 3 的 `fractions.Fraction` 和 NumPy 2.5.4。没有引用未核验文献定理，没有运行新正式搜索。

## 1. 已运行的小例与目标式

K=2,3、n=1,2,3,4、R=3，共80次 Fraction 多项式DP与独立 K^n 标签展开完全相等；其中56次使用任意非负有理权重，24次使用真实周期小盒积分，覆盖跨周期边界及长度大于一个周期。每例检查源质量和matched prior质量均为1、整数查询平移及源与标签格共同平移的不变性。另穷举K=2,3、n=2、R=2的全部72个selector，质量均精确为1。n=8,K=4,R=4的独立浮点全标签展开与DP的绝对差为2.22e−12；这项是数值交叉检查，不能代替区间认证。

另实际运行一个n=2,K=2,R=3、8×8查询格的严格有理下包：半径1/8,1/4,3/8、阈值1/2，共60/64格通过，得到周期超水平集测度下界15/16、weak指标下界15/32。截断到[-10,10)^2后，内部周期格给合法L1 weak比值至少243/640。此小例仅验证接口，不声称改进已有最佳常数。

本轮真正需探索的目标为：固定n、同一个selector s、同一组有限半径T及一个阈值λ，寻找

\[
W_T(f_s;\lambda)=\lambda\,|\{x\in[0,1)^n:\max_{r\in T}A_r f_s(x)>\lambda\}|,
\quad A_rf(x)=(2r)^{-n}\int_{x+[-r,r]^n}f.
\]

这里必须先固定s，再对各查询与半径使用该s。仅在一个n上使W_T大于某个对照，并不等于证明随n的增长阶。

## 2. 定义与精确质量

以I_k=[k/K,(k+1)/K), k=0,...,K−1划分一周期，延拓到实线。每个table b=0,...,R−1对应非负可测周期密度ρ_b，并要求

\[
\int_{I_k}\rho_b(t)\,dt=1/K\qquad\text{对全部 }b,k.
\]

对y∈R^n，用每个坐标的周期格标签定义α(y)=(α_0,...,α_{K−1})，Σα_k=n。selector是固定映射s: {α≥0,Σα=n}→{0,...,R−1}，定义

\[
f_s(y)=\sum_b\prod_{i=1}^n\rho_b(y_i)\,1_{\{s(\alpha(y))=b\}}.
\]

对每个y，指示函数恰选一个table；这里不是同时累加所有table的全部密度。给定一个有α_k次标签k的有序标签串，其在任何table b下的质量都是K^{−n}。该type共有n!/∏α_k!个有序标签串，故在任何table下的type质量相同：

\[
p_\alpha=\frac{n!}{\prod_k\alpha_k!}K^{-n},\qquad
\int_{[0,1)^n}f_s=\sum_\alpha p_\alpha=1.
\]

只需格质量一致；ρ_b在格内可以不同、可以加正背景。给定selector的matched prior是

\[
\pi_b=\sum_{\alpha:s(\alpha)=b}p_\alpha,\qquad
g_s(y)=\sum_b\pi_b\prod_i\rho_b(y_i),\quad \sum_b\pi_b=1.
\]

因此g_s的周期质量也为1。π_b是多项分布下的type质量，不是分给b的histogram数除以histogram总数。此区别在n=8,K=4时尤其重要：165个type的质量并不相同。π_b均可由分母4^8的整数精确计算。

如果selector能够在同一type内部依赖y，质量1结论一般失效。例如每个格内两个table各自占不相交的半格并保留格质量1/K，点wise挑非零table会把总质量变成2。按query选择s_x虽每个f_{s_x}分别质量1，也不能解释成一个共同输入f的最大函数。

## 3. 独立多变量多项式DP

令

\[
w_{bik}(x,r)=\int_{x_i-r}^{x_i+r}\rho_b(t)1_{\{t\bmod1\in I_k\}}dt,
\quad P_b(z)=\prod_{i=1}^n\left(\sum_k w_{bik}z_k\right).
\]

有序标签串展开及Tonelli给

\[
\int_{x+[-r,r]^n}f_s
=\sum_b\sum_{\alpha:s(\alpha)=b}[z^\alpha]P_b(z).
\]

DP从D_0(0)=1开始，依次更新D_i(α+e_k) += D_{i−1}(α)w_{bik}。非负系数无抵消。由归纳法，D_i就是前i个因子的系数；无需独立性近似或正态近似。DP最后选type，不能先对每个坐标独立选table，也不能先对每个query选type-table映射。

prior的积分则为

\[
\int_{x+[-r,r]^n}g_s
=\sum_b\pi_b\prod_i\sum_k w_{bik}.
\]

除以(2r)^n得到平均值；也可先将每个w除以2r再做DP，结果相同。中心立方体的全部坐标必须共用同一r，不能逐坐标分别最大化。有限半径的max是全尺度最大函数的下界。

`audit.py` 的 `histogram_coefficients(weights[n][K])`、`selected_numerator(weights[R][n][K], selector)`、`prior_weights` 是与主程序无依赖的checker。有理输入给精确有理输出。`direct_numerator` 单独枚举K^n个有序标签串，供小n独立核对。selector可以是tuple→int字典或固定callable。

## 4. 周期平移与L1截断转移

f_s的周期质量是1，故它在R^n上的积分无限，**f_s本身不属于L1(R^n)**。不能直接把周期比值代入欧氏弱(1,1)常数。

记有限半径集合T的最大值为R_0，取正整数q≥R_0，整数N>q，令

\[
F_N=f_s1_{[-N,N)^n},\qquad \|F_N\|_1=(2N)^n.
\]

在内部[-N+q,N−q)^n，每个r∈T的查询立方体都在截断盒内（零测度边界无影响），故A_rF_N=A_rf_s。周期超水平集E在内部整数周期格中重复，给

\[
|\{M_cF_N>\lambda\}|\ge (2N-2q)^n |E|,
\quad
\frac{\lambda|\{M_cF_N>\lambda\}|}{\|F_N\|_1}
\ge (1-q/N)^nW_T(f_s;\lambda).
\]

若C_n是欧氏中心立方体的最佳弱(1,1)常数，它对每个F_N都成立；固定n,T后令N→∞，即得严格转移

\[
C_n\ge W_T(f_s;\lambda).
\]

不需要证明全尺度周期max的上界；有限半径已足以转移下界。对一族n，N允许随n选；例如取N远大于nq即可令损失一致接近1，但任何渐近增长结论仍需独立控制W_T随n。

整数平移查询不改变A_rf_s。任意τ的全局源平移必须同时移动格标签：ρ_b^τ(t)=ρ_b(t−τ)、α^τ(y)=α(y−τ)，于是f_s^τ(y)=f_s(y−τ)，A_rf_s^τ(x+τ)=A_rf_s(x)。若只平移小盒却保持原格标签，并让质量跨格，原质量证明不自动成立。主程序的交替相位各自留在原格内，满足本证明假设。

## 5. 有限query-grid怎样给严格下界

有限点集合本身测度为零，所以λ×“格点超阈值比例”只是求积或经验统计，不能直接当作W_T的确定性下界。可以采用下面的正多项式下包。

把[0,1)^n划成互不相交的查询盒Q=∏_i[a_i,c_i]（共享边界忽略），对每个table、坐标、source格定义

\[
\underline w_{bik}(Q,r)=\min_{t\in[a_i,c_i]}w_{bik}(t,r).
\]

用这些下权重跑同一个正系数DP并除以(2r)^n，得到L(Q,r)。系数单调性保证对全部x∈Q有A_rf_s(x)≥L(Q,r)。因此

\[
W_T(f_s;\lambda)\ge
\lambda\sum_{Q:\max_{r\in T}L(Q,r)>\lambda}|Q|.
\]

每个Q可以使用不同的见证半径，因为同一个f_s的max本来逐点选半径。selector仍然必须全部共用。这里使用严格 >λ；L=λ的格不能认证到严格超水平集。

piecewise constant的ρ_b使每个移动窗口积分w(t,r)连续分段线性，折点仅在窗口端点碰到source段端点，即t=e+m±r。取查询区间端点与内部所有这些折点即可精确求min。`coordinate_cell_min` 与 `query_cell_lower` 已按此实现，使用任意长度的周期积分，不依赖主程序只适用当前长度的有限shift快捷式。它可能较保守，但给出的界方向正确。也可由有界ρ给Lipschitz邻域界；本实现直接处理分段线性核。

## 6. 学习selector、holdout及渐近目标的边界

训练可以优化同一个全局selector，也可以先随机初始化；冻结后的任何selector均满足质量与DP公式。随机初始化本身没有产生对数增益的定理。若对每个type独立随机选table，且选择概率是q_b，则E f_s=Σ_bq_bρ_b^{⊗n}；max与严格超水平集均非线性，不能把期望与weak指标无条件互换后宣称改善。即便某个有限维selector超过其matched prior，也未必超过任一固定table，更未必超过已有C_n下界。

独立holdout必须在selector、source参数、半径集合、阈值与检验模型冻结后生成。条件于训练信息，查询iid uniform时，每个严格超阈值指示是固定目标的Bernoulli观测。对固定有限的M个方向分配失败概率并union即可得到同时有效的概率界，不要求同一样本上各模型独立。若在test上又选择阈值或selector，必须另作统一覆盖或再用独立样本。有限样本CI是有置信水平的数值证据，不是无随机失败概率的确定性数值证书。

主目标若是 `C_n ≥ c sqrt(n) (log n)^a`（固定c,a>0、充分大n），需要提供一族合法来源、selector和阈值，严格给出全体充分大n的W_T下界。若目标只是加性 `c sqrt(n)+c' log n`，也必须另给相应一致下界。n=8的实验、随机selector、比matched prior好的一个阈值，均不足以推出任何一种目标。当前审计没有取得这些渐近结论。

## 7. 数学与数值证据的分层

质量、DP、prior、周期截断和query-cell下包的证明是上述自含一般结论。`AUDIT.json`中的Fraction小例是实际执行的精确证据。真实n8实验的浮点核及测试统计需另见 `INDEPENDENT_SELECTOR_CHECK.json` 与 `REVIEW.md`；浮点近似、Monte Carlo覆盖以及向外舍入认证是不同问题。独立重算一致不能把其余浮点运算自动升级成区间认证。

## 8. 实际发现的严格边界问题及修复

独立重算发现，主source每个table都具有1/2周期，而L=3/2的fixed与prior真实平均恒等1。原浮点 `>1` 把舍入误差计为成功，原naive统计失效；本目录 `BOUNDARY_DIAGNOSTIC.json`保留三case完整独立重放及原计数差异。1e−10 guard仅用于诊断，不是认证。

随后主程序冻结了 `interval_recheck.py`：真实有理参数、精确dyadic查询、每个基本运算向外nextafter、正DP及求和逐项向外，并对fixed/prior的L3/2分支使用精确[1,1]恒等式。全部三case、18model的上下成功数相同，与本审计此前独立重放的修正计数一致。`INTERVAL_CHECK.py`另用Fraction检查8个真实查询、全部7尺度的8400个kernel和336个完整gated/prior/fixed平均，全部包含通过。36个CI端点以80位Decimal独立核对，最大差约1e−15；这不认证log函数向外舍入。

修复后相对于matched prior的有限阈值差仍正，最佳fixed table仍更强。`INTERVAL_REVIEW.md`解释了局部算术认证、MC假设及渐近结论的区别。本证明没有因此取得C_n的新增长阶。
