# 固定真实构造类中的迭代保护、任意精度恢复与近零帽边界

2026-11-05。固定维数。本文只研究明确的斜格孔与指定见证规则；不把类内固定点、质量根或恢复程序解释为真实 \(C_n\) 的极值唯一性，也不假设参数依赖解析。

新附件已全文读取；本地只读版本快照见 provenance.json：main 分支，提交 205596954c5a3265532b8dfbda0354e20561dd68，读取时工作区干净。论文源和 PDF 的实际 SHA 同记录。以下是新的独立结构命题，仓库及其他代理源未改。

## 1. 类、真实算子和非极值范围

固定整数 \(n\ge2\)，令
\[
 u=t^{1/n},\quad a=(2u)^n,\quad g=(1+2u)^{-1},\quad
 e=\kappa u^3,\quad T=(1-e)I+eJ,
\]
\[
 w=ug+\Lambda u^6,\quad \rho=1-2w,\quad C_w=[-w,w]^n,
 \quad M=a\rho^n,\quad d=\det T.
 \tag{1}
\]
\(\kappa\) 在固定正紧区间内，\(\Lambda\) 在固定有界区间内。所有常数可以依这些区间和 n；先固定它们再取小 u。当前 \(M\) 是**设计质量**，规定长见证为
\[
 s=w-\|x\|_\infty,\qquad r(x)=\rho+s.
 \tag{2}
\]
实际切出的总孔质量可以小于 M；若如此，下文会显式记为 \(M'\)，并继续使用设计容量的支配，不静默重算根。

对任意 Lebesgue 可测 \(V\subset C_w\)，定义
\[
 F_w(s)=\frac a{(2w)^{n-2}}\big[(\rho+s)^n-\rho^n\big],\quad
 f=F_w'(0),\quad \delta=e/f,
\]
\[
 R_V(x)=\frac1{(2w)^{n-2}}\sum_{k\in\mathbb Z^n\setminus\{0\}}
       |(Tk+(C_w\setminus V))\cap Q(x,\rho+s)|,
\]
\[
 \mathcal B(V)=\{s=0\}\cup\{x\in C_w:R_V(x)\ge F_w(s)\}.
 \tag{3}
\]
它取真正外孔交量，不是形式核。每个输入集合的变化均不改变 \(T,w,\rho,M,F_w\) 或原规则半径。

几何清单如下。由
\(r+w=1-\|x\|_\infty\) 和
\(\|Tk\|_\infty\ge[1-(n-1)e]\|k\|_\infty\)，模至少二的指标不能进入。
对 \(k\in\{-1,0,1\}^n\)，若 \(\sum k_i\ne0\)，某坐标的晶格中心模至少为一，至多边界接触。只有平衡非零指标进入，支撑 j 为偶数且 \(j\ge2\)，\(Tk=(1-e)k\)。支撑交长不超过 e，其余坐标全宽进入。
又 \(\rho>2w\)，故本胞全部实际孔总质量进入原长见证。

所有 \(\mathcal B(V)\) 除外边界外在
\[
 s\le h_H:=F_w^{-1}(3^ne^2)\le C_ne\delta=O_n(u^4)
 \tag{4}
\]
内。以最近面符号和全局最浅轴 q 分图，写 \(d_i=s+y_i,\ y_q=0\)。
Jacobian 为一，精确纤维顶端 \(H(y)=w-\max_i y_i\)。
活动支撑前缀宽为 \((e-y_i)_+\)，与 s 无关，所以每个坏层为
\[
 0\le s\le h_V(y):=\min\{H(y),F_w^{-1}(R_V(y))\}.
 \tag{5}
\]
这里保留最浅轴不在某个目标支撑内的图、真实顶端截断及零帽。有限图分界和平局均为零测；直接用 (3) 的闭条件取 canonical 代表。帽交量对前缀端点连续，故输出为闭集，且其边界零测。

## 2. 完备空间中的唯一固定保护层

令 \(\mathfrak X_w\) 为 \(C_w\) 中全部可测集合按零测对称差取商，距离
\[
 \mathrm d(V,W)=|V\triangle W|.
 \]
这是完备空间：指标函数的 \(L^1\) Cauchy 列有 \(L^1\) 极限，取几乎处处收敛子列后该极限仍只取 0、1。

**定理 1（固定类的收缩与夹逼）。** 存在 \(u_0>0,C_0,C_1<\infty\)，只依上述固定参数，使 \(0<u<u_0\) 时
\[
 \mathrm d(\mathcal B(V),\mathcal B(W))
 \le q_u\,\mathrm d(V,W),\qquad q_u=C_0u<1/2,
 \quad |\mathcal B(\varnothing)|\le C_1au^5.
 \tag{6}
\]
因此有唯一模零固定类 \(U_*\)。取 canonical 代表后，
\(\mathcal B(U_*)=U_*\) 点态成立。
定义 \(U_0=\mathcal B(\varnothing),U_{j+1}=\mathcal B(U_j)\)，则
\[
 U_1\subset U_3\subset\cdots\subset U_*
       \subset\cdots\subset U_4\subset U_2\subset U_0,
 \tag{7}
\]
\[
 |U_j\triangle U_*|\le q_u^{j+1}|U_0|,\qquad
 |U_j\triangle U_{j+1}|\le q_u^{j+1}|U_0|.
 \tag{8}
\]
每个偶层为可给出严格见证的超解，\(\mathcal B(U_{2m})\subset U_{2m}\)。有限计算还给可观测体积区间
\[
 |U_{2m+1}|\le |U_*|\le|U_{2m}|,
 \quad |U_{2m}|-|U_{2m+1}|\le q_u^{2m+1}|U_0|.
 \tag{9}
\]

证明。在 (5) 中使用 \(F_w'\ge f\) 和真实截断的 1-Lipschitz 性。
打开 Tonelli 固定一个目标物理点；支撑大小 j 的目标，其源底面测度在最浅轴属于支撑时至多 \(C_ne^{j-1}w^{n-j}\)，否则至多 \(C_ne^jw^{n-j-1}\)。
除以 \(f(2w)^{n-2}\)，分别为
\[
 C_n\delta(e/w)^{j-2},\qquad C_n\delta(e/w)^{j-1}.
\]
有限平衡指标求和，结合 \(f\asymp u^2,\delta\asymp u,e/w=O(u^2)\)，得到 (6) 的小常数。完整名义帽同样积分给 \(|U_0|\le C_1au^5\)。
输入只通过帽体积出现，所以算子在模零空间良定义。

Banach 收缩定理给唯一固定类。若取该类任一代表 E，再置 \(U_*=\mathcal B(E)\)，则 \(U_*=E\) 几乎处处；对所有帽零测差都不改变交量，故
\(\mathcal B(U_*)=\mathcal B(E)=U_*\) 点态成立。
此 canonical 代表闭且边界零测。任何其他点态固定代表也因帽交量相同而与之逐点相同。

反单调性给有限偶／奇夹逼。特别 \(U_*\subset U_0\)，所以从 \(\varnothing\) 到固定层迭代 \(j+1\) 次给 (8) 第一式；相邻差同理由 \(|\mathcal B(\varnothing)\triangle\varnothing|=|U_0|\) 取得。点态 (7) 可由闭图高度的极限或帽交量连续性核合；高度函数的 uniform 差至多
\[
 C_n\,|V\triangle W|/[f(2w)^{n-2}],
\]
故固定 u 时奇偶极限的图一致。(9) 随即成立。证毕。

这个唯一性只针对固定 (1)–(3)。它不说所有 supersolution 相同，不说固定点包含于每个 supersolution，也不说任意真实近极值输入具有唯一结构。

## 3. 强制保留微见证 host 的版本

设 \(H\subset K\subset C_w\) 可测。H 是给定的微孔，K 包含它们所需的全部短见证；不预设微孔的帽内分布。若微孔本身已对缺集 H 有原阈值的严格短见证且见证均在 K，则以后只需保证宏孔避开 K。

在 \(\mathfrak X_K=\{V:K\subset V\subset C_w\}\) 上定义
\[
 \mathcal T_{K,H}(V)=K\cup\mathcal B(V\setminus H).
 \tag{10}
\]
它将 K 强制保留，不能仅强制 H 后迭代，否则宏孔可能重新侵入微见证 host。
映射仍反单调，并有相同收缩常数：
\[
 |\mathcal T_{K,H}(V)\triangle\mathcal T_{K,H}(W)|
 \le q_u|V\triangle W|.
 \tag{11}
\]

**定理 2（任意分布的固定 host 保护）。** (10) 有唯一模零固定保护层 \(V_{K,H}\)，并可取
\[
 V_{K,H}=K\cup\mathcal B(V_{K,H}\setminus H)
 \tag{12}
\]
的点态 canonical 代表。
编号必须从 bottom K 开始：令
\[
 W_0=\mathcal T_{K,H}(K),\qquad
 W_{j+1}=\mathcal T_{K,H}(W_j).
 \tag{13}
\]
则 \(W_0\) 是第一个 upper 层，\(W_1\) 是 lower 层；
\[
 W_1\subset W_3\subset V_{K,H}\subset W_2\subset W_0,
 \quad |\ W_j\triangle V_{K,H}\ |\le q_u^{j+1}|U_0|.
 \tag{14}
\]
所有偶层含 K 且是实际宏／微 supersolution。

证明。(11) 来自并固定 K、去固定 H 均不增加对称差，再用 (6)。
\(W_0\setminus K\subset U_0\)，所以从最小元 K 开始的误差初值仍是 \(|U_0|\)，得到 (14)。点态代表与定理 1 同理。
对偶层 W，\(\mathcal B(W\setminus H)\subset W\)。
外国实际总孔 \((C_w\setminus W)\cup H\) 恰为 \(C_w\setminus(W\setminus H)\)；
因此宏孔用同一个原长见证安全，微孔的 host 则由 \(K\subset W\) 保留。证毕。

若 K 本来已是 \(\mathcal B(K)\subset K\) 的保护层，则
\[
 |W_0\setminus K|\le q_u|H|,\qquad
 |V_{K,H}\setminus K|\le q_u|H|.
 \tag{15}
\]
这里 \(W_0=K\cup\mathcal B(K\setminus H)\) 正是一次显式保护；后续偶／奇迭代在不释放 host 的前提下给任意更高精度。

## 4. 固定设计下的有限积分恢复与真实质量

本节额外取 K、H 为有限个闭轴平行盒之并（允许任意位置及不同尺度），设 \(H\subset K\)、\(|H|<M\)，而且微孔短见证已经在 K 中验证。记
\[
 A_*=(2w)^n-|V_{K,H}|,\qquad
 G=(2w)^n+|H|-|V_{K,H}|-M.
 \tag{16}
\]
G 是这个 canonical 构造的预算，不是所有真实输入上的极值。

**定理 3（有限安全恢复）。** 给任意 \(\varepsilon>0\)，通过有限次数的 (13) 和有限精度积分，可构造宏孔 \(A_\varepsilon\subset C_w\setminus W_{2m}\) 及总实际质量
\(M'=|A_\varepsilon|+|H|\le M\)，满足
\[
 0\le M-M'\le(-G)_++\varepsilon,
 \tag{17}
\]
只要目标宏孔质量非负；若 G≥0 且 \(\varepsilon<M-|H|\)，则 \(M'\ge M-\varepsilon\)。
宏／微全部孔点都有原阈值的严格见证。

证明。选 m 使 \(q_u^{2m+1}|U_0|\le\varepsilon/8\)。
有限偶层 W 已知，(14) 给
\[
 0\le |W|-|V_{K,H}|\le\varepsilon/8.
\]
W 是闭集，所以可用区开。对它作内部有限紧盒近似，取得不交紧盒总体积 b，使
\[
 (2w)^n-|W|-\varepsilon/8\le b\le(2w)^n-|W|.
\]
切取宏孔质量
\[
 |A_\varepsilon|=\min\{M-|H|,\ (b-\varepsilon/4)_+\}.
 \tag{18}
\]
有限盒最后一个沿一轴切开即可取得此精确质量。对任意小实参数这首先是测度论构造；下节只在可计算参数与有限盒数据上证明其可认证算法，不对不可计算实数声称有读取算法。
当 G≥0 且精度充分小时非零分支有效，(17) 由上述三个误差直接得到；一般情形的多切量至多 \((-G)_++\varepsilon\)。

关键是**不重定义设计根**：长见证仍用 (2)，容量仍按设计 \(M=a\rho^n\)。
实际自身总孔质量 \(M'\le M\)，外国孔由
\((C_w\setminus W)\cup H\) 上包。
宏孔点 \(x\notin W\) 因 W 是 supersolution，有
\[
 M'+\text{外国实际孔质量}
 \le M'+(2w)^{n-2}R_{W\setminus H}(x)
 < M+ (2w)^{n-2}F_w(s)=ar(x)^n.
 \tag{19}
\]
微见证在 K 内，宏孔避开 K，所以其原严格性不变。
小量实际质量切损不能反过来塞回算子重算几何。证毕。

严格的 (19) 在 \(M'=M\) 时也成立；它来自 \(x\notin W\) 的坏集补集，而非只靠少切质量制造一个平均裕量。

有限化取 \(\Omega_N=T([-N,N]^n)\)，保留全部边界切孔。
输入质量精确为 \((2N)^n(d-M')\)，内层孔全部见证保留，边界比值误差为 \(O_n(M'/N)\)。
凸外盒之外任意中心任意半径至多半体积在分离半空间内，原阈值超过 \(1/2\) 后不进入；输入密度点进入 halo。因此得到真实全部半径的比值
\[
 \frac{|\mathcal H_n(1-t;E)|}{|E|}
 =\frac d{d-M'}+O_n(M'/N).
 \tag{20}
\]
只要 G≥0，选择 \(\varepsilon=au^{r+1}\)、\(N=\lceil u^{-(r+3)}\rceil\)，固定阶 r 的反馈深度可取任一满足 \(2m+6>r+1\) 的整数。于是
\[
 \frac{|\mathcal H_n(1-t;E)|}{|E|}
 \ge\frac d{d-M}-o_{n,r}(au^r).
 \tag{21}
\]
每个固定小实 u 都有此真实有限输入；这里没有宣称右边函数存在纯幂展开。

## 5. 不使用体积 oracle 的可认证有限过程

本节的算法陈述限定参数为有理数，或具有有效误差界的可计算实数；K、H 用有限闭轴盒给出。存在性定理仍适用于每个足够小实 u，但不声称能读取任意不可计算实数。

算法不必知道 (6) 中最优的常数：在 \(e/w\le1/2\) 时，对有限图／目标清单直接过计数可取
\(q_{\rm cert}=2^{3n}n3^n e/f\) 作为收缩上界，先认证 \(q_{\rm cert}<1/2\)。
初值可粗取 \((2w)^n\) 而非未知固定层体积。选 m 使
\(q_{\rm cert}^{2m+1}(2w)^n\le\varepsilon/8\)，是显式有限整数搜索。

所有积分都可还原为有限 chart 上连续 Lipschitz 函数的积分。对任意测试缺集，改变帽的某个前缀端点，交量之差不超过移动厚度乘完整其他横截面。支撑 j 的前缀端点 Lipschitz 常数至多
\[
 C_ne^{j-1}(2w)^{n-j}.
\]
除以 \(f(2w)^{n-2}\) 后可显式计算；min、正部和顶端截断不增加最大 Lipschitz 常数。因此每个输出图 \(h_V\) 都有一个只依当前可计算几何数据的统一 Lipschitz 上界，粗取 \(1+C_n\delta\) 即可，不依测试集边界复杂性或迭代深度。

若输入本身为已构造的有限图层，帽与该层的交量在目标 chart 中是
\[
 \int_{\Omega}
 \left[\min\{h(z),\,\min_{i\in\mathrm{support}}(P_i-z_i)\}\right]_+\,dz,
 \tag{22}
\]
再对有限符号、最浅轴求和。原帽只取临界面对端的前缀；非临界坐标全宽，所以没有其他移动的下端。
K/H 为有限物理轴盒时，这里始终使用上述 global-min 纤维
\(x_h=\sigma_h(w-s-y_h)\)，不是只沿一个物理坐标移动的纤维。
对盒 \(\prod_h[a_h,b_h]\)，置
\(\ell_h^\sigma=\min(\sigma_ha_h,\sigma_hb_h)\)、
\(u_h^\sigma=\max(\sigma_ha_h,\sigma_hb_h)\)，则其在该纤维中的区间为
\[
 \left[\max\{0,\max_h(w-y_h-u_h^\sigma)\},\
       \min\{H(y),\min_h(w-y_h-\ell_h^\sigma)\}\right].
\]
上下端均为有限仿射函数的 max/min；再与图 \(s\le h(y)\) 和帽区间相交，区间长度取正部。有限个这种区间的并／差长度对全部端点连续 Lipschitz，所以由盒清单可显式计算积分函数的 Lipschitz 常数。若另用普通轴向纤维会有投影边界跳跃，须另作分片；本算法不使用该坐标。输入图高度的数值误差经这些区间操作至多线性放大；各深度为有限值，所以可以预先把内部精度分配到所需总误差。

对一个 Lipschitz 常数 L 的函数，在边长网格 \(\eta\) 上取有理上下包，积分误差至多 \(L\operatorname{diam}(\text{网格胞})|\Omega|\)，另加点值计算误差。取足够细的确定性网格，有限终止，得到每个帽质量和每个层体积的认证上下界。
有理参数下 F 是有理系数多项式，严格递增且 \(F'\ge f_{\rm lower}>0\)；对唯一正根作有理二分隔离，达到任意精度也有限终止。可计算实参数用其有效区间包替代，正的宽度、f 和几何裕量先取有理下界。

这些是有限积分与有理根隔离的数学算法，不要求任何现成程序已经运行区间认证。本轮成功软件只核后文有限标量公式，未把浮点误差估计当作认证。

对于 (18) 的内部盒，不用未知可用体积作符号判断。枚举可用开集内闭包仍在其中的 dyadic 盒，利用上述图值／Lipschitz 区间包认证严格不相交。正安全偏移 \(\varepsilon/4\) 使所需有限盒总体积严格低于实际可用体积；内正则逼近保证枚举最终找到足量盒。最后一个盒沿一轴截取精确剩余质量。有理参数和有理目标质量时该截面坐标也有理。
具体可先得可用体积上包 \(A^+\)，其与实际可用体积之差不超过 \(\varepsilon/16\)，再累计不交盒直到体积至少 \((A^+-\varepsilon/8)_+\)；目标严格低于实际可用体积，所以有限终止。可计算目标质量即使恰等于某个前缀盒体积，也可用逐盒的 min/max 截取量直接定义切孔，避免等待等号判定。

因此“有限次数”有实质内容：

* 反馈深度由 (6)、所需体积精度和已知初始体积上界决定，不依未知 fixedpoint gap；
* 每次积分用已知 Lipschitz 模量加细网格终止；
* 若某个符号比较恰为零，不让程序无限等候其符号，而改用体积上下包与预定安全偏移；
* 微孔相对 host 的严格见证是输入证书，保护算法不读取未知最优分布。

代价可以随固定 n、阶数、当前参数及有限盒数量增长；本稿不声称高效复杂度或维数／阶数一致性。

## 6. 一个可有限夹逼的 canonical 质量根

这一节只取纯宏观 \(K=H=\varnothing\)，并允许有界 \(\Lambda\) 变化。定义
\[
 D_u(\Lambda)=(2w)^n-a\rho^n,\qquad
 G_u(\Lambda)=D_u(\Lambda)-|U_*(\Lambda)|.
 \tag{23}
\]
这是 canonical 固定保护层的预算函数，不是所有可行孔形上的上确界。

内禀最近面坐标 \(x_i=\sigma_i(w-d_i)\) 把不同 w 嵌入同一有限符号直和的 \(\mathbb R_+^n\)。所有支撑前缀仍为 \([0,(e-y_i)_+]\)，不随 w 变；非支撑积分限、归一化、F 和中分面顶端才变。固定高度窗中
\[
 |\partial_wF_w(s)|\le C_nF_w(s)/w.
\]
逐名义帽的 Fubini 边带估计给
\[
 \|\mathcal I_w\mathcal B_w(V)-\mathcal I_{w'}\mathcal B_{w'}(V')\|_1
 \le q_u\|\mathcal I_wV-\mathcal I_{w'}V'\|_1+
 C_nau^5|w-w'|/w .
\]
临界变量远离中分面；源顶端变化只在自由变量端产生同一比例边带。此估计保留全部实际多面图，也不需要输入边界规则。对固定点移项，因 \(w-w'=u^6(\Lambda-\Lambda')\)，得到
\[
 \big||U_*(\Lambda')|-|U_*(\Lambda)|\big|
 \le C_nau^{10}|\Lambda'-\Lambda|/(1-q_u).
 \tag{24}
\]
与此同时
\[
 \partial_\Lambda D_u(\Lambda)
 =na u^5g^{n-2}+O_n(au^{10}).
 \tag{25}
\]
故有显式正数 \(c_-,c_+\)，使固定有界区间中 \(\Lambda'>\Lambda\) 时
\[
 c_-au^5(\Lambda'-\Lambda)\le
 G_u(\Lambda')-G_u(\Lambda)
 \le c_+au^5(\Lambda'-\Lambda).
 \tag{26}
\]

令 \(\Lambda_0=B_n\kappa^3/(8n^2g)\)。名义二支撑积分和一次收缩给
\[
 |U_*|/a=n\Lambda_0u^5g^{n-2}+O_n(u^6),\qquad
 D_u(\Lambda)/a=n\Lambda u^5g^{n-2}+O_n(u^{10}).
\]
所以 \([\Lambda_0-1,\Lambda_0+1]\) 中有唯一根 \(\Lambda_*(u)\)。只用连续性、正斜率和夹逼，不假设其关于 u 可微或解析。定义其 canonical 比值
\[
 \mathcal Q_*(u)=\frac d{d-a[1-2w(\Lambda_*)]^n}.
 \tag{27}
\]

可通过真正有限的体积区间二分来近似此根。当前根区间 \([L,R]\) 长 \(\ell\)，在中点 z 选有限 m 并认证体积，使
\[
 G_-(z)=D_u(z)-|U_{2m}(z)|\le G_u(z)\le
 G_+(z)=D_u(z)-|U_{2m+1}(z)|,
\]
总认证区间宽度不超过 \(c_-au^5\ell/16\)。
若整个区间为正，取左半；若整个区间为负，取右半。
若区间含零，(26) 直接把根限制在
\[
 [z-\ell/16,z+\ell/16]\cap[L,R].
 \tag{28}
\]
即便中点恰等于根，此步也缩短区间，不等待其精确符号。
每次至少按固定比例缩短，达到预定宽度 \(\eta_\Lambda>0\) 有显式有限步数。

取最终上端再加 \(\eta_\Lambda\)，则 canonical 预算至少为 \(c_-au^5\eta_\Lambda\)。选偶层误差和认证积分误差小于其四分之一，即可按定理 3 切取**完整设计质量** M。也可不偏移根，而采用定理 3 少切质量；两种恢复都明确区分设计 M 与实际 M′。
由于 \(|\partial_\Lambda M|\le C_nau^6\)，取
\(\eta_\Lambda=u^{r+1}\) 并取足够大的固定反馈深度 m（例如 \(m=r+2\)），然后选 \(N=\lceil u^{-(r+8)}\rceil\)，得到实际有限输入满足
\[
 0\le\mathcal Q_*(u)-
   \frac{|\mathcal H_n(1-t;E_u)|}{|E_u|}
 \le C_{n,r}au^{r+7}+C_na/N=o_{n,r}(au^r).
 \tag{29}
\]
取上端加偏移使质量不超过根处质量，所以左侧非负；凸盒外的全半径排除给相同比值上包。
反馈次数对固定 r 可先固定，标量二分查询数至多 \(O_{n,r}(\log(1/u))\)，每个查询又由有限认证积分完成。

因此类内确有任意指定阶精度的有限实际恢复递归；它不只是把未知系数写成一个未知极限或 supremum。另一方面，(29) 只逼近明确函数 \(\mathcal Q_*\)，未证明此函数有全阶纯幂展开，更未证明它等于真实 \(C_n\)。

## 7. 二维 min 帽、矩捷径与两个完整积分

本节先取常斜率容量的精确平面算子。它保留两个排序、真实 min 帽和小帽裁剪，但不是把完整非线性 \(F_w\) 换掉后还称为原高维算子。令 \(0<\delta<1/10\)，
\[
 H_0(d)=\delta(1-d),\quad
 \mathcal M(H)(q)=\int_0^q\min\{H(d),q-d\}\,dd+
                    \int_0^1\min\{H(d),q\}\,dd,
\]
\[
 H_{j+1}(d)=\delta[(1-d)-\mathcal M(H_j)(1-d)].
 \tag{30}
\]
对应平面层在 \(X\le Y\) 中为 \(0\le X\le H(Y-X)\)，另一个排序交换两轴。
\(\mathcal M(H)(q)\) 是该层与 \([0,1]\times[0,q]\) 的真实面积；
所以 \(0\le\mathcal M(H)(q)\le q\)，递归保持
\[
 0\le H_j(d)\le\delta(1-d)\le1-d.
 \tag{31}
\]

交换 q、d 的积分次序，对任意 \(0\le H\le1-d\)，有
\[
 \int_0^1\mathcal M(H)(q)\,dq
 =\int_0^1[(2-d)H(d)-H(d)^2]\,dd.
 \tag{32}
\]
第一项内积分为 \(H(1-d)-H^2/2\)，第二项为 \(H-H^2/2\)，因此该恒等式不依小参数近似。令
\[
 I_j=2\int_0^1\mathcal M(H_j)(q)\,dq,\qquad
 |\mathcal U_j|=2\int_0^1H_j(d)\,dd.
\]
则
\[
 |\mathcal U_{j+1}|=\delta-\delta I_j.
 \tag{33}
\]
这里 \(I_2\) 不是 \(|\mathcal U_2|\)，也不能只由一个 \(I_1\) 数值还原。

可完整算出第一帽：
\[
 \mathcal M(H_0)(q)=
 \begin{cases}
 q-(1-\delta)q^2/(2\delta),&0\le q\le\delta,\\
 \delta(1-2\delta+2q-q^2)/[2(1-\delta)],&\delta\le q\le1.
 \end{cases}
 \tag{34}
\]
分界来自两条直线高度与 cap 相等，故包含 \(q\ll\delta\) 的全域。
写 \(h(r)=H_1(1-r)=\delta[r-\mathcal M(H_0)(r)]\)。
由 (32) 将两段多项式直接积分，得到精确式
\[
 I_1=\frac{\delta(\delta-1)
   (-100+95\delta-60\delta^2+26\delta^3+4\delta^4)}{60}.
 \tag{35}
\]
特别它在零点解析，不含 \(\delta\)-对数或非整数幂。

下一积分也无需嵌套根号数值积分。置
\[
 P(t)=\int_0^t h(r)\,dr,\qquad A=P(1),\qquad
 \Psi(q,m)=(1+q)\delta(q-m)-\delta^2(q-m)^2.
\]
\(h\) 严格递增且 \(0\le h'\le\delta<1\)，因此两种反变量分离：
当 \(q<h(1)\)，取唯一 \(q=h(t)\)，
\[
 \mathcal M(H_1)(q)=q^2/2+q(1-t)+P(t).
\]
当 \(q>h(1)\)，取唯一 \(q=1-t+h(t)\)，
\[
 \mathcal M(H_1)(q)=2A-P(t)+h(t)^2/2.
\]
第一分支中帽的第一排序全取三角区，第二分支中第二排序已取全层；
因而没有同时使用两个独立反根。于是精确得到
\[
 \begin{aligned}
 I_2={}&2\int_0^1
 \Psi\bigl(h(t),h(t)^2/2+h(t)(1-t)+P(t)\bigr)h'(t)\,dt\\
 &+2\int_0^1
 \Psi\bigl(1-t+h(t),2A-P(t)+h(t)^2/2\bigr)(1-h'(t))\,dt .
 \end{aligned}
 \tag{36}
\]
将 t 只在 \(\delta\) 处分段，h、P 都是 t 的多项式且系数为 \(\delta\) 的有理函数。因此 (36) 属于 \(\mathbb Q(\delta)\)；又 \(I_2=O(\delta)\)，零点可能的极点可去。所以这一具体完整积分的零点展开只能含整数幂，没有 log 或分数幂。

实际本机精确积分进一步化简为分母 30240 的十四次多项式，首项为
\[
 I_2=\frac53\delta-\frac{13}{4}\delta^2+
                \frac{33}{5}\delta^3-\frac{227}{30}\delta^4+O(\delta^5).
 \tag{37}
\]
\([\,\delta^3\,](I_2-I_1)=241/60\)，与实际第三平面层差的独立 \(L^1\) 推导回归一致。此处精确代数验证不承担任何高维几何或真实 \(C_n\) 的证明。

## 8. 近零帽的 Fibonacci 空间指数

全局参数的整数幂不意味着内部轮廓没有分数幂。以下给一个真实边界现象，先固定 \(\delta>0\)，再令空间帽宽 \(r\downarrow0\)。
写 \(h_j(r)=H_j(1-r)\)。每个 \(h_j\) 递增且 \(0\le h'_j\le\delta\)；证明是 \(\mathcal M(H_j)'(r)\) 等于水平截面长度，介于零和一之间。该长度在 \(r>0\) 严格小于一，因为全部层都包含于首层的薄臂包络。故 \(h_j\) 在内部严格递增，可逆，且 \(h_j(1)>0\)。

对足够小 r，第一排序的 cap 恰为 \(r^2/2\)。若
\(x_r=h_j^{-1}(r)\)，第二排序恰为
\[
 r-\int_0^{x_r}[r-h_j(s)]\,ds.
\]
因此
\[
 h_{j+1}(r)=
 \delta\left[\int_0^{x_r}(r-h_j(s))\,ds-r^2/2\right].
 \tag{38}
\]
从 (34) 有
\(h_1(r)=(1-\delta)r^2/2\)（在线性容量模型的 \(0\le r\le\delta\) 内甚至精确）。
若
\[
 h_j(r)\sim A_jr^{\beta_j},\quad \beta_j>1,
\]
则 \(x_r\sim(r/A_j)^{1/\beta_j}\)，实际积分给
\[
 \int_0^{x_r}(r-h_j(s))\,ds
 \sim\frac{\beta_j}{\beta_j+1}A_j^{-1/\beta_j}
                       r^{1+1/\beta_j}.
\]
\(1+1/\beta_j<2\)，所以 (38) 中 \(r^2/2\) 是更高阶。由此
\[
 \beta_1=2,\quad
 \beta_{j+1}=1+\frac1{\beta_j},\qquad
 A_1=(1-\delta)/2,\quad
 A_{j+1}=\delta\frac{\beta_j}{\beta_j+1}A_j^{-1/\beta_j}.
 \tag{39}
\]
若 \(F_1=F_2=1\) 为 Fibonacci 数列，则
\[
 \beta_j=F_{j+2}/F_{j+1}
       =2,\ 3/2,\ 5/3,\ 8/5,\ 13/8,\ldots .
 \tag{40}
\]
这不是对数值曲线拟合的空间指数。

真实二维宏孔的非线性容量也保留这些 leading 空间指数：其高度递归将 (30) 最后一步改为 \(\Phi^{-1}(q-\mathcal M(H_j)(q))\)，其中
\(\Phi(h)=F_w(eh)/e^2,\ (\Phi^{-1})'(0)=\delta\)。
当空间 r→0 时只用这个严格正导数，故 (38) 的右边改为一个在零点与 \(\delta z\) 等价的函数；(39) 的首项递推不变。这里仍只取二维孤立短角，没有把多维联合高度套成一个单角。

若随后再看 \(\delta\downarrow0\)，系数 \(A_j\) 的 leading \(\delta\) 指数满足
\[
 a_1=0,\qquad a_{j+1}=1-a_j/\beta_j,
\]
前六项为 \(0,1,1/3,4/5,1/2,9/13\)。这是空间 tail 系数的参数依赖；其有效空间范围也随 \(\delta\) 缩小，不能将这些数直接当作总面积的 \(\delta\) 分数幂。
正如 (35)–(37) 所示，具体积分可经反变量消去根号而完全没有这类项。后续是否仍能这样闭合，需要另行证明；收缩性本身既不排除 log，也不保证纯幂或解析。

## 9. 本轮新增、计算证据与边界

| 对象 | 本轮结果 | 准确类别 | 尚未供应的内容 |
|---|---|---|---|
| 固定实际 \(\mathcal B_w\) | 唯一模零固定保护层、偶奇误差和严格安全恢复 | 显式真实构造类定理 | 不是 \(C_n\) 极值唯一性 |
| 固定微孔及 host | \(K\cup\mathcal B(V\setminus H)\) 的同常数收缩 | 任意分布的类内保护 | 不优化所有微孔分布 |
| canonical 质量根 | 有限认证积分、无符号停滞二分、任意阶精度恢复 | 非 oracle 递归及实际有限输入 | 未证该根函数全阶展开或等于 \(C_n\) |
| 真实二维空间帽 | Fibonacci 分数空间指数 | 固定参数下的局部几何 | 不推出参数展开非纯幂 |
| 常斜率 min 平面 | \(I_1,I_2\) 完整有理／多项式积分 | 指定算子内精确计算 | 不推广为全层或非线性实际 B 全阶 |

实际本机材料在 r067-computation/。本作者已全文读取 plane_exact.wl/log、tail_exponents.wl/log、i2_common.wl、i2_branches.wl/log、i2_fraction.py/json 及成功核验日志。其完整实际 SHA 清单另存本目录 computation_manifest.json。
Wolfram 15 的 I2 窄核用 15 个实际子核分别计算四分支的十二个多项式矩和三个接点；接点残差为零。独立 Python Fraction 在 \(\delta=1/20,1/40,1/80\) 三点的有理积分与完整多项式完全吻合。它们只核有限积分代数；区间认证恢复算法尚未由本轮程序实际执行，未声称已运行机器区间证明。

本轮不计算或宣布真实 \(C_n\) 的任意阶系数，不假设固定层随 u 解析，不将类内收缩交换为全输入上确界。若要进一步得到 canonical 类的系数递归，必须补充每个有限积分所处的函数类及端点闭合；若要得到 \(C_n\) 的全阶规律，还需独立的全部输入上界与恢复序列匹配。


