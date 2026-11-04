# R046　一般维数八层共同结算（明确以 R045 联合价格为条件）

2026-11-05。完整条件证明，供独立审查。主仓库和其他代理文件只读。

**当前状态。** R043、R044 已有独立通过及正式归档；本稿实际读取 R045 最终冻结稿，SHA **8f0ad8b01bb443789b33f55188c7a04d000a3ceaee3fd61d0a83eb458c9bf8ea**。本稿不预支它尚未完成的独审。下文证明的是精确的 \((\mathrm H_{45})\Rightarrow\) 八阶余项结论；只有该同层接口独立通过并与本账本合核后，才能去条件化。

## 1. 目标、分支及唯一新增条件

固定 \(n\ge2\)，置
\[
 a=v^n=2^nt,\quad A_v=(1+v)^n,\quad
 B_n=\binom n2,\quad k_n=n^2/27,\quad
 K_n=B_nk_n,\quad \Gamma_n=(10/9)k_n,
\]
\[
 L_n(t)=\frac{A_v}{A_v-a},\qquad d_n=n+26/9,\qquad
 P_7(t)=L_n(t)+K_nav^6-d_nK_nav^7.
 \tag{1}
\]
所求为
\[
 C_n(1-t)\le P_7(t)+O_n(av^8).
 \tag{2}
\]
再与已审 R041 的实际有限族下界组合，可得
\[
 0<c_n av^8\le C_n(1-t)-P_7(t)\le C_n' av^8
 \quad(0<t<t_n).
 \tag{3}
\]
这里 \(C_n'\) 只表示一个有限维数常数；(3) 不识别八阶精确系数或其极限。

任意输入 \(E\) 有 \(V=|E|\in(0,\infty)\)。从全部正半径严格 halo 的外部部分取紧集 \(K\)，记 \(M=|K|\)、\(T=V+M\)，沿用原有限严格见证选族 \(Q_i,B_i,r_i\) 及 \(S=\sum r_i^n\)。只须研究任意满足
\[
 v\downarrow0,\qquad
 1+M/V\ge L_n(t)+K_nav^6-\beta_0av^7,\quad
 \beta_0=d_nK_n+1
 \tag{4}
\]
的输入序列。低于此分支的比值已经低于目标上包。所有 little-o 都沿任意这样的输入序列，不是只对固定输入成立。

在 §2 一次新建同一 \(I\) 的八层族。第八层的实际量记为
\[
 R_i=(g_{i,8}/a)^{1/n},\quad
 J_i=(R_i+|W_i|^{1/n})^n-(R_i+\mu_i^{1/n})^n,
\]
\[
 H_2=\sum_{q(e)=2}|X_i\cap X_j|-\sum_iJ_i,\qquad
 H_8=\sum_{i<j}|X_i\cap X_j|-\sum_iJ_i.
 \tag{5}
\]
源和 pair 都是该层全部外部实际对象，含边界、异常端及退化交长；每个 \(J_i\) 只减一次。

**唯一新增局部条件 \((\mathrm H_{45})\)。** 对上述新第八层及本稿随后建立的同层粗先验，R045 冻结稿的 (32′) 成立：正常集合 \(\mathcal N\)、实际槽及正常 clean 占用须与 R043 完全相同。记
\[
 \mathcal C_8=\frac14\sum_{\text{实际槽 }s}R_i^2|F_s|,\qquad
 D_{\mathcal N}=\frac14\sum_sR_i^2|F_s\setminus O_s|\ge0.
 \tag{6}
\]
其全部有界实际层 \(Y_{\rm all}\) 位于一个最后选定的固定大管 \(\mathcal C_i(K_{\rm all})\)，令 \(Z_8=Z(K_{\rm all})\)。存在固定 \(c_R,c_{\rm slot}>0\)，对任意预先指定 \(\chi>0\)，最后取小 \(v\) 时
\[
 \boxed{\ 
 H_2+c_RZ_8+c_{\rm slot}v^6D_{\mathcal N}
 \le\left[\frac{k_n}{(1+v)^2}v^6+\Gamma_nv^7\right]\mathcal C_8
       +\chi\mathcal F_8+C_\chi v^8S.\ }
 \tag{H45}
\]
这是同一次完整 CDF 盈余结算的联合价格，不是把 R042 减去全体 \(J\) 的 packet 价与另一份全盈余价相加。\(\chi\) 只改变小参数起效阈值；\(c_R,c_{\rm slot}\) 固定正值先保留。

条件包含每个固定 \(n\ge2\) 和新固定第八层的适用性；\(n=2\) 的桥接在 §4 明列。常数可依固定 \(n,A_8\) 及已固定局部阈值，不可依输入数目、半径尺度跨度或该序列未知的收敛率。

## 2. 只选一个 I，一次重建八层

先固定纯维数回收价格 \(c_{\rm rec}\)、边界价格 \(C_{\rm b}\) 及正预留，使原图选择给
\[
 C_*P_\partial\le c_{\rm rec}P_I/4,\qquad
 2nvS_I\le c_{\rm rec}P_I/4,\qquad C_*\ge C_{\rm b}+2.
 \tag{7}
\]
依既有顺序固定全图参数，用原 \(P,S\) 选择同一标记集 \(I\)。它不依赖源盒、CDF 或 R045 的长角阈值。所有层使用
\[
 K'=K\setminus\bigcup_{i\in I}Q_i,\qquad M'=|K'|,\qquad M_I=M-M'.
\]

在这个 \(I,K'\) 上，先固定八个有限深度所需的目标和中心邻域常数，统一为 \(A_8\)，再一次选择适用于所有八层及下述有限辅助工具的正强 CDF 斜率。最后才令 \(v\to0\)。边界未删除源保持原私有 \(P_i\)，深处后层取前层的中心邻域；所有源点保持原指定严格见证。不从已有七层中另挑一族，也不将旧第七层换名为第八层。

按原有限递推记
\[
 \mathcal D_h=\sum_i g_{i,h}-M',\quad
 \mathfrak d_h^0=M'-\sum_i m^0_{i,h},\quad
 \Delta_h^{\rm mass}=\sum_i(m^0_{i,h}-\mu_{i,h}),\quad 1\le h\le8.
\]
删除质量和 CDF 端移始终分开。精确恒等式和递推为
\[
 \sum_i(g_{i,h}-\mu_{i,h})
 =\mathcal D_h+\mathfrak d_h^0+\Delta_h^{\rm mass},
 \qquad \mathcal D_h\ge f_h\ge0,
 \tag{8}
\]
\[
 \mathfrak d_1^0\le f_0,\quad f_0=A_0vM',\qquad
 \mathfrak d_h^0\le f_{h-1}+aP_\partial\quad(2\le h\le8).
 \tag{9}
\]
小源删除不藏在 \(\mathfrak d_h^0\) 中。外部为空时所有外部源量取零，不定义空 CDF；内部精确计数仍给下文同方向不等式。

## 3. 新八层族上的先验、窗口与两个旧精度前缀

原保留正量粗递推逐层适用于任意预先固定有限深度。以共同左侧
\[
 Z_0=a[T-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+2nvS_I]
\]
表示，有固定正系数 \(c_0,b_0\) 及非负 \(\mathcal J_h\)，满足
\[
 Z_0\ge c_0f_h-b_0vf_{h-1}
       +a\mathcal J_h+\mathcal D_h/4
       -Ca(v^6S+vP_\partial).
 \tag{10}
\]
\(\mathcal J_h\) 保留源过量、形状、原免费面及每个预先固定稳健裕量的免费面正量。这个粗递推没有使用 (H45)。

分别对同一族前缀 \(1\le h\le m\)，\(m=6,7,8\)，使用 \((b_0v/c_0)^{m-h}\) 加权。初端除以 \(a\) 为 \(O(v^{m+1}S)\)。原六阶粗上界及 (4) 给 \(|T-A_vM/a|=O(v^6S)\)，共同图条件吸收 (10) 的边界价格，于是对每个 \(m=6,7,8\)：
\[
 P_I+\mathcal D_m/a+
 v\sum_i(\varepsilon_{i,m}+\mathfrak q_{i,m})r_i^n+
 \mathcal F_m+c_\tau\mathcal F_{\tau,m}=O(v^6S),
 \tag{11}
\]
\[
 \mathcal D_{m-1}/a+\mathfrak d_m^0/a+
 S_I+S(\partial I)=O(v^5S),\qquad
 \Delta_m^{\rm mass}/a=O(v^6S).
 \tag{12}
\]
同层宽度平方先验给三种盒相对宽误差的加权 \(L^1\) 为 \(O(v^{5/2}S)\)，端移总权为 \(O(v^6S)\)，原缺量满足 \(|\mathscr D_i|\le C\varepsilon_i ar_i^n\)。这些正是新第八层 R045、R043、R044 所需的定量粗先验。每个固定稳健面窗口的预算也在 (11) 中，而非仅对最初面邻居有效。

原输入到 LT 的证明、坏根构造和完整固定物理窗口现在使用这些**新层**实际盒。原全图与半盒未变；保留原中心严格分离，排除微小漏盒和远大漏盒的有限图球／扩窗证明对固定 \(A_8\) 仍成立。窗口匹配只须误差最终小于固定几何阈值，不要求其为 \(O(v)\)。坏根权仍为 \(O(v^{9/2}S)\)，固定窗口实际根传输为 \(O(v^6S)\)，全部清单数和重复次数统一有界。不能拿旧层的好坏根或 helper 清单直接套在新层上。

基本质量同时给
\[
 V\ge(1-a)S,\quad M\le aS,\qquad
 M/V=a[1+O(v)],\quad T/S=1+O(v).
 \tag{13}
\]
仅使用最后一式的领先精度，不提取未知的一阶 \(T/S\) 系数。

### 第六前缀

对当前同 \(I\)、同一前六层和统一 \(A_8\)/CDF，重放已审六阶局部价格及领先 incidence。它们只用固定有限源常数和 (11)–(12)，所以
\[
 H_6\le[K_n+o(1)]v^6S.
 \tag{14}
\]
小源和环境在此仅需粗 \(O(v^7S)=o(v^6S)\) 精度，完整源缺量继续保留，得到下文账本的
\[
 E_6=[K_n+o(1)]v^6S.
 \tag{15}
\]

### 第七前缀

同一新八层族的前七层重放 R038/R040 的**证明**。先保留第七层有符号 \(H_7-c\mathcal F_7\)，前六层用粗价和 (9) 抵消；(4)、(13)、(14) 及第七层领先价给
\[
 H_7/(v^6S)\to K_n,\qquad \mathcal F_7=o(v^6S).
 \tag{16}
\]
这先于任何七阶新价格，故无循环。随后在这同一第七层重做未加权实际原缺集轮廓、真实轴帽反馈、退化软权与管外高支撑支付；原证明只用固定 \(A_8\)、同层原指定见证、上述六阶饱和和原完整清单。所有目标质量在同一见证内相加，源盈余仅计一次。二维使用原唯一角槽／实际根 incidence，高维使用 adaptive-clean 四份分配。因此已审七阶局部价和精删除步骤供应
\[
 E_7=Q(v)T+e_7,\qquad e_7=o(v^7S),
 \quad
 Q(v)=K_nv^6\left[(1+v)^{-4}+\frac{10}{9}v(1+v)^{-2}\right].
 \tag{17}
\]
这里的对象是当前族第七前缀，不是旧第七层的对象。统一 CDF 斜率可以预先取得更小的固定正值，所有原尾斜率、宽度及误差常数仍允许依赖这个固定值。没有把 \(E_7\) 简化成未知 \(Cv^6S\)，也没有借全输入七阶定理反推某个实际源层的能量。

前五层只使用原粗并集接口，取 \(E_h=C_0v^6S\)，\(h\le5\)。各层所需基线系数最终统一增大为一个固定 \(b\)。

## 4. 第八层联合局部价及二维桥接

固定 R045 的共同正常集合、其小管、全域非正常强制性常数，然后依其顺序取 \(T,m,C_{\rm soft}\)、共同实际合并层和最后大管 \(K_{\rm all}\)。R043 的槽几何只需要 R045 (6) 的正常窗后果，故可用这个**相同**正常集合和相同 clean 分配。正常集合的中心 \(n/[3(1+v)]\) 不被换成另一个集合。

取
\[
 0<\eta<\frac14\min(c_R,c_{\rm slot}).
 \tag{18}
\]
在新第八层原地应用已审 R043/R044：好所有者高支撑由实际有限面对路径和原中心分离给 gross \(O(v^8S)\)；坏根权由真实缺槽量满足
\[
 S({\rm Bad})\le C[D_{\mathcal N}+v^{9/2}S].
 \tag{19}
\]
R044 的一般 \(q\) 矩形仍取同一层原有效盒和实际最大质量根所有者。目标截宽、目标缺量例外、面对／端移例外及总相交度与层号无关；固定管后取长尾阈值，不随 \(v\) 增长。其短窗三临界 \(C_Kv^7S({\rm Bad})\) 和长窗净余 \(C_\eta v^{15/2}S({\rm Bad})\) 由 (19) 在小 \(v\) 下被任意固定 \(\eta v^6D_{\mathcal N}\) 吸收，全部 \(q\ge4\) 与边界也包含在内。因此
\[
 \sum_{q(e)\ge3}I_e
 \le\eta Z_8+\eta v^6D_{\mathcal N}+C_\eta v^8S.
 \tag{20}
\]
这只花 R045 同一次余量结算留下的固定小份额；没有另减任何 \(J_i\)。

将 (20) 与 (H45) 相加，留下的 \(Z_8,D_{\mathcal N}\) 系数仍为正，故
\[
 H_8\le\left[\frac{k_n}{(1+v)^2}v^6+\Gamma_nv^7\right]\mathcal C_8
                 +\chi\mathcal F_8+Cv^8S.
 \tag{21}
\]

### \(n=2\) 不能越界引用高维源稿

二维没有 \(q\ge3\)，(20) 左侧为零，无需调用只为 \(n\ge3\) 陈述的 R043/R044。R045 条件必须包含其 §9 的二维独立重放：横向积为空、\(b_e=1\)，只存在同一轴对的四个实际角槽；同见证轴帽和合并目标质量、完整 \(J\) 的正剩余、长条容量及全域标量折扣均逐式成立。实际角位分配使用已审 R027 的二维唯一角槽和实际根权传输，不引用高维 helper 声明来跨维数。

二维全槽容量恰为 \(\mathcal C_8=\sum_iR_i^2\)，无横向误差；以下质量步骤同样给所需分母 \((1+v)^{-2}\)。因此在明确的二维 (H45) 下，(21) 对每个固定 \(n\ge2\) 成立。仅有旧二维 \(o(av^7)\) 结论本身不足以给本稿 \(O(av^8)\)；本稿实际使用的是新共同二轴价格的二维版本。

### 实际全槽容量上界

每轴对有四个符号槽。由新层宽度先验，逐真实乘积求和给
\[
 \mathcal C_8
 =\sum_i\sum_{q<l}R_i^2\prod_{h\ne q,l}(R_i+w_{ih})
 =B_n(1+v)^{n-2}\sum_iR_i^n+O(v^{7/2}S).
 \tag{22}
\]
精确质量及 (11) 给
\[
 \sum_iR_i^n=(M'+\mathcal D_8)/a
 \le T/A_v+O(v^6S).
 \tag{23}
\]
这里 \(M'/a\le M/a\le T/A_v+O(v^6S)\) 只用原粗临界比较。故
\[
 \mathcal C_8\le \frac{B_n}{(1+v)^2}T+O(v^3S).
 \tag{24}
\]
这控制的是完整实际槽容量，不把仅关于 \(W_\omega\) 的不等式反向当作容量上界。二维 (22) 的空积式精确成立。乘 (21) 的六阶主价，容量余项为 \(O(v^9S)\)，于是
\[
 H_8\le Q(v)T+\chi\mathcal F_8+Cv^8S.
 \tag{25}
\]

## 5. 只回收一次免费面，并同时支付删除和环境

用原与全部 \(X_{i,8}\)、标记部分分离的实际面条，在第八层并集上恢复一次，正体积至少
\[
 c_{\rm face}\mathcal F_8,\qquad c_{\rm face}>0.
 \tag{26}
\]
面条自身可重叠，但已按固定重数取得这个系数；不再另造第二份同样收入。

先在 (H45) 取 \(\chi<c_{\rm face}/8\)。对第八层深处小源，原扩展盒 \(\mathcal R_i\) 的入向扩张为 \(\gamma^\pm_{il}\le2\Gamma_{ip}\)，所有原指定小半径至多为 \(A_8vr_i\)。盒内 dyadic 估计保留原系数 \(a\)：
\[
 \delta^{\rm mass}_{i,\rm in}\le Ca|\mathcal R_i\setminus T_{i,8}|,
 \quad
 |\mathcal R_i\setminus T_{i,8}|
 \le C\varepsilon_i\mu_i+
       C(vr_i)^{n-1}\sum_{l,\pm}\gamma^\pm_{il}.
 \tag{27}
\]
同一固定 dyadic 网格的最大高缺量格互不交，覆盖所有被计小源；盒缺量为零时给零小源。此常数不依原半径层数，未将 \(a\) 粗化成 \(v\)。乘最终账本系数 \(v/a\) 后，
\[
 (v/a)\sum_{\rm deep}\delta^{\rm mass}_{i,\rm in}
 \le C a v\sum_i\varepsilon_i r_i^n+C v^n\mathcal F_8
 \le C v^{n+6}S+C v^n\mathcal F_8.
 \tag{28}
\]
第一项对 \(n=2\) 恰为八阶，对更高维更小。

盒外尾由原入向邻居的同一严格见证、与自身完整目标排他的小源尾给
\[
 (v/a)\sum_{\rm deep}\delta^{\rm mass}_{i,\rm out}
       \le C v\mathcal F_8.
 \tag{29}
\]
原端点以源点逼近和半空间质量连续性处理，无须统一严格 margin。面尺度可比及次数有界。

边界不用深面邻居：
\[
 (v/a)\sum_{\partial I}\delta_i^{\rm mass}
 \le C avS(\partial I)=O(v^{n+6}S).
 \tag{30}
\]
共同实际环境满足
\[
 \sum_i u_{i,8}\le Ca(\mathcal F_8+P_\partial).
 \tag{31}
\]
其中 \(aP_\partial=O(v^{n+6}S)\) 由 (7)、(11) 控制；固定突出边界价格 \(C_{\rm b}P_\partial\) 仍留给图条件，不能混入八阶小量。

完整系数先固定后，最后取 \(v\) 小，使
\[
 Cv\mathcal F_8\le c_{\rm face}\mathcal F_8/8,\quad
 Cv^n\mathcal F_8\le c_{\rm face}\mathcal F_8/8,\quad
 Ca\mathcal F_8\le c_{\rm face}\mathcal F_8/8.
 \tag{32}
\]
连同 \(\chi\mathcal F_8<c_{\rm face}\mathcal F_8/8\)，总共最多花一半收入，仍留正余量。这里分别支付 R045 长角面成本、盒外删除、盒内扩张和环境，没有让它们各自使用完整 (26)。\(\mathcal F_8=0\) 时所有这些面成本也为零。

在恢复后的真实并集上，Bonferroni 减完整 \(H_8\)，再用
\[
 (R_i+\mu_i^{1/n})^n
 \ge A_vg_i/a-C(v/a)(g_i-\mu_i)
 \tag{33}
\]
及 (8)。扣除实际重复质量 \(\mathcal D_8\) 和 (31)，重复净系数
\((A_v-Cv-a)/a\ge1/(2a)\)。用一次面收入吸收上述支出后，保留 \(\mathfrak d_8^0\) 不估掉，得到第八层实际并集接口，其误差可取
\[
 E_8=Q(v)T+e_8,\qquad 0\le e_8\le C_8v^8S.
 \tag{34}
\]
没有要求先证明 \(\Delta_8^{\rm mass}/a=O(v^7S)\) 或 \(\mathcal F_8=O(v^7S)\)。后者只需粗 \(O(v^6S)\)；盒外的 \(Cv\mathcal F_8\) 由正恢复吸收，不能直接丢入八阶。

## 6. 内部恒等式、边界预留与八层精确抵消

对每层分别加入同一标记半盒和内部突出集。共同真实输入体积只计一次，原内部精确式仍为
\[
 (1-a)S_I+M_I-A_vM_I/a
 =(1-a)(S_I-M_I/a)-(A_v-1)M_I/a.
 \tag{35}
\]
内部重数／共享质量给 \(c_{\rm rec}P_I\)，内部体积价格不超过 \(2nvS_I\)。以 \(c_{\rm mass}=1/2\)，八层分别得到同一左侧
\[
 T\ge A_vM/a+c_{\rm rec}P_I-C_{\rm b}P_\partial-2nvS_I
 +c_{\rm mass}\mathcal D_h/a-b(v/a)\mathfrak d_h^0-E_h .
 \tag{36}
\]
所有源常数造成的 \(CvP_\partial\) 在固定全部参数后由纯维数边界预留吸收。图 \(I\) 不须为 R045 重新选择。

各层不等式是同一体积 \(T\) 的独立有效下界；随后作正权平均，不是把八个不同物理并集相加成一个并集。每层面收入只在该层使用一次，故第七前缀的恢复不占第八层的 (26)。

置
\[
 \tau=bv/c_{\rm mass},\quad \Omega_8=\sum_{j=0}^7\tau^j,\quad
 Z=a[T-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+2nvS_I].
 \tag{37}
\]
用 (9)、\(\mathcal D_h\ge f_h\)，(36) 给
\[
 Z\ge c_{\rm mass}f_h-bvf_{h-1}-aE_h-abvP_\partial.
 \tag{38}
\]
乘 \(\tau^{8-h}\) 求和的精确恒等式为
\[
 \boxed{\ 
 \sum_{h=1}^8\tau^{8-h}(c_{\rm mass}f_h-bvf_{h-1})
 =c_{\rm mass}f_8-bv\tau^7f_0.\ }
 \tag{39}
\]
原源缺量被精确抵消，未将 \((v/a)\mathfrak d_h^0\) 先粗估成误差。初端除以 \(a\Omega_8\) 恰为
\[
 \frac{A_0b^8v^9(M'/a)}
      {c_{\rm mass}^7\Omega_8}=O(v^9S),
 \tag{40}
\]
其中 \(M'/a\le S\)。

设 \(\bar E=\Omega_8^{-1}\sum_h\tau^{8-h}E_h\)。必须先精确减同一 \(Q(v)T\)：
\[
 \bar E-Q(v)T
 =\frac{e_8+\tau e_7+\tau^2[E_6-Q(v)T]
       +\sum_{h=1}^5\tau^{8-h}[E_h-Q(v)T]}{\Omega_8}.
 \tag{41}
\]
(13)、(15)、(17) 给
\[
 E_6-Q(v)T=o(v^6S)+O(v^7S).
 \tag{42}
\]
于是 (41) 的四类项依次为 \(O(v^8S)\)、\(o(v^8S)\)、\(o(v^8S)\)、\(O(v^9S)\)，故
\[
 \bar E\le Q(v)T+C_8v^8S+o(v^8S).
 \tag{43}
\]
第七层与第六层分别保留正确七阶价和领先常数，才能抵消权归一化的低阶影响；只知它们是未知 \(Cv^6S\) 不足以推出 (43)。不需要 \(b\) 的精确数值。

图条件 (7) 吸收 \((C_{\rm b}+bv)P_\partial+2nvS_I\)，剩余内部量和 \(f_8\) 非负可丢。由 (38)–(43)，
\[
 \boxed{\ T\ge A_vM/a-Q(v)T-C_8v^8S-o(v^8S).\ }
 \tag{44}
\]
这一步保留了全部 internal、boundary 和源缺量原账本，没有把固定边界价格或七阶缺源成本隐入八阶误差。

## 7. 实际 M 的分母与全部输入上确界

令 (44) 的归一化余项为 \(r\)，则 \(|r|\le C v^8S+o(v^8S)\)。直接移实际 \(M\)，而不是先用未知 \(S/V\) 一阶项，得
\[
 1+M/V
 \le\frac{A_v}{A_v-a(1+Q(v))}
       +Cav^8\,S/V+o(av^8).
 \tag{45}
\]
由 \(S/V\le(1-a)^{-1}\)，后一项为 \(O(av^8)\)。精确差为
\[
 \frac{A_v}{A_v-a(1+Q)}-\frac{A_v}{A_v-a}
 =\frac{aA_vQ}{(A_v-a)(A_v-a-aQ)}.
 \tag{46}
\]
因 \(a=O(v^2)\)、\(Q=O(v^6)\)，替换第二分母产生的误差 \(O(a^2Q^2)\) 晚于 \(av^8\)。且
\[
 Q(v)=K_nv^6-\frac{26}{9}K_nv^7+O(v^8),\qquad
 \frac{A_v}{(A_v-a)^2}=1-nv+O_n(v^2).
 \tag{47}
\]
所以对任意 (4) 序列
\[
 1+M/V\le
 L_n(t)+K_nav^6-d_nK_nav^7+C_8'av^8+o(av^8).
 \tag{48}
\]
这里不需要求八阶精确系数；二维 \(a=v^2\) 对分母的影响正落在允许的八阶上包中。

从逐序列到全部输入：若 (2) 失败，则有 \(t_j\downarrow0\) 和有限正测度输入，使归一化余量
\[
 \frac{|\mathcal H_n(1-t_j;E_j)|/|E_j|-P_7(t_j)}{a_jv_j^8}
 \longrightarrow+\infty.
 \tag{49}
\]
可先按上确界选输入，再对外部严格 halo 作紧集内逼近，两个比值损失各取 \(o(a_jv_j^8)\)。输入密度点几乎处处在严格 halo，外部 halo 可测有限且有内正则性，故实际比值为 \(1+|{\cal H}\setminus E|/V\)。有限严格见证选择仍保留原全部正半径，不截尺度跨度。所得紧集序列最终满足固定分支 (4)，与 (48) 矛盾。

局部八阶常数只依固定先验常数和固定几何阈值；第六／七层 little-o 可随序列有不同模量，其在 (41) 的归一化贡献仍趋零。因此反证没有把只对固定输入成立的小量交换到上确界外，也没有预设一个未证明的统一收敛速度。由此在 (H45) 下得到 (2)。

## 8. 有限下界给正 Theta，而非精确八阶

已审 R041 对每个充分小 \(t\) 给真实有限输入，全部正半径严格比值满足
\[
 P_7(t)+J_nav^8+O_n(av^9),
 \qquad
 J_n=\frac{(n-1)n^3(2609+1971n+81n^2)}{26244}
             +\frac8{27}\mathbf1_{\{n=2\}}>0.
 \tag{50}
\]
故 \(C_n(1-t)-P_7(t)\ge(J_n/2)av^8\) 对所有充分小 \(t\) 成立。与条件上界 (2) 合并即 (3)，亦可记
\[
 C_n(1-t)=L_n(t)+K_nav^6-d_nK_nav^7+\Theta_n(av^8).
 \tag{51}
\]
本稿没有证明 \(J_n\) 是最优八阶系数，也没有证明该归一化余量有极限。若 R045 尚未独审通过，(2)、(3)、(51) 都保持上述条件状态，不能在摘要中无条件宣布。

## 9. 条件核合清单、常数顺序和计算

去条件化前须独立核对：

1. R045 最终冻结 SHA 的 (32′) 确由同一次完整 \(J\) 供应 \(Z_8\) 和 \(D_{\mathcal N}\)，且其两管顺序无 \(T,K_{\rm all}\) 循环。
2. 正常集合、actual root、槽占用、最后大管都与新第八层 R043/R044 相同；(20) 不从另一份 \(J\) 重新购买高支撑成本。
3. R045 的 \(n=2\) 独立桥接成立；不能仅引用 R042 高维声明或旧二维七阶 little-o。
4. 一次新建八层及其六、七前缀确供应 (11)–(17)，前缀使用原指定见证和统一新 CDF，不借旧层标签。
5. 一次 \(\mathcal F_8\) 物理恢复确按 §5 同时支付联合价、盒外删除、盒内扩张与环境；内部和固定边界仍按 §6 结算。

常数顺序为：固定 \(n,N=8\) 和纯维数回收／边界预留、原图参数及同一 \(I\)；固定全部八层有限深度和统一 \(A_8\)/强 CDF；在本族取得先验与固定窗口工具；固定正常窗、小管及全域非正常折扣；取 R045 的 \(T,m,C_{\rm soft},Y_{\rm all},K_{\rm all}\)；取 (18) 的高支撑份额及 \(\chi<c_{\rm face}/8\)；最后取小 \(v\)。辅助固定稳健裕量可在其允许范围内取得更小值，不依 \(v\) 变化。不允许增长深度或将常数解释为维数一致。

计算由 Terra 按分工完成，已交付 r046-computation/eight_layers.wl/log，复用 R038 的真实分母一般维小核。它只核精确望远镜、减 \(\Omega_8QT\) 后权误差和最终标量展开；源层、面收入、严格水平、n2桥及 (H45) 仍由上述解析证明／待审接口承担。没有为重复符号批次新增负载。

依赖实际读取版本：

- R043 正式 proofs/043-normal-slots.md：**75d0f89e79e8de4d00a6d3926532357f08a60a1ed5cbbc20b6abfd83e82ab4ac**。
- R044 正式 proofs/044-high-support-net.md：**ef3a6df6ee4581e990fd55fbe76c2b03f43b2e5bdbd0cb27f7eba18231cf0cc4**。
- R045 本次条件供货：**8f0ad8b01bb443789b33f55188c7a04d000a3ceaee3fd61d0a83eb458c9bf8ea**。
- R038 正式 proofs/038-general-seven-layer.md：**0c20c209c2e9c3ad3925fd0f78870e6d8af2edeff97261d799f6890844f5d2ce**，以及已独审 R040 的现行 seventh.tex。
- R041 正式 proofs/041-eighth-feedback-lower.md：**bba23ac1dba9cfe6c091426595fa162bebe3c585d3afce6d9d2968f497d11845**。

**最终交付口径：完整证明 \((\mathrm H_{45})\Rightarrow\)(2)、(3)。R045 未独审通过前条件不取消；本稿自身也须独立复核，不能仅因末步算术通过就宣布一般维八阶余项。**
