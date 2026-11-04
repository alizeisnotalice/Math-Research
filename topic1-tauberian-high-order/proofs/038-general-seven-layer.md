# R038：一般维七层账本的条件性下一阶总装

日期：2026-11-05。完整条件候选，待独立审查。仓库及所有旧冻结稿只读。

本稿以已审 R30 的七层结构为参照，重新核对一般 \(n\ge2\) 的同层对象、小源系数和真实分母。R36、R37 的各自独审状态与本总装独立；以下仍以新第七层的完整局部价格为明确条件，待同层接口与本账本合审后再去条件化，不在本稿宣布一般维真实七阶定理。

## 1. 任意序列及精确末层条件

固定整数 \(n\ge2\)，记
\[
 H=2^n,\quad a=v^n=2^nt,\quad A_v=(1+v)^n,\quad
 L_n(t)=\frac{A_v}{A_v-a},\quad K_n=\binom n2\frac{n^2}{27}.
\]
任意有限正测度输入 \(E\)，\(V=|E|\)，外部严格 halo 中紧集 \(K\) 的质量为 \(M\)，写 \(T=V+M\)。原有限严格见证选族及半盒仍为 \(Q_i,B_i\)，\(S=\sum_i r_i^n\)。只研究任意满足
\[
 v\to0,\qquad 1+M/V\ge L_n(t)+K_n av^6-Bav^7,
 \qquad B<\infty\text{ 固定} \tag{1}
\]
的输入序列。全部 little-o 都沿任意这样的序列成立，可依赖该序列的收敛模量，不能只对固定 \(E\) 陈述后再取上确界。

在 §2 一次新建七层族。其第七层定义实际
\[
 R_i=(g_{i,7}/a)^{1/n},\quad X_i=x_i+W_{i,7}+[-R_i/2,R_i/2]^n,
\]
\[
 J_{i,7}=(R_i+|W_{i,7}|^{1/n})^n-(R_i+\mu_{i,7}^{1/n})^n,
 \quad H_7=\sum_{i<j}|X_i\cap X_j|-\sum_iJ_{i,7}.
\]
和取全部外部 pair，包含边界及异常端点，盈余为完整实际 \(J\)。

**唯一外接的局部条件 \(\mathrm H_7\)。** 对这同一新第七层及任意 (1) 序列，
\[
 \boxed{H_7\le Q(v)T+o(v^7S),\qquad
 Q(v)=K_n v^6\left[(1+v)^{-4}+\frac{10}{9}v(1+v)^{-2}\right].} \tag{H7}
\]
当 \(n=2\) 时，已审 R28 的真实局部价格与 R27 实际半径 incidence，经 R30 的新固定第七层核合，供应此形式；其 \(K_2=4/27\) 正是旧二维接口。当 \(n\ge3\) 时，R36 软容量结合 R37 真实软价格供应此形式。两份 R37 独审已分别签署同一新固定第七层的前置轮廓和完整价格；本稿只记录供货来源，保留 (H7) 为待本次总审核合的条件。

供货必须覆盖**预先固定七个深度后的新末层**，所有几何常数固定后才令 \(v\to0\)。只在旧六层对象上证明价格、不确认其对任意预先固定有限末深度（至少本次七层）及新的共同 CDF 斜率适用，不足以满足 (H7)。§9 列出逐项核合口径。

**条件总装定理。** 若 (H7) 成立，则
\[
 \limsup_{t\downarrow0}
 \frac{C_n(1-t)-L_n(t)-K_nav^6}{av^7}
 \le-\left(n+\frac{26}{9}\right)K_n. \tag{2}
\]
此处严格是一条条件性全部输入结论。

## 2. 先选同一 I，再一次构造七层

先固定纯维数内部及边界价格
\[
 c_{\rm rec}=1/(2H),\quad C_{\rm b}=2C_{\rm pr}(n)+2,
 \quad C_*=C_{\rm b}+2,
\]
\(C_{\rm pr}\) 是原突出几何常数，与 Tauberian 系数区分。原全图参数 \(K_0,\eta_{\rm w},D\) 依既有次序固定。标记 \(I\) 只用原 \(P,S\) 和完整图，取第一个满足加权边界条件的层，故
\[
 C_*P_\partial\le c_{\rm rec}P_I/4,
 \quad 2nvS_I\le c_{\rm rec}P_I/4. \tag{3}
\]
它不依赖源盒、深度、CDF 斜率或局部价格阈值。

同一 \(K'=K\setminus\bigcup_IQ_i\)、\(M'=|K'|\)、\(M_I=M-M'\) 用于全部七层。逐层固定深目标系数和下一中心邻域系数，一次执行 R14 有限递推至七层；统一为 \(A_7\) 后再选正强 CDF 斜率，最后取小 \(v\)。边界未删除源始终是原私有 \(P_i\)，深处在 \(h\ge2\) 用前一中心邻域。所有点保持原指定严格见证。

写
\[
 \mathcal D_h=\sum_i g_{i,h}-M',\quad
 \mathfrak d_h^0=M'-\sum_i m^0_{i,h},\quad
 \Delta_h=\sum_i(m^0_{i,h}-\mu_{i,h}).
\]
\(\Delta_h\) 是真实小源删除质量，不是 CDF 端移。精确有
\[
 \sum_i(g_{i,h}-\mu_{i,h})=\mathcal D_h+\mathfrak d_h^0+\Delta_h, \tag{4}
\]
\[
 \mathcal D_h\ge f_h\ge0,\quad
 \mathfrak d_1^0\le A_0vM',\quad
 \mathfrak d_h^0\le f_{h-1}+aP_\partial\quad(2\le h\le7), \tag{5}
\]
其中 \(f_h=M'-\sum_i|K'\cap Q(x_i,L_hvr_i)|\)，\(f_0=A_0vM'\)。小源删除不藏入 \(\mathfrak d_h^0\)。

## 3. 同一七层族的两个前缀先验

R17 保留正量粗递推对全部七层同时成立：
\[
 Z_0\ge c_0f_h-b_0vf_{h-1}+a\mathcal J_h+\mathcal D_h/4
                 -Ca(v^6S+vP_\partial). \tag{6}
\]
\(\mathcal J_h\) 含固定正比例的源过量、形状余量、原免费面量和每个固定稳健裕量的免费面量。\(Z_0\) 是同一 \(a[T-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+2nvS_I]\)。

分别对 \(1\le h\le m\) 用 \((b_0v/c_0)^{m-h}\) 加权，\(m=6,7\)。这两个前缀针对同一 \(I\) 及实际同一前六层，不重新选更有利的集合。初端除以 \(a\) 分别为 \(O(v^7S)\)、\(O(v^8S)\)。共同图条件 (3) 吸收价格后，(1) 与已证粗上界给 \(|T-A_vM/a|=O(v^6S)\)，从而对 \(m=6,7\) 有
\[
 P_I+\mathcal D_m/a+
 v\sum_i(\varepsilon_{i,m}+q_{i,m})r_i^n
 +\mathcal F_m+c_\tau\mathcal F_{\tau,m}=O(v^6S), \tag{7}
\]
\[
 \mathcal D_{m-1}/a+\mathfrak d_m^0/a+S_I+S(\partial I)=O(v^5S),
 \quad \Delta_m/a=O(v^6S). \tag{8}
\]
同层宽度、端移及稳健实际根差因此可分别用于第六、第七层。R19 的原输入 LT 及六阶 incidence，以及已审六阶真实源价格，对每个这样的固定末层给
\[
 H_m\le[K_n+o(1)]v^6S,\qquad m=6,7. \tag{9}
\]
它仅需原有限固定源常数，不要求 (H7)。旧证明的固定份额和阈值按先固定后 \(v\to0\)、最后份额趋零执行。

补充一个有用但不可过度使用的改进：同一第六前缀的 \(\mathcal D_6/a=O(v^6S)\) 及 (5) 实际还给 \(\mathfrak d_7^0/a=O(v^6S)\)。它的成本仍可正好达到七阶 \(O(v^7S)\)，故最终仍须精确抵消，不能扔进 little-o。

基本容量还给
\[
 M/V=a[1+O(v)],\qquad S/V=1+O(v),\qquad S/T=1+O(v). \tag{10}
\]
这里只作领先比较，不提取未知的一阶 \(S/T\) 系数。

## 4. 先证明末层免费面饱和，不用 H7

在第七层实际 \(\bigcup X_i\) 上恢复原免费面条。它们避开全部 \(X_i\) 及标记部分，有界重数给正增量 \(c\mathcal F_7\)。共同环境成本至多 \(Ca(\mathcal F_7+P_\partial)\)，粗先验使其为 \(O(v^{n+6}S)\)；小源原成本是 \(O(v^7S)\)。固定正比例的 \(\mathcal F_7\) 可保留。

此时精确保留有符号 \(H_7\)。前六层使用粗 \(Cv^6S\)，末层误差取 \(H_7-c\mathcal F_7+O(v^7S)\)，正权同左侧七层相加合法。源缺量按 (5) 抵消，早层总权 \(O(v)\)、初端 \(O(v^8S)\)，图项按 (3) 吸收，得到
\[
 T\ge A_vM/a-\frac{H_7-c\mathcal F_7}{1+O(v)}-O(v^7S). \tag{11}
\]
(1) 另一方面给
\[
 A_vM/a-T\ge(A_v-a)(K_n-Bv)v^6V=[K_n+O(v)]v^6V.
\]
用 \(V\ge(1-a)S\) 及 (9)，(11) 推出
\[
 \boxed{H_7/(v^6S)\to K_n,\qquad \mathcal F_7=o(v^6S).} \tag{12}
\]
这是 R37 在新末层使用的六阶饱和与面严格小量，完全在 (H7) 之前建立；没有循环。

## 5. 一般维真实小源删除的七阶精度

以下只讨论第七层。深处扩展盒 \(\mathcal R_i\) 各边是固定倍 \(vr_i\)，侧扩张量为入向 \(\gamma_{il}^{\pm}\le2\Gamma_{ip}\)。令 \(D_i=\mathcal R_i\setminus T_{i,7}\)，实际小源分盒内、盒外两部分。

### 盒内保留 a 的 dyadic 引理

若盒 \(\mathcal R\) 各边在 \([cL,CL]\)，\(T\subset\mathcal R\setminus E\)，每个被计 \(x\in T\) 有原严格见证半径 \(r\le AL\)，则
\[
 |S_{\rm small}|\le C_{n,c,C,A}a|\mathcal R\setminus T|. \tag{13}
\]
简证：同时二分各边得到一个固定可数 dyadic 网格。对每个原见证取含 \(x\) 的内格 \(P\subset\mathcal R\cap Q(x,r)\)，使 \(|P|\ge\beta r^n\)，其中可取
\(\beta=\min\{(c/(2C))^n,(c/A)^n\}>0\)。若 \(r\) 大于原盒最大边，直接取根盒；否则选最大边刚小于 \(r\) 的格层。排除统一零测网格边界后，严格见证及 \(T\subset E^c\) 给
\(|P\cap T|<a r^n\le(a/\beta)|P|\)。在所有 \(D=\mathcal R\setminus T\) 占比严格大于 \(1-a/\beta\) 的格子中，按包含关系取最大者；每格到根盒只有有限个祖格，它们不交至边界且覆盖全部被计点。令 \(\eta=a/\beta<1/2\)，其并 \(H_D\) 满足
\(|H_D\setminus D|\le\eta|H_D|\le\eta|D|/(1-\eta)\)，即 (13)。\(|D|=0\) 亦给零质量。常数独立输入、半径层及点。

在本层取 \(L=vr_i\)、\(T=T_{i,7}\cap\mathcal R_i\)，保留 (13) 的 \(a\)，不要把它粗化为 \(v\)。原端实际缺量和扩张体积给
\[
 \delta_{i,\rm in}/a\le C|D_i|,
 \quad |D_i|\le C\varepsilon_i\mu_i+
             C(vr_i)^{n-1}\sum_{l,\pm}\gamma_{il}^{\pm}.
\]
因此
\[
 \sum_{\rm deep}\delta_{i,\rm in}/a
 \le C\left[a\sum_i\varepsilon_i r_i^n+v^{n-1}\mathcal F_7\right]. \tag{14}
\]

### 盒外尾及边界

原入向邻居的同一个严格见证、与其自身目标排他的小源尾，给
\[
 \delta_{i,\rm out}/a\le Cr_i^{n-1}\sum_{l,\pm}\gamma_{il}^{\pm},
 \qquad\sum_{\rm deep}\delta_{i,\rm out}/a\le C\mathcal F_7. \tag{15}
\]
端点通过源点逼近和半空间质量连续性取极限，不要求严格余量统一。相邻尺度可比和有限面重数使入向 \(\gamma\) 总量确由原 \(\mathcal F_7\) 支付。

边界不使用这个深处论证。原局部小源估计给
\[
 \delta_i/a\le Ca r_i^n,
 \quad \sum_{i\in\partial I}\delta_i/a
       \le CaS(\partial I)=O(v^{n+5}S). \tag{16}
\]
局部输入零测时小源为空；不留退化分支。

由 (7)、(12)，(14) 的首项为 \(O(v^{n+5}S)\)，次项为 \(o(v^{n+5}S)\)；(15) 是 \(o(v^6S)\)。所有 \(n\ge2\) 都有
\[
 \boxed{\Delta_7/a=o(v^6S),\qquad (v/a)\Delta_7=o(v^7S).} \tag{17}
\]
共同环境另外满足
\[
 \sum_i u_{i,7}\le Ca(\mathcal F_7+P_\partial)
       =O(v^{n+6}S)=o(v^7S). \tag{18}
\]
特别 \(n\ge3\) 至少为九阶，二维也为八阶。固定边界突出价格仍保留，不混入这些小量。

## 6. 末层与第六层实际并集接口

末层取 \(U_7=\bigcup X_{i,7}\)，不加 §4 为饱和而恢复的面条。Bonferroni 及完整鲜盈余给
\[
 |U_7|\ge\sum_i(R_i+\mu_{i,7}^{1/n})^n-H_7.
\]
用标量基线
\[
 (R_i+\mu_i^{1/n})^n\ge A_vg_i/a-C(v/a)(g_i-\mu_i),
\]
代入 (4)，再扣物理 \(\mathcal D_7+\sum_i u_{i,7}\)。重复量净系数 \((A_v-Cv-a)/a\ge1/(2a)\)，删除和环境按 (17)、(18) 为 little-o。因此在 (H7) 条件下
\[
 |U_7|-\mathcal D_7-\sum_i u_{i,7}
 \ge A_vM'/a+\mathcal D_7/(2a)-b(v/a)\mathfrak d_7^0-E_7,
 \quad E_7=Q(v)T+e_7,\ e_7=o(v^7S). \tag{19}
\]
可增大 \(e_7\) 为非负误差。\(\mathfrak d_7^0\) 保留到望远镜。

同一七层族第六前缀有 (7)--(9)。其小源原成本只需 \(O(v^7S)=o(v^6S)\)，故已证领先价格给同样并集接口而
\[
 E_6=[K_n+o(1)]v^6S. \tag{20}
\]
这里没有调用第六层七阶价。前五层取 \(E_h=C_0v^6S\)。全部层共同采用一个固定足够大的 \(b\)。

## 7. 内部恢复、七层抵消及权归一化

对每层并集加入同一标记原半盒和内部突出集。真实并集上界只计输入一次，为
\(V+M'+\mathcal D_h+\sum_i u_{i,h}+aS_I\)。内部重数及共享质量给
\[
 S_I-M_I/a\ge\tfrac12\sum_I s_ir_i^n-P_\partial/2,
 \quad M_I/a\le S_I,
\]
并有精确恒等式
\[
 (1-a)S_I+M_I-A_vM_I/a
 =(1-a)(S_I-M_I/a)-(A_v-1)M_I/a.
\]
纯维数突出价格与它共同得到每层同一左侧
\[
 T\ge A_vM/a+c_{\rm rec}P_I-C_{\rm b}P_\partial-2nvS_I
      +c_D\mathcal D_h/a-b(v/a)\mathfrak d_h^0-E_h,
 \quad c_D=1/2. \tag{21}
\]
源常数造成的 \(CvP_\partial\) 在所有参数固定后由纯维数预留吸收。

置
\[
 Z=a[T-A_vM/a-c_{\rm rec}P_I+C_{\rm b}P_\partial+2nvS_I],
 \quad\tau=bv/c_D,\quad \Omega=1+\tau+\cdots+\tau^6.
\]
由 (5)、(21)
\[
 Z\ge c_Df_h-bvf_{h-1}-aE_h-abvP_\partial.
\]
按 \(\tau^{7-h}\) 加权，精确有
\[
 \sum_{h=1}^7\tau^{7-h}(c_Df_h-bvf_{h-1})
       =c_Df_7-bv\tau^6f_0. \tag{22}
\]
末层 \(-bvf_6\) 与第六层正重复量完全消去；初端除 \(a\Omega\) 为
\[
 \frac{bv\tau^6f_0}{a\Omega}
   =\frac{A_0b^7}{c_D^6\Omega}v^8M'/a=O(v^8S). \tag{23}
\]

设 \(\bar E=\Omega^{-1}\sum_{h=1}^7\tau^{7-h}E_h\)。由 (10)、(20)
\(E_6-Q(v)T=o(v^6S)\)。准确比较为
\[
 \begin{aligned}
 \bar E-Q(v)T=\Omega^{-1}\big[&e_7+
       \tau\{E_6-Q(v)T\}\\
       &+\sum_{h=1}^5\tau^{7-h}\{E_h-Q(v)T\}\big]
       =o(v^7S).
 \end{aligned} \tag{24}
\]
前五层总权 \(O(v^2)\)，成本 \(O(v^8S)\)；第六层的一阶领先价格抵消第七层除以 \(1+\tau+O(v^2)\) 的一阶变化。不能把第六层只记作未知 \(O(v^7S)\)，也不需要 \(b\) 的精确值。

除以 \(a\Omega\)，用 \(f_7\ge0\)、\(bv\le1\) 及同一 (3) 吸收边界和内部体积价格，仍余非负 \(c_{\rm rec}P_I/2\)。从而
\[
 \boxed{T\ge A_vM/a-Q(v)T-o(v^7S).} \tag{25}
\]
零外部时取零源而不定义空 CDF，内部计数自身给同方向；没有遗漏此分支。

## 8. 全维真实分母与七阶系数

\(S/V\le(1-a)^{-1}\) 把 (25) 的误差换成 \(o(v^7V)\)。令 \(x=M/V\)，则
\[
 (1+Q(v))(1+x)\ge(A_v/a)x-o(v^7).
\]
分母 \(A_v-a(1+Q(v))\to1\) 正且远离零，直接移实际 \(M\) 得
\[
 1+x\le\frac{A_v}{A_v-a(1+Q(v))}+o(av^7). \tag{26}
\]
精确差式为
\[
 \frac{A_v}{A_v-a(1+Q)}-\frac{A_v}{A_v-a}
 =\frac{aA_vQ}{(A_v-a)(A_v-a-aQ)}. \tag{27}
\]
因 \(a=v^n=O(v^2)\)、\(Q=O(v^6)\)，替换第二分母的误差 \(O(a^2Q^2)\) 是 \(o(av^7)\)。只须一般维一次项：
\[
 Q(v)=K_nv^6-\tfrac{26}{9}K_nv^7+O(v^8),\qquad
 \frac{A_v}{(A_v-a)^2}=1-nv+O(v^2). \tag{28}
\]
因此
\[
 \boxed{1+M/V\le L_n(t)+K_nav^6
       -\left(n+\frac{26}{9}\right)K_nav^7+o(av^7).} \tag{29}
\]
二维代入得到 \(-176/243\)，与 R30 已审结算相同。真实 \((1+v)^n\) 分母的 \(-nK_n\) 没有被二维的 \(-2K_n\) 取代；没有用未知一阶的 \(S/V\) 直接读系数。

若 (2) 失败，沿违反序列按上确界选有限输入、再对外部严格 halo 取紧集，各损失一个足够小的比值 \(\epsilon av^7\)。所得序列满足 (1) 的某个固定 \(B\)，与 (29) 矛盾。严格 halo 的可测有限性及输入几乎处处进入由已证性质保证；紧集不需是 \(E^c\) 中的开集。全部原半径层和尺度跨度均保留，没有新增等尺度／周期假设。

## 9. 供货核合与常数顺序

(H7) 必须逐项交付如下对象与量词：

1. 同一新七层递推的 \(I,K',F_{i,7},G_{i,7},W_{i,7}\)、原指定见证及实际 \(R_i\)，不借用旧六层更有利对象。
2. 全部外部交对减完整 \(J_{i,7}\)，包括边界、混合高支撑及退化尾；若先给深处软权价，须同时给出其他交损的合法支付，再由 R36 在同层给环境价格。
3. (12) 的六阶饱和及免费面严格小量由本稿先给，(7)、(8)、宽度／端移及未删除源排他亦已在同层建立。不能用 (H7) 反推这些前置量。
4. 若 R37 调用 R34 的二轴管外支付，需确认其原证明只用已固定有限定位常数、同层真实源和上述六阶饱和，故可在任意预先固定末深度重放。只写“把 6 改成 7”不算供货；本稿条件定理不预支此核合已经完成。
5. 所有 \(o(v^7S)\) 沿任意 (1) 输入序列成立；先固定有限末深度及全部源常数，再统一 CDF 斜率，随后各局部阈值按其证明次序固定，最后 \(v\to0\)。

总装的常数次序为：纯维数回收与边界预留、图参数 → 只由原 \(P,S\) 选同一 \(I\) → 七层全部常数与 \(A_7\) → 正强 CDF 斜率 → 六／七前缀先验 → 先证明领先饱和 → dyadic 精删除量 → 假设 (H7) → 七层精确消去及真实分母。源常数只影响最后起效阈值，不反馈到纯维数边界价格。

本稿外接而未自行重证的接口只有新第七层的 (H7) 及上述对象适用性。R36 与 R37 原稿各自已有正式独审；本条件稿仍等待对本次同层供货、一般维小源及最终总装的独立核合。一般维七阶真实结论须在该核合后才能去条件化。

依赖版本：R36 `adaptive_soft_slots_r036.md`，SHA256 `c6b761427b26c1ab88d2ecf17f47ef856cb7c9565d53dc980bcb78fbd6e75c3d`；R37 `general_soft_price.md`，SHA256 `74c2cedf52300d117ac8f8e9c2e06a3ec95cd9fd89d5d5ea9b4b3f78a6c6d853`。已实际全文读取 `r037-review-rigidity/audit.md` 与 `r037-review-price/audit.md`，后者 §6 在新固定末层逐项重建 R31、R34；这不涵盖随 \(v\) 增长的层数。二维供货与七层账本参照已审 R30 SHA256 `79378aa279dbb7d90ef26ceb7597b2b141bad7a756c074702687c38fd08fcd93`。

## 10. 实际符号核验

已实际读取 Terra 的 `r038-computation/general_seven.wl/log`，Mathematica 15.0 通式核验：七层望远镜残差为零；\(f_0=A_0vM'\) 的归一化初端恰为 \(A_0b^7M'v^8/(ac_D^6)\)（再除 \(\Omega\)）；前五层权从 \(v^2\) 起；\(T/S=1+dv\) 模型中权平均减 \(Q(v)T\) 到七阶为零。

最终精确差式残差为零。为处理全部固定 \(n\ge2\)，用 \(a=v^2w\)、\(0<w\le1\) 保留一般维分母，一次项结果为 \(\{K_n,-(n+26/9)K_n\}\)，比较返回 True。实际 \(a=v^n\) 对应 \(w=v^{n-2}\in(0,1]\)，故没有用有限维枚举代替通式。数值／符号证据不承担 dyadic、严格水平、源层或本稿外接 (H7) 的证明。
