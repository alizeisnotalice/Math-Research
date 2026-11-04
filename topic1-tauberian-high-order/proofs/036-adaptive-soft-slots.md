# R036：退化角的自适应软权容量与主价格

日期：2026-11-05。完整新几何候选，待独立审查；七阶局部价格衔接仍为条件命题。只写独立 scratch。

几何部分精确依赖 R35 的实际 probe/helper 及跨锚点清单论证。本稿已完整读取 R35 独审报告，冻结 SHA `17cd4a7a82a58b273c536c5af0e080265232c469966f0013da9149d1b8953b22` 已判 PROOF-VALID。本轮新增的每边自适应阈值和软权抵消仍待本稿独审。本稿不假设“每固定阈值”自动具有阈值趋零的一致性；改用每边自身的真实交长设碰撞阈值。

## 1. 目标与已证输入

固定 \(n\ge3\)，\(a=v^n\)、\(A_v=(1+v)^n\)、\(V=|E|\)、\(M=|K|\)、\(S=\sum_i r_i^n\)，在原临界粗先验的同一固定末层工作。\(X_i\) 为实际 CDF 盒，\(R_i=(g_i/a)^{1/n}\) 为实际目标根；不改原源及指定严格见证。

取所有恰两临界轴、真实正交叠的外部无向边 \(e=\{i,j\}\)。也可只取任一子集，例如全部深处双角。两临界交长为 \(L_{e,1},L_{e,2}>0\)，其余轴交投影为 \(T_e\)、\(b_e=|T_e|\)。以下固定
\[
 R_e=\max(R_i,R_j),\quad
 p_e=L_{e,1}/(v^3R_e),\quad q_e=L_{e,2}/(v^3R_e),
 \quad m_e=\min(p_e,q_e).
\]
任意固定 \(A\ge1\) 定义
\[
 \omega_e=\min(1,A m_e),\qquad W_\omega=\sum_e b_eR_e^2\omega_e.
\]

**几何命题。** 在上述已有实际源与 R35 的已审几何输入下，
\[
 \boxed{W_\omega\le
       \frac{\binom n2}{(1+v)^2}(V+M)+O_{n,B,A}(v^3S).} \tag{1}
\]
因此乘任何 \(O(v^6)\) 局部价格后，容量误差是 \(O(v^9S)=o(v^7S)\)。这是全体正双角的带权界，允许 \(m_e\to0\)，不通过丢掉一个 \(o(S)\) 的角集合来损失六阶精度。

直接使用 R29 的真实量化结果：
\[
 |\mathcal Q\setminus U|=O((v^{n+1}+v^5)S),\quad
 \sum_{k<l}|X_k\cap X_l|=O(v^6S),\quad S(Z_v)=O(v^{9/2}S), \tag{2}
\]
以及固定窗口实际根传输
\[
 \sum_{i\in G_v}\sum_{k\in\mathcal I_i(R_*)}
      |R_i-R_k|r_i^{n-1}=O(v^6S),\qquad R_*=20. \tag{3}
\]
\(G_v=Z_v^c\)。好根窗口的相关指标数量有界、尺度一致趋一、全部源深处。相交 \(X\) 的度有纯维数上界。LT 的未定速率只用于固定几何清单，不作为面积误差。

## 2. 每边自身的 clean 阈值

对两端均好的边，以固定较小编号端点 \(i\) 为锚，沿临界两轴反射为右上布局。令
\[
 \ell_e=\min(L_{e,1},L_{e,2}),\qquad
 \lambda_e=\frac{\ell_e}{1000r_i}>0. \tag{4}
\]
R14 临界坐标界给 \(\ell_e\le C_Avr_i\)，因此对足够小 \(v\)，**所有边同时**满足 \(\lambda_e<\epsilon^2\)，其中 \(\epsilon\) 是 R35 固定 probe 的绝对常数。没有对 \(\lambda_e\) 设正下界。

使用 R35 的同两个实际 probe：
\[
 P_B=[1/2+\epsilon,1/2+2\epsilon]\times
                  [1/2-2\epsilon,1/2-\epsilon],\quad P_C=P_B^{\rm transpose},
\]
按 \(r_i\) 缩放并移到实际中心。其平面位置在 \(Q_i\) 内。所有根、轴对、符号总 probe 数统一为 \(8\binom n2\)。非覆盖集 \(V_{i,P}\) 与 R35 相同：probe 至少一半未被实际 \(U_y\) 覆盖。

碰撞集改成 \(C_i^{qr}(\lambda_e)\)：根的固定清单中某个不同胞对二维交面积至少为 \(\lambda_e r_i^2\)，且两胞 transverse 投影均含该 fibre。取
\[
 \Gamma_e=T_e\setminus
       (V_{i,P_B}\cup V_{i,P_C}\cup C_i^{qr}(\lambda_e)). \tag{5}
\]
这些都是有限 Borel 判据，边阈值是该有限边的常数，不造成选择可测性问题。

## 3. 条件 helper 几何如何适用于任意退化边

R35 的 helper 选择只用两件事：固定清单归一化误差小于 \(\epsilon/100\)，以及碰撞阈值小于 \(\epsilon^2\)。它未用阈值的下界。

具体地，probe 被实际盒覆盖时，覆盖者若不在正确 helper 列，则与一端胞的二维交面积严格大于 \(2\epsilon^2r_i^2\)；同列两个胞则交面积大于固定比例 \(r_i^2\)。这些都触发 \(C_i^{qr}(\lambda_e)\)。因此 (5) 上仍有唯一的两个实际 helper，并按有限清单指标形成 Borel 分片；四胞都确实包含该 transverse 点。此处使用 R35 的有限平面证明，而不是仅引用其固定 \(\tau\) 容量不等式。

需要重新核对的是**不同阈值下的角面单射**。若两条 clean 正双角占同一个实际胞角面点，R35 的跨锚点清单论证先迫使它们属于同一四胞的两个相反对角。不同列身份会形成固定大面积碰撞，已被任一边自身的 \(\lambda_e<\epsilon^2\) 排除。

取第一条边为 \(AD\)，令其锚尺度为 \(r_i\)。另一条 \(BC\) 只需交长为正，不要求与第一条有共同的 resolved 下界。实际端点有精确恒等式
\[
 L_{AD,x}+L_{BC,x}
 =(A_x^+-B_x^-)+(C_x^+-D_x^-).
\]
因 \(L_{AD,x}\ge\ell_e\)、\(L_{BC,x}>0\)，至少一条相邻面侵入量大于 \(\ell_e/2\)。其另一平面轴交长大于 \(r_i/2\)，故二维交面积大于
\[
 \ell_er_i/4>\lambda_er_i^2.
\]
所有四胞在第一清单并出现于该 fibre，矛盾。所以每个实际角面点至多收到一条 adaptive-clean 边；**不需要共同的正阈值**。

每条片仍分给四个实际角面位置，用同一个 transverse 点，不把片迁移到容量之外。完整原族和所有每胞基线项都保留。

## 4. 软权精确抵消碰撞收费

关键比率是
\[
 \boxed{\frac{\omega_e}{\lambda_e}
 \le\frac{A\ell_e/(v^3R_e)}{\ell_e/(1000r_i)}
 =\frac{1000A r_i}{v^3R_e}\le\frac{C_nA}{v^3}.} \tag{6}
\]
它对任意小正 \(\ell_e\) 都成立，饱和分支 \(\omega_e=1\) 也由 \(\omega_e\le A m_e\) 涵盖。

不能先对每边套一个常数可能不同的 \(O(1/\lambda_e)\) 再宣称一致；应先打开碰撞并取并的求和上界。设某清单 pair \(k,l\) 的平面交面积为 \(A_{kl}^{qr}\ge\lambda_er_i^2\)。轴平行性给体积精确分解。因此
\[
 \begin{aligned}
 R_e^2\omega_e
 |\pi_{\widehat{qr}}(X_k\cap X_l)|
 &=\frac{R_e^2\omega_e}{A_{kl}^{qr}}|X_k\cap X_l|\\
 &\le C\frac{\omega_e}{\lambda_e}|X_k\cap X_l|
 \le\frac{C_nA}{v^3}|X_k\cap X_l|.
 \end{aligned} \tag{7}
\]
每锚点只有固定数相交边；每个全局 pair 被固定窗口清单重复的根数也有统一上界。将 (7) 对边、清单 pair 和轴对求和，使用 (2)，得到
\[
 \sum_{e\text{ 两端好}}R_e^2\omega_e
      |T_e\cap C_i^{qr}(\lambda_e)|=O_{n,B,A}(v^3S). \tag{8}
\]
这里抵消发生在求和之前，故不出现退化角尾端的 \(1/m_e\) 损失。

其余费用更小：非覆盖 probe 由 \(\omega_e\le1\) 与 R29 的 Fubini 给 \(O((v^{n+1}+v^5)S)\)；坏端点边给 \(O(v^{9/2}S)\)；实际根权转移给 \(O(v^6S)\)。最后一项不受小阈值影响：\(\omega_e\le1\)，每片面积至多为固定倍 \(r_i^{n-2}\)，helper 分片互不交，故 (3) 逐式支付，全部重数如 R35 固定窗口论证。

在 clean 片上，每个实际角面点至多使用一次且 \(0\le\omega_e\le1\)，因此
\[
 4\sum_eR_e^2\omega_e|\Gamma_e|
 \le\sum_{c,qr,\pm,\pm}R_c^2|F_{c,qr}|+O(v^6S).
\]
加入所有已付费用，得到带权四份容量
\[
 \boxed{4W_\omega\le
 \sum_{c,qr,\pm,\pm}R_c^2|F_{c,qr}|+O_{n,B,A}(v^3S).} \tag{9}
\]

源宽 \(L^1\) 误差为 \(O(v^{7/2}S)\)，精确质量给
\[
 \sum_iR_i^n=(M'+\mathcal D)/a\le(V+M)/A_v+O(v^6S).
\]
所以与 R35 同一乘积计算将 (9) 化为 (1)。二维可直接用已审 R27 的 raw 实际根 incidence，因 \(W_\omega\le W\) 得同一方向、更小误差；无需把高维 helper 再用于二维。

## 5. 主净收益的逐边校准

这一节是严格标量结论，不代替实际角成本。设 \(s\in[1/2,2]\)，定义
\[
 \phi_{s,n}(p,q)=pq\left(1-\frac{s(p+q)}n\right),
 \quad k_{n,s}=\frac{n^2}{27s^2}.
\]
若固定 \(A\ge\max\{1,27/(2n)\}\)，则对全部 \(p,q\ge0\)
\[
 \boxed{[\phi_{s,n}(p,q)]_+
       \le k_{n,s}\min(1,A\min(p,q)).} \tag{10}
\]
证明：令 \(x=sp/n,y=sq/n\)。归一化左侧是
\(f=\max\{27xy(1-x-y),0\}\)。若 \(x+y\le1\)，三数 AM--GM 给 \(f\le1\)，否则 \(f=0\)。又
\[
 xy(1-x-y)\le xy(1-y)\le x/4,
\]
对称交换给 \(f\le(27/4)\min(x,y)\)。于是
\[
 f\le\min\{1,\tfrac{27s}{4n}\min(p,q)\}
 \le\min\{1,A\min(p,q)\}.
\]
证毕。

故对于实际边的明确主价多项式之和，有同时结算
\[
 \sum_e b_ev^6R_e^2\phi_{s,n}(p_e,q_e)
 \le k_{n,s}v^6W_\omega. \tag{11}
\]
这不是先把 \(m_e\) 小的边按 \(o(S)\) 删除；它逐边控制其全部正收益，并保留原始全部基线。

若 \(s=1+v\)，(1)、(11) 给
\[
 \sum_e b_ev^6R_e^2\phi_{1+v,n}(p_e,q_e)
 \le\mathfrak k_n\frac{v^6}{(1+v)^4}(V+M)+O_A(v^9S),
 \quad\mathfrak k_n=\binom n2\frac{n^2}{27}. \tag{12}
\]
它只对显示的主价多项式作结论，不假设真实目标退款能无条件塞进这个多项式。

### 有界 gross 小角分支

固定 \(T<\infty\)。若 \(p,q\le T\)、\(m<1/A\)，则不计任何源成本或反馈退款的 gross 也满足
\[
 pq\le Tm=(T/A)\omega.
\]
令 \(k_0=\inf_{s\in[1/2,2]}k_{n,s}=n^2/108\)。可选固定 \(A\ge4T/k_0\)，使 gross 至多为 \(k_0\omega/4\)。这给局部真实证明一个有余量的小角分支。\(m\ge1/A\) 时 \(\omega=1\)，任何已经证明的 resolved 精确系数保持原值，不因软化降低。

本段不处理 \(\max(p,q)>T\) 的大交长，更不擅自给其目标退款硬乘 \(m\)。该部分必须由局部真实成本另行支付。软容量不能代替这条几何论证。

## 6. 与拟议七阶局部价格的准确条件衔接

若局部实际源证明对同一边集、同一固定大 \(A\) 给出
\[
 H\le \frac{k_n}{(1+v)^2}v^6W_\omega
       +\frac{10}{9}k_nv^7W_\omega+o(v^7S),
 \qquad k_n=n^2/27, \tag{P}
\]
其中 \(H\) 是实际全部交损减完整实际源盈余，那么代入 (1) 严格得到
\[
 H\le\mathfrak k_n v^6(V+M)
       \left[(1+v)^{-4}+\frac{10}{9}v(1+v)^{-2}\right]
       +o(v^7S). \tag{13}
\]
容量误差乘主价是 \(O(v^9S)\)，可忽略于七阶。括号的一次展开为
\[
 1-\frac{26}{9}v+O(v^2).
\]
这是**以实际局部假设 (P) 为条件**的结算；本稿未证明 (P)，不以标量模型代替 R37 的目标退款、大角例外及余缺量处理，也未完成全输入的七层 \(E+M\) 组装。因此不由此宣布真实 Tauberian 下一系数是 \(-(n+26/9)\mathfrak k_n\)。

## 7. 可审范围与本机证据

本稿使用 R35 的已审 actual helper／清单证明，并逐项给出从固定阈值改为每边阈值后的单射加强。R29 的 void、固定窗口 actual root 相干性及碰撞计费亦已独审；本轮新自适应几何结论须另经独审。raw 无权 Slot 仍不是本稿结论。

实际读取 `r036-computation/soft_price.wl/log`。Mathematica 15.0 对 \(x\le y\) 的 \(f>\min(1,27\min(x,y)/4)\) 违反集返回 False，对称性覆盖全域；独立 \(f>1\)、线性界和四分之一界的违反集也均为 False。没有有限配置枚举或旧优化大表。软件另核了层积分 \(\int_{v^3}^1v^3/\theta\,d\theta=3v^3\log(1/v)\)，但本稿不使用 layer-cake，因此不以该积分替代几何的一致阈值量词。

严格证明是每边正阈值、\(\omega_e/\lambda_e\) 的求和前抵消、实际角面点单射、真实 root 权传输，以及 (10) 的标量推导；计算只核所列代数。
