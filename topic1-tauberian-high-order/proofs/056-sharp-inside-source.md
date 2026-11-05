# R056：二维真实 inside 小源的领先常数一

2026-11-05。新引理及同层应用，待独立严格审查。只写本独占 scratch，不改旧源或仓库。不据此宣布二维八阶最优系数。

## 1. 目标和结论

固定二维 \(a=v^2=4t\)。在已审 R049 同一新八层竞争族中，记深处未删源 \(T_i\)、有效大见证源 \(F_i\)、原扩展盒 \(\mathcal R_i\)、正常实际层 \(Y_{i,N}\)，及
\[
 \delta_{i,\rm in}
 =|\{x\in T_i\cap\mathcal R_i:r(x)<\tau_i v r_i\}|,
 \qquad 0<\tau_i\le A.
\]
每个 \(r(x)\) 都保留原严格见证；它可跨越源原盒、宏观正常栏和 \(\mathcal R_i\) 的边界，不额外假设见证包含于某个缺集。

本稿证明的候选精确回收式为
\[
 \boxed{\quad
 \sum_{i\ {\rm deep}}\delta_{i,\rm in}
 \le\frac{a}{1-a}\sum_{i\ {\rm deep}}|Y_{i,N}\setminus T_i|
                  +C a^2v^7S.\quad}
\tag{1}
\]
因此
\[
 \frac va\sum_{i\ {\rm deep}}\delta_{i,\rm in}
 \le [1+O(a)]v\sum_i|Y_{i,N}\setminus T_i|+O(av^8S).
\tag{2}
\]
这里第二项在二维为 \(O(v^{10}S)=o(v^8S)\)。再由同层正常窗和 R049 的参数稳定性，
\[
 \sum_i|Y_{i,N}|\le[2K_2+O(v)]av^5S,
 \qquad K_2=4/27,
\tag{3}
\]
故
\[
 \boxed{\quad
 \frac va\sum_{i\ {\rm deep}}\delta_{i,\rm in}
 \le\frac8{27}v^8S+O(v^9S).\quad}
\tag{4}
\]
这是一个实际小源上界，系数只支付一次删除费用，不把正常层的完整 J 当作额外收入。它与 R054 的新增微孔下界机制同阶；它既不单独完成正常 packet 的八阶上价格，也不排除其他八阶优化机制。

## 2. 有限值原半径的真正全点选择

以下 elementary 选择会在局部引理中用到。设有界集 \(A\) 上给定一个取有限正值
\(\rho_1>\cdots>\rho_m>0\) 的原半径函数 \(r(x)\)。从**整个点集**作如下操作，而不是先固定一个有限 probe／subcover 的中心集：

依次处理 \(A_j=\{x\in A:r(x)=\rho_j\}\)。在该层尚未被之前全部闭方体覆盖的点集中，选任一点 \(x\)，加入原方体 \(Q(x,\rho_j)\)，重复到本层没有未覆盖点。

任何新点不在此前的闭方体内，所以任意两所选中心满足
\[
 \|x_k-x_l\|_\infty>\max(r_k,r_l).
\tag{5}
\]
不能无限选点：所有半径至少 \(\rho_m\)，这些中心位于有界集，其半径 \(\rho_m/2\) 的盒内部两两不交，故选中数有一个有限 packing 上界。若一层还存在未覆盖点，就可继续选择；超过该上界矛盾。因此每层经过有限次选择后确实穷尽，得到覆盖**全部 A**的有限族。

此证明没有从“某个有限子覆盖的中心已经被覆盖”推出其原开集仍覆盖 A。这里每一步都从尚未覆盖的 A 本身取点，半径原值不变，有限终止给真正全点覆盖。无需统一严格余量、原子化质量误差或改变任何半径。

本项目的实际原函数确有上述有限正值性质：正式 main 的 `lem:finitegreedy` 和 global 开头先在紧 halo 集上选有限严格半径列表，再取列表中实际可用的最大值。所有新 N8 的源／删除集合只限制这个同一原函数，不另换半径。本稿的局部源集不必紧或开；有界、有限值即足够。即使最小半径随 v 或输入趋零，每次固定输入仍为正，packing 只用于有限终止，不进入价格常数。

## 3. 同向锚定方块保持强分离

固定符号 \(\sigma=(\sigma_1,\sigma_2)\in\{\pm1\}^2\)。对满足 (5) 的中心方体定义边长为 r 的锚定方块
\[
 P(x,r;\sigma)=\prod_{h=1}^2
 \begin{cases}[x_h,x_h+r],&\sigma_h=+1,\\[x_h-r,x_h],&\sigma_h=-1.\end{cases}
\tag{6}
\]
它包含于原 \(Q(x,r)\)，且面积为 \(r^2\)。同一符号下这些 P 两两不交：若中心分离由第 h 轴实现，则该轴中心距大于两个半径最大值，两个同向长度分别为 r、r′ 的区间也不交。不要求半径相等。

注意这不是把不交的原半盒平移后仍直接断言不交；它使用了更强的实际中心分离 (5)，逐坐标重新证明。

## 4. 一个独立的可测局部盒引理

设轴平行矩形 \(\mathcal R\) 各边介于 \(cL,CL\)，\(c,C,A>0\) 固定。设
\(T\subset\mathcal R\setminus E\) 可测，\(Z\subset\mathcal R\) 可测，且 Z 位于四个角各边长 \(\epsilon L\) 的角块之并。写
\[
 D=\mathcal R\setminus T,
 \quad z=|Z\setminus T|,
 \quad d=|D\setminus Z|,
 \quad |Z|\le bL^2,\quad d\le bL^2.
\tag{7}
\]
考虑任意可测 \(\mathcal S\subset T\)，其原严格见证取有限个正值、满足 \(r(x)\le AL\)。令 \(a=4t\to0\)，参数满足
\[
 b\to0,\qquad \epsilon\to0,
 \qquad\sqrt b/\epsilon\to0.
\tag{8}
\]
则全部常数固定、参数充分小时
\[
 \boxed{\quad |\mathcal S|\le\frac{a}{1-a}z+C_{c,C,A}a d.\quad}
\tag{9}
\]
Z 不需要包含于 D；因此允许 Z 中一部分本身属于小源 T。原输入在盒外可任意，不作“盒外没有输入”或“见证全在盒内”的假设。

### 4.1 原半径在小缺量盒中自动很小

由实际严格见证和 \(T\subset E^c\)，对每个原点 x 有
\[
 |Q(x,r)\cap T|<a r^2.
\tag{10}
\]
设 \(\ell_{\min}\) 为 \(\mathcal R\) 的最短边。若 \(r\ge\ell_{\min}\)，中心在盒内给
\(|Q(x,r)\cap\mathcal R|\ge\ell_{\min}^2\ge c^2L^2\)。而 \(|D|\le2bL^2\)、\(r\le AL\)，故
\(|Q\cap T|\ge(c^2-2b)L^2>aA^2L^2\)，小参数下违反 (10)。因此 \(r<\ell_{\min}\)。此时每轴裁入盒内的交长至少为 r，于是
\[
 r^2-|D|\le|Q\cap T|<a r^2,
 \qquad r\le\sqrt{2b/(1-a)}\,L=o(\epsilon L).
\tag{11}
\]
这一推断处理了跨盒边界的实际大 r 分支，没有先假设所有小源微见证已经处于 Z 或 D 中。

### 4.2 角内点：精确常数一

将四个角块扩大到边长 \(2\epsilon L\)，记其并为 V。对 \(\mathcal S\cap V\) 每个角分别按 §2 从全部该角源点选有限下降族，并将 §3 的符号固定为向盒内的两个方向。

由于中心距本角端点至多 \(2\epsilon L\)，且所有半径为 \(o(\epsilon L)\)，每个锚定 P 确实包含于 \(\mathcal R\)，并仍落在该角边长 \(3\epsilon L\) 的块内。四个这样的块两两分离，因为 \(3\epsilon L\ll\ell_{\min}\)。所以不同角的 P 也物理不交；不是只在各角各自套容量后重复使用 D。

同角 P 由 (5)–(6) 不交，且
\[
 |P\cap D|\ge(1-a)r^2
\]
来自 \(P\subset Q\cap\mathcal R\) 及 (10)。所选原 Q **覆盖整个角内源点集**，故
\[
 |\mathcal S\cap V|\le\sum_k|Q_k\cap T|
 <a\sum_k r_k^2
 \le\frac{a}{1-a}|D|
 =\frac{a}{1-a}(z+d).
\tag{12}
\]
价格常数正好为 a，不以一个固定覆盖重数代替，也无质量逼近损失。

### 4.3 角外点：四个向盒中心的锚定组

对 \(x\in\mathcal S\setminus V\)，由 (11)，其原 Q 和 Z 不交（Z 在角块 \(\epsilon L\) 内，x 在所有 \(2\epsilon L\) 角块外，到 Z 的距离至少 \(\epsilon L\)）。按 x 在两个盒中点坐标的哪一侧，把角外源点分成四组；平局按固定编号。每组的符号固定为两个方向都朝盒中心。

由 (11)，全部 r 最终小于最短边的一半。朝盒中心方向的可用轴长至少为对应盒边长的一半，所以每个该符号的锚定 P 都在盒内。对每组分别用 §2 从该组**全部源点**作有限下降选择，得到原 Q 覆盖全组；同组 P 按 §3 不交，并且 \(P\subset Q\) 保证避开 Z。于是对每组以 \(D\setminus Z\) 代 D 重复 (12)，跨四组仅取四倍，得
\[
 |\mathcal S\setminus V|\le\frac{4a}{1-a}|D\setminus Z|
                          =\frac{4a}{1-a}d.
\tag{13}
\]
合并 (12)、(13)，实际得到 \(|\mathcal S|\le a z/(1-a)+5a d/(1-a)\)，从而得 (9)。四组之间 P 可以相交，但这里只在 residual 上收固定四倍，并未让正常 z 支付四次。所有集合可测；有限源点选择只产生有限实际方体，没有不可测的逐点全局选择，也没有盒反射或网格边界噪声。

## 5. 在实际 N8 宏盒上分好坏

下文只取深处源。定义实际非正常残量
\[
 q_i=|\mathscr D_i\setminus Y_{i,N}|
               +|\mathcal R_i\setminus(x_i+W_i^0)|,
 \qquad \mathscr D_i=(x_i+W_i^0)\setminus F_i.
\tag{14}
\]
R049 已审 raw_N 和原入向扩张体积给
\[
 \sum_iq_i\le(a/v)R_N+Cv\mathcal F=O(av^7S)
 \qquad(n=2).
\tag{15}
\]
这里使用 \(R_N,\mathcal F=O(v^8S)\)，且 \(v\mathcal F=O(v^9S)=O(av^7S)\)。没有把全部正常层重新视作未用 raw。

令坏宏盒满足 \(q_i>av^5r_i^2\)，其原权
\[
 \sum_{i\ {\rm bad}}r_i^2=O(v^2S).
\tag{16}
\]
所有宏盒都有 \(|Y_{i,N}|\le Ca v^5r_i^2\)。对坏宏盒用原粗 dyadic，而非试图保留领先常数：
\[
 \sum_{\rm bad}\delta_{i,\rm in}
 \le Ca\sum_{\rm bad}(|Y_{i,N}\setminus T_i|+q_i)
 =O(a^2v^7S).
\tag{17}
\]
坏盒处的正常层只乘其 \(O(v^2S)\) 原权；残余 q 直接求和。费用晚于八阶。

### 好盒的真实端点和角位置

好盒上 \(q_i\le av^5r_i^2\)，所以
\[
 |\mathscr D_i|\le|Y_{i,N}|+q_i\le Ca v^5r_i^2.
\]
实际质量 \(\mu_i\asymp ar_i^2\) 及 \(\varepsilon_i\mu_i\le|\mathscr D_i|\) 给逐盒 \(\varepsilon_i\le Cv^5\)。强 CDF 端移界遂给
\(\Delta_i^{\rm end}\le C\varepsilon_i v r_i\le Cv^6r_i\)。

原盒每边 \(\ge cvr_i\)，故扩张体积控制每侧扩张长度：由 (14)
\[
 \text{每侧扩张}\le Cq_i/(vr_i)\le Cv^6r_i.
\]
因此正常层相对于有效源盒的两个近面深 \(\le K_Nv^3r_i\)，转到扩展盒角点后只新增 \(O(v^6r_i)\)。好盒上的 \(Y_{i,N}\) 确位于四个实际 \(\mathcal R_i\) 角的边长 \(C_Nv^3r_i\) 块内，不把有效端点、原端点或扩展端点静默视为相等。

在局部引理取
\[
 L_i=vr_i,\quad Z=Y_{i,N},\quad T=T_i\cap\mathcal R_i,
 \quad b=Cv^5,\quad \epsilon=C_Nv^2.
\]
实际 \(D=\mathcal R_i\setminus T_i\) 满足 \(|D\setminus Z|\le q_i\)；\(|Z|\le bL_i^2\)。好盒的 \(q_i\le av^5r_i^2=v^5L_i^2\)，满足 (7)。并且 \(\sqrt b/\epsilon=O(v^{1/2})\to0\)，全部固定源常数仅改变起效阈值。原小半径 \(\le Avr_i\) 和有限值原见证函数保持不变。

故每个好盒确有
\[
 \delta_{i,\rm in}
 \le\frac a{1-a}|Y_{i,N}\setminus T_i|+Ca q_i.
\tag{18}
\]
对好盒求和、加 (17)，用 (15)，即得 (1)。没有假设任何实际 Q 位于正常栏，也没有以额外局部输入截断损失替换其原严格平均。

## 6. 正常层总体积的领先系数

二维 \(b_e=1\)。正常无向边以实际 \(p,q,z\) 和 \(R_e=\max(R_i,R_j)\) 表示；源 \(i\to j\) 的两个排序标称面积之和为
\[
 \frac{P_{ij}Q_{ij}(P_{ij}+Q_{ij})}{4aR_i},
 \quad P_{ij}=(p+z)v^3R_e,\quad Q_{ij}=(q+z)v^3R_e.
\tag{19}
\]
实际目标缺量、共同侵蚀、正部及盒纤维截断只减小层，因此可用此上包；反向端为 \(p-z,q-z,R_j\)。正常区未触 clip，且两根比 \(1+O(v^3)\)。

在标称 \(p=q=2/[3(1+v)],z=0\)，两个端点 (19) 除 \(av^5R_e^2\) 的系数为
\[
 \frac8{27(1+v)^3}=\frac{2K_2}{(1+v)^3}.
\]
在固定正常紧窗内该多项式和根比是一致 Lipschitz 的。R049 正常强制性给
\(\sum_{\mathcal N}d_e^2R_e^2=O(v^2S)\)，有限度和 Cauchy 因而给
\(\sum_{\mathcal N}d_eR_e^2=O(vS)\)。所以
\[
 \sum_i|Y_{i,N}|
 \le2K_2av^5\sum_{e\in\mathcal N}R_e^2+O(av^6S).
\tag{20}
\]
此处同源同角正常目标唯一，其他角作并至多减小体积；所有正常边均保留，不按异常端丢弃物理层。

二维实际四槽版本 R045 (31) 给
\[
 W_N=\sum_{e\in\mathcal N}R_e^2
 \le\mathcal C+O(v^3S),\quad
 \mathcal C=\sum_iR_i^2\le S.
\tag{21}
\]
这里调用的是已审二维四角／实际根桥，不跨维引用 R043 的高维 helper 定理，也不把软权容量上界反向当作实际全容量界。代入 (20) 得 (3)，再结合 (2) 与 \(a=v^2\) 得 (4)。

边界小源、环境及深处盒外尾由 R054 §2 的独立内部正预留细化，分别为 \(O(v^{10}S)\)、\(O(v^{10}S)\)、\(O(v^9S)\)。因此本稿的 inside 常数不会再带一个未计的同阶边界项。

## 7. 真实基线中删除质量的系数也为一加小量

式 (4) 尚不能直接放到旧的 \(C(v/a)\delta\) 后面，否则其领先常数会被 C 放大。这里另给一个精确同层基线接口。

固定任何实际源，令 \(g_i=aR_i^2\)、\(m_i^0=|T_i|\)、\(\mu_i=|F_i|\)、\(\delta_i=m_i^0-\mu_i\)。实际包含和质量下界给
\[
 (1-Av)g_i\le\mu_i\le m_i^0\le g_i.
\]
记 \(f_i(m)=(R_i+\sqrt m)^2\)。在整个真实区间 \([\mu_i,g_i]\)，
\[
 \frac{f_i'(m)}{v/a}=v+\frac{vR_i}{\sqrt m}
 \le v+\frac1{\sqrt{1-Av}}
 \le1+(1+2A)v
\tag{22}
\]
（最后取 \(Av\le1/2\)）。所以
\[
 f_i(\mu_i)\ge A_vg_i/a
 -[1+O_A(v)]\frac va(g_i-m_i^0)
 -[1+O_A(v)]\frac va\delta_i.
\tag{23}
\]
求和后 \(\sum(g_i-m_i^0)=\mathcal D_h+\mathfrak d_h^0\)。原真实并集的环境／重复量再扣 \(\mathcal D_h\)，其净正系数仍为
\((A_v-a-[1+O(v)]v)/a\ge1/(2a)\)。\(\mathfrak d_h^0\) 可使用统一固定 b 作原八层精确望远镜，但 \(\delta_i\) 单独保留 (23) 的 \(1+O(v)\)，不被该 b 放大。

更明确地，二维实际乘积给 \(|X_i|\ge f_i(|W_i|)=f_i(\mu_i)+J_i\)，所以 Bonferroni 原样给 \(|\bigcup_iX_i|\ge\sum_i f_i(\mu_i)-H_h\)，其中 \(H_h=\sum I_e-\sum J_i\)。将 (23) 放到这同一个式子后，才扣真实重复量和环境质量。这里只是把 \(f_i(g_i)-f_i(\mu_i)\) 的一次真实积分拆成 \(g_i-m_i^0\) 与 \(\delta_i\) 两份，没有改变源、目标或原见证；J 仍基于同一有效 \(F_i\)，不再添加一份 J、也不在 H 中重复扣删除质量。

因此本稿 inside 费用进入第八层真实体积账本时恰至多为
\[
 (1+O(v))\frac va\sum\delta_{i,\rm in}
 \le (2K_2+O(v))v^8S.
\]
其余缺源量继续望远镜，盒外／边界删除和环境仍是晚阶项。这提供所需的 sharp 系数接口；正常 packet 本身八阶源成本尚须独立供应，不能据此完成总最优系数结论。

## 8. 计算、量词及待审范围

此证明的核心是从整个有限值源点集作全点下降选择、实际 inward P 的不交与盒包含、好盒强制所有原半径远小于角窗、远角按四个向盒中心的符号组只向 residual 收费。它不是只在模型盒内部取平均，也不依赖“微见证全在 D 中”或不可检查的近似覆盖。

所有常数先固定：维数二、原 N8 源常数与强 CDF、正常窗和 K_N，再取 v 小。坏盒阈值 \(av^5r_i^2\) 是显示的参数函数；其费用由 (15) 逐项支付，不偷带 LT 速率。最小原半径只用于每个固定输入的有限终止，价格与它、原有限列表长度、输入尺度跨度均无关。所有点集可测，中点平局按固定编号；本局部 sharp 引理没有网格边界问题。坏宏盒仅调用已审的粗 dyadic 界，其系数只乘晚阶费用。

已向 Terra 请求并实际读取本机最小新核 `r056-computation/sharp_local.wl/log`。Mathematica 15 ARM：同向 1D 区间正交叠与强分离的违反集 False；正常双端系数残差 0，值为 \(8/[27(1+v)^3]\)；fresh 导数残差 0，m 下界及 \((1-x)^{-1/2}\le1+2x\)、\(0\le x\le1/2\) 的违反集均 False；最终坏宏盒费用的显示幂次残差为 0。首次脚本曾混淆坏盒权与 normal_bad 成本，非零幂次输出留在 first-invalid 日志，不作为通过证据；改正后的成功日志才用于这里。没有重跑旧 R041/R046 等批次。这些输出不替代全点覆盖、实际锚盒包含／不交、跨盒见证或同层残量的解析证明。

新结论为 (1)–(4) 的实际 inside 小源精确领先上界及 (22)–(23) 的真实基线系数接口，待独立审查。它不识别正常 packet 自身八阶系数，不主张 \(J_2+8/27\) 已是全局最优，也不回改 R054 的停止范围。
