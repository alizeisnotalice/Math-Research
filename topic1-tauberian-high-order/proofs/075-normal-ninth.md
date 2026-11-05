# R075：真实正常场的第三次反馈与九阶有效源泛函

2026-11-05。固定每个整数 \(n\ge2\)。只写独占 scratch，不修改旧源或主仓库。本稿不识别 \(C_n\) 的九阶系数。

## 1. 实际范围与新结论

沿任意真实输入及严格外部 halo 紧集序列
\[
v\downarrow0,\quad a=v^n,\quad 1+M/V\ge P_{8,n}(t)-Bav^9,\quad B<\infty. \tag{1}
\]
使用同一次实际 N9 构造、原有限正值 Borel 半径、源／目标、CDF、正常集合 N 和全部实际几何。记
\[
D_i=B_i^0\setminus F_i,\quad D_i'=B_i'\setminus F_i,\quad
\lambda_v=s_vv/a,\quad s_v=(1+v)^{n-1},\quad S=\sum r_i^n.
\]
原始／有效盒分别为 \(B_i^0,B_i'\)。每个模型用有效源 chart、原始目标帽、实际 root、真实上下 clip；不换成独立二维 pair。

实际读取基础为 R063 foundation 702a70cdd735883d5adbe3d22da5acbe950e665c3023457a76416dcc5dd4f572、sixth 3a5d38fe5a00392af687c09ad59fcf072c552cbe418504af12fc38dc93e3cb2c、eighth 16401f6d8b9688d94ae115b2d5df1974dc36accbc404e102675e479de1e11a33，新 N9 的 R064 80f92dc770da78284a56e56a331b82112ea5912e892a456e968078e412a25699，以及 R071 fixed_tube_recentered_settlement.md SHA a014ae167a2d83a30fb2b455d827a2b13253fce5033f9223ca9c674956e6d5f2。

这里使用 R071 的实际正预算 (20)，不只用其全局上包：存在已经固定的同族 packet \(A_*=N\cup B_*\)，使
\[
\mathsf M_*+\mathcal R_{A_*}^-+v^6D_N+\mathcal L_N=O(v^9S). \tag{2}
\]
\(B_*\) 的 T、m、侵蚀及 CDF 均固定。任意 Y 强 CDF、全源原见证、排他和 whole-chart 必要性也按同一族供应。不把状态标签作为证明。

**主命题。** 对正常 packet，或任意先固定的有限上窗 packet \(\mathscr P\supset N\)，两端深处且 \(0<L_1,L_2\le Tv^3R_e\)，定义下节完整算子 \(B_{\mathscr P}\) 和
\[
Y_{0,\mathscr P}=B_{\mathscr P}(0),\qquad
Y_{j+1,\mathscr P}=B_{\mathscr P}(Y_{j,\mathscr P}).
\]
则
\[
\|D-Y_{2,\mathscr P}\|_1=O_T(av^8S),\quad
\|B_{\mathscr P}(D)-Y_{3,\mathscr P}\|_1=O_T(av^9S).       \tag{3}
\]
对任意 \(0\le\ell_i\le C_\ell v/a\)，
\[
\left|\sum_i\ell_i|B_{\mathscr P}(D)_i|
-\sum_i\ell_i|Y_{3,\mathscr P,i}|\right|
\le C_{T,\ell}v^{10}S=o(v^9S).                            \tag{4}
\]
在任何固定窗 Tonelli 退款核中，实际 D 可换成 \(Y_{2,\mathscr P}\)，误差同为 \(O_T(v^{10}S)\)。

若额外坏 packet \(B=\mathscr P\setminus N\) 位于固定 normal 窗外且满足固定 \(m\le p_e,q_e\le T\)、m>0，则
\[
\sum_i\ell_i(|Y_{3,\mathscr P,i}|-|Y_{3,N,i}|)
=\sum_{i,\mathrm{charts}}\ell_i\int h_{B,i}^0(y)\,dy
+O_{T,m,\ell}(v^{10}S).                                  \tag{5}
\]
此为 normal 与完整 bounded 非正常 packet 的加法，不把 normal 不同轴对拆成独立 pair。没有 m 下界时仍有 (3)–(4)，但不声称 (5) 的坏名义成本不被九阶退款抵消。

## 2. 算子保留全部实际几何

在有效盒 \(B_i'\) 先选最近面符号 \(\sigma\)，再选全局最浅轴 q，\(y_q=0\)。令
\[
x_h=f_{ih}^{\sigma_h}-\sigma_h(s+y_h),\quad
\Omega_{i,q,\sigma}=\prod_{h\ne q}[0,w'_{ih}/2],\quad
H_i(y)=\min_h(w'_{ih}/2-y_h).
\]
Jacobian 一；这些图覆盖整个有效盒至零测分界。最浅轴不在某 target pair 的图也保留。

对朝向匹配的 \(i\to j\)，
\[
c_{ij}=\operatorname{clip}((R_i-R_j)/2;-\max L_h,\max L_h),\quad
P_{ij,h}=(L_h+c_{ij})_+.
\]
原始目标帽 \(\Pi_{ij}(y)\) 的两临界深度为
\[
\min(w_{jl}^0,(P_{ij,l}-y_l)_+),\quad
\min(w_{jk}^0,(P_{ij,k}-y_k)_+),
\]
其余取原始整宽。定义
\[
h^0_{\mathscr P,i}(y)=\frac{\sum_j|\Pi_{ij}(y)|}{naR_i^{n-1}},
\quad K_{\mathscr P,i}(E;y)=\frac{\sum_j|E_j\cap\Pi_{ij}(y)|}{naR_i^{n-1}},
\]
\[
B_{\mathscr P}(E)_i=\bigcup_{q,\sigma}
\{x(s,y):0\le s<\min(H_i(y),[h^0_{\mathscr P,i}-K_{\mathscr P,i}(E;y)]_+)\}.
                                                               \tag{6}
\]
全部实际目标先在同一高度相加，随后只取一次正部和一次 source clip。测试 E 不改 root、盒、原见证、目标清单或 CDF。每个测试集合在原始盒内，输出在有效盒内；所有求和有限，可用 Tonelli。

打开积分，固定目标物理点：q 在其 pair 内时，可覆盖它的源底面积为 \(O_T((v^3R_i)(vR_i)^{n-2})\)，除容量导数得 \(C_Tv\)；q 在 pair 外时得 \(C_Tv^3\)。有限逆度、正部及 min 的 1-Lipschitz 性给
\[
\|B_{\mathscr P}(E)-B_{\mathscr P}(F)\|_1
\le\Lambda_T\|E-F\|_1,\quad \Lambda_T\le C_Tv<1/2,\quad
|Y_{0,\mathscr P}|_1\le C_Tav^5S.                         \tag{7}
\]
二维后一类为空。薄帽和零高度也满足同式。

原每点严格见证同时捕获完整自身 G、真正浅面邻轴帽和全部 packet 目标，先排 full cap，再比较 partial cap／容量斜率，给实际保证层
\[
B^-_{\mathscr P}(D)\subset B_{\mathscr P}(D)\cap D',\quad
|B_{\mathscr P}(D)\setminus B^-_{\mathscr P}(D)|_1\le C_Tav^9S. \tag{8}
\]
这里实际供应的是更强的侵蚀积分预算。令共同非负侵蚀为 \(e_{\mathscr P,i}(y)\)，保证层的高度是在同一个 source clip 内再扣该侵蚀。原必要性构造给
\[
\mathscr E_{\mathscr P}:=
\sum_{i,q,\sigma}\int_{\{h^0_{\mathscr P,i}>0\}}e_{\mathscr P,i}(y)\,dy
\le C_Tav^9S.                                             \tag{8a}
\]
它由正帽 base 测度及 \(\sum R_i^{n-1}e_i(T)=O_T(v^8S)\) 相乘得到，不是由联合层已经裁剪后的实际损失反推。对正常高度使用同一个侵蚀，记保证子层 \(B^-_{N,\mathscr P}(D)\)。正常正高度的 base 包含于 \(\{h^0_{\mathscr P}>0\}\)，每条 fiber 的高度损失至多该侵蚀，所以
\[
|B_N(D)\setminus B^-_{N,\mathscr P}(D)|_1
\le\mathscr E_{\mathscr P}\le C_Tav^9S.                    \tag{8b}
\]
即使 full 层先被顶端 clip 饱和，这个子层估计仍成立；不能仅以 full 实际 loss 的大小推出它。上述均为物理误差，乘源导数后为十阶；不删除整个三近面。只取固定 T，不用 fixed-window 未知阈值代入增长参数。

## 3. 新整源 raw，而非旧正常核心 raw

对 (2) 的同一个侵蚀，非负实际高度为 \(h_{A_*}=h_N+h_{B_*}\)。正部／clip 给常数一
\[
\mathcal R_{N,*}^-:=\frac va|D\setminus B^-_{N,*}(D)|_1
\le\mathcal R_{A_*}^-+s_v^{-1}\mathsf M_*=O(v^9S).         \tag{9}
\]
侵蚀只保证实际缺源，不改变未侵蚀正常算子 \(B_N\)。于是
\[
|D\setminus B_N(D)|_1\le(a/v)\mathcal R_{N,*}^-=O(av^8S).
\]
反方向由 (8) 为 \(O(av^9S)\)，所以
\[
\boxed{\rho_N:=\|D-B_N(D)\|_1=O(av^8S).}                 \tag{10}
\]
旧核心 \(Y^{old}\) 的 raw 仍可为 \(O(v^8S)\)；没有把它改名为 (9)。多出的一阶来自真实 raw_A 与 nominal M 转移。

正常层逐源名义体积与旧 raw 给
\[
|D|_1=O(av^5S).                                         \tag{11}
\]
这是正常 mandatory 主量，不能误说成 \(O(av^8S)\)。

对其他任意固定上窗 packet \(\mathscr P\supset N\)，取它的共同侵蚀。即使额外帽任意薄，也有 \(B^-_{N,\mathscr P}\subset B^-_{\mathscr P}\)，正常未侵蚀模型仍是同一 \(B_N(D)\)。由 (8b)、(10) 及 \(B^-_{N,\mathscr P}\subset B^-_{\mathscr P}\)，
\[
\mathcal R_{\mathscr P}^-\le(v/a)|D\setminus B_N(D)|_1+C_Tv^{10}S=O_T(v^9S).
\]
再用 (8)，得
\[
\rho_{\mathscr P}:=\|D-B_{\mathscr P}(D)\|_1=O_T(av^8S).  \tag{12}
\]
这把 (3) 推到所有后来固定的 bounded 上窗。名义坏成本的简单加法仍需第 6 节的 m 下界。

## 4. 有率递推与有限第三反馈

固定 packet，写 \(B=B_{\mathscr P},\Lambda=\Lambda_T,\rho=\rho_{\mathscr P},M_0=|D|_1\)。直接三角不等式给
\[
\|D-Y_0\|\le\rho+\Lambda M_0,\quad
\|D-Y_1\|\le\rho+\Lambda(\rho+\Lambda M_0),
\]
\[
\|D-Y_2\|\le\rho(1+\Lambda+\Lambda^2)+\Lambda^3M_0=O_T(av^8S). \tag{13}
\]
再作用 B，
\[
\|B(D)-Y_3\|\le\Lambda\|D-Y_2\|=O_T(av^9S).
\]
乘最大源体积权 \(C_\ell v/a\) 就是 (4)，不作逐胞估计。

保留 fiber clip 的非线性退款
\[
\Phi_\ell(E)=\sum_i\ell_i(|Y_{0,i}|-|B(E)_i|)
\]
的 L1 常数至多 \(C_{T,\ell}v^2/a\)，实际 Tonelli 退款核也有此上界。因此 D→Y2 的费用为 \((C_Tv^2/a)O_T(av^8S)=O_T(v^{10}S)\)。一般权 (4) 只是模型积分比较，不声称任意比 fresh 导数大的权也是 J 下界；真实 J 仍用 \(\lambda_v\)。

反单调性还给
\[
Y_1\subset Y_3\subset Y_2\subset Y_0,\quad
|Y_2\setminus Y_3|_1\le\Lambda^3|Y_0|_1=O_T(av^8S).       \tag{14}
\]
所以正常第三反馈净价
\[
\mathsf P_N^{[9]}:=I_N-\lambda_v|Y_{3,N}|
=\mathsf P_N+\lambda_v|Y_{2,N}\setminus Y_{3,N}|
\]
新增一个非负 \(O(v^9S)\) 项，不能自动丢入 o9。下一次 Y3 与 Y4 的成本差才从十阶起。

## 5. 一次 CDF 的实际九阶源泛函

对 (8) 的真实保证层使用同一次任意 Y 强 CDF，系数一主层与固定 raw／端移同时保留。完整导数 epsilon 费为
\[
C_T\lambda_v\sum_i\varepsilon_i|Y_{0,\mathscr P,i}|
\le C_Tv^6\sum_i\varepsilon_iR_i^n=O_T(v^{11}S).
\]
侵蚀成本十阶。由 (3)–(4)，
\[
\boxed{\sum_iJ_i\ge\lambda_v|Y_{3,\mathscr P}|_1
+c\,\mathcal R_{\mathscr P}^-+c\,\mathfrak E-C_Tv^{10}S.}  \tag{15}
\]
主项是明确的有限 box/cap 算子从空集开始应用四次，不是用原极大算子或未知 Tauberian 上确界重新定义未知系数。没有把全部 D 换成 Y3 后丢掉九阶 raw。

还有一个等价的几何 fixed-layer：模零测集合族在对称差 L1 距离下完备，B 是严格压缩，故给定有限几何有唯一 \(Z_{\mathscr P}=B(Z_{\mathscr P})\)。这不是极值输入唯一性。由压缩余项
\[
\|Z_{\mathscr P}-Y_{3,\mathscr P}\|_1
\le\frac{\Lambda^4}{1-\Lambda}|Y_{0,\mathscr P}|_1=O_T(av^9S),
\]
(15) 主成本可换成 \(\lambda_v|Z_{\mathscr P}|\)，费用仍十阶。有限第三反馈已足够，不需要无限迭代作为证明前提。

## 6. bounded 非正常 packet 的联合加权层

另要求 \(B=\mathscr P\setminus N\) 在固定 normal 窗外且 \(m\le p,q\le T\)，m>0。完整坏名义成本为
\[
\mathsf M_B=\lambda_v\sum_{i,\mathrm{charts}}\int h_{B,i}^0(y)\,dy.
\]
未侵蚀真实增量精确满足
\[
\Delta_{\mathscr P}:=\lambda_v|B_{\mathscr P}(D)\setminus B_N(D)|_1
=\mathsf M_B-\mathfrak F_B-\mathfrak L_B.                 \tag{16}
\]
非负共同高度和真实上截断保留；固定窗给 \(\mathfrak L_B\le C_Tv^3\mathsf M_B\)。名义三近面项在 M 内，不被这项误差删除。

退款目标缺集现在分解为
\[
D_j=(D_j\cap B_N(D)_j)\,\dot\cup\,(D_j\setminus B_N(D)_j).
\]
前者被完整正常模型逐源体积 \(Cav^5R_j^n\) 控制；后者的全族量由 (10) 为 \(O(av^8S)\)，而非旧核心 raw 的 \(O(av^7S)\)。两端完整坏名义成本至少 \(c_mv^6R_e^n\)，来自真实原整宽下界和根差正部凸性，不要求坏边根差正常或横权有下界。因此同一个 Tonelli 核给
\[
\mathfrak F_B\le C_{T,m}v\mathsf M_B+C_TvR_N^{un},\quad
R_N^{un}=(v/a)|D\setminus B_N(D)|_1=O(v^9S).              \tag{17}
\]
实际增量在 D 内的部分属于 D\B_N，D 外的部分由 (8) 控制，故
\[
\Delta_{\mathscr P}\le s_vR_N^{un}+C_Tv^{10}S.
\]
与 (16)–(17) 合并，固定 T,m 后小 v 吸收系数，得到
\[
\boxed{\mathsf M_B=O_{T,m}(v^9S),\quad
\mathfrak F_B=O_{T,m}(v^{10}S),\quad
\mathfrak L_B=O_{T,m}(v^{12}S).}                          \tag{18}
\]
这对任意后来固定的 bounded-away-from-zero packet 成立，不只 R071 的 \(B_*\)，未反用 (5) 证明 M 先验。

一般权 \(0\le\ell_i\le C_\ell v/a\) 与 \(\lambda_v\) 的比一致有界；逐源增量的名义、退款和 clip 分解同样成立，故退款／clip 的加权和仍十阶。不要求权有正下界。对 normal 和 joint 算子分别用 (4)，便得 (5)。两套 Y3 不必逐点嵌套；加法来自真实 B(D) 的嵌套和 o9 比较，不能假定迭代本身嵌套。

真实源权取 \(\lambda_v\)，(15) 与 (5) 给
\[
\boxed{\sum J_i\ge\lambda_v|Y_{3,N}|_1+\mathsf M_B
+c\mathcal R_{\mathscr P}^-+c\mathfrak E-C_{T,m}v^{10}S.} \tag{19}
\]
只用一次完整 J。对应部分净价保留为
\[
I_{N\cup B}-\sum J_i\le \mathsf P_N^{[9]}+I_B-\mathsf M_B
-c\mathcal R_{\mathscr P}^- -c\mathfrak E+C_{T,m}v^{10}S.  \tag{20}
\]
非正常 M、raw、端移仍可能在九阶贡献，未从主阶删除。

## 7. RMS、罕见源及跨 pair 的范围

R071 的 moving-normal 亏损为
\[
\mathcal L_N=v^6\sum[(1-\bar b_e)+(x_e-r_v)^2+
(y_e-r_v)^2+\zeta_e^2]R_e^n=O(v^9S).
\]
没有把它改成逐边 \(O(v^{3/2})\)，也没有改善旧固定中心 GN。由旧正常模型定量上包和 (14)，仍可保留
\[
\mathsf P_N^{[9]}\le q_NW_N-c\mathcal L_N+Cv^9S.
\]
这个粗上包不评定九阶系数；九阶有效泛函是未逐边展开的 \(\mathsf P_N^{[9]}\)。

抽象 normal 参数场若在 \(v^3\) 权上有固定偏差，二次罚和三次 Taylor 余项都可能在 \(v^9S\)。这仅说明 RMS 不足以删除三次项，不声称该抽象场已被几何实现。本稿的全族 L1/operator 推导不使用这种删除，允许全部实际稀有源。

q 不在 pair 的图、free-axis source clip、两个不同根帽、normal 全部目标的共同正部及 Y0/Y1/Y2/Y3 中的跨 pair 反馈，都在 (6) 内。没有写成独立二维 pair 成本之和。这与 R072 signed-shear 的物理比较相容；其 cross-pair 和 determinant 可能是独立机制，本稿不由 pair 系数推断九阶极值几何。

## 8. 实际计算与当前交付

已实际运行本机 third_feedback.wl/log，Wolfram 15 ARM，PID86683。取 \(\Lambda=Cv,\rho=Rav^8,M_0=Mav^5\)，程序打印 (13) 的有限递推：源成本和退款上界均从 v10 起，常数不同，不断言两者相等。Y2−Y3 成本从 v9 起，Y3−Y4 从 v10 起。软件只核代数阶数；真实 raw 改善、原见证、盒帽及同一次 J 由正文证明承担。

**交付是 general-fixed-n 的局部九阶源泛函与第三反馈替换引理：任意 bounded 上窗的 (3)–(4)，以及有固定 m 下界的联合加权 (5)/(19)。**

未计算第三模型的九阶极大值，未闭 long/nearzero、inside、原环境／端移、非正常 nominal 或全局前缀的九阶系数，不识别真实 C_n 九阶系数。固定新 N≥9 只有实际供应 (2) 与完整同层 G/C/B/L/Q 后才可迁移；进行中的 R076 分工不是前提。冻结后交非作者独立复核。

