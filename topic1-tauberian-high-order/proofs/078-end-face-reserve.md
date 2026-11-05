# R078：同一次九层结算中的端移与免费面正量

2026-11-05。固定每个整数 \(n\ge2\)。本稿只写独占 scratch，不修改 R071、R075 或正式正文。
这是供非作者审核的完整新候选；以下证明不以旧审核标签代替实际几何。

## 1. 范围与结论

令 \(a=v^n,\ A_v=(1+v)^n,\ \mathcal T=V+M,\ S=\sum_i r_i^n\)。
固定 \(B<\infty\)，沿任意真实输入和严格外部 halo 紧集序列，假设
\[
 v\downarrow0,\qquad 1+M/V\ge P_{8,n}(t)-Bav^9.                 \tag{1}
\]
\(P_{8,n}\) 为 R071 的已确定八阶多项式，包括 \(n=2\) 的真实 inside 增量。
先固定 \(0<\tau<1/100\)（或有限多个这样的值），再一次构造实际 N9。
原有限正值 Borel 半径、原贪心族、标记集 I、九层目标／源、三端盒、CDF、
normal 集合与 clean 槽，全部取这一族。没有重新选源或再使用一次输入。

未标层号的对象取第九层。记 \(B_i^0,B_i',B_i^m\) 为原始、有效、中点源盒，
\[
 D_i=B_i^0\setminus F_i,\quad D_i'=B_i'\setminus F_i,\quad
 \lambda_v=s_vv/a,\quad s_v=(1+v)^{n-1},
\]
\[
 \mathfrak E=\sum_i R_i^{n-1}\Delta_i,\qquad
 \mathcal F=\sum_{i\ {\rm 深处},q,\pm}\Gamma_{ip}r_i^{n-1}.
\]
\(\Delta_i\) 为全部 CDF 端移之和；普通面对 \(p\) 是原真正面邻居。
\(\mathcal F_\tau\) 为 §4 所定义的全部固定裕量稳健面对的有向和。
两种面量均非负，不要求它们的条互不相交。

R071 的固定 packet \(A=N\cup B_*\)、侵蚀及 normal 整管 K 在全文不变。
令 \(\mathsf M=\lambda_v\sum\int h_{B_*}^0\)、
\(\mathcal R_A^-=(v/a)\sum|D_i\setminus B^-_{A,i}(D)|\)。
槽欠量 \(D_N\) 和移动 normal 亏损 \(\mathcal L_N\) 均按 R071 的实际定义。

**命题。** 可按 R071 的非循环固定参数次序，并保留同一族，使
\[
 \boxed{\mathsf M+\mathcal R_A^-+v^6D_N+\mathcal L_N+
        \mathfrak E+\mathcal F+\mathcal F_\tau
       =O_{n,B,\tau}(v^9S).}                              \tag{2}
\]
结论是原先八阶数值预算的端移／面量升级，不识别九阶系数。
允许有限多个固定 \(\tau\)，常数可依这些值；不允许 \(\tau=\tau(v)\)。
二维所有高支撑求和为空，证明仍适用；二维 inside 仍只有已供应的八阶主项加 O9。

## 2. 本族的前置数值与一次 CDF 的端移余量

本次 N9 的粗 bootstrap 已经给出
\[
 R_*+\mathfrak E+\mathcal F+\mathcal F_{1/128}
       +P_I+P_\partial=O(v^8S),\quad D_N=O(v^2S),          \tag{3}
\]
以及本族宽度 \(\sum\xi_i^2R_i^n=O(v^7S)\)、
\(\sum\varepsilon_iR_i^n=O(v^5S)\)。
这些由同族六阶、七／八阶粗前缀先得，不由 (2) 反供。
对一般固定 \(\tau\)，六阶面预算即可用于几何身份；本稿不需要预设其八阶界。
此外，R071 原证明在不改善端移／面量时，已经独立给出同一 packet 的
\[
 \mathsf M+\mathcal R_A^-+v^6D_N+\mathcal L_N=O(v^9S).
                                                               \tag{3a}
\]
这是 R075 第三反馈的数值前提。本稿允许引用这个已证数值；
并不将旧并集或旧 J 的正量再加进新收入。
下面升级所用的端移、面量仍只知道 (3) 的旧界，没有反用待证 (2)。

下面重开 R071 使用的同一份 CDF，不把另一条源价格相加。
设 \(\mu_i=|F_i|\)，\(Y_i\subset D_i'\) 任意可测。
强 CDF 构造的实际估计为
\[
 |F_i\setminus B_i'|\le\lambda_0(vr_i)^{n-1}\Delta_i,\quad
 |B_i^m|-|B_i'|\ge c_0(vr_i)^{n-1}\Delta_i,
\]
\[
 |B_i^0\setminus B_i'|\le C_0(vr_i)^{n-1}\Delta_i.
\]
有限深度／宽度常数先固定，随后取 \(\lambda_0<c_0/2\)。
精确恒等式
\[
 |B_i^m|-\mu_i-|Y_i|
 =|(B_i'\setminus F_i)\setminus Y_i|
   +|B_i^m|-|B_i'|-|F_i\setminus B_i'|
\]
同时留出固定 \(\kappa>0\)：
\[
 |B_i^m|-\mu_i\ge |Y_i|
   +\kappa\{|D_i\setminus Y_i|+(vR_i)^{n-1}\Delta_i\}.   \tag{4}
\]
原盒环与端移同阶，故 \(\kappa\) 只需缩小固定倍数，不依后来 T、m。

对 \(f_i(z)=(R_i+z^{1/n})^n\)，\(g_i=aR_i^n\)，有
\[
 f_i'(g_i)=\lambda_v,\qquad f_i'(z)\ge c\,v/a
\]
于全部实际积分区间；锐主段另用
\(f_i'(z)\ge(1-C\varepsilon_i)\lambda_v\)。
将 (4) 中长度 \(|Y_i|\) 的主段与剩余段分开积分，得到
\[
 \sum_iJ_i\ge\lambda_v\sum_i|Y_i|
   +c_R\frac va\sum_i|D_i\setminus Y_i|
   +c_E\mathfrak E
   -C\lambda_v\sum_i\varepsilon_i|Y_i|,\qquad c_R,c_E>0. \tag{5}
\]
端移系数的归一化恰为
\((v/a)(vR_i)^{n-1}\Delta_i=R_i^{n-1}\Delta_i\)，因为 \(a=v^n\)。
即使单个 \(\varepsilon_i\) 大，主段锐下界变负仍合法；只用全族加权界结算。
所有 \(J_i\) 是这一次源盒扩张的真实盈余。

取同一个实际保证层 \(Y=B_A^-(D)\)。
R075 v2 的 positive nominal base 积分预算给
\[
 \sum_{i,q,\sigma}\int_{\{h^0_{A,i}>0\}}e_{A,i}(y)\,dy
       \le C_Tav^9S.                                    \tag{6}
\]
这是对积分侵蚀的预算，不是 full 已 clip 后的损失。
normal 正 base 包含于 joint 正 base，所以同一侵蚀下 normal 子层亦损失至多 (6)。
乘真实源导数，侵蚀费为 \(C_Tv^{10}S\)；导数 epsilon 费为 \(C_Tv^{11}S\)。
保留 \(c_E\mathfrak E\) 不需删除整个三近面，也没有额外的端移主费用。
R075 的第三反馈可在 (5) 同一次使用，给
\[
 \sum J_i\ge \lambda_v|Y_{3,N}|+\mathsf M
   +c_R\mathcal R_A^-+c_E\mathfrak E-C_{T,m}v^{10}S.     \tag{7}
\]
固定 m 下的坏 nominal 加法已包含共同正部、两根帽、off-min 图和真实 clip。
本稿不会把 (7) 与 R071 的另一份完整 J 再加一次。

## 3. 联合价格中保端移，吸收仍用同一份额

先于 T、m，固定窗外 coercivity 余量、CDF 的 \(c_R,c_E\)、
槽严格折扣系数和 normal 强凹罚。各系数可先缩小固定倍数，以留足有限份额。
对实际 normal 第三模型，令 \(W_N=\sum_{e\in N}b_eR_e^2\) 为实际 normal 权，
其中 \(b_e\) 是该边实际横投影量；它不是完整槽容量 \(\mathcal C\)。
\[
 \mathsf P_N^{[9]}=I_N-\lambda_v|Y_{3,N}|
       \le q_NW_N-c\mathcal L_N+Cv^9S.                  \tag{8}
\]
这里只在此后作一次 soft 槽转换以保留 \(D_N\)，未先用完整 \(\mathcal C\) 替换 \(W_N\)。
仅用已有实际积分的 O9 上包，不展开所有稀有边的三阶 Taylor。
\(Y_2\setminus Y_3\) 的源成本可为九阶，已包含在 \(Cv^9S\)，未当作 o9 删除。

使用 (7) 重开 bounded 非正常强制性和一次 soft 槽转换。
其余尚待吸收的费用逐项为
\[
 \theta_{T,m}(v)\mathsf M,\quad
 \theta_{T,m}(v)\le C(Tm^{-3}v+T^2v^3+v),
\]
\[
 I_{\rm long}\le(C_L/T)(\mathcal R_A^-+s_v^{-1}\mathsf M)
            +Cv(R_*+\mathcal F+\mathfrak E),
\]
\[
 I_{\rm small}\le CTm\,v^6D_N+C_{T,m}v^9S,\quad
 I_{\ge3}\le\eta(\mathcal R_A^-+s_v^{-1}\mathsf M)
                              +C_{\eta,K}v^9S.          \tag{9}
\]
长条的 K 是仅含 normal 整模型的固定管，独立 T、m；
条在整个 K 管外，原始目标帽的原／有效环另外收费 \(Cv\mathfrak E\)。
不得将这个环直接记到 Z。高支撑只取同一 Z 的一个固定小份额。
其证明仍用 (3) 的旧数值，未先用 (2)。
normal 实际端锚费也至多 \(Cv\mathfrak E\)；ordinary 面费为 \(Cv\mathcal F\)。
width 费 \(O(v^{19/2}S)\)、必要性侵蚀费 O10、soft 费 \(C_mv^9S\)
以及其余 normal O9 合入固定常数。

现在先取固定大 T，使长条首项只耗 raw、nominal 预留的小份；
再固定小 m，使近零首项只耗槽预留的小份；
再固定小 \(\eta\)，最后取 v 小吸收 \(\theta_{T,m}(v)\mathsf M\)。
保留的端移系数与 T、m 无关，所有 \(C_{T,m}v\mathfrak E\) 最后再吸入其一半。
所得同一次完整净交价 \(H=\sum_{i<j}|X_i\cap X_j|-\sum_iJ_i\) 满足
\[
 \boxed{H+c_M\mathsf M+c_R'\mathcal R_A^-+c_Dv^6D_N+
              c_L\mathcal L_N+c_E'\mathfrak E
       \le q_N\mathcal C+C_1v\mathcal F+Cv^9S.}           \tag{10}
\]
全部显示系数固定正。没有为端移或高支撑重新减一次 J。
实际完整槽容量和本族宽度给
\[
 q_N\mathcal C\le Q_{8,N}(v)\mathcal T+O(v^{19/2}S+v^{12}S). \tag{11}
\]
这里 \(Q_{8,N}\) 是 R071 §9 的 normal 八阶中心，不含 n2 inside。
系数 \(q_N>0\) 整体乘实际容量；不能按其负八阶项逐项反向替换。

## 4. 真实免费条：全外部盒之外的一次并集

普通面条用本族第九层 \(X_i=B_i^m+[-R_i/2,R_i/2]^n\) 的两个中点面对端。
在原真正面邻的固定切向块上作开条 \(Z_{ip}^0\)。
有向面对容量给两向 \(g\ge0\) 及其半和 \(\Gamma_{ip}\ge0\)；
条体积为固定倍数的 \(\Gamma_{ip}r_i^{n-1}\)。

稳健面对的定义保持 R017 的实际范围：i 深处，p 外部，
\(7r_i/8\le r_p\le9r_i/8\)，相应原半盒法向缝
\(0\le b\le r_i/8\)，各切向原投影交长至少 \(\tau r_i\)。
两符号均计。其 CDF 中点面对端的 \(\Gamma_{ip}\) 加权有向和定义
\(\mathcal F_\tau\)。
先在切向投影交内取边长 \(\tau r_i/2\)、距双方切向边至少 \(\tau r_i/4\) 的面块，
再在两个中点面对端间作同块上的开条 \(Z_{ip}^\tau\)。

为说明没有漏第三盒，归一化 \(x_i=0,r_i=1\)。
原半盒间缝上的条在 \(A_i^0\) 内。任何穿条的原半盒 k 由中心严格分离
与尺度几何先有 \(1/2<r_k<3/2\)；其切向不能完成分离，故须沿法向位于 i 右侧。
它还与 p 的切向内部相交，法向不能以边长大于 \(1/2\) 塞入宽至多 \(1/8\) 的缝，
也不能从 p 右侧穿回。因此原条无半盒穿过。
原空隙量便给 \(c_\tau b\le e_i^0\le Dv\)，中心严格分离还给
\(|r_i-r_p|<2b\)。所以 \(b+|r_i-r_p|=O_\tau(v)\)。
固定 \(\tau\) 后，源、根和 CDF 的 O(v) 位移都小于切向固定裕量；
同一实际见证完整捕获真正相邻源的横截面，强 CDF 两向容量仍给中点面对顺序。

若 \(z\in Z_{ip}^\tau\cap X_k\)，
\(z_q=1/2+O_\tau(v)\) 和
\(X_k\subset Q(x_k,(1/2+A_9v)r_k)\) 先迫 \(r_k=1+O_\tau(v)\)。
切向距不够中心分离，所以 k 必在法向右侧。
其原半盒切向与 p 有正交长；于是 k 若在 p 左，近单位边长不能塞入 O(v) 缝，
若在 p 右，\(X_k\) 只能从约 \(3/2\) 的端向左移动 O(v)，不能达到 z。
这同时排除了全部外部 k，不只 normal 邻居。
普通条取固定面块时同证。因此
\[
 (Z^0\cup Z^\tau)\cap\bigcup_{i\notin I}X_i=\varnothing. \tag{12}
\]
所有条附着于深处源的 \(Q_i\)，而深处 \(Q_i\) 避开每个标记 \(Q_j\)，
故它们还避开内部半盒及标记空隙。

每源候选数有固定上界，且每条在 \(A_i^0\) 内；
半盒 packing 给固定有限重数。因此分别有
\[
 |Z^0|\ge\alpha_0\mathcal F,\qquad
 |Z^\tau|\ge\alpha_\tau\mathcal F_\tau,\quad
 \alpha_0,\alpha_\tau>0.
\]
两类可以重叠。只令 \(Z=Z^0\cup Z^\tau\)，则
\[
 |Z|\ge\max(|Z^0|,|Z^\tau|)
       \ge\tfrac12\alpha_0\mathcal F+
                    \tfrac12\alpha_\tau\mathcal F_\tau. \tag{13}
\]
有限多个 \(\tau\) 改用相应有限正权平均。未给任何一类第二份收入。

## 5. 增广并集的实际环境只花一次输入

对每个深处源，令
\[
 \widehat X_i=X_i\cup\bigcup_{p}(Z_{ip}^0\cup Z_{ip}^\tau),
\]
边界源仍令 \(\widehat X_i=X_i\)。记 \(U^+=\bigcup\widehat X_i\)。
所有附着点距 \(x_i\) 至多 \((1/2+C_\tau v)r_i\)，
每个 \(x\in F_i\) 距 \(x_i\) 为 \(O(A_9vr_i)\)，
而 \(R_i=(1+O(A_9v))r_i\)。
因此固定 \(\tau\) 后统一小 v，整个 \(\widehat X_i\) 对每个同一原源点均位于
\(Q(x,R_i)\)。没有为不同条换用不同源见证或另选输入。

令 \(u_i=|(\widehat X_i\setminus E)\setminus G_i|\)。
完整 \(G_i\) 与这部分排他，原严格见证给
\[
 u_i\le Ca r_i^{n-1}[r(x)-R_i]\quad(x\in F_i).
\]
选原普通面邻的源点趋原面对端，有实际
\(r(x)-R_i\le2(g_{i\to p}+d)\)，\(d\downarrow0\)；
固定的整个 \(u_i\) 因而满足
\[
 \sum_{\rm 深处}u_i\le C_\tau a\mathcal F,\qquad
 \sum_{\partial I}u_i\le CaP_\partial.                  \tag{14}
\]
新增多少稳健条都不按条数再付质量。可仍用一个原普通面方向控制整个 \(u_i\)，
所以环境费不需要 \(\mathcal F_\tau\) 的八阶先验。
面对端点未取到时用有效源点逼近，半空间质量连续性只将严格界变成非严格界。

由 (12)–(13) 和 Bonferroni，
\[
 |U^+|\ge\sum|X_i|-\sum_{i<j}|X_i\cap X_j|
       +\tfrac12\alpha_0\mathcal F+
                    \tfrac12\alpha_\tau\mathcal F_\tau. \tag{15}
\]
而同一个并集有
\[
 |U^+\setminus E|\le\sum_i g_i+\sum_i u_i
                         =M'+\mathcal D_9+\sum_i u_i.   \tag{16}
\]
给它加入原标记半盒和原标记空隙时，面条完全避开标记区域，
所以原内部回收恒等式及两次纯维数突出价格原封保留。
上界仍是 \(V+M'+\mathcal D_9+\sum u_i+aS_I\)，V 只出现一次。
normal／joint 保证缺源层是 (10) 内唯一 J 的证明工具；
本节条是避开全部 \(X_i\) 的真实空间增量，二者没有相加第二份 J。
粗 bootstrap 曾恢复面的数值结论 (3) 仅作先验，本节没有把旧并集的体积再加入。

## 6. 一次支付小源、环境及面锚费用后的终端正量

在 (15) 使用 \(\sum|X_i|=\) 基线 \(+\sum J_i+\) 非负形状，
并以 (10) 结算唯一 \(H=\sum I-\sum J\)。
普通面收入的费用统一列为
\[
 C_1v\mathcal F\quad\hbox{(联合价与长角／面锚)},\quad
 C_2v\mathcal F\quad\hbox{(盒外小源)},\quad
 C_3v^n\mathcal F\quad\hbox{(盒内扩张)},\quad
 C_4a\mathcal F\quad\hbox{(整个环境)}.                  \tag{17}
\]
它们共同由 (13) 的普通半份支付，不能各耗整份收入。
因 \(a=v^n\le v\)，固定全部参数后小 v 使总系数不超过 \(\alpha_0/4\)，
仍余 \((\alpha_0/4)\mathcal F+(\alpha_\tau/2)\mathcal F_\tau\)。
所有 \(C_{T,m,\tau}v\mathfrak E\) 已在 §3 由 \(c_E/2\) 支付，仍余固定正端移。

实际小源仍逐点使用原指定半径：盒内 dyadic 保留 a 因子，
盒外入向尾使用全部 \(T_i\) 而不只 \(F_i\)。
本次 n2 的 inside 完整费用是
\[
 (8/27)v^8\mathcal T+O(v^9S),
\]
由本族全点原半径和实际四角窗给出；n≥3 的绝对盒内项从 \(v^{n+6}S\) 起，
故至少九阶。保留这些 O9 数值费用不使用待证 (2)。
边界／环境价格仍由原纯维数图预留吸收；
其附加系数 \(C_{T,m,\tau}vP_\partial\) 最后小 v 后使用原独立一份预留。
在整个实际 \(\mu/g\in[1-A_9v,1]\) 区间，删除质量的导数为
\((1+O(v))v/a\)；未删 \(d_9^0\) 按原递推单独处理，不给 inside 主费错误乘固定 b。

令 \(Q_8=Q_{8,N}+(8/27)\mathbf1_{n=2}v^8\)，
记固定正系数之和
\[
 \mathcal L_{78}=
 c_M'\mathsf M+c_R''\mathcal R_A^-+c_D'v^6D_N+c_L'\mathcal L_N+
 c_E''\mathfrak E+c_F\mathcal F+c_\tau\mathcal F_\tau\ge0.
\]
以上对同一个真实增广并集和原内部恒等式给本次第九层公共账本
\[
 \begin{split}
 \mathcal T\ge {}&A_vM/a+c_{\rm rec}P_I-C_bP_\partial-2nvS_I\\
 &+c_{\rm mass}\mathcal D_9/a-b(v/a)d_9^0
        -Q_8\mathcal T-Cv^9S+\mathcal L_{78}.
 \end{split}                                             \tag{18}
\]
这里 \(\mathcal D_9\) 的系数可先缩小固定倍数，未花掉所有重复量正份额。
原 \(1-a\) 内部恒等式、标记集 I 和纯维数边界价格没有改动。

## 7. 同 I 九层抵消及实际质量反向界

第 1–8 层仍用本族原账本，不要求它们也有 \(\mathcal L_{78}\)。
第八、第七前缀是粗 \(Q(v)\mathcal T+O(v^8S)\)；
第六前缀是 \((K_n+o(1))v^6S\)，前五层是 \(O(v^6S)\)。
全部与 (18) 有同一 I、原半径、公共左侧和固定 \(c_{\rm mass},b\)。
用 \(\mathcal D_h\ge f_h\)、
\(d_h^0\le f_{h-1}+aP_\partial\)，令
\[
 \theta=bv/c_{\rm mass},\quad
 w_h=\theta^{9-h},\quad \Omega_9=\sum_{j=0}^8\theta^j.
\]
精确望远镜为
\[
 \sum_{h=1}^9w_h(c_{\rm mass}f_h-bvf_{h-1})
     =c_{\rm mass}f_9-bv\theta^8f_0.                    \tag{19}
\]
由 \(f_0=A_0vM'\le Ca vS\)，除 a 后初端为 O10。
先减同一个 \(\Omega_9Q_8\mathcal T\) 再除 \(\Omega_9\)，不能分开粗估权分母。
末层误差 O9，第八层 O8 乘 v 为 O9，第七层 O8 乘 \(v^2\) 为 O10，
第六层 \(o(v^6)\) 乘 \(v^3\) 为 o9；早五层最多 O10。
边界与内部价格由同一个原图条件吸收，\(\Omega_9=1+O(v)\)，
故末层新正量至少保留固定一半：
\[
 (1+Q_8)\mathcal T-A_vM/a
       \ge \tfrac12\mathcal L_{78}-Cv^9S-o(v^9S).       \tag{20}
\]

另一方面，定义已有解析八阶上包
\[
 \mathcal R_{Q_8}(v)=\frac{A_v}{A_v-a(1+Q_8)}
                          =P_{8,n}(t)+O_n(av^9).
\]
它的定义是明确有理式，不是未知真实 \(C_n\)。
对实际 \(x=M/V\)，恒等式为
\[
 (1+Q_8)\mathcal T-A_vM/a
 =\frac{A_v-a(1+Q_8)}a V
             [\mathcal R_{Q_8}-(1+x)].                 \tag{21}
\]
由 (1) 得其上界 \(C_Bv^9V\le C_B'v^9S\)。
最后一步只用粗临界质量比较 \(V\asymp S\)，不提取未知 \(S/V\) 的下一阶项。
将其与 (20) 合并，各非负量分别为 O9，即证 (2)。
同样由 (20) 丢正量可得原 \(P_8+O(av^9)\) 上界；本稿没有改变它。

量词保持序列口径：固定 n、B、有限个 \(\tau\) 及原有限工具之后，
常数不依输入数、最小半径或尺度跨度。LT 仅保证固定几何身份最终成立，
无速率模量只影响序列阈值，不能把它冒称为另一个幂余项。
若统一 (2) 失败，可取归一化正量趋无穷的 (1) 序列，
在该序列重建上述同 N9 预算，(20)–(21) 的固定常数给矛盾。
外部为空时源量均零，原内部计数处理，不定义空 CDF。

## 8. 可用后果与仍须保留的九阶量

固定参数下，原来出现的 \(Cv\mathfrak E,Cv\mathcal F,Ca\mathcal F\)
现在分别为 O10、O10、\(O(v^{n+9}S)\)。
这只是这些特定费用升级；侵蚀中的其他 O9 数值、inside 和 normal 标量余项未自动升级。
旧核心 raw \(R_*\) 仍只为 O8，\(CvR_*\) 仍可能九阶；不更名成新 raw。

由原 R071 已有 \(v^6D_N=O(v^9S)\)，真正缺槽坏根权
\[
 S({\rm Bad})\le C(D_N+v^{9/2}S)=O(v^3S).
\]
在同一 normal 整管和 R044 的固定高支撑矩形中，坏 owner 短窗 q3 的 gross
\(C_Kv^7S({\rm Bad})=O(v^{10}S)\)；长窗 Young 固定小 \(\eta\) 余项
\(C_\eta v^{15/2}S({\rm Bad})=O_\eta(v^{21/2}S)\)。
这里仍须保留 \(\eta Z(K)\)，不能把带净价的句子冒充 gross。
q≥4 的坏短窗至少 O11，坏长窗更晚；目标例外另按实际缺量结算，不由此删除。

对好 owner，实际固定面对路径的结论本来是
\[
 I_e\le C[v^{3q}r_i^n+v^{q-1}\mathcal E_i r_i^{n-1}],\qquad q\ge3.
\]
在预定有限稳健面列表包含路径所用的 \(\tau_*=1/128\) 时，
新增 (2) 只改善其 gap／端移第二项；第一项 \(v^{3q}\) 没有改善。
q3 的 \(v^9S\) 真实项仍在，可与 signed triangle 机制同阶。
本稿绝不宣称全部 high support 为 o9，也未给完整九阶价格或其极大值。

## 9. 实际读取、计算及冻结状态

实际使用的数学源：

* R063 foundation proof_record.md，SHA 702a70cdd735883d5adbe3d22da5acbe950e665c3023457a76416dcc5dd4f572，G3–G5、G7–G12；
* R063 eighth reconstruction.md，SHA 16401f6d8b9688d94ae115b2d5df1974dc36accbc404e102675e479de1e11a33，§3 强 CDF 及同族正量结算；
* R017 reviewed，SHA d045835884aa6c3c8169e231d67b743174c34b170094f47761b7869081b8cd37，§1 实际稳健面条；
* R071 fixed_tube_recentered_settlement.md，SHA a014ae167a2d83a30fb2b455d827a2b13253fce5033f9223ca9c674956e6d5f2，实际 N9、固定 K／参数吸收、真实 M 及同族前缀；
* R075 v2，SHA 658fc3de8b8971b290d758606c8b29ecddfd3bda16bb12a2f827575cf6458830，positivebase 侵蚀及一次 J 的第三反馈；
* R062 quantitative_complement.md，SHA 96f536af48cb942fe2e42f7043ebc3539bbc6a63737eb8fb4ef45ecfe814754f，同份 CDF 的端移余量、所有固定参数费用；
* R061 high_support_rate.md，SHA a037675512160f4717f49f37f466f0579fb386e8cdc67078fa1a4a34b8b9c663，仅重开好／坏 owner 显式实际公式，不照搬其旧 N8 身份；
* 正式 global 的 free-face-strips、common-defect、interior-recovery 原证明，实读 R068 project/sections/global.tex 对应段；它们的有限深度延长已在本次 N9 执行。

本机 Mathematica 15 的 end_face_ledger.wl/log 已实际成功运行，只窄核 (19)、(21)、质量归一化和显示幂次；
条避全部 X／标记区域、共同原见证环境、一次 CDF 份额由 §§2–7 的集合证明承担。
九层抵消、真实反向比值及端移归一化残差均为 0，全部精确断言 True。
先试用的 Python 环境缺少 SymPy，未产生计算证据；随后改用已有 Wolfram 内核，未安装软件。
计算记录与 SHA 另存 manifest，未复跑旧批次。

**当前状态：完整实际同族候选 (2)；冻结后交非作者独立审核。
不修改既有冻结源，不把作者自核视为终审，不向正文预供未经独审的新结论。**
