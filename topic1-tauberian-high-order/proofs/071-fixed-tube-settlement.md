# R071：固定正常整管与移动最大点的九阶余项候选

2026-11-05。非作者独审前的完整数学候选。本文件只写独占 scratch，不改论文、前轮冻结源或他人材料。

## 0. 目标、版本和使用的基础

固定整数 \(n\ge2\)，\(a=v^n=2^nt,\ v=2t^{1/n}\)，
\[
 A_v=(1+v)^n,\quad K_n=\binom n2n^2/27,\quad
 P_{8,n}=L_n+K_nav^6-(n+26/9)K_nav^7+\widehat J_nav^8 .
\]
本候选的目标是
\[
 \boxed{C_n(1-t)=P_{8,n}(t)+O_n(av^9).}                 \tag{T}
\]
不识别九阶系数，也不宣称全阶纯幂形式。

实际全文读取并使用本轮解析重建：

* R063 foundation proof_record.md，SHA 702a70cdd735883d5adbe3d22da5acbe950e665c3023457a76416dcc5dd4f572；
* R063 sixth rebuild.md，SHA 3a5d38fe5a00392af687c09ad59fcf072c552cbe418504af12fc38dc93e3cb2c；
* R063 eighth reconstruction.md，SHA 16401f6d8b9688d94ae115b2d5df1974dc36accbc404e102675e479de1e11a33；
* 新 N9 的实际构造 power_remainder_r064.md，实读 SHA 80f92dc770da78284a56e56a331b82112ea5912e892a456e968078e412a25699；
* 统一局部证明 R062 quantitative_complement.md，SHA 96f536af48cb942fe2e42f7043ebc3539bbc6a63737eb8fb4ef45ecfe814754f。

这里使用的是 R064 §1 的九步实际构造、粗先验和粗稳定性，不以其主结论反供本证明。
所需对象及来源在下节明确列出。R059 宽度结论不直接改名用于 N9，而在 §2 重做所需论证。
R062 的增长窗率也不是本证明的输入；本稿重开其有限参数推导保留正项。
旧 AI 裁决不是任何一个数学前提。

## 1. 同一新 N9 族和旧数值预算

取任意真实输入及严格外部 halo 紧集序列，令 \(V=|E|>0,M=|K|\)，假设
\[
 1+M/V\ge P_{8,n}(t)-Bav^9,\qquad B<\infty .             \tag{1}
\]
它最终进入固定的 \(P_7-B_0av^8\) 竞争分支。以下全部原半径、I、九层、CDF 及正常集合一次固定。

新九步构造为：原有限正值 Borel r(x) 不变；在原图选同一 I 后，
第 h 步先 \(H_h>L_{h-1}+D+1\)，取深目标 \(G_{i,h}\)，深处未删除源
\(T_{i,h}=N_{i,h-1}\)，边界仍取原 \(P_i\)，仅删除原小半径，
再用收缩族唯一覆盖定位选择 Lh。九步后统一 A9，继而统一 CDF。
由实际原捕获和线性半径控制，每个 \(x\in F_{i,h}\) 的原见证满足
\[
 G_{i,h}\subset Q(x,r(x)),\quad R_{i,h}<r(x)\le(1+A_9v)r_i;
 \quad F_j\cap G_i=\varnothing\ (j\ne i).
\]
盒宽可比于 \(vR_i\)，原始／有效／中点端盒、真面邻、全横捕获、
实际 clean 四槽、全部尺度固定窗和有限相交度均在本层重建。
这些具体步骤由上述 foundation G10、sixth §4/§14、R064 §1 展开；
不能将旧 N8 末层改名。以下未标层号的量取本次第九层。

记 \(\mathcal T=V+M,\ S=\sum r_i^n,\ g_i=aR_i^n,\mu_i=|F_i|\)，
\[
 D_i^0=(x_i+W_i^0)\setminus F_i,\quad
 D_i'=(x_i+W_i')\setminus F_i,\quad
 \lambda_v=s_vv/a,\quad s_v=(1+v)^{n-1}.
\]
旧正常实际保证层 \(Y_i^{old}\)、原正常集合 N、旧管和旧共同层均为本次 N9 的固定对象。
粗 bootstrap 在该族重放后给固定数值
\[
 R_*:=\frac va\sum|D_i^0\setminus Y_i^{old}|=O(v^8S),
 \quad F+F_{1/128}+\mathfrak E+P_I+P_\partial=O(v^8S),
 \quad D_N=O(v^2S),                                      \tag{2}
\]
以及 \(\sum\varepsilon_iR_i^n=O(v^5S)\)，旧正常层逐源体积 \(O(av^5R_i^n)\)。
这些先验不含本稿目标，也不由 \(P_8\) 的未知全局上界导出：
在本 N9 的 m=7,8,9 末层重放粗联合价和一次面恢复，先得 \(P_7+O(av^8)\)，
再在同次证明保留正份额并用真实比值反向界，即得 (2)。
普通／稳健面条重叠只取下界半和，没有重复面收入。

对九层有真实质量和公共账本
\[
 \mathcal D_h\ge f_h\ge0,\quad
 d_h^0\le f_{h-1}+aP_\partial,\quad f_0=A_0vM',
\]
\[
 \sum_i(g_{i,h}-\mu_{i,h})=\mathcal D_h+d_h^0+\Delta_h^{mass},
\]
\[
 \mathcal T\ge A_vM/a+c_{rec}P_I-C_bP_\partial-2nvS_I
 +\frac{\mathcal D_h}{2a}-b(v/a)d_h^0-E_h.                \tag{3}
\]
第七、第八前缀各可取 \(E_h=Q\mathcal T+O(v^8S)\)，
\(Q=K_n[v^6(1+v)^{-4}+(10/9)v^7(1+v)^{-2}]\)；
第六前缀为 \((K_n+o(1))v^6S\)，早五层粗界 \(Cv^6S\)。
全部层同一个 I、原半径和公共左侧。纯维数图条件留有固定边界吸收预留。

## 2. 在这个 N9 上重新得到宽度 L1 \(O(v^{7/2}S)\)

写 \(b_i=|W_i|^{1/n}\)、\(q_i=\sum_l\log^2(w_{il}/b_i)\)，
\(\mathcal A_i=|X_i|-(R_i+|W_i|^{1/n})^n\)。
实际恒等式是
\[
 |X_i|=(R_i+\mu_i^{1/n})^n+J_i+\mathcal A_i .
\]
令 \(s_l=\log(w_{il}/b_i)\)，\(\sum s_l=0\)，
归一化盒宽在固定正紧区间。对
\(\sum_l\log(1+(b_i/R_i)e^{s_l})\)，沿线段的 Hessian 下界为 cv，
零点的一阶项相消，故
\[
 \mathcal A_i\ge cvq_iR_i^n.                              \tag{4}
\]
它独立于 J；在原粗账本里保留末层 \(\sum\mathcal A_i\) 不重复盈余。

在 (3) 把重复量系数 \(1/2\) 拆为两份 \(1/4\)，一份用 \(f_h\) 消项，
另一份保留。取 \(\tau=4bv,\ w_h=\tau^{9-h}\)，九层正权平均的初端为 \(O(v^{10}S)\)。
末层保留形状，另有
\[
 \Omega^{-1}\sum w_h\mathcal D_h/(4a)
 \ge c(\mathcal D_9/a+v\mathcal D_8/a).
\]
同中心 Q 的加权误差是 \(O(v^8S)\)：第九层八阶，第八层误差乘 v，
第七及更早粗价最多从 \(v^2v^6\) 起。真实竞争比值也只需
\(1+M/V\ge P_7-B_0av^8\)，所以按真实 M 的反向恒等式得
\[
 v\sum q_iR_i^n+\mathcal D_9/a+v\mathcal D_8/a=O(v^8S).
                                                               \tag{5}
\]
此为本 N9 的证明，没有引用 N8 的同名正量。

设 \(\eta_i=1-\mu_i/g_i\)。旧 dyadic/盒外/边界删除公式与 (2) 给
\(\Delta_9^{mass}/a=O(v^7S)\)（最坏 n=2），从而
\[
 \sum\eta_iR_i^n
 =(\mathcal D_9+d_9^0+\Delta_9^{mass})/a=O(v^7S);
\]
这里 \(d_9^0\le\mathcal D_8+aP_\partial\)，故前层重复量不能省。
又 \(0\le\eta_i\le1\)，其平方和同阶。
由
\[
 \varepsilon_i\mu_i\le|D_i^0|
 \le Cav^5R_i^n+|D_i^0\setminus Y_i^{old}|
\]
及 raw 归一化高度一致有界，得 \(\sum\varepsilon_i^2R_i^n=O(v^7S)\)。
端移 \(\rho_i=\Delta_i/(vR_i)\) 一致有界，\(\sum\rho_i^2R_i^n=O(v^7S)\)。
由 \(|W_i|/g_i=(1+\varepsilon_i)(1-\eta_i)\) 和三端盒边长差不超过 \(\Delta_i\)，
同一个最大相对宽差 \(\xi_i\) 满足
\[
 \xi_i^2\le C(q_i+\varepsilon_i^2+\eta_i^2+\rho_i^2).
\]
因此
\[
 \sum\xi_i^2R_i^n=O(v^7S),\qquad X_\xi:=\sum\xi_iR_i^n=O(v^{7/2}S). \tag{6}
\]
所有异常源保留；只用加权平方和，没有逐源小偏差。

## 3. 真正移动的 normal 最大点

原正常 whole-chart 模型为 \(Y_{2,N}\)，定义
\[
 \mathsf P_N=\sum_{e\in N}I_e-\lambda_v\sum|Y_{2,N,i}|,
 \quad W_N=\sum_{e\in N}b_eR_e^2 .
\]
它没有扣全部 J。实际正常积分保持两个不同根帽及三近面名义层。
固定正常窗的全部非宽度误差为 \(O(v^9S)\)：
第二次实际反馈替换、顶端 clip、非中央片／其他角、三近面退款、
端锚退款 \(Cv\mathfrak E\)、根比系数、双帽统一余项。
双帽误差为 \(C(v^7z^2+v^8|z|+v^9)R_e^n\)，用旧 RMS 得九阶。
宽度名义误差由 (6) 为 \(Cv^6X_\xi=O(v^{19/2}S)\)。
这些是实际积分的有率误差，不是把旧 little-o 改名。

归一化为 \(x=(1+v)p/n,y=(1+v)q/n,\zeta=(1+v)z/n\)，
令单边标量 \(b=\bar b_e:=\min(1,\widehat b_e)\)，固定正常窗内 \(b\ge1-\nu\)；
物理横投影量 \(b_e\) 仍只用于 \(W_N\)，不与 \(\bar b_e\) 混用。
定义
\[
 Q_4(x,y)=xy(x^2+3xy+y^2)/3,\quad
 Q_5(x,y)=xy(x+y)(5x^2+29xy+5y^2)/24,
\]
\[
 \Phi_v(x,y)=xy(1-x-y)+\frac{2v}{1+v}Q_4
 -\frac{2v^2}{(1+v)^2}Q_5-\frac{n(n-2)v^2}{1+v}x^2y^2.    \tag{7}
\]
零参数 Hessian 在 \((1/3,1/3)\) 为
\(\left(\begin{smallmatrix}-2/3&-1/3\\-1/3&-2/3\end{smallmatrix}\right)\)，
特征值 \(-1,-1/3\)。先取固定小正常窗，Hessian 在其内统一负定；
固定二次强制性又使所有最大点在该窗内部且距基点 \(O(v)\)。
隐函数定理及对称性给唯一最大点 \((r_v,r_v)\)，
\[
 r_v=1/3+(20/81)v+O(v^2).
\]
所以亏损必须围绕 \(r_v\)，不能围绕固定 \(1/3\) 强行改进旧 \(G_N\)。

带 b、\(\zeta\) 的单边泛函为
\[
 \Phi_v(x,y)-(1-b)xy-3(x+y)\zeta^2.
\]
它也等于 \(b\Phi_v(x,y)-(1-b)[xy(x+y)-\text{扰动}]-3(x+y)\zeta^2\)。
方括号在固定小正常窗内有固定正下界。因此统一有
\[
 F_v(b,x,y,\zeta)\le b\Phi_v(r_v,r_v)
 -c[(1-b)+(x-r_v)^2+(y-r_v)^2+\zeta^2].                  \tag{8}
\]
令
\[
 \mathcal L_N=v^6\sum_{e\in N}
 [(1-\bar b_e)+(x_e-r_v)^2+(y_e-r_v)^2+\zeta_e^2]R_e^n\ge0.
\]
只在真实最大点作二阶 Taylor，三阶余项在固定紧窗上一致 \(O(v^3)\)。
两阶值与已有 \(c_{8,n}\) 配方相同，所以实际每权峰值
\[
 q_{\rm peak}(v)=\frac{n^2v^6}{(1+v)^2}\Phi_v(r_v,r_v)
 =q_N(v)+O(v^9),
\]
\[
 q_N=\frac{kv^6}{(1+v)^2}+\frac{\Gamma v^7}{(1+v)^3}
       +kc_{8,n}v^8 .
\]
横权截幅和实际根宽误差分别由 (6) 支付。因而
\[
 \boxed{\mathsf P_N\le q_NW_N-c\mathcal L_N+Cv^9S.}       \tag{9}
\]
下面其实只需其上包；保留 \(\mathcal L_N\) 另给真实 recentered 稳定性。

## 4. 一次固定的整个 normal 管

固定正常窗后，对每个 normal 有向目标，原始帽参数
\(P_{ij,l},P_{ij,k}\le C_Nv^3R_i\)，正常名义总高度 \(\le C_Nv^4R_i\)。
这些常数与后选 T,m 无关。

在有效源 global-min chart 中，若 \(x(s,y)\in B_N(0)\)，
至少一个 normal 目标的帽质量为正，故其两个临界轴
\(y_l,y_k\le C_Nv^3R_i\)，并且 \(s\le C_Nv^4R_i\)。
这两轴的物理有效面深均不超过 \(C_N(v^3+v^4)R_i\)。
此论证也包括最浅轴不在该 pair 的整张图；没有删除三近面。

因此可一次选择有限固定 \(K\ge K_N^{old}\)，只依正常窗和 A9，使
\[
 B_N(D)_i\subset B_N(0)_i\subset\mathcal C_i(K)
 :=\{x\in x_i+W_i':\text{至少两轴距某有效面}\le Kv^3r_i\}. \tag{10}
\]
它包含的是整个正常名义模型层 \(B_N(0)\)，不仅旧核心层 \(Y_i^{old}\)；
K 因而连测试缺集 D 都不依赖。
任意后来共享侵蚀的 \(B_N^-\) 也在此管内。
K 不依 T、m、输入数、尺度跨度或迭代次数。

对任意 packet \(A=N\cup B\)，同 chart 上 \(h_A=h_N+h_B\ge h_N\)，
用同一个非负侵蚀和同一个真实上截断，令 \(B_N^-\subset B_A^-\subset D'\)。
记
\[
 \mathsf M=\lambda_v\sum\int h_B^0,\qquad
 \mathcal R_A^-=\frac va\sum|D_i^0\setminus B_{A,i}^-|.
\]
正部／clip 的单调 1-Lipschitz 性给
\[
 \frac va|B_A^-\setminus B_N^-|_1\le s_v^{-1}\mathsf M .
\]
在管外把缺集分为 \(D'\setminus B_A^-\) 和 \(B_A^-\setminus\mathcal C(K)\)，
后者因 (10) 包含于 \(B_A^-\setminus B_N^-\)。故得到精确几何式
\[
 \boxed{Z(K):=\frac va\sum|D_i'\setminus\mathcal C_i(K)|
          \le\mathcal R_A^-+s_v^{-1}\mathsf M.}           \tag{11}
\]
这里没有侵蚀误差项，也没有把旧 \(R_*\) 改名为新 raw。

## 5. 二轴联合价中的固定正份额

以下所有常数先于 T,m 固定。全域窗外标量强制性
\(G\ge c_\nu(b+D)\) 来自固定正常窗，与 T,m 无关。
任意 Y 的强 CDF raw／端移系数也只依固定 A9。
由它们分别预留足够小的固定 \(8c_B,8c_R>0\)，以及由角权折扣给固定 \(8c_D>0\)。
取 \(8c_B<1\)，三近面名义量始终非负，可以全部保留在左侧的 \(8c_B\mathsf M\) 中。

重开 R062 的证明，而不是引用其已经上包掉正量的 \(O(v^{65/8})\) 式。
对任意先固定的有限 T,m，实际层增量为
\(\Delta_A=\mathsf M-\mathfrak F_B-\mathfrak L_B\)，且
\[
 \mathfrak F_B\le CTm^{-3}v\,\mathsf M+CTvR_*,
 \qquad \mathfrak L_B\le CT^2v^3\mathsf M.                \tag{12}
\]
必须保留 \(\mathsf M\) 因子；不先以 \(O(v^8S)\) 把这两项变成绝对误差。
根比的 \(Cv\mathsf M\) 亦可在此吸收。

whole-height 必要性侵蚀成本 \(CT^5v^{10}S\)，固定正常反馈替换 \(Cv^9S\)，
导数 epsilon 费 \(CT^3v^{11}S\)，整宽费 \(CT^3v^{19/2}S\) 由 (6)；
固定 T,m 后它们合计 \(C_{T,m}v^9S\)。
本段没有让 T,m 随 v 增长，所以不存在未知 fixed-window 阈值的非法代入。

用 (9)，先对 bounded 非正常类作强制性，再以 actual soft 槽作一次转换，
而把长角和近零 gross 尚保留在右侧，得到
\[
 \begin{split}
 H_2+8c_B\mathsf M+8c_R\mathcal R_A^-+8c_Dv^6D_N
       +c\mathcal L_N
 \le q_N\mathcal C+
       \theta_{T,m}(v)\mathsf M+
       I_{\rm long}+I_{\rm small}+C_{T,m}v^9S,\\
 \theta_{T,m}(v)\le C(Tm^{-3}v+T^2v^3+v).
 \end{split}                                             \tag{13}
\]
完整 J 只在联合真实层中使用一次。
soft 费为 \(C(1+C_s)v^3S\)，取 \(C_s=C_0/m\) 后乘六阶价为 \(C_m v^9S\)。
其 normal clean 和 helper 阈值未改变。正槽系数来自
\(q_N-(1-\gamma_\nu)A_6\ge c v^6\)，先于 T,m 固定；
没有另加正常模型单独转换的第二份槽收入。
单独的 CDF 端移正份额可以丢弃；本稿不要求它改善。

## 6. 长条真在整个固定管外

长 q2 边的大根 owner 记
\(P=\max L/r_i,\ Q=\min L/r_i,\ H_0=\min(Q,PQ/v^2)\)。
原实际缺源条的深度为长轴 \([c_2P,2c_2P]r_i\)，
短轴 \([c_3H_0,2c_3H_0]r_i\)，其余轴取到两有效面均 \(>Kv^3r_i\) 的核心。
K 是 (10) 的固定常数；固定小 v 后该核心仍至少保留每宽的一半。

选择固定 \(T\ge T_0(K)\) 足够大，使长轴近面深 \(>Kv^3r_i\)；
相反面仍距至少 \(3/4\) 有效宽，亦 \(>Kv^3r_i\)。
所以条内至多短轴可接近任何一面，整条避开整个 \(\mathcal C_i(K)\)。
不是仅避开列出的旧正常角层。

原面对端距离、真实目标截宽和同一个原见证证明该整条缺源。
其成本／gross 比至少
\[
 c\min(1/v,P/v^3)\ge cT
\]
（固定 T 后最终 \(1/v\ge T\)）。
有限相交度统一控制条的重叠。目标 raw 例外只用旧层，
费用为 \(CvR_*\)，面对／端移例外为 \(Cv(F+\mathfrak E)\)。
更明确地，目标原缺集先分为旧正常层、有效盒内其余 raw 和原／有效壳。
壳不属于 \(Z(K)\)，不把它记进 Z；其费用有独立界
\[
 (v^2/a)\sum|B_j^0\setminus B_j'|
 \le Cv\sum R_j^{n-1}\Delta_j=Cv\mathfrak E=O(v^9S).
\]
剩余 target raw 以 \(CvR_*\) 上包。也可把壳包含在旧 \(R_*\) 中，
但不能把目标原盒管外缺量直接等同有效盒的 Z。
故在这个同一 K 上有
\[
 I_{\rm long}\le(C_L/T)Z(K)+Cv(R_*+F+\mathfrak E)
 \le(C_L/T)(\mathcal R_A^-+\mathsf M)+Cv^9S.              \tag{14}
\]
\(C_L\) 不依 T,m；核心比例常数不依 K，K 只影响最后固定小 v 阈值。
原正常层在目标帽中的占比 \(C(v^3/P+v^2)\) 一致于任意小 Q，
所以目标例外收费没有隐含短轴正下界。

对近零 bounded 边，\(pq\le Tm\,\omega_e\)，实际槽给
\[
 I_{\rm small}\le CTm\,v^6D_N+C_Tv^9S.                  \tag{15}
\]
不以 \(D_N=O(v^2S)\) 把首项提前变成八阶误差；
它将由保留的槽正份额直接支付。

## 7. 高支撑只花一次剩余份额

固定 K 后重开原实际高支撑矩形证明，取任意随后固定的小 \(\eta>0\)。
好 owner 的稳健面路径由 (2) 给 \(O(v^9S)\)。
坏 owner 集满足 \(S(\mathcal B)\le C(D_N+v^{9/2}S)\)；
目标 epsilon 的局部 raw 分解及面／端移例外从十阶起。
坏长窗的共同 Young 价为
\[
 \eta Z(K)+C_\eta\sum_{q=3}^n v^{q(2q-1)/(q-1)}S(\mathcal B).
\]
矩形在整个固定 K 管外，所有 q 共用度数 D 和一份 \(\eta Z\)。
最早 q3 指数 \(15/2\)，用本次旧稳定数值 \(S(\mathcal B)=O(v^2S)\)
后，长窗净余从 \(19/2\) 阶起。坏短窗从九阶起。
这里只把已经证明的旧坏根数值作为上界，不再消费一份槽收入。
边界高支撑为 \(O(v^{19/2}S)\)，故
\[
 I_{\ge3}\le\eta Z(K)+C_{\eta,K}v^9S.                    \tag{16}
\]
二维为空和。这是重开固定份额净价后的式子，不是把旧 gross \(43/5\) 改成九阶。
只有这一固定份额在下面并入完整 H；没有再减第二份 J。

## 8. 参数次序与直接闭合

先固定原 N9、normal、CDF、normal whole-model 管 K，以及 (13) 的正系数。
依次选择：

1. 固定大 \(T\ge T_0(K)\)，使 \(C_L/T\le\min(c_B,c_R)\)；
2. 固定小 \(m>0\)，使 \(CTm\le c_D\)，并 \(m<c_N/2\)；
3. 固定 \(\eta\le\min(c_B,c_R)\)；
4. 最后取 v 小，使 \(\theta_{T,m}(v)\le c_B\)，且所有原几何、
   高支撑固定-\(\eta\) 和真实帽／clip 条件同时生效。

coercivity、CDF 和槽的正系数均先于 T,m。
K 只含固定 normal whole-model，不是覆盖整个新 packet 的大管；
因而第一步没有 \(K(T,m)\) 循环。m 只影响固定常数和最终阈值。

把 (14)–(16) 代入 (13)，M、raw 和槽分别最多支付
\(3c_B\mathsf M,\ 2c_R\mathcal R_A^-,\ c_Dv^6D_N\)；
原预留有八份，仍留固定正份额。因此
\[
 \boxed{H+c_B\mathsf M+c_R\mathcal R_A^-+c_Dv^6D_N
                  +c\mathcal L_N
       \le q_N\mathcal C+C_nv^9S.}                       \tag{17}
\]
所有数值预算均来自同一次实际族；重叠条由共同度数和同一次 Z 分份。
旧 \(R_*=O(v^8S)\) 只出现在乘 v 的费用中，未声称它改善。
旧固定中心 \(G_N\) 也未声称改善。

这一步比有限 bootstrap 更直接：没有把 M 的小系数费用先粗估为八阶再逐次改善。
若只想 \(9-\epsilon\)，它是 (17) 的弱推论；不需要无限迭代或增加递归深度。

## 9. 实际质量、inside 与九层总装

整体小 v 时 \(q_N>0\)，完整实际槽乘积和、(6) 及质量恒等式给
\[
 q_N\mathcal C\le Q_{8,N}(v)\mathcal T+O(v^{19/2}S+v^{12}S),
\]
\[
 Q_{8,N}=K_n[v^6(1+v)^{-4}
 +(10/9)v^7(1+v)^{-5}+c_{8,n}v^8(1+v)^{-2}].
\]
不能把软接触上界反向当槽容量，也不逐项按负八阶系数拆乘。

对终端精确并集可以只用 \(\bigcup X_i\)，不用再次加入免费面条；
粗稳定性阶段的一次面恢复只用于已经证明的数值界 (2)。
Bonferroni 使用完整 H。基线真实导数在整个 \(\mu/g\in[1-A_9v,1]\) 区间为
\([1+O(v)]v/a\)，删除质量与未删除 d0 分段。
深处盒外费用 \(CvF=O(v^9S)\)，环境 \(Ca(F+P_\partial)\) 更晚。
n≥3 的盒内／边界绝对项 \(v^{n+6}S\) 从九阶起。

二维在本 N9 原有限值半径上重放全点入向方块局部费用：
坏宏盒权 \(O(v^2S)\)，好宏盒 raw／扩张和旧正常层四角窗给
\[
 (1+O(v))(v/a)\sum\delta_{in}
 \le(8/27+O(v))v^8S.
\]
这是实际全点覆盖；旧参数 RMS 给其领先正常体积误差九阶。
其边界／环境十阶，盒外九阶。令
\(\kappa_n=(8/27)\mathbf1_{n=2}\)，则 inside 成本
\(\kappa_nv^8\mathcal T+O(v^9S)\)。

故本次第九层可取
\[
 E_9=Q_8\mathcal T+O(v^9S),\quad Q_8=Q_{8,N}+\kappa_nv^8,
\]
并同时保留 (17) 的正剩余。第八、第七前缀仍用各自的粗 Q，
与 \(Q_8\) 差 \(O(v^8)\)；第六层 leading，前五粗价。
权 \(\tau^{9-h}\) 的精确消项为
\[
 \sum_{h=1}^9\tau^{9-h}(c_D'f_h-bvf_{h-1})
       =c_D'f_9-bv\tau^8f_0,\quad \tau=bv/c_D'.
\]
初端十阶，整权平均先减同一 \(Q_8\mathcal T\) 后四类误差为
\[
 O(v^9S),\quad O(v^9S),\quad O(v^{10}S),
       \quad o(v^9S)+O(v^{10}S).
\]
内部精确式的 \(1-a\) 和纯维数边界预留不变，图条件吸收后
\[
 \mathcal T\ge A_vM/a-Q_8\mathcal T-C_nv^9S.              \tag{18}
\]
直接移动实际 M，得
\[
 1+M/V\le\frac{A_v}{A_v-a(1+Q_8)}+C_nav^9
         =P_{8,n}(t)+O_n(av^9).                          \tag{19}
\]
真实有理分母的后三阶 Taylor 误差为九阶；n2 分母和 inside 两笔 \(8/27\) 分开保留。
外部为空时源量取零，不定义空 CDF，原内部计数给同方向。

若全输入 (19) 失败，先选超额／\(av^9\) 趋无穷的实际输入序列，
再以内正则性取外部严格 halo 紧集，损失小于 \(Vav^9\)。
该序列最终满足 (1)，如 B=1；同一次 N9 工具及本稿固定 T,m、\(\eta\) 后，
(19) 的固定常数矛盾。只需固定精度 LT 最终成立，无新的 LT 速率。
每个小实参数已有真实有限 lower \(P_8-C_nav^9\)，遂与上界合并得到 (T)。

同样在 (1) 竞争族中，保留 (17) 的正量再用真实比值反向界，得到
\[
 \mathsf M+\mathcal R_A^-+v^6D_N+\mathcal L_N=O_{n,B}(v^9S),
 \qquad Z(K)=O_{n,B}(v^9S).
 \tag{20}
\]
这是固定测试 packet 的非正常剩余与 moving-normal 亏损。
不包含旧 \(R_*\)，不包含固定中心 \(G_N\)，不宣称每一 normal 边逐点正常。

## 10. 实际计算和当前状态

本机实际运行 bootstrap_scalar.wl，Wolfram 15 ARM，PID 79446。
日志给 normal 基点 Hessian 特征值 \(-1,-1/3\)，移动梯度一阶残差 0，
峰值与 qN 二阶残差 0；原保守有限 bootstrap 两支指数违反集均 False。
这些只核标量代数。管包含、无侵蚀 Z 比较、whole-pipe 长条和一次预算吸收
由 §§4–8 的实际集合与容量证明承担，不由指数 LP 承担。

本轮使用了近极值稳定性技能整理非负亏损和等号中心，但没有借用其示例作外部定理。
没有新增文献结论，也没有修改仓库或发布。

**状态：完整新候选 (T) 及 (20) 已给出上述实际证明，待非作者对冻结字节独立敌意审核。
尤其须核固定 K、whole-pipe 长条和 cB/cR/cD 的先后分份；作者不以本句代替独审。**
