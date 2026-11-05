# 固定有限深度的全点供货与误差传播

2026-11-05。本稿只写独占 scratch。目标是把已成立的全点构造整理为一个有准确量词的可复用接口，避免每增加一层都隐含改变见证、CDF 或正常集合。它不识别九阶系数，不把递归深度本身当作更高阶几何。

## 1. 量词、输入和对象

固定整数 \(n\ge2,N\ge9\)。下面所有常数允许依赖 (n,N)、固定竞争常数及预定有限工具，**不要求关于 N 一致**，不依输入数、最小原半径或尺度跨度。N 在 \(v\to0\) 前固定。

记
\[
 a=v^n=2^nt,\quad A_v=(1+v)^n,\quad V=|E|>0,\quad
 K\subset\mathcal H_n(1-t;E)\setminus E\text{ 紧},\quad M=|K|,\quad
 \mathcal T=V+M.
\]
原有限降序贪心给有限值 Borel 严格见证 \(r(x)>0\) 和 \(Q_i=Q(x_i,r_i)\)，满足
\[
 \|x_i-x_j\|_\infty>\max(r_i,r_j),\quad
 |Q_i\setminus E|<ar_i^n,\quad S=\sum_i r_i^n.
 \tag{1}
\]
每个 x 被某个 \(Q_j\ni x,r_j\ge r(x)\) 捕获。私有簇
\(P_i=(K\cap Q_i)\setminus\bigcup_{j\ne i}Q_j\) 上 \(r(x)\le r_i\)。原半盒不交给
\[
 M/a\le S\le V/(1-a).
 \tag{2}
\]
以下原图、规则性参数和内部／外部指标，使用正文已证明的定义；尤其标记集 I 包含全部非规则指标。外部规则簇满足私有质量下界与 \(P_i\subset Q(x_i,cDv r_i)\)。这些规则性条件来自原输入的固定图选层，不是新增输入假设。

先于深度取纯维数价格
\[
 c_{\rm rec}=1/(2^{n+1}),\quad C_b=2C_2(n)+2,\quad C_*=C_b+2,
 \quad C_2(n)=2n4^n3^{n-1},
\]
并在**原 (Q,P,S) 图**选择同一个 I，使
\[
 C_*P_\partial\le c_{\rm rec}P_I/4,\qquad
 2nvS_I\le c_{\rm rec}P_I/4.
 \tag{3}
\]
这里 \(K'=K\setminus\bigcup_{i\in I}Q_i,M'=|K'|\)，\(M_I=M-M'\)。
图只看原数据；后来的 N、\(H_h,L_h,A_N\) 与 CDF 不反馈改变它。原图选层的固定预留可先多留一份，以吸收本文最后出现的 \(C_NvP_\partial\)，因为对每个固定 N 最终 \(C_Nv\le1\)。这不要求把 \(C_*\) 改成依赖 N 的价格。

记
\[
 L_n=A_v/(A_v-a),\quad K_n=\binom n2n^2/27,\quad
 P_{7,n}=L_n+K_nav^6-(n+26/9)K_nav^7,\quad
 P_{8,n}=P_{7,n}+\widehat J_nav^8,
\]
\[
 \widehat J_n=\frac{(n-1)n^3(2609+1971n+81n^2)}{26244}
                  +\frac{16}{27}\mathbf1_{n=2}.
\]
为供给价格与稳定性，分两个实际竞争分支：
\[
 \text{(C8)}\quad1+M/V\ge P_{7,n}(t)-Bav^8,
\qquad
 \text{(C9)}\quad1+M/V\ge P_{8,n}(t)-Bav^9.
 \tag{4}
\]
B 固定有限；C9 最终蕴含 C8（取一个固定更大的常数）。构造本身只用原规则族；下文数值结论在相应分支中给出。

## 2. 供货定理

**定理 A（固定 N 的实际供货）。** 在上述规则下，可先固定 I，再按每步 \(H_h\to L_h\) 执行全部 \(h=1,\ldots,N\)，最后统一 \(A_N\) 和一个强 CDF 斜率。所有原指定见证保持不变。对每个 h 得
\[
 G_{i,h},\ T_{i,h},\ F_{i,h},\ N_{i,h},\quad
 g_{i,h}=|G_{i,h}|,\ m^0_{i,h}=|T_{i,h}|,\ \mu_{i,h}=|F_{i,h}|,
\]
\[
 R_{i,h}=(g_{i,h}/a)^{1/n},\quad
 \mathcal D_h=\sum_i g_{i,h}-M',\quad
 d_h^0=M'-\sum_i m^0_{i,h},\quad
 \delta_h=\sum_i(m^0_{i,h}-\mu_{i,h}).
\]
它们满足：

1. 全点完整目标捕获、源／异目标排他和严格原半径：
\[
 F_{i,h}\subset T_{i,h}\subset G_{i,h}\subset K'\subset E^c,
\quad T_{i,h}\cap G_{j,h}=\varnothing\ (i\ne j),
\]
\[
 G_{i,h}\subset Q(x,r(x)),\quad
 R_{i,h}<r(x)\le(1+A_Nv)r_i\quad(x\in F_{i,h}),
\quad(1-A_Nv)ar_i^n\le\mu_{i,h}\le g_{i,h}\le ar_i^n.
 \tag{5}
\]
边界源仍在原 \(P_i\)，有更强 \(r(x)\le r_i\)。
2. 目标覆盖 (K')，且有精确递归与质量式
\[
 \mathcal D_h\ge f_h\ge0,\quad
 d_1^0\le f_0=A_0vM',\quad
 d_h^0\le f_{h-1}+aP_\partial\ (h\ge2),
\]
\[
 f_h=M'-\sum_i|N_{i,h}|,\qquad
 \sum_i(g_{i,h}-\mu_{i,h})=\mathcal D_h+d_h^0+\delta_h.
 \tag{6}
\]
3. 每层都有自身的三端盒、面邻和全点小源回收；每个本族前缀末层 \(\ell=6,\ldots,N\) 在粗六阶先验后另有本层正常集合、actual clean/helper 与槽。§4 标明各数值结论所需的最低前缀长度；O8 raw、准确 inside 等仅在 \(\ell\ge8\) 调用，不对早层越界供货。不会混用另一深度的盒、根、保证层或半径。
4. C8 下，每个 \(\ell=6,\ldots,N\) 有粗六阶正先验；每个 \(\ell=7,\ldots,N\) 可取同一中心价
\[
 E_\ell=Q(v)\mathcal T+O_N(v^8S),\qquad
 Q(v)=K_nv^6[(1+v)^{-4}+(10/9)v(1+v)^{-2}].
 \tag{7}
\]
每个 \(\ell\ge8\) 有同层八阶正预算及完整宽度
\[
 \mathfrak G_\ell+\mathsf M^0_{B,\ell}+R^0_{\rm all,\ell}+R^0_{N,\ell}
 +\mathfrak E_\ell+\mathcal F_\ell+\mathcal F_{1/128,\ell}
 +P_I+P_\partial+v^6D_{\mathcal N_\ell}=O_N(v^8S),
 \tag{8}
\]
\[
 \sum_i\xi_{i,\ell}^2R_{i,\ell}^n=O_N(v^7S),\quad
 \sum_i\xi_{i,\ell}R_{i,\ell}^n=O_N(v^{7/2}S).
 \tag{9}
\]
上标0标记粗 bootstrap 的固定旧 packet 与旧核心 raw，避免与后选 packet 混淆。
5. C9 下，每个 \(\ell=9,\ldots,N\) 的**本族本层**可取
\[
 E_\ell=Q_8(v)\mathcal T+O_N(v^9S),
 \tag{10}
\]
并可在同一层保留固定 packet \(A=\mathcal N_\ell\cup B\) 的正预算
\[
 \mathsf M_{B,\ell}+\mathcal R^-_{A,\ell}
 +v^6D_{\mathcal N_\ell}+\mathcal L_{\mathcal N_\ell}
 =O_{N,B}(v^9S),\qquad Z_\ell(K_N)=O_{N,B}(v^9S).
 \tag{11}
\]
\(K_N\) 是包含整个正常名义模型的固定管，允许依固定 N，独立于后选 T、\(\mu_*\)。旧 \(R^0_{N,\ell}\) 仍只断言 O8。

这里
\[
 Q_8=K_n\left[\frac{v^6}{(1+v)^4}
 +\frac{10v^7}{9(1+v)^5}+\frac{c_{8,n}v^8}{(1+v)^2}\right]
 +\frac8{27}\mathbf1_{n=2}v^8,
\quad c_{8,n}=-13/18-n(n-2)/3+400/243.
\]
定理不供未知九阶价格，也不将当前 inside 的 O9 余项升级为 o9。

## 3. 任意固定 N 的全点归纳

取 \(L_0=C_0>cD+1\)。对任意预定固定 \(H_*>0\)，原小空隙定位先把每个 (K') 点送进一个外部 \(0.9Q_j\)。当 \(H_*v<1/10\) 时，该点仍在 \(Q(x_j,(1-H_*v)r_j)\)。返回指标不可能标记，否则该点不在 (K')。

共同缩半径后仍中心分离；其新3/4盒相遇的原半盒只有固定个尺度可比者。每个被缩去的部分体积 \(O(H_*v)r_i^n\)，所以新空隙为 \(O_{n,D,H_*}(v)\)。原私有定位证明适用于整个收缩族：唯一覆盖点距该中心至多 \(C(n,D,H_*)vr_i\)。这是一个关于任意先固定 \(H_*\) 的陈述，未给 N 增长的统一界。

已有 \(L_{h-1}\) 后，先取
\[
 H_h>L_{h-1}+D+1,\qquad
 G_{i,h}=K'\cap Q(x_i,(1-H_hv)r_i).
\]
首层取 \(T_{i,1}=P_i\)；其余层取 \(T_{i,h}=N_{i,h-1}\)，边界改回 \(P_i\)。删除原 \(r(x)<2L_{h-1}vr_i\) 的点得 \(F_{i,h}\)。局部输入
\(E\cap Q(x_i,3L_{h-1}vr_i)\) 的严格外部 halo 估计给删除量
\(C_{n,L_{h-1}}av^nr_i^n\)；输入零测时这个小半径集为空。原半径函数从未重选。

保留点的见证包含整个 \(N_{i,h-1}\supset P_i\)，因而 \(r(x)>(1-Dv)r_i\)。贪心捕获给外部 \(Q_j\ni x,r_j\ge r(x)\)；若 \(r_j>r_i\)，原线性半径引理给 \(r_j-r_i\le Ce_j^0r_j=O(vr_j)\)，故仍有统一 \(r(x)\le(1+C_{n,N}v)r_i\)。\(H_h>L_{h-1}+D+1\) 使整个 \(G_{i,h}\) 进入同一见证。

另一方面，对 \(x\in N_{i,h-1},j\ne i\)，
\[
 \|x-x_j\|_\infty>\max(r_i,r_j)-L_{h-1}vr_i
 \ge(1-L_{h-1}v)r_j.
\]
它大于第 j 个深目标的半径，故未删除源与所有异目标排他。边界 \(P_i\subset N_{i,h-1}\)，另有 \(r(x)\le r_i\)，两个方向均合法。

**此后**才取
\[
 L_h>\max\{L_{h-1}+1,C(n,D,H_h)\},\quad
 N_{i,h}=K'\cap Q(x_i,L_hvr_i).
\]
最终小 v 时中心邻域两两不交，所有未删除源也在自身深目标内。中心邻域之外每点至少双重覆盖，故 \(\mathcal D_h\ge f_h\)。边界改源的额外缺量
\[
 |N_{i,h-1}|-|P_i|\le |K\cap Q_i|-|P_i|\le as_ir_i^n
\]
求和至多 \(aP_\partial\)。删除质量只计入 \(\delta_h\)，不混进 \(d_h^0\)。由此 (6) 精确成立。

这一步只依以前的 \(L_{h-1}\)。有限 N 次后，取所有局部源常数的最大值为 \(A_N\)，再选强 CDF 斜率 \(\lambda_0\)，最后统一小 v。所有目标与原半径的严格条件仍在；零外部直接用原内部计数，不定义空 CDF。

**纯维数边界价格不反馈深度。** 每层边界有效源均在原 \(P_i\)，投影端点仍在原私有定位区间，且 \(R_{i,h}\le r_i\)。所以每层同样有
\[
 |X_{i,h}\setminus B_i|\le C_2(n)e_i^\beta r_i^n.
\]
因此图 (3) 的价格只依 n。深度只缩小最终阈值，没有图/I/CDF 选择循环。

## 4. 每层实际结构清单及深度依赖核查

以下几何构造在本族每层重新定义；全尺度清单／actual槽只在已有粗六阶先验的末层 \(\ell\ge6\) 调用。凡使用O8 raw、坏owner权O2或准确inside的数值段，明确限定 \(\ell\ge8\)；早五层只请求粗局部价，不请求这些精度。允许有限多个拟用层的阈值取最小。单独构造的深度 \(\ell\) 族，通常不是 N 族的前缀：其 CDF 可以不同。本定理只使用 N 族已经构造出的真实前缀，统一 \(A_N,\lambda_0\) 不变。

### G：原见证、原目标与面对几何

(5) 给同一原点见证同时捕获完整自身 \(G_i\subset E^c\) 及其他真实目标 \(F_j\)，后者彼此不交且与自身 \(G_i\) 排他。根比 \(R_i/r_i=1+O_N(v)\)；相交根可比，原半盒 packing 给固定相交度。

非临界轴捕获目标的整个**原始**投影宽，临界点簇同序。q≤1 的正体积交由原端—有效端双向尾斜率与完整自身目标外容量排除。真正面邻居从轴向5/8覆盖及固定切向裕量取得；第三盒不能塞进窄面缝。原半径对所有源点捕获邻居的完整横截面，因而
\[
 0\le g_{i\to p}\le2\Gamma_{ip}\le C_NvR_i,\qquad
 \sum_{i,p}R_i^{n-1}g_{i\to p}\le C_N\mathcal F_h.
\]
面逆度固定；不要求两膨胀 X 预先相交。每个推导仅用 (5)、固定源窗和原分离，未用 h=9。

### C：强 CDF、任意后选层、成本系数一

统一 \(A_N\) 后选择 \(\lambda_i=\lambda_0(vr_i)^{n-1}\)。连续坐标分布减 \(\lambda_is\) 的极小点给有效端，尾质量≤\(\lambda_i\Delta_i\)；有效盒各边≥\(c_{n,A_N}vr_i\)，中点壳体积≥\(c(vR_i)^{n-1}\Delta_i\)。先取小 \(\lambda_0\)，并支配**预定固定**原半径窗的端外延伸，故对任意后来选出的可测层 \(Y_i\subset D'_i=(x_i+W'_i)\setminus F_i\)，
\[
 |W_i|-\mu_i\ge |Y_i|+
 c_N[|D_i^0\setminus Y_i|+(vR_i)^{n-1}\Delta_i].
 \tag{12}
\]
精确相减的两项是 \(D'_i\setminus Y_i\) 与中点壳减源尾；原盒剩余只多端移控制的外环。

在完整质量区间，\(f_i(z)=(R_i+z^{1/n})^n\) 同时有
\((1-C_N\varepsilon_i)(1+v)^{n-1}v/a\) 和 \(c_Nv/a\) 的导数下界。先把积分长度分为 Y 与剩余，故层主系数一、raw 与端移正份额可同留；不许对另一价格再减全 J。以后选的 T、\(\mu_*\)、保证层不重选 CDF。包含整个正常模型的固定管只需最终 \(K_Nv^2\ll1\)，不是一次新的端外延伸要求。

### B：粗价、保留正量和正常模型的实际锚

原粗并集证明保留形状、重复量、面与稳健面；小半径删除和环境分别只用预留份额。每个拟用末层先在本族 (6) 上作粗权平均，得到六阶先验后，才重放本层 coarse joint。这一顺序在 §5 写清，不用尚未知的高精度价反供先验。

旧正常保证层 \(Y^{old}_{i,h}\subset D'_{i,h}\) 按本层有效源排序、原始目标帽及固定横向核心定义。高度≤\(Cv^4R_i\)，伴随短轴≤\(Cv^3R_i\)，逐源体积≤\(Cav^5R_i^n\)。在同一 root/clip/原始目标的 whole-min chart 上，只增加非负质量、不减旧侵蚀，精确
\[
 Y^{old}_{i,h}\subset B_{\mathcal N_h}(D_h^0)_i.
 \tag{13}
\]
全局最浅轴 q 不在 pair 的三近面图保留。

以下数值取 \(\ell\ge8\)。固定正常算子反单调，Tonelli 的 L1 范数≤\(C_Nv\)，pair 外部分≤\(C_Nv^3\)。本层正预算与 (13) 给 \(\rho_h=\|D_h^0-B_{\mathcal N_h}(D_h^0)\|_1=O_N(av^7S)\)，两次迭代的原误差分别为 \(O_N(av^7S)\)、\(O_N(av^8S)\)。每次迭代固定几何；不将正常模型价冒充减全 J 的价。

**共同侵蚀必须在整个正名义 base 上计费。** 对本层任意先固定的 bounded joint packet \(A\supset\mathcal N_h\)，记 \(\Gamma_i(A;\sigma,q)=\{y:h^0_A(y)>0\}\)，并在 normal/full 层使用同一个非负 \(e_i(T)\)。原实际两分支必要性给
\[
 e_i(T)=C_N[T^4v^8R_i+\Delta_i+
          \max_p(g_{i\to p}+\eta_p/\beta_p)],\qquad
 \sum_iR_i^{n-1}e_i(T)\le C_N(1+T^4)v^8S,
\]
其中 \(\eta_p=|D_p^0\setminus Y^{old}_{p,h}|\)、\(\beta_p\asymp(vR_p)^{n-1}\)。当最浅轴属于目标 pair 时，正名义 base 的测度至多 \(C_NTv^{n+1}R_i^{n-1}\)；不在 pair 时至多 \(C_NT^2v^{n+3}R_i^{n-1}\)。有限度与 chart 数固定，故先固定 \(T\ge1\)、最终 \(Tv^2\le c_0\) 时有真正的积分预算
\[
 \boxed{\sum_{i,\sigma,q}\int_{\Gamma_i(A;\sigma,q)} e_i(T)\,dy
                 \le C_NT^5av^9S.}
 \tag{13a}
\]
它由上述 base 测度乘 \(\sum R_i^{n-1}e_i\) 直接得到，先于任何真实 fibre 上截断。由于 \(\Gamma_i(\mathcal N_h)\subset\Gamma_i(A)\)，normal 与 full 的侵蚀损失**各自**不超过(13a)，即使 full 层已饱和而其 clipped loss 为零。不得由 full clipped loss 上界推 normal loss；这里只使用正名义 base 上的共同积分。因此 §6(23) 所需的 normal 侵蚀物理量 \(O_{N,T}(av^9S)\) 已有独立来源。该数值使用本层O8 raw，适用 \(h\ge8\)，不对早层作同阶声明。

### L：全尺度清单、根传输及实际槽

原 Q 图有限邻域增长只用原半盒不交和规则邻居尺度上界，故本 I 的扩大坏根权仍为 \(O(v^{9/2}S)\)。每个好根固定窗的**全部尺度**极限为单位铺砌：遗漏微盒/大盒先由原分离锁定尺度与中心，再与已覆盖单位胞的距离≤1/2矛盾。匹配只需最终固定精度；不把未知收敛模量当作幂率。

本层粗 Hölder 给 pair gross \(O_N(v^6S)\)。原并集与单盒基线给 void
\(O_N((v^{n+1}+v^5)S)\)。固定模型窗内随机面路径选取有切向密度控制，沿真实稳健面对精确求和根差及端移，固定反向重数给
\[
 \sum_{i\in G_v}\sum_{j\in\mathcal I_i(R)}
 |R_{i,h}-R_{j,h}|r_i^{n-1}=O_{N,R}(v^6S).
 \tag{14}
\]
固定 probe/清单20的 actual helper 按本层端点定义。错误列有固定大碰撞；正确列唯一；同槽相反对角端点恒等式触发该边自己的 \(\min L/(1000r_i)\) 碰撞，不要求匹配误差小于 minL。横片先取并再作实际槽单射，故
\[
 W_{\mathcal N_h}=\mathcal C_h-D_{\mathcal N_h}+O_N(v^3S),\quad
 D_{\mathcal N_h}\ge0,\quad
 \mathcal C_h=\sum_i\sum_{l<k}R_{i,h}^2\prod_{j\ne l,k}(R_{i,h}+w_{ij,h}).
 \tag{15}
\]
软权 \(\omega=\min(1,C_s\min(p,q))\) 的唯一增长比是 \(\omega/\lambda_e\le C_NC_s/v^3\)，从而同一分配同时给 \(W_{C_s}\le\mathcal C_h+C_N(1+C_s)v^3S\)；clean 身份阈值先于 \(C_s,T,\mu_*\)。二维用真实四槽与面对序，不借高维坏根结论。

### Q：边界、任意薄长角和 inside

边界的原始突出估计在每层成立。以下O7/O9数值取 \(\ell\ge8\)：边界收费使用正向事实 \(\varepsilon_i\mu_i=|W_i|-\mu_i\le|D_i^0|\) 及 \(\mu_i\asymp aR_i^n\)，不倒置一个粗上界。旧层体积及 raw 分解给本层邻边界 epsilon 和为 O7，边界 q2 为 O9。非临界目标原整宽捕获、实际长角截宽及 \(C(v^3/P+v^2)\) 正常帽比例都一致于任意小短轴 Q。

全部 q 的真实源矩形使用同一层原见证及共同度数。固定管后第2至q临界深带远离两面，因此至多第一短轴近面，能对同一管外 Z 分一个固定小份额。二维高支撑为空；n≥3 才定义坏所有者，其权由本层槽欠量控制。

局部小源始终是原有限值半径的整个点集：盒内用 dyadic 含 a 的回收系数，盒外由本层真正面邻控制；零输入为空。二维准确角窗中，近角定向贪心块不交，远角 residual 只计一次，因此在 \(\ell\ge8\) 的O8预算下
\[
 (1+O_N(v))(v/a)\delta_{\ell,\rm in}
 \le(8/27+O_N(v))v^8S\quad(n=2).
 \tag{16}
\]
边界和环境十阶、盒外九阶；n≥3 的绝对盒内／环境项从 \(v^{n+6}S\) 起。此为当前精度，**不是任意阶精度**。全区间导数为 \((1+O_N(v))v/a\)，删除与未删除 \(d_h^0\) 分段；不把固定 b 错乘到 inside 主费用上。

## 5. 同一个 N 族的前缀 bootstrap 与宽度

固定一个 \(\ell\in\{6,\ldots,N\}\)。从第一层至第 \(\ell\) 层，几何对象、CDF 和 I 都不换；只更改加权长度。

原保留粗预算在每层给同一左侧的
\[
 Z_0\ge c_0 f_h-b_0vf_{h-1}+a\mathcal J_h+\mathcal D_h/4
       -C_Na(v^6S+vP_\partial),
 \tag{17}
\]
其中 \(\mathcal J_h\) 包含固定正比例的
\(v\sum_i(\varepsilon_{i,h}+q_{i,h})r_i^n\)、面与稳健面。用
\((b_0v/c_0)^{\ell-h}\) 加权，中间项精确消去；初端除 a 为
\(O_N(v^{\ell+1}S)\)。C8 蕴含六阶临界，真实比值反向界与粗全输入上界给
\(\mathcal T-A_vM/a=O_N(v^6S)\) 双侧界。故在本层得到六阶正先验、
\(\mathcal D_{\ell-1}/a+d_\ell^0/a+S_I+S(\partial I)=O_N(v^5S)\)，以及
\(\delta_\ell/a=O_N(v^6S)\)。这先供给 §4 的清单、面路径与 coarse joint。

每个 \(\ell\ge7\) 的 local joint 只用这些先验：正常退款由固定份额的原平方罚与 Cv8 支付；窗外完整强制性保名义成本和角折扣；长条先取固定窗、高支撑只用固定小份额。面条只恢复一次，普通与稳健面并集取两个下界的半和。于是 (7) 在本层成立。

对每个 \(\ell\ge8\)，在第 \(\ell\) 层重开同一次证明保留其所有正份额。以本前缀权平均，末层权一，第六层领先误差为
\(v^{\ell-6}o(v^6S)=o(v^\ell S)\)，早五层及初端为 O\((v^{\ell+1}S)\)。因 \(\ell\ge8\)，均属于 O8；前面粗中心 Q 的误差也属于 O8。C8 与
\[
 (1+Q)\mathcal T-A_vM/a
 =\frac{A_v-a(1+Q)}aV
  \left[\frac{A_v}{A_v-a(1+Q)}-(1+M/V)\right]
 \tag{18}
\]
将这一个实际公共左侧压为 O8。因而取得本层 (8)，不借精确八阶上界。

为得完整宽度，实际恒等式
\[
 |X_i|=(R_i+\mu_i^{1/n})^n+J_i+\mathcal A_i,\qquad
 \mathcal A_i\ge c_Nv q_iR_i^n
\]
的形状量独立于 J。把公共重复量系数1/2拆成两份1/4，只用一份消项，另一份在前缀末端保留，得到
\[
 v\sum_iq_{i,\ell}R_{i,\ell}^n+\mathcal D_\ell/a
       +v\mathcal D_{\ell-1}/a=O_N(v^8S).
 \tag{19}
\]
原 dyadic/盒外/边界删除式及 (8) 给 \(\delta_\ell/a=O_N(v^7S)\)；由 (6)，源质量缺陷 \(\eta_i=1-\mu_i/g_i\) 的加权和为 O7。原 raw 与逐源旧层体积给 \(\sum\varepsilon_i^2R_i^n=O_N(v^7S)\)。端移比 \(\rho_i=\Delta_i/(vR_i)\) 一致有界且其和为 \(\mathfrak E_\ell/v=O_N(v^7S)\)，故平方和同阶。三盒宽度满足
\[
 \xi_i^2\le C_N(q_i+\varepsilon_i^2+\eta_i^2+\rho_i^2),
\]
得 (9)。这保留全部异常源；只证明 q 并不足以控制同一个完整 \(\xi\)。

## 6. C9 下各真实前缀的九阶正剩余

对任一固定 \(\ell\ge9\)，§5 在本层已给 (8)、(9)。正常峰值积分的非宽度误差是 O9，宽度费为 \(v^6\sum\xi_iR_i^n=O_N(v^{19/2}S)\)。固定小正常窗内的真正最大点移至 \(r_v=1/3+(20/81)v+O(v^2)\)；仅在这个标量最大点 Taylor，保留移动平方罚 \(\mathcal L_{\mathcal N_\ell}\)。此阶段不含未知九阶系数。

选择包含整个 \(B_{\mathcal N_\ell}(0)\) 的固定管 \(K_N\)：其两个底临界参数≤\(C_Nv^3R_i\)，高度≤\(C_Nv^4R_i\)，包括最浅轴不在 pair 的图。对同一侵蚀的联合层，有**常数一**关系
\[
 Z_\ell(K_N)\le\mathcal R^-_{A,\ell}+s_v^{-1}\mathsf M_{B,\ell},
 \qquad s_v=(1+v)^{n-1}.
 \tag{20}
\]
coercivity、CDF raw 与槽的固定正系数都先于 T、\(\mu_*\)。完整名义成本的退款/clip保为小系数乘 \(\mathsf M_B\)，而不先粗估为O8；长条在整个固定管外，代价≤\((C_L/T)Z\)；近零 gross≤\(CT\mu_*v^6D_{\mathcal N}\)；高支撑只用同一 Z 的固定小 \(\eta\) 份额。

次序是固定正常窗/CDF/\(K_N\)/正系数，再固定大 T、小 \(\mu_*\)、小 \(\eta\)，最后取小 v。共同 J 只用一次，得到本层完整
\[
 H_\ell+c_B\mathsf M_{B,\ell}+c_R\mathcal R^-_{A,\ell}
       +c_Dv^6D_{\mathcal N_\ell}+c\mathcal L_{\mathcal N_\ell}
 \le q_N\mathcal C_\ell+C_Nv^9S.
 \tag{21}
\]
目标原／有效壳另付 \(Cv\mathfrak E\)，不是误记入 Z。所有 q 共同度数只分一份管外量。当前 inside (16) 与准确真实导数给 (10)，并同时保留 (21) 的正项。

在本族长度 \(\ell\) 前缀使用共同中心 \(Q_8\)：h≥9用(10)，h=8,7用(7)，h=6用本层领先价，早五粗。它们的加权误差分别为
\[
 O_N(v^9S),\quad O_N(v^{\ell}S),\quad O_N(v^{\ell+1}S),\quad
 o(v^\ell S),\quad O_N(v^{\ell+1}S),
\]
初端也是 O\((v^{\ell+1}S)\)。所以都属于O9。注意 h≥9 的局部式(21)是从 (8)、(9) 推出的，**没有先用本段的全局结算**；同时重放有限多个末层不构成归纳循环。

令 \(D_v=A_v-a(1+Q_8)>c'>0\)、\(R=1+M/V\)。真实 M 的反向式准确为
\[
 (1+Q_8)\mathcal T-A_vM/a
   =(D_v/a)V[A_v/D_v-R].
 \tag{18'}
\]
并且
\[
 \frac{A_v}{D_v}-L_n
  =\frac{aA_vQ_8}{(A_v-a)(A_v-a-aQ_8)}.
\]
将已显示的 \(Q_8\) 代入，以 \(a=v^2w,0\le w\le1\) 在统一正分母上作二阶 Taylor，再代 \(w=v^{n-2}\)，得
\(A_v/D_v=P_{8,n}+O_n(av^9)\)。二维分母的 \(8/27\) 与 inside 的 \(8/27\) 各保留一次；这是明示有理式的有限展开，不使用 \(C_n\) 上界。
C9 及 (18') 因而使公共左侧至多 \(C_{N,B}v^9V=O_{N,B}(v^9S)\)，最后一步用竞争分支的 \(V\asymp S\)。遂保留本末层的 (11)，再由 (20) 得管外 Z。
全程使用同一个原 M，不能换成 \(M'\)。末层权除以 \(\Omega_\ell=1+O_N(v)\) 仍有固定正系数。

**给后续正常反馈的直接推论。** 固定同一 packet 的共同侵蚀后，
\[
 \mathcal R^-_{\mathcal N,\ell}
 \le\mathcal R^-_{A,\ell}+s_v^{-1}\mathsf M_{B,\ell}=O_N(v^9S).
 \tag{22}
\]
固定正常必要性侵蚀的物理损失 \(O_N(av^9S)\) 更晚，因此
\[
 \|D_\ell^0-B_{\mathcal N_\ell}(D_\ell^0)\|_1=O_N(av^8S).
 \tag{23}
\]
这是**新 whole-normal raw** 的结论，不是改善旧核心 \(R^0_{N,\ell}\)。它与任意Y的(12)、固定算子及全部实际 chart 同对象，可供第三次反馈；本文不预先签任何新的九阶积分系数。

## 7. 精确共同账本和一般误差传播

每层使用同一次内部回收，公共左侧输入只计一次。不同层的正量通过权平均保留，不把各层体积当成互不相交的新增收入。写
\[
 B(v)=A_vM/a+c_{\rm rec}P_I-C_bP_\partial-2nvS_I.
\]
存在固定 \(c,b>0\)，可按需要把重复量再拆份，使所有层同时有
\[
 \mathcal T\ge B(v)+c\mathcal D_h/a-bv d_h^0/a-E_h+\mathcal P_h,
 \qquad\mathcal P_h\ge0.
 \tag{24}
\]
例如不另保重复量取 \(c=1/2\)；若总系数1/2拆成两份，则消项取 \(c=1/4\)，另一份实际 \(\mathcal D_h/(4a)\) 放入 \(\mathcal P_h\)。形状也放在这同一个正余量内，不能既被忽略又在另一价格中重复加入。

设
\[
 \tau=bv/c,\quad w_h=\tau^{N-h},\quad
 \Omega_N=\sum_{h=1}^Nw_h.
\]
精确恒等式为
\[
 \sum_{h=1}^Nw_h(cf_h-bvf_{h-1})
      =cf_N-bv\tau^{N-1}f_0.
 \tag{25}
\]
由 \(d_h^0\le f_{h-1}+aP_\partial\) 和 \(\mathcal D_h\ge f_h\)，可得
\[
 \mathcal T\ge B(v)-\Omega_N^{-1}\sum_hw_hE_h
 -\frac{bv\tau^{N-1}f_0}{a\Omega_N}-bvP_\partial
 +\Omega_N^{-1}\sum_hw_h\mathcal P_h.
 \tag{26}
\]
未写出的 \(cf_N/(a\Omega_N)\) 与 \(c(\mathcal D_h-f_h)/a\) 均非负。首层并无边界误差，故 \(-bvP_\partial\) 只是安全上包。内部必须保留精确式
\[
 (1-a)S_I+M_I-A_vM_I/a
 =(1-a)(S_I-M_I/a)-(A_v-1)M_I/a.
\]
图(3)和事先预留吸收全部边界／内部项，不以 \(aS=M\) 代换。

对**任意共同中心** q(v)，准确中心化是
\[
 \Omega_N^{-1}\sum_hw_hE_h-q\mathcal T
 =\Omega_N^{-1}\sum_hw_h(E_h-q\mathcal T).
 \tag{27}
\]
若各本层实际上包满足
\(E_h\le q\mathcal T+\eta_h(v)S\)，则
\[
 (1+q)\mathcal T\ge A_vM/a-
 S\Omega_N^{-1}\sum_h\tau^{N-h}\eta_h(v)-C_Nv^{N+1}S
 +\Omega_N^{-1}\sum_hw_h\mathcal P_h.
 \tag{28}
\]
初端来自 \(f_0=A_0vM'\)：
\[
 \frac{bv\tau^{N-1}f_0}{a\Omega_N}
 =\frac{A_0b^N}{c^{N-1}\Omega_N}v^{N+1}(M'/a)
 =O_N(v^{N+1}S).
 \tag{29}
\]
给定 \(E_h\) 的误差有符号时，(27)按原符号保留；用 \(|\eta_h|\) 只是安全界。

一般地，\(\eta_h=O(v^{s_h})\) 的贡献为
\(O_N(v^{s_h+N-h})\)，little-o 同样准确平移；固定有限和不要求各层之间有同一个收敛模量。真正全输入结论仍用失败序列选择固定精度工具，只把起效阈值留给该序列。

最后用原 \(S/V\le(1-a)^{-1}\)，对小v的统一正分母直接移实际 M。若 q 是实际共同价且分母 \(A_v-a(1+q)>c'>0\)，(28)给归一化价误差
\[
 O_N\!\left(a\sum_hv^{N-h}|\eta_h(v)|+av^{N+1}\right),
 \tag{30}
\]
或相应 little-o。此处 q 的有理分母到所需目标多项式的截断误差必须另外实际证明；(30)不是未知系数的生成器。

## 8. N=10 的明确接口及任意固定阶的边界

设 R075 或以后局部工作在本族第10层供应某个共同中心 \(q_9(v)=Q_8(v)+O(v^9)\) 的末层完整误差 \(\eta_{10}\)，其中**已经包含 inside、环境、边界和正常补集**。则同一 N10 中可用：

| 原本层 | 已供精度 | 权 | 对共同中心的贡献 |
|---|---:|---:|---:|
| 10 | \(\eta_{10}\) | 1 | \(\eta_{10}\) |
| 9 | O9 | O(v) | O10 |
| 8 | O8 | O(v²) | O10 |
| 7 | O8 | O(v³) | O11 |
| 6 | \(o(v^6)+O(v^7)\) | O(v⁴) | o10 + O11 |
| 1…5 | O6 | 至少O(v⁵) | O11 |
| 初端 | — | — | O11 |

第六层仅定性领先价现在严格晚于九阶；没有把它的未知模量改成一个更快幂率。如果末层完整 \(\eta_{10}=o(v^9)\)，则账本误差为 o9；若末层供 O\((v^{9+\alpha})\)、\(0<\alpha\le1\)，账本供同率（\(\alpha=1\) 时包括 o10 项）。**当前 (16) 自身只有 O9，不自动满足这一末层新条件。** R075 仅正常九阶积分即使闭合，也不能省略该条件。

对一般固定目标代价 \(v^r\)，选择 N>r 可把第六层及旧7/8层、早端移到更晚；但仍需最后有限个近末前缀满足
\[
 \eta_h(v)=o(v^{r-(N-h)})
 \quad\text{或严格更高的已证 O 幂率}.
 \tag{31}
\]
只让 N 增大、却给倒数第二层一直保留 O9，会留下 O10 地板；它不推出任意固定阶。所需高精度前缀数有限，全部必须在同一个 N 族本层实供。N 不得随 v 增长，固定窗口的未知紧性模量也不因此获得幂率。

## 9. 实际依据和验证边界

本轮实际逐段读了正式正文的 global/sharp/eighth/quantitative，特别是固定深度收缩定位、边界纯维数C2、共同质量式、coarse bootstrap、全点 inside、固定整管/精确Z及正份额吸收；并回读 R064 新N9供货和 R071 全证。引用这些展开的局部证明只因它们的假设已在上文逐层重新满足，不以历史审裁决为前提。

本轮已实际运行 `finite_depth_fraction.py`，原 stdout 保存在 `finite_depth_fraction.stdout`。标准库 Fraction 在 N=9、10、13 的三个有理实例上核得望远镜、共同中心化和精确初端残差全为零；保留末两层重复量的不等式有严格正余量。N10 阶表实际返回 `{10,10,11,10,11,11}`，其中第6层的10只表示 little-o 的平移指数，不将它当作有率 O 的新常数。没有补造时间或 PID。

另实际读取并复用既有 `r064-computation/nine_layers.wl/log`：Wolfram 15 的一般N重索引内项残差为0，精确初端目标式残差为0。历史日志还含65/8专用量，本文不把那些旧值当作N10或本稿的新证据，也未无意义重跑。

`source_and_computation_manifest.json` 记录本轮实际读取的正文、两份供货主源与全部上述脚本／日志的精确路径、字节及SHA256。有限有理样本不证明任意N；一般恒等式由(25)的内项系数 `c*tau-b*v=0` 解析消去。计算更不认证 G/C/B/L/Q 的几何，这些由 §§3–6 的全点与实际集合论证承担。

本稿是可复用供货定理与误差接口整理。它实际新增了对任何固定N及所有本族前缀的量词闭合、九阶正剩余接口和一般中心化传播规则；不把旧算法重写称作新的渐近突破。
