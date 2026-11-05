# R052：全局最浅面坐标中的实际帽反单调迭代

2026-11-05。固定 \(n\ge2\)。条件性局部引理，供独立审核。仓库、论文及前轮冻结源只读。

本稿的新增结论是：在原实际盒、实际目标帽和原源缺集上构造一个全族 \(L^1\) Lipschitz 常数为 \(O(v)\) 的反单调算子。若 R051 交付下述确切的同见证高度必要性，则 R049 真实余量使第一次迭代足以计算八阶退款，第二次迭代足以计算八阶强制源成本，误差均为 \(O(v^9S)\)。这里不识别八阶系数，不把安全构造的 \(B_2\) 条件当作原极大算子的必要条件。

## 1. 同层输入与唯一新增几何条件

在 R049 的同一新固定八层族、同一固定正常集合 \(\mathcal N\) 上工作。所有实际源 \(F_i\subset G_i\)、质量根 \(R_i\)、原指定严格见证、原始盒 \(W_i^0\)、有效盒 \(W_i'\)、中点盒、原始目标缺集
\[
 \mathscr D_i=(x_i+W_i^0)\setminus F_i
\]
均不改变。记 \(a=v^n\)、\(S=\sum_i r_i^n\)。输出层在**有效源盒**中，目标缺量始终在**原始目标盒**中；不可将两种盒或两种坐标静默交换。

已审 R049 的实际稳定性给
\[
 \sum_i|\mathscr D_i|\le C av^5S,\qquad
 r:=\sum_i|\mathscr D_i\setminus Y_{i,N}|\le C av^7S,
\tag{1}
\]
其中 \(Y_{i,N}\) 是 R045/R049 共用同一侵蚀和核心的正常实际层。第一式也直接来自原粗 epsilon 先验。强 CDF、盒宽 \(\asymp vR_i\)、实际相交度有界及相交根比仍保持。

可只取正常 packet \(\mathscr P=\mathcal N\)。更一般地，也允许任意一个**预先固定、两端均深处**的有界二轴 packet \(\mathscr P\supset\mathcal N\)，其中
\[
 0<L_{ij,1},L_{ij,2}\le T v^3R_e,\qquad R_e=\max(R_i,R_j),\qquad T<\infty.
\tag{2}
\]
它可以含多个同源同符号目标；不要求这些目标的 source 禁层分开，也不要求交长有统一正下界。每个无向 pair 在两个有向端点各提供实际帽，仍用原有限度。不得令 \(T=T(v)\) 增长。

对有向 \(i\to j\)，临界轴为 \(l,m\)，朝向符号为 \(\tau_{ij,l},\tau_{ij,m}\)。定义与 R045 完全相同的中点帽参数
\[
 c_{ij}=\operatorname{clip}\left((R_i-R_j)/2;\,
             -\max(L_{ij,1},L_{ij,2}),\max(L_{ij,1},L_{ij,2})\right),
\]
\[
 P_{ij,l}=(L_{ij,1}+c_{ij})_+,\qquad
 P_{ij,m}=(L_{ij,2}+c_{ij})_+.
\tag{3}
\]
故 \(0\le P_{ij,h}\le C_Tv^3R_i\)。正常边根差未截断；有界补集不预设未截断。目标帽总取原始面对端、其余轴原始整宽。

**唯一新增几何条件 \(\mathrm{Nec}_{51}\)。** 对下面第 2–3 节定义的模型高度，存在非负实际误差函数 \(\epsilon_{i,q,\sigma}(y)\)，使同一个原指定严格见证推出
\[
 x=\Phi_{i,q,\sigma}(s,y)\in F_i
 \quad\Longrightarrow\quad
 s\ge h^0_{i,q,\sigma}(y)
       -\mathcal K_{i,q,\sigma}(\mathscr D)(y)
       -\epsilon_{i,q,\sigma}(y)
\tag{N1}
\]
（只在模型高度可能为正且 \(s\le C_Tv^4R_i\) 的域需要证明；更高的 \(s\) 自动满足此下界），并有
\[
 e:=\sum_{i,q,\sigma}
       \int_{\{h^0_{i,q,\sigma}>0\}}
                 \epsilon_{i,q,\sigma}(y)\,dy
       \le C_T av^9S.
\tag{N2}
\]
这必须来自所有 packet 目标真实质量在**同一原见证**内先相加、完整自身目标排他及真正最浅面邻居的实际轴帽。它不能由周期安全测试、理想单胞算子、只对某个指定角点的估计或三面核心剪裁代替。

R051 冻结源 scratch/r051-refined-axial-cap/refined_axial_cap.md，SHA256 **37bbb0ec63fac1dce33fb9883e14499488baef490341df32a6fdfb0da427c979**，以正常层的轴帽占比 \(O(v^2)\)、只扣真实 raw-neighbor 余量的机制供给 (N1)–(N2)：其 (22)、(24)–(27) 给所需的实际高度和物理侵蚀。已实际读取该冻结源，但其独审尚未完成，本文保持此条件。第 4 节的旧层包含也与其 (30) 逐坐标相符；不额外假设原缺集已经等于某个模型层。

## 2. 精确的整个有效源盒坐标

写有效源盒第 \(h\) 轴端点为 \(a'_{ih},b'_{ih}\)，全宽 \(w'_{ih}\)。忽略有限个最近面／最浅轴平局超平面。对每点先取每轴的最近面符号
\(\sigma\in\{-1,+1\}^n\)，再取全局最近的轴 \(q\)。令
\[
 f_{ih}^{+}=b'_{ih},\quad f_{ih}^{-}=a'_{ih},\qquad
 y_q=0,\quad
 x_h=f_{ih}^{\sigma_h}-\sigma_h(s+y_h).
\tag{4}
\]
其中 \(s\) 是所有 \(2n\) 有效面中的最小距离，\(y_h\ge0\)。精确 base 与纤维高度为
\[
 \Omega_{i,q,\sigma}=\prod_{h\ne q}[0,w'_{ih}/2],\qquad
 H_{i,q}(y)=\min\{w'_{iq}/2,\ \min_{h\ne q}(w'_{ih}/2-y_h)\}.
\tag{5}
\]
\(0\le s\le H_{i,q}(y)\)，映射 \(\Phi_{i,q,\sigma}\) 的绝对 Jacobian 为一。全部 \(q,\sigma\) 图覆盖整个有效盒至零测边界；不存在被删除的三近面区域。

在 \(s=O_T(v^4R_i)\) 时，一个 normal／有界 packet 目标只有当其两临界朝向与 \(\sigma\) 对应符号相同时才提供本模型帽。相反面距至少为半个源宽 \(\asymp vR_i\)，大于 \(O_T(v^3R_i)\) 的标称帽。这个选择不删掉正模型帽。

若 \(q\) 不在某个目标的两临界轴内，该目标仍完整参与同一图上的高度，不能仅保留“最浅轴属于本 pair”的两个排序。正是这些图保留了三近面八阶源成本。

## 3. 实际目标帽、完整高度与反单调算子

对图 \((i,q,\sigma)\) 中方向匹配的 packet 目标 \(j\)，定义原始目标盒内的帽 \(\mathcal P_{ij}(y)\)：第 \(l,m\) 轴从其朝源原始面向内的深度分别为
\[
 \min\{w^0_{jl},(P_{ij,l}-y_l)_+\},\qquad
 \min\{w^0_{jm},(P_{ij,m}-y_m)_+\},
\tag{6}
\]
其余轴取原始目标整宽。实际存在的截宽始终保留；固定 \(T\)、小 \(v\) 后这些上截断自动不活跃，但证明不需要先删去它们。

定义高度量
\[
 h^0_{i,q,\sigma}(y)=
   \frac1{naR_i^{n-1}}\sum_{j\in\mathscr P(i,\sigma)}
                                |\mathcal P_{ij}(y)|,
\]
\[
 \mathcal K_{i,q,\sigma}(E)(y)=
   \frac1{naR_i^{n-1}}\sum_{j\in\mathscr P(i,\sigma)}
                                |E_j\cap\mathcal P_{ij}(y)|.
\tag{7}
\]
输入 \(E=(E_j)\) 是各原始目标盒中的任意可测集合族。改变测试 \(E\) 时不重算实际 \(G_i,R_i,c_{ij}\)、盒宽、目标清单或 chart；这些全部冻结。各目标质量先相加，再取一个高度正部。实际 \(E=\mathscr D\) 时，
\[
 h^0-\mathcal K(\mathscr D)
   =\frac1{naR_i^{n-1}}\sum_j|F_j\cap\mathcal P_{ij}(y)|\ge0.
\tag{8}
\]
尚未以 (8) 断言源点必须缺失；那一必要性恰是 (N1) 所列真实见证条件。

令
\[
 B(E)_i=
  \bigcup_{q,\sigma}
   \{\Phi_{i,q,\sigma}(s,y):
       0\le s\le
       \min(H_{i,q}(y),[h^0_{i,q,\sigma}(y)
                       -\mathcal K_{i,q,\sigma}(E)(y)]_+)\}.
\tag{9}
\]
各图几乎处处分开，(9) 是真正的物理集合，并映回原始盒中的有效子集。有限 Borel／Lebesgue 可测集及 Tonelli 足以保证可测性，无新的可测选择。

逐点高度至多 \(C_Tv^4R_i\)，因帽体积至多
\(C_T(vR_i)^{n-2}(v^3R_i)^2\)。同时
\[
 E\subset F\quad\Longrightarrow\quad B(F)\subset B(E).
\tag{10}
\]
正部及真实纤维上截断均保留，零高度附近不要求相对误差或统一正下界。

## 4. 旧正常实际层的精确包含

取 R045 的同一 \(Y_{i,N}\)，保持其原有效源坐标、共同 \(E_i(T_{45})\) 和核心。其一个排序的点深为 \(s,s+y\)，其他坐标到有效盒两面均大于 \(H_{T_{45}}v^3R_i\)。该层高度为 \(O(v^4R_i)\)，所以对小 \(v\)，全局最浅轴恰为原 pair 的浅轴；两参与轴的最近面符号恰为原有向角，其他轴符号随点确定。

因此 (4) 把这个点送入原来的同一排序，原实际帽
\(P_{ij}(Q_{ij}-y)_+\) 及目标 raw 缺量积分完全相同。新高度 (7) 包含全部原正常目标，其余目标帽的实际 \(F_j\) 质量非负；旧侵蚀 \(E_i(T_{45})\) 也非负。故
\[
 \boxed{Y_{i,N}\subset B(\mathscr D)_i\subset B(0)_i.}
\tag{11}
\]
没有把图改成原始源盒后省略端移：**源图一直是有效盒**，target 缺集一直是原始盒。(11) 依赖 \(\mathscr P\supset\mathcal N\)，不要求旧正常层与新增有界层互不相交。

由 (N1) 定义实际保证层 \(B^-(\mathscr D)\)：在 (9) 高度中再减
\(\epsilon_{i,q,\sigma}(y)\) 后取同一正部及上截断。则
\[
 B^-(\mathscr D)_i\subset D'_i,\qquad
 \sum_i|B(\mathscr D)_i\setminus B^-(\mathscr D)_i|\le e.
\tag{12}
\]
所以实际缺集满足
\[
 \|\mathscr D-B(\mathscr D)\|_1\le r+e,\qquad
 \sum_i|\mathscr D_i\setminus B^-(\mathscr D)_i|\le r+e.
\tag{13}
\]
这里 \(\|E-F\|_1=\sum_i|E_i\triangle F_i|\)。这是由原实际必要性和实际旧层余量导出的近固定点；没有假设实际缺集属于一个安全测试类。

## 5. 全族实际 \(L^1\) 小算子常数

存在 \(C_T<\infty\)，只依赖固定窗口及原源常数，使
\[
 \boxed{\|B(E)-B(F)\|_1\le\lambda\|E-F\|_1,\qquad
                         \lambda=C_Tv.}
\tag{14}
\]
证明先用
\[
 |\min(H,[u]_+)-\min(H,[w]_+)|\le|u-w|.
\]
对每张图的 \(y\) 积分，目标差只涉及
\(|(E_j\triangle F_j)\cap\mathcal P_{ij}(y)|\)。打开 Tonelli，固定目标物理点 \(z\)，估计可把它包含在帽内的源 base 测度。

若 \(q\in\{l,m\}\)，其中一轴 \(y_q=0\)，另一轴的有效 \(y\) 范围至多 \(C_Tv^3R_i\)；其余 \(n-2\) 轴范围至多 \(C(vR_i)\)。除以 (7) 的 \(naR_i^{n-1}\)，得
\[
 \frac{C_T(v^3R_i)(vR_i)^{n-2}}
            {aR_i^{n-1}}\le C_Tv.
\tag{15}
\]
若 \(q\notin\{l,m\}\)，两轴范围均至多 \(C_Tv^3R_i\)，其余只有 \(n-3\) 轴：
\[
 \frac{C_T(v^3R_i)^2(vR_i)^{n-3}}
            {aR_i^{n-1}}\le C_Tv^3.
\tag{16}
\]
二维第二种情形为空。源／目标实际根仅使用原相交尺度可比；每个 \(R_i\) 的幂在自己的分母中抵消，没有将不同实际根直接设为相等。全族每个目标被有界个相关源／图重复，故 (15)–(16) 相加即得 (14)。

这是一条对**任意输入缺集**的全族物理 \(L^1\) 界，不只对逐边有界的一阶场成立。正部、源纤维上截断及任意薄目标帽都已包含在其常数中。

另令 \(Y_0=B(0)\)。直接积分名义完整帽给
\[
 |Y_0|_1\le C_Tav^5S.
\tag{17}
\]
逐目标，在 \(q\in\{l,m\}\) 图的量级为
\((vR_i)^{n-2}(vR_j)^{n-2}(v^3R_i)^3/(aR_i^{n-1})
 =O_T(av^5R_i^n)\)；
在 \(q\notin\{l,m\}\) 图为 \(O_T(av^7R_i^n)\)。后者的源成本是八阶，**此处只是给大小上界，没有在模型中删去它**。真实源纤维上截断只减少 (17)。

## 6. 两次迭代与实际缺集的误差

定义
\[
 Y_0=B(0),\qquad Y_1=B(Y_0),\qquad Y_2=B(Y_1).
\tag{18}
\]
反单调性给
\[
 Y_1\subset Y_2\subset Y_0,\qquad
 \|Y_2-Y_1\|_1
 \le\lambda^2|Y_0|_1=O_T(av^7S).
\tag{19}
\]
第二式先用 \(\|Y_0-Y_1\|_1\le\lambda|Y_0|_1\)，再用一次 (14)。

记 \(M_0=\|\mathscr D\|_1\le Cav^5S\)、\(\rho=r+e=O_T(av^7S)\)。由实际近固定点 (13) 和小算子界，
\[
 \|\mathscr D-Y_0\|_1
 \le \rho+\lambda M_0=O_T(av^6S),
\tag{20}
\]
\[
 \boxed{\|\mathscr D-Y_1\|_1
 \le \rho+\lambda(\rho+\lambda M_0)
       =O_T(av^7S).}
\tag{21}
\]
因此
\[
 \boxed{\|B(\mathscr D)-Y_2\|_1
       \le\lambda\|\mathscr D-Y_1\|_1
       =O_T(av^8S).}
\tag{22}
\]
所有等号后的常数可依固定 \(T\)，不依原有限族数目、尺度跨度或 rescaled 一阶参数场的逐点上界。尚未对这些集合或全源剖面作 Taylor 展开。

## 7. 八阶退款与强制源成本的准确替换

任取实际源体积权 \(0\le\ell_i\le C v/a\)，包括完整鲜导数
\(\ell_i=(1+v)^{n-1}v/a\)。由 (22)，
\[
 \left|\sum_i\ell_i|B(\mathscr D)_i|
           -\sum_i\ell_i|Y_{2,i}|\right|
       \le C_Tv^9S=o(v^8S).
\tag{23}
\]
这是强制模型层的**第二次迭代**替换；不能用
\(\|\mathscr D-Y_2\|_1=O(av^7S)\) 直接丢弃整个实际缺集，后者乘源导数仍可能贡献八阶正余量。

若将模型退款写成保留真实纤维截断的非线性量
\[
 \Phi_\ell(E)=\sum_i\ell_i[|Y_{0,i}|-|B(E)_i|],
\]
则 (14) 给
\[
 |\Phi_\ell(E)-\Phi_\ell(F)|
       \le C_T(v^2/a)\|E-F\|_1.
\]
从 (21) 得
\[
 \boxed{|\Phi_\ell(\mathscr D)-\Phi_\ell(Y_1)|
                  \le C_Tv^9S=o(v^8S).}
\tag{24}
\]
同样适用于 R042/R045 的任何实际帽 Tonelli 退款核，其全族核范数为 \(C_Tv^2/a\)。因此退款只需真实第一次迭代 \(Y_1\)，不再把目标正常层缺量整体扣成一个无结构 epsilon。

在 (N1)–(N2) 已供应的条件下，强 CDF 可对真实保证层
\(B^-(\mathscr D)\subset D'\) 使用一次。原鲜导数损失至多
\[
 C\frac va\sum_i\varepsilon_i|Y_{0,i}|
 \le C_Tv^6\sum_i\varepsilon_iR_i^n=O_T(v^{11}S).
\]
保证层与模型层之差的源成本为 \(O((v/a)e)=O_T(v^{10}S)\)。由 (23) 可以得到实际源成本下界
\[
 \boxed{\displaystyle
 \sum_iJ_i\ge (1+v)^{n-1}\frac va\sum_i|Y_{2,i}|
       +c\frac va\sum_i|\mathscr D_i\setminus B^-(\mathscr D)_i|
       +c\mathfrak E-O_T(v^9S).}
\tag{25}
\]
这是同一次完整 CDF 结算，raw 与端移由原强 CDF 同时分配。packet 若只取正常集合，本式只给正常模型层和未用真实余量；不能与另一条减全体 \(J\) 的补集价格直接相加。若要同全有界二轴结算，必须在最初的 (7) 中就取相应同一有界 packet，使所有目标在同一高度先相加。

## 8. clip、三近面和参数场的适用边界

1. **零高度。** (9)、(14)、(19)–(25) 始终保留正部和真实纤维上截断。没有用除以高度的相对误差，也没有在帽趋零处省略 clip。写成 \(Y_1=Y_0-\mathcal K Y_0\) 时必须理解为 (9) 的实际高度运算；本稿不将集合减法或忽略 clip 的形式表达认作该等式。

2. **三近面。** 所有 \(q\notin\{l,m\}\) 图、以及全局最浅轴下多个目标高度相加的区域均保留。这些正好可能携带八阶主成本。本文只证明把实际缺集反馈换成 \(Y_1\) 的误差晚于八阶；没有计算这些图的八阶系数或对其使用旧核心剪裁。

3. **实际根和宽度。** 每个目标帽保留其原整横宽，每个源保留自身 \(R_i^{n-1}\) 分母和有效宽。根差截幅只按 (3) 原定义使用。没有用等尺度模型替换，再把一个未经控制的根差误差宣布为 \(o(v^8S)\)。

4. **非逐边参数。** R049 的正常一阶场只有加权 \(L^2\) 控制，本稿的算子小常数对整个固定 packet 窗一致，所以不需要先把这些场截成逐点 \(O(v)\)。真正评估 \(\sum|Y_2|\) 的第二变分仍须合法处理参数尾、槽损失、真实三面几何及所有有界非正常项；本文不预支那一优化。

## 9. 实际计算与当前状态

Terra 在本机 Mathematica 15.0 实际运行 r052-computation/operator.wl/log：clip 的单调与 1-Lipschitz 违反集均为 False；代入
\(\lambda=Cv,M_0=Aav^5,\rho=Rav^7+Eav^9\)
的递推后，(23) 和 (24) 的最低可能阶均为 \(v^9\)，\(\lambda^2M_0=A C^2av^7\)。这只核 clip 和有限递推，不替代实际帽 Tonelli、坐标分割、原见证必要性或 CDF 物理分配。

依赖：R049 已审实际稳定性；R045 的同一正常层及原帽参数；R042 的强 CDF 及实际退款核。唯一尚须新供货并独立核合的几何接口是 (N1)–(N2)。本稿不因 R051 的口头预告就去条件化。

**当前结论：\(\mathrm{Nec}_{51}\) 加同层既有输入推出 (14)、(21)–(25) 的实际第二反馈引理。** 没有宣布 \(J_n\) 是真实八阶上系数，也没有给出全输入精确八阶展开。
