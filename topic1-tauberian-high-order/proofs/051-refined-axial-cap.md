# R051：实际薄轴帽的斜率保留与全局最浅面禁层

2026-11-05。固定 \(n\ge2\)，\(a=v^n\)。本稿改进 R037/R042 中用整个目标缺量平移轴帽的误差；不识别八阶系数，不对未固定的大角窗作统一断言。仓库、旧冻结源及他人文件均只读。

## 1. 同一实际层和精确范围

采用 R049 同一新第八层的实际对象，或任何满足下列明确条件的有限层。所有常数仅依赖固定维数、原选择常数和下文预先固定的管宽、角窗，不依赖输入数目及尺度跨度。指标只指真正带有 \(F_i,G_i,R_i\) 的实际外部源；标记集的无源根不被补造为这种指标。

源 \(F_i\subset G_i\subset E^c\)，质量根 \(R_i=(|G_i|/a)^{1/n}\)，原始、有效、中点源盒分别为
\[
 B_i^0=x_i+W_i^0,\qquad B_i'=x_i+W_i',\qquad B_i=x_i+W_i,
 \qquad B_i'\subset B_i\subset B_i^0.
\]
有效端总移距为 \(\Delta_i=\Delta_i^{\rm end}\)，不是删除质量。原始盒 \(B_i^0\subset Q(x_i,Avr_i)\)。三个盒的每条全边长均在 \([c_wvR_i,C_wvR_i]\) 内，\(|F_i|\ge c\,aR_i^n\)。原指定的严格见证对每个 \(x\in F_i\) 满足
\[
 G_i\subset Q(x,r(x)),\qquad
 R_i<r(x)\le(1+Av)r_i,\qquad R_i/r_i=1+O_A(v).
 \tag{1}
\]
对 \(j\ne i\)，\(F_j\cap G_i=\varnothing\)，各 \(F_j\) 两两不交。本文始终使用 (1) 中这个见证，不替换或优化其半径。

取深处源。每个有向轴面有原 R014 的真正单临界面邻居 \(h\)，其根与 \(R_i\) 可比，逆邻接次数统一有界。原有向端余量 \(g_{i\to h}\ge0\) 以目标有效端和源原始端定义，满足
\[
 g_{i\to h}\le 2\Gamma_{ih}\le C_AvR_i,\qquad
 \sum_{i\ {\rm deep}}R_i^{n-1}\max_h g_{i\to h}\le C_A\mathcal F.
 \tag{2}
\]
这里的整横截面性质可直接核验。R014 真正面邻居引理在归一化 \(r_i=1,x_i=0\) 后，给每个切向轴 \(|x_{h,\ell}|\le1/2+1/128\) 和 \(r_h/r_i=1+O(v)\)。因此任意 \(x\in B_i^0,z\in B_h^0\) 满足
\[
 |z_\ell-x_\ell|
 \le(1/2+1/128)r_i+Av(r_i+r_h)<R_i
       \qquad(\ell\ne q).
 \tag{2a}
\]
故邻居的原始整横截面对所有源点、所有 \(r\ge R_i\) 均精确进入见证，甚至保留固定物理裕量。这不依赖源点离其他面的距离；有效端错移也不会改变两个原始盒的包含。不需要额外 \(F_\tau\) 或切向 \(\Delta\) 支出。更一般地，原临界判据保证：任意二临界目标的每个非临界轴，两个完整原始源区间的任意坐标对相距小于 \(\min(R_i,R_j)\)。这些是原实际面容量和原临界盒的性质，不是新作的近晶格假设。

固定 R045/R049 实际正常层
\[
 Y_{h,N}\subset B_h'\setminus F_h,\qquad
 \mathscr D_h=B_h^0\setminus F_h,\qquad
 m_h=|\mathscr D_h\setminus Y_{h,N}|.
 \tag{3}
\]
这里 \(Y_{h,N}\) 正是 R049 量 \(R_N\) 使用的层，不在本稿重新选择。它包含于固定两轴管
\[
 \mathcal C_h(K_N)
 =\{z\in B_h':\text{至少两个不同坐标各距其有效盒某一端不超过 }
                      K_Nv^3R_h\}.
 \tag{4}
\]
采用已审 R049 的定量输入
\[
 R_N=\frac va\sum_hm_h=O(v^8S),\qquad
 \mathfrak E=\sum_hR_h^{n-1}\Delta_h=O(v^8S),\qquad
 \mathcal F=O(v^8S),\qquad S=\sum_h r_h^n.
 \tag{H49}
\]
以下先给不代入 (H49) 的误差式。因此也可把 (H49) 保留为显式条件。正常参数的逐边刚性、未加权轮廓接近及任何 LT 收敛速率都不是本稿假设。

主结论是：任意预先固定 \(T<\infty\) 的二临界目标，可在全局最浅有效面坐标中同时形成一个实际联合禁层，其帽高度误差总成本不超过
\[
 C_{T,K_N}v^2[v^8S+R_N+\mathfrak E+\mathcal F]
       =O_{T,K_N}(v^{10}S)=o(v^8S).
 \tag{5}
\]
包括源有三个或更多坐标近面的区域。这里只改进轴帽/端移/容量曲率的误差，不删除联合禁层本身的三近面八阶主量。

## 2. 任意深度薄帽：正规缺层只降低斜率

在 \(B_h^0\) 的任意有向轴 \(q\) 上，从其原始端向内取深度 \(d_*= \min(d,w_{hq}^0)\) 的整横截面帽 \(C_{h,q}(d)\)，其中 \(d\ge0\)，
\[
 A_{h,q}=\prod_{\ell\ne q}w_{h\ell}^0.
\]
对包含轴 \(q\) 的管内轴对，忽略其 \(q\) 约束，另一轴在原始整宽中的占比至多 \(2K_Nv^2/c_w\)。对不包含 \(q\) 的轴对，两条轴各给该占比。有效面可以相对原始面移动；限制区间的总长度仍不超过 \(2K_Nv^3R_h\)。有限并估计故给
\[
 |Y_{h,N}\cap C_{h,q}(d)|
 \le\theta_N(v)A_{h,q}d_*,
 \qquad
 \theta_N(v)=\frac{2(n-1)K_N}{c_w}v^2
       +\frac{4\binom{n-1}{2}K_N^2}{c_w^2}v^4.
 \tag{6}
\]
此式对每个 \(d\ge0\) 同时成立，包括 \(d\downarrow0\)；没有用除以 \(d\) 的误差。\(n=2\) 的第二项为空。由 (3)，取
\[
 \beta_{h,q}=(1-\theta_N(v))A_{h,q}\ge A_{h,q}/2
\]
并令 \(v\) 足够小，有精确集合下界
\[
 |F_h\cap C_{h,q}(d)|
 \ge [\beta_{h,q}\min(d,w_{hq}^0)-m_h]_+.
 \tag{7}
\]
在 (7) 中正常部分只改变 \(\beta_{h,q}\)；平移帽的唯一质量是 \(m_h\)。没有再次扣除整个 \(|\mathscr D_h|\)。

## 3. 原始全部见证半径上的转折点

反射坐标，使源 \(i\) 的选定有效浅面为第 \(q\) 轴右端。\(x\) 距该有效面为 \(s\)。写源原始、有效右端为 \(b_i^0,b_i'\)，目标原始、有效左端为 \(a_h^0,a_h'\)。定义
\[
 \delta_i^+=b_i^0-b_i',\quad
 \delta_h^-=a_h'-a_h^0,\quad
 a_h'-b_i^0=R_i+g_{i\to h}.
\]
实际从 \(x\) 到目标原始左端的距离等于
\[
 a_h^0-x_q=R_i+s+g_{i\to h}+\delta_i^+-\delta_h^-.
 \tag{8}
\]
这里 \(\delta_i^+\le\Delta_i\)，目标端移 \(\delta_h^-\ge0\) 只改善捕获。所有端余量的定义与原 R014 一致，不把中点隙或有效对有效距离换作 \(g\)。

对原半径写 \(\ell=r(x)-R_i>0\)。由 (8) 及原非临界整宽捕获，保证进入 \(Q(x,r(x))\) 的原始轴帽深度至少为
\[
 d=(\ell-s-g_{i\to h}-\Delta_i)_+.
 \tag{9}
\]
这只是下界，未把实际穿透与 (9) 认作相等。

若 \(d\ge w_{hq}^0\)，整个 \(F_h\) 与完整 \(G_i\) 同时进入原严格见证；但 \(|F_h|\ge c aR_i^n\)，而 (1) 的剩余容量
\[
 C_i(\ell)=a[(R_i+\ell)^n-R_i^n]\le C_AavR_i^n
 \tag{10}
\]
更小。这对足够小 \(v\) 矛盾。因此所有实际有效见证均有 \(d<w_{hq}^0\)。这里先用真实整目标排除了 full-cap 饱和，未从可能为零的 (7) 逆推出它。

省略方向下标，置
\[
 b=s+g_{i\to h}+\Delta_i,\qquad z=b+m_h/\beta_h.
\]
在整个原半径区间 \(0\le\ell\le L_i\le C_AvR_i\) 上，
\[
 C_i'(\ell)=na(R_i+\ell)^{n-1}<\beta_h,\qquad
 [\beta_h(\ell-b)_+-m_h]_+=\beta_h(\ell-z)_+.
 \tag{11}
\]
第一式因 \(\beta_h\asymp(vR_i)^{n-1}\)，容量导数与它之比为 \(O_A(v)\)。于是容量减 (7) 的下包，在转折点 \(z\) 前递增、之后递减，给
\[
 C_i(\ell)-|F_h\cap Q(x,r(x))|\le C_i(z).
 \tag{12}
\]
即使 \(z>L_i\)，也由 \(C_i\) 的单调性成立。更直接地，\(z\ge\ell\) 时比较 \(C_i(\ell)\le C_i(z)\)；\(z<\ell\) 时只在 \([z,\ell]\subset[0,L_i]\) 积分 (11)。全部使用原指定半径；没有假定 \(r(x)=R_i+s\)。

粗宽界、(2) 及 \(m_h\le|B_h^0|\) 还给 \(z\le C_{A,T}vR_i\)（以下只用 \(s\le C_Tv^4R_i\)）。因此
\[
 C_i(z)\le C_i(s)+C_AaR_i^{n-1}
                   (g_{i\to h}+\Delta_i+m_h/\beta_h).
 \tag{13}
\]
将正规缺量并入斜率后，并无新的 \(v^2s\) 半径损失；这正是只扣全 \(|\mathscr D_h|\) 的旧估计所不能给出的结论。

## 4. 全局最浅有效面坐标和同源全部有界帽

对每个 \(B_i'\) 作如下几乎处处唯一的有限分片。各坐标先选最近的有效面，记符号为 \(\sigma_\ell\)；平分面及等深并列按固定顺序分配。再选全局最浅轴 \(q\)。用各自有效面的向内深度记
\[
 d_q=s,\qquad d_\ell=s+y_\ell\quad(\ell\ne q),\qquad
 y_q=0,\quad y_\ell\ge0,\quad
 0\le s\le d_i(y):=\min_\ell(w'_{i\ell}/2-y_\ell).
 \tag{14}
\]
只保留 \(d_i(y)\ge0\) 的底面。映射 \((s,y_{\ell\ne q})\mapsto x\) 的 Jacobian 绝对值为 1，且 (14) 是精确的真实盒域；不把其他坐标剪成核心。

固定 \(T<\infty\)，选取任意一组实际二临界邻居，其两临界交长均不超过 \(Tv^3R_e\)，\(R_e=\max(R_i,R_j)\)。允许其中有正常、非正常及任意小的正交长。每一 chart 只保留临界方向与两个所选最近面符号一致的目标。对于不一致的目标，在 \(s\le C_Tv^4R_i\) 的本高度窗内，朝向对面所需源深至少 \(w'_{i\ell}/2\)，远大于 \(O_T(v^3R_i)\)，所以它没有本章所定义的正标称帽；忽略它不损害下界。

对所选 \(j\) 的两临界轴 \(k,l\)，固定
\[
 c_{ij}=\operatorname{clip}\left((R_i-R_j)/2;\
           -\max(L_k,L_l),\,\max(L_k,L_l)\right),\qquad
 P_{ij,k}(y)=(L_k+c_{ij}-y_k)_+,
 \quad P_{ij,l}(y)=(L_l+c_{ij}-y_l)_+.
 \tag{15}
\]
\(\Pi_{ij}(y)\) 是锚在目标 \(B_j^0\) 原始面对端的两个此前缀区间，乘其余轴的原始整宽。所有正帽宽 \(O_T(v^3R_i)\) 均小于目标相应全宽。若下截幅生效，两帽宽全为零；若上截幅生效，它只缩小由真实根差可保证的帽。因此不把截幅后根差当作实际根差。

置
\[
 \mathcal S_i(y)=\sum_j|F_j\cap\Pi_{ij}(y)|
   =\sum_j\bigl(|\Pi_{ij}(y)|
                  -|\mathscr D_j\cap\Pi_{ij}(y)|\bigr),
 \qquad
 h_i(\mathscr D;y)=\frac{\mathcal S_i(y)}{naR_i^{n-1}}.
 \tag{16}
\]
每项均为实际非负量；不另扣整个目标缺量。有限度及整宽界给 \(0\le h_i\le C_Tv^4R_i\)。同一实际目标在该 chart 只有一套帽，源上的 chart 分片也只有一份体积。

对任意 \(x\in F_i\cap B_i'\) 且 \(s\le C_Tv^4R_i\)，使用其原半径 \(\ell=r(x)-R_i\)。两中点面对端间距恰为 \((R_i+R_j)/2-L_k\)。从中点源端换成有效源端只多至多 \(\Delta_i\)；换到目标原始端只减小端距。因此，若 \(\ell\le s\)，实际帽相对 \(\Pi_{ij}(y)\) 至多丢失沿两轴厚为 \((s-\ell)+\Delta_i\) 的边带；若 \(\ell\ge s\)，至多丢失厚为 \(\Delta_i\) 的边带。

所有所选边带的总质量上包为
\[
 D_{i,T}\bigl((s-\ell)_++\Delta_i\bigr),\qquad
 0\le D_{i,T}\le C_Tv^{n+1}R_i^{n-1}
                \le\tfrac12naR_i^{n-1}
 \tag{17}
\]
（最后令 \(v\) 足够小）。边带式对任意厚度仍有效：每个矩形截宽之损至多为该厚度乘另一短宽及全部横宽，不要求实际丢失条彼此不交。

若 \(\ell\le s\)，所有角目标实际质量与完整 \(G_i\) 在同一见证中相加，得
\[
 \mathcal S_i(y)\le C_i(\ell)+D_{i,T}(s-\ell+\Delta_i)
                  \le C_i(s)+D_{i,T}\Delta_i.
 \tag{18}
\]
最后一步由 \(C_i'(t)\ge naR_i^{n-1}\ge D_{i,T}\) 积分得到。若 \(\ell\ge s\)，取全局最浅面 \(q,\sigma_q\) 的真正单临界邻居 \(h\)。它不同于每个二临界 \(j\)；其实际质量、各个 \(F_j\) 与自身 \(G_i\) 全部排他。故必须在同一个见证中先把轴帽与所有角帽相加，再用 (12)–(13)，得到
\[
 \mathcal S_i(y)\le C_i(s)+
 C_AaR_i^{n-1}(g_{i\to h}+\Delta_i+m_h/\beta_h)
                     +D_{i,T}\Delta_i.
 \tag{19}
\]
这没有向每个目标分别分配一份完整轴容量，也没有将 \(q\) 预设成某个目标的临界轴。

在当前高度窗，容量曲率满足
\[
 0\le\frac{C_i(s)}{naR_i^{n-1}}-s\le C_Tv^8R_i.
 \tag{20}
\]
于是以长度单位定义
\[
 e_i=v^8R_i+\Delta_i+
       \max_{h\ {\rm 真正有向面邻}}
           (g_{i\to h}+m_h/\beta_h),
 \tag{21}
\]
对所有上述实际源点同时有
\[
 \boxed{\quad s\ge h_i(\mathscr D;y)-C_T e_i.\quad}
 \tag{22}
\]
证明覆盖三个、四个乃至所有源坐标同时靠近有效面的情形；所有目标在 (16) 中相加后才产生一个源高度。

## 5. 物理误差、鲜成本与不重复计费

面邻逆度和 \(R_h/R_i\asymp1\) 给
\[
 \sum_{i\ {\rm deep}}R_i^{n-1}
       \max_h\frac{m_h}{\beta_h}
 \le C v^{1-n}\sum_hm_h=C R_N.
 \tag{23}
\]
由 (2)、(21)，
\[
 \sum_{i\ {\rm deep}}R_i^{n-1}e_i
 \le C[v^8S+R_N+\mathfrak E+\mathcal F].
 \tag{24}
\]
无须逐源 \(m_h/\beta_h=o(v^6R_i)\)；所有异常源已由此加权和覆盖。

对固定 chart，只有某个正标称帽存在时，(22) 的高度侵蚀才有损失。若最浅轴 \(q\) 属于该目标轴对，另一临界底变量占宽 \(O_T(v^3R_i)\)，其余 \(n-2\) 个底变量各占 \(O(vR_i)\)，故底面测度不超过
\[
 C_Tv^{n+1}R_i^{n-1}.
 \tag{25}
\]
若 \(q\) 不属于该目标轴对，则两条临界底变量各占 \(O_T(v^3R_i)\)，其余 \(n-3\) 个底变量各占 \(O(vR_i)\)，所得 \(C_Tv^{n+3}R_i^{n-1}\) 更小；二维此情形为空。有限 chart 数、有限邻居度可吸入常数。这是对误差支撑的估计，不是丢掉主层的三近面区域。

正部、以及真实域截断 \(h\mapsto\min(d_i(y),h_+)\)，均为 1-Lipschitz。因此定义未侵蚀的联合层
\[
 \mathcal B_i(\mathscr D)
 =\bigcup_{\rm charts}
    \{x(s,y):0\le s<\min(d_i(y),h_i(\mathscr D;y))\}\subset B_i',
 \tag{26}
\]
以及将高度换成 \((h_i-C_Te_i)_+\) 的侵蚀层 \(\widetilde{\mathcal B}_i\)，有
\[
 \widetilde{\mathcal B}_i\subset\mathscr D_i,\qquad
 \sum_i|\mathcal B_i(\mathscr D)\setminus\mathscr D_i|
 \le C_{T,K_N}v^{n+1}[v^8S+R_N+\mathfrak E+\mathcal F].
 \tag{27}
\]
所有集合的零测边界可按固定规则处理。代入 (H49) 后，右侧是 \(O(av^9S)\)，特别是 \(o(av^7S)\)。

实际鲜导数在相关质量区间上下均与 \(v/a\) 可比；上界由 \(|F_i|\asymp aR_i^n\) 直接得到。将 (27) 的物理体积损失乘任一所需的完整鲜导数，损失不超过 (5)。没有把 (6) 的 \(1-\theta_N\) 乘到主角成本上：它只用于维持 (11) 的严格下降斜率。

本稿的 \(m_h\) 是已有 (H49) 的数值数据，不是新收入。若下一步以 \(\widetilde{\mathcal B}\) 替换旧 \(Y_N\) 进行强 CDF 结算，仍只使用同一份 \(J_i\) 一次；不能把旧层成本与新层成本相加，也不能再以全部 \(|\mathscr D_h|\) 平移同一轴帽。目标缺量只在 (16) 的真实帽中扣一次。

## 6. 旧正常层的精确锚定及反单调接口

本节令所选有界目标至少包含 R049 的全部正常集合 \(\mathcal N\)，并固定 \(T\) 覆盖其短角。可以只选正常目标，也可以加上其他固定有界目标；相关边的两端均在本次允许的深处族内。所有几何数据、原始帽、根、clip、chart 与目标清单在本节固定，不随测试缺集改变。未使用的源令 \(\mathcal B_i=\widetilde{\mathcal B}_i=\varnothing\)；其旧 \(Y_{i,N}\) 也为空。如此 (30)–(34) 可在整个实际外部源族求和，仍保留 (H49) 对这些源的 raw 缺量，不因定义空算子而删去它们。

R049 的 \(Y_{i,N}\) 按 R045 §4 定义：源坐标是 \(B_i'\) 的有效面对深度；在原有向二轴组的两个排序中，先对正常目标的实际帽质量相加，减去非负的旧高度侵蚀 \(C R_iE_i(T)\)，再取正部；其余坐标取固定横向核心。旧层中
\[
 s=O_N(v^4R_i),\qquad s+y=O_N(v^3R_i),
 \qquad\text{其余坐标距有效两面}>H_Tv^3R_i.
 \tag{28}
\]
故对足够小 \(v\)，全局最浅轴正是该排序中的第一轴，两个近面符号也是 (14) 的最近面符号。转换到 (14) 后，此二轴目标的帽仍逐字为原始目标上的
\[
 [0,P_{ij}]\times[0,(Q_{ij}-y)_+]\times W_{j,\perp}^0.
 \tag{29}
\]
源的其他坐标改写成 \(s+y_\ell\) 不改变目标的原始整横截面。因此新 (16) 至少包含原组全部非负实际帽，而 (26) 没有扣旧侵蚀。原点本就在 \(B_i'\) 内，亦满足真实域截断。这逐点证明
\[
 \boxed{\,Y_{i,N}\subset\mathcal B_i(\mathscr D)\,}
 \quad\text{（零测边界除外）}.
 \tag{30}
\]
没有比较新旧侵蚀常数，也没有将源 chart 从有效盒移回原始盒。若改用原始源 chart，(29) 就不再是无误差的恒等；本稿不作该更改。

为说明其算子意义，对任意测试缺集族
\(D=(D_j)\)、\(D_j\subset B_j^0\)，以
\[
 h_i(D;y)=\frac1{naR_i^{n-1}}
       \sum_j\bigl(|\Pi_{ij}(y)|-|D_j\cap\Pi_{ij}(y)|\bigr)
 \tag{31}
\]
代替 (16)，定义同一 (26) 的 \(\mathcal B(D)\)。每项非负，几何域固定，故
\[
 D\subset D'\ \Longrightarrow\ \mathcal B(D')\subset\mathcal B(D).
 \tag{32}
\]
clip 和根在这里不依赖测试缺集质量；测试算子不重新定义原 \(G_i\) 或实际根。输出仍是有效盒的子集，输入则允许原始盒中的全部缺集。更换 chart 没有破坏反单调性。

由 (27)、(30)，实际缺集满足两个独立接口
\[
 \sum_i|\mathcal B_i(\mathscr D)\setminus\mathscr D_i|
       =O(av^9S)=o(av^7S),\qquad
 \sum_i|Y_{i,N}\setminus\mathcal B_i(\mathscr D)|=0.
 \tag{33}
\]
以及
\[
 \sum_i|\mathscr D_i\setminus\mathcal B_i(\mathscr D)|
 \le\sum_i|\mathscr D_i\setminus Y_{i,N}|=O(av^7S),\qquad
 \sum_i|\mathscr D_i\triangle\mathcal B_i(\mathscr D)|=O(av^7S).
 \tag{34}
\]
这里只给出真实算子的锚定与近固定点误差；第二次反馈迭代的定量 Lipschitz 性、三近面完整积分、参数二次变分及八阶上界仍由各自后续命题证明。

## 7. 常数次序、边界和实际核验

原有限层、CDF 斜率、定位常数与真正面邻居先固定；再固定 R045 的正常窗、小管 \(K_N\)、旧正常层及共同核心；本稿的 \(T\) 是任意随后预先固定的有限值。最后取小 \(v\)，同时满足 \(\theta_N\le1/2\)、(11)、(17)、短帽小于真实目标宽及 (28) 的 chart 识别。没有 \(T(v)\)，没有对一串趋大的窗声称统一 \(o(v^8)\)。

所有结论局限于真正面邻居仍属同一实际外部族的深处源，以及所列固定有界二临界目标。原外部边界、无限大角尾、高支撑交量及其他层不借本稿自动消失；(H49) 和 R050 等接口如需使用，必须在同层另外接入。本稿没有删去任何原盒基线。固定数量的源有向面和所有目标的有限逆度已在 (23)、(25) 中计入。

本轮实际读取的主要源：

* R014 正式工作源 global.tex 的真正面邻居、全部源点面对容量及临界轴判据；尤其有向余量是目标有效端减源原始端再减 \(R_i\)，本稿 (8) 保留此定义；
* R037 general_soft_price.md §2 的原指定见证、原始整横截面帽及两半径分支；
* R042 normal_packet_and_tail_obstruction.md，SHA dd291e7bc4294483571eaa8dcfb85e8e409266fbfb8eb8fb59ebcea6b9e2dcff，§3；
* R045 joint_two_axis_price.md，SHA 8f0ad8b01bb443789b33f55188c7a04d000a3ceaee3fd61d0a83eb458c9bf8ea，§4 的同一实际 \(Y_N\)；
* R049 competitive_eighth_stability.md，SHA 5c4d7b371a0c33c39449afd686de599c49570aaa639d4f7db7d3c1f166a040f4，及其正式独审；本稿只需其中 (15) 的 raw 正余量、端移和普通面预算。

Terra 已在本机 Mathematica 15 ARM 实际运行，文件位于 r051-computation/：

* axis_cap.wl / axis_cap.log：正部转折恒等的违例集为 False；归一化 \(n=2,3\) 净帽不等式的违例集均为 False；记录一般幂函数的导数与低半径差因式。一般 \(n\) 的证明是本文 (11)–(12)、(18) 的解析积分，不将有限维 Reduce 当作全维证明。
* axis_fraction.py / axis_fraction.json：8 组精确有理边界全部通过，包括零深度、转折点等号及 \(z>L\)。
* README.md 明列这些标量证书不证明实际几何、不证明相交度和不重复计费。本文已实际读取脚本和成功输出，未重跑旧价格多项式。

状态：供独立审查的新局部证明。证明对象是 (5)、(22)、(27)、(30)–(34)，不是 \(C_n\) 的精确八阶系数。
