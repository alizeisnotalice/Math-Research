# R029：实际 void、固定窗口半径传输与切片碰撞

日期：2026-11-05。完整新引理候选，待独立审查。只写独立 scratch。

本稿证明两个全输入新接口：实际 CDF 盒并集在原见证覆盖中的小 void，以及任意固定物理窗口内 actual roots 的加权 \(O(v^6S)\) 相干性。另给出覆盖 probe 和切片碰撞的定量收费。**不声称高维 raw Slot 或 resolved Slot 已闭合**；最后精确列出 helper 分配与退化角的剩余步骤，不把未定量 LT 面积误差当作 \(o(v)\)。

## 1. 输入与记号

固定 \(n\ge2\)、\(a=v^n\)、\(A_v=(1+v)^n\)，\(V=|E|\)、\(M=|K|\)。原有限输入及外部紧集在临界分支
\[
 M/V\ge a/(A_v-a)-Bav^6,\qquad B<\infty\text{ 固定}.
\]
采用已审 R17/R21 的同一固定末层粗先验，不物理删除源或输入。原族为 \(Q_i=Q(x_i,r_i)\)，半盒 \(B_i=Q(x_i,r_i/2)\)，\(S=\sum_i r_i^n\)。标记集 \(I\) 及末层实际 CDF 盒为
\[
 R_i=(g_i/a)^{1/n},\quad X_i=x_i+W_i+[-R_i/2,R_i/2]^n,
 \quad U=\bigcup_{i\notin I}X_i,
 \quad \mathcal Q=\bigcup_i Q_i.
\]
写 \(M'=|K\setminus\bigcup_IQ_i|\)、\(M_I=M-M'\)、\(\mathcal D=\sum g_i-M'\)。明确使用的已有量为
\[
 \Xi=V+M-A_vM/a=O(v^6S),\quad S\asymp V,\quad
 S(I)=O(v^5S),\quad M_I/a\le S(I), \tag{1}
\]
\[
 \mathcal D/a=O(v^6S),\quad
 \sum_i(g_i-\mu_i)/a=O(v^5S),\quad
 \sum_i\varepsilon_i r_i^n=O(v^5S),\quad\mathcal F=O(v^6S). \tag{2}
\]
\(\Xi\) 的双侧控制分别来自临界假设及 R14 粗全局界。原半盒不交，原中心强分离。所有 \(X_i\subset Q_i\)，相交邻居数有界且尺度比为 \(1+O(v)\)。R19 的全源粗 Hölder 界还给
\[
 \sum_{i<j}|X_i\cap X_j|=O(v^6S). \tag{3}
\]
每个固定稳健面裕量 \(\tau>0\)，R17 给全部该类面边
\[
 \sum_{(p,q)\in\mathscr E_\tau}|R_p-R_q|r_p^{n-1}=O_\tau(v^6S). \tag{4}
\]
R19 给坏根 \(Z_v\)，\(S(Z_v)=O(v^{9/2}S)\)；其补集 \(G_v\) 在任意固定物理窗口一致匹配完整单位立方铺砌，相关指标最终都在外部深处。只用 LT 的 \(\delta_R(v)\to0\)，不要求速率。

## 2. 实际全域 void 预算

**定理1。** 令 \(\mathcal H=\mathcal Q\setminus U\)。则
\[
 \boxed{|\mathcal H|\le C_{n,B}(v^{n+1}+v^5)S.} \tag{5}
\]
尤其 \(n=3\) 为 \(O(v^4S)\)，\(n\ge4\) 为 \(O(v^5S)\)，均为 \(o(vS)\)。这个集合是真实原覆盖中未被实际盒占据的部分，不是某个模型铺砌的对称差。

证明：由单盒基线体积、(2) 及二项 Bonferroni，
\[
 \begin{aligned}
 |U|&\ge\sum_i|X_i|-\sum_{i<j}|X_i\cap X_j|\\
 &\ge A_v\sum_i g_i/a-C(v/a)\sum_i(g_i-\mu_i)-O(v^6S)\\
 &\ge A_vM'/a-O(v^6S)\\
 &\ge V+M-O(v^5S).
 \end{aligned} \tag{6}
\]
最后用 \(A_vM/a=V+M-\Xi\) 及 \(M_I/a\le S(I)\)。原严格见证给 \(|\mathcal Q|\le V+aS\)，又 \(U\subset\mathcal Q\)，所以
\[
 |\mathcal H|=|\mathcal Q|-|U|\le aS-M+O(v^5S). \tag{7}
\]
基本容量 \(S\le V/(1-a)\) 和临界假设给
\[
 0\le S-M/a
 \le V\left[\frac1{1-a}-\frac1{A_v-a}+Bv^6\right]
 =O(vS),
\]
因为括号内两倒数之差为
\((A_v-1)/[(1-a)(A_v-a)]\)。故 \(aS-M=O(v^{n+1}S)\)，代回即得 (5)。证毕。

### 真实 probe 的 Fubini 付费

固定轴对 \(q,r\)。对根 \(i\) 选有限个平面矩形 probe \(P_i\)，每个面积至少为 \(c r_i^2\)，其两坐标均在 \(Q_i\) 投影内。令 \(Y_i\subset\pi_{\widehat{qr}}X_i\)，若 \(y\in Y_i\)，要求
\[
 |P_i\setminus U_y|\ge c r_i^2/2,
 \qquad U_y=\{z\in\mathbb R^2:(z,y)\in U\}.
\]
则
\[
 \boxed{\sum_i r_i^2|Y_i|\le C_{n,c}|\mathcal H|.} \tag{8}
\]
证明：\(P_i\times Y_i\subset Q_i\)。逐 fibre 的空面积积分至少为 \(c r_i^2|Y_i|/2\)；原 Besicovitch 族 \(Q_i\) 的重数至多为 \(H=2^n\)，每根 probe 数固定，故对 \(i\) 求和最多把 \(\mathcal H\) 重复纯维数次。有限盒保证所有切片可测，Fubini 合法。证毕。

这只支付“确实有一块固定二维面积未覆盖”的 fibre。尚不能仅从一个指定 helper 不出现，就断言整个 probe 空；其他胞可能覆盖它，必须另做实际局部选择。

## 3. 单位铺砌内存在统一稳健面路径

**引理1。** 固定 \(R\ge1\)。存在只依赖 \(n,R\) 的 \(N_R<\infty\)、\(\tau_R>0\)，使任意轴平行单位立方铺砌中，中心在 \(Q(0,R)\) 内的两胞可用至多 \(N_R\) 条相邻面边相连。沿路的每个面接触片面积至少为 \(\tau_R\)，故每条切向交边长也至少为 \(\tau_R\)。路径可选在 \(Q(0,R+2)\) 所遇的胞中。

证明：在两个指定胞中心各取全边长 \(1/4\) 的内盒，独立均匀取端点 \(A,B\)。连接线段始终在固定凸窗口内。遇到该窗口的单位胞数有上界，例如取 \(N_R=\lceil(2R+6)^n\rceil\) 即可：中心位于固定扩张窗口，内部不交的单位胞按体积计数。

对某个面片 \(F\subset\{x_l=c\}\)，条件化两个端点的第 \(l\) 坐标。若线段经过该平面，交点参数 \(\theta\in[0,1]\) 已确定，切向点为
\((1-\theta)A_{\widehat l}+\theta B_{\widehat l}\)。两个切向端点仍独立均匀。卷积密度不超过
\[
 \min\left\{\frac{4^{n-1}}{(1-\theta)^{n-1}},
                 \frac{4^{n-1}}{\theta^{n-1}}\right\}\le8^{n-1};
\]
\(\theta=0,1\) 直接用另一个端点密度。因此穿过面积为 \(|F|\) 的面片的概率至多 \(8^{n-1}|F|\)。

窗口内至多有 \(N_R^2\) 个两胞面接触片。取
\[
 \tau_R=[4\cdot8^{n-1}N_R^2]^{-1}.
\]
所有面积小于 \(\tau_R\) 的面片被穿过的总概率至多 \(1/4\)。同一密度论证说明线段碰余维至少二的胞边界的概率为零。故存在一条线段只穿过面积至少 \(\tau_R\) 的面片；它给出所需路径。每个单位胞凸，线段不能离开后重入，路径长度至多为所遇胞数。面片各边不大于一，乘积至少为 \(\tau_R\) 则各边均至少为 \(\tau_R\)。证毕。

本引理容许 T 接缝、多面片及宏观错列；未用面面对齐或整数格点标签。

## 4. 实际固定窗口的 root 相干性

对 \(i\in G_v\) 定义完整物理清单
\[
 \mathcal I_i(R)=\{k:B_k\cap Q(x_i,Rr_i)\ne\varnothing\}.
\]
R19 保证固定 \(R\) 下这份清单数量有界、所有半径与 \(r_i\) 比值趋一，全部相关源有效。实际根仍是 \(R_k=(g_k/a)^{1/n}\)。

**定理2。** 对每个固定 \(R<\infty\)，
\[
 \boxed{\sum_{i\in G_v}\sum_{k\in\mathcal I_i(R)}
       |R_i-R_k|r_i^{n-1}=O_{n,B,R}(v^6S).} \tag{9}
\]
不存在 \(\delta_R(v)\) 的加法面积误差。

证明：按 \(r_i\) 归一化，用 LT 的全空间单位铺砌。清单内胞及所需固定扩张窗口里的胞全部有原族对应。引理1以固定 \(\tau_R\) 给任一模型根到模型 \(k\) 的有限稳健面路径。对足够小 \(v\)，一致匹配误差小于一个依赖 \(\tau_R\) 的固定阈值，所以实际相邻原半盒每个切向交长至少为 \(\tau_Rr_i/2\)，法向缝趋零，半径比趋一。原中心分离使法向缝非负。于是它们是 R17 命题2中固定裕量（可取 \(\tau_R/4\)）的稳健候选；全部在深处。沿实际路径精确望远镜，
\[
 |R_i-R_k|\le\sum_{\{p,q\}\text{ 在路径上}}|R_p-R_q|.
\]
各 \(r_p/r_i\to1\)，所以乘以 \(r_i^{n-1}\) 可换为相应 \(r_p^{n-1}\)，只有固定有界因子。

每根有固定数量目标和固定长度路径。一个实际面边能被计入的根中心都在其固定倍尺度的窗口内，且这些根的半径与它可比；原半盒不交给根数的纯 \(n,R\) 上界。因此所有路径对同一实际面边的重复次数有固定上界。对根、目标及路径求和后用 (4)，即得 (9)。证毕。

特别，即便 helper 与原源在某个 transverse 方向靠近边界而不是直接的稳健面邻居，只要它在同一固定物理窗口，也能用这样的路径支付根权转移。此处只构造固定阈值路径，再收费真实 root 差；没有把模型中心误差乘上面积。

## 5. 切片碰撞的严格收费

固定轴对 \(q,r\)，把 \(\pi_{qr}\) 用于二维平面投影，\(\pi_{\widehat{qr}}\) 用于其余坐标。取同一个固定窗口 \(R\)。对 \(i\in G_v\)，在 \(\mathcal I_i(R)\) 中取所有不同胞对 \(k,l\)。定义
\[
 A_{kl}^{qr}=|\pi_{qr}X_k\cap\pi_{qr}X_l|,
\]
及坏 fibre 集
\[
 C_i(\lambda)=\pi_{\widehat{qr}}X_i\cap
 \bigcup_{\substack{k,l\in\mathcal I_i(R),\ k\ne l\\
                    A_{kl}^{qr}\ge\lambda r_i^2}}
       (\pi_{\widehat{qr}}X_k\cap\pi_{\widehat{qr}}X_l),\quad\lambda>0.
\]

**引理2。** 有
\[
 \boxed{\sum_{i\in G_v}r_i^2|C_i(\lambda)|
 \le\frac{C_{n,R}}{\lambda}\sum_{k<l}|X_k\cap X_l|
 =O_{n,B,R}(v^6S/\lambda).} \tag{10}
\]
证明：轴平行盒的交体积精确分解为二维交面积乘 transverse 投影交量。对每个被计 pair，\(r_i^2\) 乘该 fibre 量不超过 \(|X_k\cap X_l|/\lambda\)。固定物理窗口内的尺度可比及半盒 packing，又使一个 pair 只会被纯 \(n,R\) 个根重复计入。取并的测度用求和上界，再用 (3)。证毕。

例如对固定 \(\tau>0\) 取 \(\lambda=c\tau v^3\)，坏 fibre 容量是 \(O_{\tau,R}(v^3S)=o(vS)\)。这是 full-dimensional 真实交体积的收费；并未声称全部混合 q3 交量是 \(o(v^7S)\)。

### resolved 双对角的确定代数

设某同一 fibre 中已经选出四胞 \(A,B,C,D\)，依次位于左下、右下、左上、右上，且面对端依该次序（近单位布局保证无包含投影）。记 \(x\) 方向端点为 \(a_x=A_x^+,b_x=B_x^-,c_x=C_x^+,d_x=D_x^-\)。两条对角的 \(x\) 交长若都至少为 \(\ell>0\)，则
\[
 (a_x-d_x)+(c_x-b_x)=(a_x-b_x)+(c_x-d_x)\ge2\ell.
\]
因此至少一个水平相邻面的侵入量为 \(\ge\ell\)。\(y\) 方向同理。若该相邻面对的另一平面轴交长至少为 \(cr_i\)，则其二维交面积至少为 \(c\ell r_i\)，所以本 fibre 被 (10) 的 \(\lambda=c\ell/r_i\) 收费。

当两个窄交长均为 \(\ge\tau v^3r_i\) 时，这给前述 \(O(v^3S/\tau)\) 容量量级。它允许小的 spurious planar face overlap，但不能让两条 resolved 对角无偿同时占据一个已确定四胞片。

## 6. 对 Slot 缺口的准确缩小

现在可严格支付三种误差：

* 若给定 fibre 的 helper 缺失确实造成一个固定面积 probe 未被覆盖，用 (8) 和 (5) 收费，损失 \(O((v^{n+1}+v^5)S)=o(vS)\)。
* 若一个已识别四胞片有两条 resolved 活跃对角，或同一个局部二维位置有明显的双重占据，用 (10) 收费，损失 \(O_{\tau,R}(v^3S)=o(vS)\)。
* 对任何已给定、位于固定窗口且局部重复有界的片分配，根权 \(R_e^2\) 向 helper 的 \(R_k^2\) 转移，由 (9) 收费为 \(O_R(v^6S)\)：每片面积至多 \(Cr_i^{n-2}\)，\(|R_e^2-R_k^2|\le Cr_i|R_e-R_k|\)，再求和。

**尚未完成的步骤。** 必须从这些实际盒切片中构造一个可测的四胞片分配，使每个 actual corner-face 位置最多用一次。尤其不能把“一个指定 helper 不存在”直接换成“probe 有固定面积 void”；可能有另一胞覆盖该 probe，需证明其应进入同一候选角列，或将它合法计入碰撞收费。不同锚点给同一角面片的选择也需证明相容。只有完成这些选择和单射，才能写出 resolved Slot；本稿没有这样做。

即使 resolved Slot 完成，\(\min(L_1,L_2)/(v^3R_e)\to0\) 的退化边仍有 raw 接触权；不能由其交体积小直接删除该权。这一尾部可能需要 R20 的完整 transverse 成本留下的正罚，或一个连续权的有效 Slot。R31 对混合三轴真实源额外成本的工作是另一接口，本稿不代其证明或预支其结论。

因此本稿的无条件新结论是 (5)、(8)、(9)、(10) 及确定四胞的 resolved 碰撞代数；高维 raw Slot、完整 resolved Slot 和下一尺度的真实局部价格均不在已证范围内。

## 7. 本机核验范围

实际读取 Terra 单次 Mathematica 15.0 ARM 运行的 `r029-computation/face_collision.wl/log`：上面的对角交长／面侵入恒等式残差为零；在 \(\ell>0\) 和两条对角重叠均至少 \(\ell\) 时，两个相邻面侵入都小于 \(\ell\) 的违反集合为 False。没有配置枚举、15窗试探或重跑旧大计算。

随机线段密度界、全输入 void、实际路径权重和 Fubini 收费均给出了解析证明；计算只核这一个通用线性代数核心。冻结之后若要处理 helper 片选择，将另开新命题，不静默扩张本稿结论。
