# R035：resolved 双角的实际可测四份容量

日期：2026-11-05。完整候选，待独立数学审查。仅写独立 scratch，不改 R29 或仓库。

本稿以已独审 R29 的真实 void、固定窗口 root 传输和切片碰撞为准确输入，完成 resolved 角的可测 helper 选择、角片单射及所有重数。只处理两个交长均不小于 \(\tau v^3R_e\) 的边，\(\tau>0\) 固定；**不声称 raw 退化 Slot 已证明**，不使用 R31/R34 的真实混合角成本。

## 1. 命题及输入

固定 \(n\ge3\)、\(a=v^n\)、\(A_v=(1+v)^n\)。使用同一临界粗先验末层，\(V=|E|\)、\(M=|K|\)、\(S=\sum_i r_i^n\)。实际 CDF 盒、目标根为
\[
 X_i=x_i+W_i+[-R_i/2,R_i/2]^n,
 \quad R_i=(g_i/a)^{1/n},\quad U=\bigcup_{i\notin I}X_i,
 \quad\mathcal Q=\bigcup_i Q_i.
\]
令 \(c_i\) 为 \(X_i\) 的实际中心。原强分离、\(X_i\subset Q_i\)、相交度有界及 \(R_i/r_i=1+O(v)\) 均保持。

R19 的同一坏根 \(Z_v\) 满足 \(S(Z_v)=O(v^{9/2}S)\)，好根集为 \(G_v=Z_v^c\)。对每个固定窗口，好根的全部半盒清单一致匹配单位铺砌，相关源全部深处。这里仍不假设 LT 有 \(o(v)\) 速率。

\(\mathscr A\) 表示恰两临界轴、真实交体积正的无向 pair。给每边任选
\[
 \min(R_i,R_j)\le R_e\le\max(R_i,R_j).
\]
两临界交长为 \(L_1,L_2\)，其余轴投影交为 \(T_e\)，\(b_e=|T_e|\)。对 \(0<\tau\le1\) 取
\[
 \mathscr A_\tau=\{e\in\mathscr A:\min(L_1,L_2)\ge\tau v^3R_e\}.
\]
不对 \(\max L_l/(v^3R_e)\) 加上界。

胞 \(i\) 的轴对 \(q<r\) 的四个实际角面位置记为 \((i,q,r,\sigma_q,\sigma_r)\)，\(\sigma_q,\sigma_r\in\{\pm\}\)。其横截面容量为
\[
 F_{i,qr}=\pi_{\widehat{qr}}X_i,
 \qquad |F_{i,qr}|=\prod_{l\ne q,r}(R_i+w_{il}),
\]
其中 \(w_{il}\) 是 \(W_i\) 的边长。

**主命题。** 对每个固定 \(\tau>0\)，有实际可测片分配及下列准确容量界：
\[
 \boxed{4\sum_{e\in\mathscr A_\tau}b_eR_e^2
 \le\sum_{i\notin I}\sum_{q<r}\sum_{\sigma_q,\sigma_r}
     R_i^2|F_{i,qr}|+O_{n,B,\tau}(v^3S).} \tag{1}
\]
从而
\[
 \boxed{\sum_{e\in\mathscr A_\tau}b_eR_e^2
 \le\frac{\binom n2}{(1+v)^2}(V+M)+O_{n,B,\tau}(v^3S).} \tag{2}
\]
误差对固定 \(\tau\) 是 \(o(vS)\)。常数中的 \(B\) 来自固定临界下界 \(M/V\ge a/(A_v-a)-Bav^6\)。

明确使用 R29 的三个量化接口（固定窗口取足够大的纯常数，例如 \(R_*=20\)）：
\[
 |\mathcal Q\setminus U|=O((v^{n+1}+v^5)S), \tag{3}
\]
\[
 \sum_{i\in G_v}\sum_{k\in\mathcal I_i(R_*)}
       |R_i-R_k|r_i^{n-1}=O(v^6S),
 \quad \mathcal I_i(R_*)=\{k:B_k\cap Q(x_i,R_*r_i)\ne\varnothing\}, \tag{4}
\]
以及对固定轴对，定义 \(C_i^{qr}(\lambda)\) 为同一清单中某个不同胞对二维投影交面积至少为 \(\lambda r_i^2\)、同时其 transverse 投影包含该 fibre 的集合，则
\[
 \sum_{i\in G_v}r_i^2|C_i^{qr}(\lambda)|=O(v^6S/\lambda). \tag{5}
\]
清单量有固定上界，所有胞尺度与根尺度比一致趋一。原 \(Q_i\) 重数至多为 \(H=2^n\)。

## 2. 固定实际 probe 与坏 fibre

固定一个很小的绝对数 \(\epsilon>0\)，例如 \(10^{-6}\)。对一条两端均在 \(G_v\) 的边，固定其较小编号端点 \(i\) 为锚点。取两临界轴 \(q,r\)，反射它们，使另一端 \(j\) 在右上。以实际中心 \(c_i\) 为平面原点、\(r_i\) 为尺度；其他坐标不反射。

该平面中 \(X_i\) 是中心零、边长 \(1+O(v)\) 的矩形，\(X_j\) 中心为 \((1,1)+O(v)\)、边长也为 \(1+O(v)\)。放置两个固定 probe
\[
 P_B=[1/2+\epsilon,1/2+2\epsilon]
       \times[1/2-2\epsilon,1/2-\epsilon],\qquad P_C=P_B^{\rm transpose}. \tag{6}
\]
将它们按 \(r_i\) 缩放并移回 \(c_i\)。其面积为 \(\epsilon^2r_i^2\)，对小 \(v\) 全在 \(Q_i\) 的平面投影内。所有轴对与四个符号下，每根总 probe 数恰为 \(8\binom n2\)，是统一维数常数。

对每个 probe 记
\[
 V_{i,P}=\{y\in\pi_{\widehat{qr}}X_i:
              |P\setminus U_y|\ge|P|/2\}.
\]
R29 的 probe 收费（固定上述统一 probe 数）及 (3) 给
\[
 \sum_{i,P}r_i^2|V_{i,P}|
       =O((v^{n+1}+v^5)S). \tag{7}
\]
置 \(\lambda=\tau v^3/1000\)。对锚点为 \(i\) 的该边，删去 transverse 集
\[
 T_e\cap\left(V_{i,P_B}\cup V_{i,P_C}\cup C_i^{qr}(\lambda)\right),
\]
剩余集记为 \(\Gamma_e\)。因为每锚点边数有界、\(R_e/r_i\) 可比，由 (5)、(7)
\[
 \sum_{e\text{ 两端好}}R_e^2|T_e\setminus\Gamma_e|
 \le C\left(v^{n+1}+v^5+\frac{v^3}{\tau}\right)S. \tag{8}
\]
其中 (5) 对有限个轴对求和。边的一端若不在 \(G_v\)，其全部接触权不超过 \(C S(Z_v)=O(v^{9/2}S)\)；只用相交度有界及尺度可比，不物理删基线。

所有这些集均为可测，实际为有限盒的有限 Boolean 组合：有限 transverse 端点分区上 \(U_y\) 中的矩形清单固定，probe 的被覆盖面积固定，故 \(V_{i,P}\) 是 Borel 集。忽略有限个 transverse 边界超平面只损零测集。

## 3. clean fibre 的 helper 存在、唯一及可测性

固定 \(y\in\Gamma_e\)。两端胞 \(i,j\) 在该 fibre 中出现。probe 至少一半被 \(U_y\) 覆盖，因此存在另一个实际胞矩形覆盖 probe 的正面积。

先核完整候选清单：任何这样的 \(X_k\) 都在点 \((z,y)\in Q_i\) 与 \(Q_i\) 相遇，故 \(Q_k\cap Q_i\ne\varnothing\)。好根距标记集的图距离趋于无穷，所以 \(k\) 规则且深处；规则邻尺度比固定有界，中心与半盒位于根的一个固定物理窗口。因此 \(k\in\mathcal I_i(R_*)\)。LT 给这份清单中全部矩形的归一化边长为 \(1+o(1)\)，一致成立。

以下有限平面事实把“probe 被别的胞覆盖”转成正确 helper，而不把它直接当 void。令共同边长／中心误差为 \(0\le\delta\le\epsilon/100\)。覆盖 \(P_B\) 中一点的候选胞中心 \((b_x,b_y)\) 满足
\[
 \epsilon-\delta\le b_x\le1+2\epsilon+\delta,
 \quad-2\epsilon-\delta\le b_y\le1-\epsilon+\delta.
\]
所有盒均含其中心的边长 \(1-\delta\) 核心；右上端点胞的中心距 \((1,1)\) 也至多 \(\delta\)。若 \(b_x\le1-4\epsilon\)，候选与左下端点的二维交面积至少
\[
 (4\epsilon-\delta)(\epsilon-2\delta)>2\epsilon^2.
\]
若 \(b_y\ge4\epsilon\)，候选与右上端点的二维交面积至少
\[
 (\epsilon-3\delta)(4\epsilon-2\delta)>2\epsilon^2.
\]
这些下界只用最坏的核心区间交长。对小 \(v\)，\(\lambda<\epsilon^2\)，且两端与候选都在清单并包含 \(y\)；所以以上两种情况都违背 \(y\notin C_i^{qr}(\lambda)\)。因此实际 helper 的平面中心满足
\[
 |b_x-1|\le4\epsilon,\qquad |b_y|\le4\epsilon. \tag{9}
\]
其不等于端点，因为两个 probe 都在端点投影之外并留有固定裕量。交换两轴得到另一个 helper，中心在 \((0,1)\) 的 \(4\epsilon\) 邻域。两个 helper 也不相同。

同一 helper 列不可能有两个不同胞同时出现：它们的平面中心相距至多 \(8\epsilon\)，边长为 \(1+o(1)\)，所以二维交面积大于 \(r_i^2/2\)，远大于 \(\lambda r_i^2\)。这又会触发 \(C_i^{qr}(\lambda)\)。故每个 helper 唯一。

可测选择现在是有限清单选择：helper 是唯一的 \(k\in\mathcal I_i(R_*)\)，其 transverse 投影含 \(y\)，且其平面投影与相应 probe 有正面积交。所有条件都是有限 Borel 判据。于是 \(\Gamma_e\) 可按两个 helper 指标分成有限个 Borel 片，不涉及不可测的逐点选择。

## 4. 四份实际角面片与跨锚点单射

在每个 Borel 片上写四胞为
\[
 \begin{matrix}C&D=j\\ A=i&B\end{matrix}.
\]
它们在实际平面中的中心分别位于 \((0,0),(1,0),(0,1),(1,1)\) 的固定小邻域，边长为 \(1+o(1)\)。把同一 transverse 点 \(y\) 分别送到 A 的右上、D 的左下、B 的左上、C 的右下角面位置。这四个位置都是实际 \(X\) 的角面；四胞都包含该 transverse 点，所以并未将某片送入 helper 容量之外。

**单射。** 假设同一胞的同一轴对、同一符号角面位置及同一 \(y\)，被两条 clean resolved 边使用。以该共享胞为参考，两组四胞的另三个平面位置都落在相同三个近单位列中；共享角符号固定了该参考胞在四胞中的位置。实际半径及尺度在各自完整窗口比一致趋一。由 (9) 平移到共享胞并改用其尺度后，另一组三中心与同三个整数位置的误差可统一取为小于 \(100\epsilon\)；因此两个同列矩形的平面中心差小于 \(200\epsilon\) 倍胞边长。

第一锚点的固定清单包含第二组全部胞。理由是：每组平面坐标距各自锚点为固定倍尺度；所有胞的 transverse 投影都包含同一个 \(y\)，故 transverse 中心距也至多固定倍尺度；两组通过共享胞比较尺度。每组胞中心距自身锚点在平面上至多 \((2+o(1))r_i\)，在 transverse 上至多 \((1+o(1))r_i\)；通过共享胞后，第二组中心到第一锚点至多 \((6+o(1))r_i\)，其原半盒半径为 \((1/2+o(1))r_i\)。所以 \(R_*=20\) 足以涵盖这些胞及其半盒。它们的源均有效且尺度与第一锚点可比。

如果两组的同列指标不同，则这两个胞在该 fibre 的二维投影交面积大于第一锚尺度的固定比例，例如 \(r_i^2/4\)。它们都在第一清单中并出现于 \(y\)，与第一边 clean 矛盾。因此两组四胞指标完全相同。两条不同边只能是这四胞的两条相反对角。

再用实际端点而不是模型端点。记四胞相关 \(x\) 面端为 \(a_x=A_x^+,b_x=B_x^-,c_x=C_x^+,d_x=D_x^-\)。近单位列保证相反端点有序且不存在包含投影，故
\[
 L_{AD,x}+L_{BC,x}=(a_x-b_x)+(c_x-d_x). \tag{10}
\]
若两条对角都 resolved，左边至少为 \(\tau v^3\) 乘两条边的实际根之和，故对小 \(v\) 至少为 \((3/2)\tau v^3r_i\)。于是至少一个右侧面侵入量为 \(\ge(3/4)\tau v^3r_i\)。这个相邻面的另一平面轴交长大于 \(r_i/2\)，所以它的二维交面积大于 \(\tau v^3r_i^2/4\)，严格大于 \(\lambda r_i^2\)。四胞在该 fibre 都出现，故触发第一锚的 \(C_i^{qr}(\lambda)\)，矛盾。

所以每个实际角面位置在几乎每个 transverse 点最多使用一次。此论证比较了不同锚点的完整清单和真实面积，并非只证明某个预先给定四胞内两对角不共存。

## 5. 重数、根权传输与容量求和

以 helper 指标分片后，每条 clean 边的同一 transverse 片被分配四次，且角面单射给实际根权容量
\[
 \sum_{\text{全部分配片到胞 }c}R_c^2|\text{片}|
 \le\sum_{c\notin I}\sum_{q<r}\sum_{\sigma_q,\sigma_r}R_c^2|F_{c,qr}|. \tag{11}
\]
实际边权 \(R_e\) 在端点区间内，所以对每个分配胞 \(k\in\mathcal I_i(R_*)\)，
\[
 |R_e^2-R_k^2|\le Cr_i(|R_i-R_j|+|R_i-R_k|).
\]
每片面积至多为 \(Cr_i^{n-2}\)。每锚点只有固定数相交边，清单有限且 helper 分片是该边 \(\Gamma_e\) 的互不交划分；因此所有权重误差至多为
\[
 C\sum_{i\in G_v}\sum_{k\in\mathcal I_i(R_*)}|R_i-R_k|r_i^{n-1}
 =O(v^6S), \tag{12}
\]
其中端点 \(j\) 也在清单内。这里不把每个片的误差重复乘以其分片数；按同一 helper 求并后面积仍不超过该边的总 transverse 面积。R29 已证明固定窗口实际路径／根重复次数有界，(4) 精确支付 (12)。

结合 (11)、(12)
\[
 4\sum_{e\text{ 两端好}}R_e^2|\Gamma_e|
 \le\sum_{c,qr,\pm,\pm}R_c^2|F_{c,qr}|+O(v^6S).
\]
再加回 (8) 及坏端点的 \(O(v^{9/2}S)\) 接触权，得到
\[
 4\sum_{e\in\mathscr A_\tau}b_eR_e^2
 \le\mathrm{Cap}+C\left(\frac{v^3}{\tau}+v^{n+1}+v^5+v^{9/2}+v^6\right)S.
\]
对 \(n\ge3\) 及固定 \(\tau\) 即是 (1)。这包括非覆盖 probe、伪面碰撞、跨锚点选择和权重转移的全部费用。

## 6. 容量到实际环境体积

R17 的末层源宽控制给
\[
 \sum_i r_i^n\max_l|w_{il}-vR_i|/r_i=O(v^{7/2}S).
\]
因此逐个乘积展开，每个固定轴对满足
\[
 \sum_iR_i^2\prod_{l\ne q,r}(R_i+w_{il})
 =(1+v)^{n-2}\sum_iR_i^n+O(v^{7/2}S).
\]
同一末层的精确质量及目标重复恒等式给
\[
 \sum_iR_i^n=(M'+\mathcal D)/a
 \le M/a+O(v^6S)
 \le(V+M)/A_v+O(v^6S).
\]
最后一方向使用已证 \(V+M-A_vM/a=O(v^6S)\)。将它们代入 (1)，四个符号消去系数四，得到 (2)。没有近似 \(aS=M\)，也没有把原 \(r_e\) 替作实际 \(R_e\)。

## 7. 范围及后续剩余

结论只覆盖 \(\min(L_1,L_2)\ge\tau v^3R_e\)。这些交长可以很大，若某边平面交面积达到 \(\lambda r_i^2\)，其 transverse 容量直接进入碰撞收费，并不造成遗漏。

阈值 \(\tau\) 固定后，误差为 \(o(vS)\)。让 \(\tau\downarrow0\) 不能直接得到 raw Slot：存在最小交长随 \(v\) 退化的边，其 raw 接触权不由交体积控制。本稿不将两极限交换，也不宣布所有退化边已被删除。完整 transverse 源成本的正罚或连续权容量仍是另一命题。

实际三轴混合成本、七阶目标退价及 Tauberian 下一系数不是本稿输入或结论。这里完成的只是 resolved 双角的实际四份几何容量。

## 8. 本机核验

实际读取 `r035-computation/helper_probe.wl/log`。Mathematica 15.0 在 \(0<\epsilon<1/100,0\le\delta\le\epsilon/100\) 下核得两项 probe 核心面积下界均严格大于 \(2\epsilon^2\)，违反集为 False；没有模型枚举。实际中心区间、覆盖清单、Borel helper、跨锚点单射及全部费用由以上证明承担。

四胞 (10) 的线性恒等式及“双对角 resolved 迫面侵入”的本机证书已在 R29 实跑，本稿复用其解析恒等式，不重跑旧计算。所有原冻稿保持不变。
