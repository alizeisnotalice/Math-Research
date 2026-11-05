# R054：二维小半径源的八阶贡献与真实两尺度孔候选

2026-11-05。新解析证明稿，待非作者独立审查；只写本独占 scratch。原源、论文和仓库未改。本稿不处理正常 packet 的精确八阶积分。

## 1. 结论和停止边界

目前二维 inside dyadic 成本只估到 \(O(v^8S)\)；本轮保留这一同阶费用，并构造真实两尺度新增八阶下界。下面给出真实两尺度有限输入构造，其微孔有自己的短严格见证，宏孔有保留共同实际目标质量的长严格见证；不是把标量盒模型的亏损当成真实极大函数现象。

若以下新几何和预算经独审通过，构造给出
\[
 C_2(1-t)\ge P_7(t)
       +\left(J_2+\frac8{27}\right)av^8-O(av^9),
 \qquad J_2=\frac{15694}{6561},
\tag{1}
\]
其中 \(a=v^2=4t\)、\(v=2u\)、\(u=\sqrt t\)，故候选新系数为
\(17638/6561\)。这里仅是实际有限族的一侧下界，不是最优系数声明。已发表 R041 从未声称其 \(J_2\) 最优，因而不需修正其旧结论。

另有一条独立的正预算细化：在八阶竞争输入中保留原图吸收后剩余的内部 \(P_I\) 正项，可得 \(P_I+P_\partial=O(v^8S)\)。所以二维边界小源和全部环境额外质量都为 \(O(v^{10}S)\)，深处盒外小源为 \(O(v^9S)\)。真正未被消去的八阶项只在 inside 的正常缺层部分。

## 2. 边界与环境可由原内部预留解决

采用 R049 同一新八层族及竞争分支。其八层共同账本在原图两项吸收之后，仍保留固定正的内部收入：按原
\[
 C_*P_\partial\le c_{\rm rec}P_I/4,
 \qquad 2nvS_I\le c_{\rm rec}P_I/4,
\]
并在全部常数固定后令 \(v\) 小，\((C_{\rm b}+bv)P_\partial\) 也在预定边界份额内。因此可以不丢掉该项，而把 R049 (12) 加强为
\[
 (1+Q)T-A_vM/a
 \ge c_P P_I+c_4\mathcal L-Cv^8S-o(v^8S),\qquad c_P>0.
\tag{2}
\]
这来自原真实内部恒等式的 \((1-a)\) 系数和共享质量／突出集回收，不从完整 \(J\) 再抽任何份额。正权平均后内部收入不被末层权削弱，因为每层都有同一 \(c_{\rm rec}P_I\)。

八阶竞争反向比值恒等式正如 R049 §5 给公共左侧 \(\le C_Bv^8S\)。所以
\[
 P_I+P_\partial=O(v^8S).
\tag{3}
\]
同时已有 \(\mathcal F=O(v^8S)\)。原 global 的精确边界小源式比 R046 粗写法更强：
\[
 (v/a)\sum_{\partial I}\delta_i\le C aP_\partial.
\]
二维因此为 \(O(v^{10}S)\)。环境式
\(\sum u_i\le Ca(\mathcal F+P_\partial)\) 同为 \(O(v^{10}S)\)。深处盒外成本 \(Cv\mathcal F=O(v^9S)\)，盒内扩张成本 \(Cv^2\mathcal F=O(v^{10}S)\)。这些结论保留各点原小半径见证，无需消除边界点本身。

对 inside 原格子计数可写更精确的缺集分解
\[
 \frac va\sum\delta_{i,\rm in}
 \le Cv\sum_i|Y_{i,N}\setminus T_i|
       +Cv\sum_i|\mathscr D_i\setminus Y_{i,N}|
       +Cv^2\mathcal F.
\tag{4}
\]
第二项由 R049 raw_N 为 \(O(av^8S)=O(v^{10}S)\)。第一项只知
\(|Y_N|=O(av^5S)\)，二维恰为 \(O(v^8S)\)。它是已经支付宏观禁层中的实际输入区域，不能把整个 \(|Y_N|\) 再作为未用 raw 价格收费。下面的微孔机制正使用这个区域。

## 3. 宏观基础对象：只取已审 R041 的二维集合

以下 \(u\to0\)，\(g=(1+2u)^{-1}\)、\(a=4u^2\)。固定 \(\kappa\) 在一个包含 \(16/3\) 的正紧区间；最终取 R041 的移动剪切
\[
 \kappa=\frac{16}{3}g^2+\frac{640}{81}g^3u.
\tag{5}
\]
令
\[
 T=\begin{pmatrix}1&e\\e&1\end{pmatrix},\quad
 e=\kappa u^3,\quad d=1-e^2,\quad
 C_w=[-w,w]^2,
\]
\[
 w=ug+\Lambda u^6,\quad \rho=1-2w=g-2\Lambda u^6,
 \quad M=a\rho^2.
\tag{6}
\]
\(M\) 从此固定，是每个宏观周期孔中**宏孔加微孔的总实际质量**。

对 \(x\in C_w\) 写 \(s=w-\|x\|_\infty\)、\(r(x)=\rho+s\)，并令
\[
 F(s)=a[(\rho+s)^2-\rho^2],\qquad f=F'(0)=2a\rho.
\tag{7}
\]
对删除集 \(V\subset C_w\) 定义 R041 的坏集算子
\[
 B(V)=\{s=0\}\cup
 \left\{x\in C_w:
       \sum_{k\in\mathbb Z^2\setminus\{0\}}
       |(Tk+(C_w\setminus V))\cap Q(x,\rho+s)|\ge F(s)\right\}.
\tag{8}
\]
令 \(U_0=B(\varnothing),U_1=B(U_0),U_2=B(U_1)\)。已审 R041 给
\[
 U_1\subset U_2\subset U_0,
 \qquad B(U_2)\subset U_2.
\tag{9}
\]
其证明和体积展开在所有有界 \(\Lambda\) 上统一，故本稿改变 \(\Lambda\) 的 \(u^2\) 项后仍可使用。

所有正交量只可能来自 \(k=(1,-1),(-1,1)\)；\(\|k\|_\infty\ge2\) 或不平衡符号指标不进入原见证。这仅用 \(r+w=1-\|x\|_\infty\)、\(e\ll w\) 和 \(T\) 的显式坐标。源点的两个相对角各有两种全局最浅轴排序，Jacobian 均为一；在排序坐标 \((s,y)\) 中，目标原始交帽恰为 \([0,e]\times[0,(e-y)_+]\) 或交换两轴，独立于 \(s\)。这对任意目标子集成立，不要求目标形状对称。

## 4. 在 U2 内真实造微孔

\(U_2\) 的坏层厚度在上述各 chart 中为连续高度
\(\min\{F^{-1}(R(y)),w-y\}\)，另有零测边界。故 \(U_2\) 是闭可测集，边界零测，且 \(|\operatorname{int}U_2|=|U_2|\)。其正体积主量为
\[
 |U_2|=\frac{a\kappa^3}{16g}u^5[1+O(u)].
\tag{10}
\]

固定每个小 \(u\)，取微观阈值 \(t'=t(1-u)<t\)、\(v'=2\sqrt{t'}\)、\(a'=4t'\)。微格间距记 \(h\)，孔半宽
\[
 w_\mu=\frac{h v'}{2(1+v')}.
\]
只保留微格中心 \(q\in h\mathbb Z^2\) 中满足 \(Q(q,h)\subset\operatorname{int}U_2\) 的那些，并置
\[
 H=\bigcup_q(q+[-w_\mu,w_\mu]^2).
\tag{11}
\]
集合有限且各微孔不交。对任意微孔点 \(z\) 指定
\[
 r_\mu(z)=h-w_\mu-\|z-q\|_\infty.
\tag{12}
\]
有 \(Q(z,r_\mu)\subset Q(q,h)\subset\operatorname{int}U_2\)。该见证含本微孔全部质量，和任何其他微孔无正体积交。其缺陷率不超过
\((w_\mu/(h-2w_\mu))^2=t'<t\)，故是原阈值 \(1-t\) 的严格见证。

体积分数为 \(\beta=a'/(1+v')^2=a[1+O(u)]\)。随着固定 \(u\) 下 \(h\to0\)，内部有效微格胞总体积趋于 \(|U_2|\)。可明确选第一个满足 \(h\le u^{20}\) 且有效胞体积至少 \((1-u)|U_2|\) 的 dyadic 间距；该选择终止。因此
\[
 |H|=a|U_2|[1+O(u)],\qquad |H|\le C a|U_2|.
\tag{13}
\]
微孔接口由 Astra 另独立写在 `scratch/r054-microcell-interface/`；本稿只需上述总体质量和所有短严格见证的 host 包含，不需假设每个宏观帽中的微孔质量恰等于密度乘帽体积。

## 5. 固定共同 M 后的宏观保护：不循环选孔

固定已经造好的 \(H\subset U_2\)，定义
\[
 B_H(V)=\{s=0\}\cup
 \left\{x\in C_w:
  \sum_{k\ne0}|(Tk+(C_w\setminus V))\cap Q(x,\rho+s)|
  +\sum_{k\ne0}|(Tk+H)\cap Q(x,\rho+s)|\ge F(s)\right\}.
\tag{14}
\]
这里 \(M,\rho,F\) 均保持 (6)–(7)，不由宏孔剩余质量重算。算子对 \(V\) 反单调。置
\[
 U=U_2\cup B_H(U_2).
\tag{15}
\]
因为 \(H\subset U_2\)，在计算 \(B_H(U_2)\) 时目标总集
\((C_w\setminus U_2)\cup H\subset C_w\)，两部分互不交，所以 \(B_H(U_2)\subset U_0\)。并且
\[
 B_H(U)\subset B_H(U_2)\subset U,
 \quad H\subset U_2\subset U\subset U_0.
\tag{16}
\]
这是一遍显式有限集合运算；既没有求未知孔形固定点，也没有先按宏孔质量设计后再直接加微孔。

### 新保护层的全局 L1 上界

在第 3 节的四个最浅轴 chart 中，\(B(U_2)\) 与 \(B_H(U_2)\) 的目标帽质量差为 \(|H\cap\mathcal P(y)|\)，且只在 \(0\le y\le e\) 可能为正。\(F'\ge f\)，逆函数以及真实纤维上截断均为 \(1/f\)-Lipschitz，所以
\[
 |B_H(U_2)\setminus B(U_2)|
 \le\frac1f\sum_{\text{四图}}\int_0^e
                          |H\cap\mathcal P(y)|\,dy
 \le\frac{4e}{f}|H|=O(u|H|).
\tag{17}
\]
这里可以用任意微孔分布；没有把微孔分布均匀性代替实际帽质量。由 \(B(U_2)\subset U_2\)，得
\[
 |U|\le |U_2|+C u|H|.
\tag{18}
\]
四图及 clip 已覆盖所有实际源点；源深 \(s\) 大于 \(C e^2/f=O(u^4)\) 时满外包帽也不足达到 \(F(s)\)，不会有未处理的深源坏层。

## 6. 预算的精确八阶改动

二维 R041 系数为
\[
 \Lambda_0=\frac{\kappa^3}{32g},\quad
 \Lambda_1=-\frac{5\kappa^4}{768g^2},\quad
 \Lambda_2=\frac{13\kappa^5}{8192g^3}.
\]
本稿选择
\[
 \boxed{\Lambda=\Lambda_0+\Lambda_1u+\Lambda_2u^2
                         -a\Lambda_0+c u^3,}
\tag{19}
\]
其中固定安全常数 \(c\) 最后在统一余项界之后取得足够大。

记 \(D=(2w)^2-M\)。可直接精确展开
\[
 D/a=2\Lambda u^5+\Lambda^2u^{10}(1-4u^2).
\tag{20}
\]
而已审 R041 的真实 \(U_2\) 体积在这个有界参数族上统一给
\[
 |U_2|/a=2\Lambda_0u^5+2\Lambda_1u^6
                          +2\Lambda_2u^7+O(u^8).
\tag{21}
\]
所以
\[
 (D-|U_2|)/a=-2a\Lambda_0u^5+2c u^8+O(u^8).
\]
另一方面 (13)、(21) 给
\[
 |H|/a=|U_2|[1+O(u)]
       =2a\Lambda_0u^5+O(u^8),
 \qquad u|H|/a=O(u^8).
\tag{22}
\]
因此由 (18)
\[
 \frac{D+|H|-|U|}{a}\ge 2c u^8-O(u^8)>0
\tag{23}
\]
（取固定 \(c\) 大，再取小 \(u\)）。这里安全项首次影响 \(av^9\)，不会改变八阶系数。所有余项常数在固定 \(\kappa,\Lambda\) 紧范围上统一；可先容纳 (19) 前四项再加一，选择该范围上的误差界，最后要求 \(c u^3\le1\)，没有安全常数循环。

由 (23)，在开集 \(C_w\setminus U\) 内先取紧子集再连续分位切片，可取质量恰为 \(M-|H|>0\) 的紧宏孔 \(A\)。可同时保留原点（增加一个零测点即可）。令
\[
 P=A\cup H,\qquad |P|=M.
\tag{24}
\]
两部分物理不交；总质量不是 \(M+|H|\)。

## 7. 所有宏／微点的实际原严格见证

对 \(x\in A\)，保留 \(r(x)=\rho+s\)。它含同周期整个 \(P\)，因 \(\rho>2w\)。所有其他周期孔满足
\(P\subset(C_w\setminus U)\cup H\)，故由 \(x\notin U\) 和 (16)，
\[
 \left|\bigcup_{k\in\mathbb Z^2}(Tk+P)\cap Q(x,r(x))\right|
 \le M+R_U(x)+R_H(x)<M+F(s)=a r(x)^2.
\tag{25}
\]
因此这些宏孔点确在阈值 \(1-t\) 的严格 halo，自己的完整目标质量始终是 \(M\)。

对 \(z\in H\)，其原短见证 (12) 完全位于 \(\operatorname{int}U_2\subset U\)，该区域没有任何宏孔 \(A\)。它至多包含本微孔、对其他微孔无正体积交，因此第 4 节的严格水平计算仍是实际输入的计算。后加宏保护层不破坏微见证；没有先造一个局部孤立微孔后忽略宏孔伸入它的见证。

上述实际有限下界只使用宏孔和微孔各自的原严格见证，不要求它们等于 R049 特定新八层贪心族的源或删除集合。更强的 inside 删除质量非 \(o(v^8S)\) 识别，须另作完整的同族定位；本稿不声明这一结论。

## 8. 有限测度、全部半径及系数

取 \(N=\lceil u^{-10}\rceil\)、凸外盒 \(\Omega_N=T([-N,N]^2)\)，定义
\[
 E_u=\Omega_N\setminus\bigcup_{k\in\mathbb Z^2}(Tk+P).
\tag{26}
\]
周期积分给 \(|E_u|=(2N)^2(d-M)>0\)。\(\|T^{-1}\|_{\infty\to\infty}<2\)，宏见证半径 \(<1\)、微见证更小，所以距外盒固定周期层以内的全部孔点保留上述严格见证。丢失边界孔总量 \(O(aN)\)，实际比值的损失为 \(O(a/N)\)。外盒外任一点的任意中心方体至多一半体积落在包含 \(\Omega_N\) 的分离半空间内，因此阈值 \(1-t>1/2\) 排除所有外部正半径见证。输入密度点几乎处处进入 halo。所以实际有限比值被两侧夹为
\[
 \frac{|\mathcal H_2(1-t;E_u)|}{|E_u|}
       =\frac d{d-M}+O(a/N).
\tag{27}
\]
这是全部正半径的真实比值；不是只估计所指定的宏／微见证子族。

与 R041 参数相同的旧 \(\Lambda_{\rm old}=\Lambda_0+\Lambda_1u+\Lambda_2u^2\) 相比，(19) 的主改动为 \(-a\Lambda_0\)。安全 \(cu^3\) 首次影响 \(av^9\)。相应
\[
 \Delta\rho=2a\Lambda_0u^6+O(u^9),
 \qquad\Delta M=4a^2g\Lambda_0u^6+O(av^9),
\]
并用真实比值差恒等式
\[
 \frac d{d-M_{\rm new}}-\frac d{d-M_{\rm old}}
 =\frac{d(M_{\rm new}-M_{\rm old})}
         {(d-M_{\rm new})(d-M_{\rm old})}.
\]
除 \(av^8\) 后极限为
\[
 \lim_{u\to0}\frac{g\Lambda_0}{16}
 =\frac{(16/3)^3}{512}=\frac8{27}=2K_2.
\tag{28}
\]
R041 的实际旧族展开和 \(a/N=O(a u^{10})\) 遂给 (1)。此结论待新两尺度几何与预算独审，不能仅由 (28) 的标量计算签署。

## 9. 计算和独审待办

已向 Terra 明确提交 (19) 的精确有理 \(\kappa,\Lambda\)、真实 \(d/(d-M)\) 差及候选分数；只运行这个最小新批次，未重跑 R041。我实际读取 `r054-computation/two_scales.wl/log` 与 `fraction_series.py/json`：Mathematica 15 ARM（主 PID 13682，退出 0）给三项 \(\Lambda\) 因子差均为 0，实际比值差至 \(u^{10}\) 为 \(8192u^{10}/27\)，除 \(av^8=1024u^{10}\) 为 \(8/27\)；安全 \(cu^3\) 对 \(u^{10}\) 系数无影响。独立 Python Fraction 形式幂级数对 \(c=0,1,7\) 全部精确吻合，旧／新有理系数分别为 \(15694/6561\)、\(17638/6561\)。这是标量和系数证据，不能代替本稿的新坏预算、宏／微严格见证或有限化。

必要独审点：U2 内微格所有短严格见证与总体密度；BH 的任意实际目标帽及全最浅 chart 的 \(O(u)\) 上界；共同 M 下的可测精确选孔；\(\Lambda\) 改动后 R041 统一余项和所有 small-u 守卫；全部半径的有限化；原内部 \(P_I\) 预留细化。尤其不能把微孔质量直接加到旧 M 后仍沿用旧根，也不能把已经付正常禁层的完整 J 再计为可用正收入。
