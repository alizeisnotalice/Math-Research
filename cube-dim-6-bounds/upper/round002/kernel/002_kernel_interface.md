# 第002轮：六节点截时核、sharp tube 与 face quotient 审计

本轮只建立连续认证的分析接口；没有 cone 覆盖、\(g_0\)、\(I_0\) 或 \(C_6\)
新上界。

## 六节点截时表达

在任一满维 Fourier cone 中，写
\(\xi=\sum_{j=1}^6u_jr_j\)、\(u_j\ge0\)、
\(\psi(\xi)=\sum q_ju_j\)，并设
\(b_j=x\cdot r_j/q_j\)、\(\alpha=|\det(r_1,\ldots,r_6)|/\prod q_j\)。
Fourier 反演给
\[
p_t(x)=(2\pi)^{-6}\sum_R\alpha_R\prod_{j=1}^6(t-i b_{R,j})^{-1}.
\]
完整实和的截时势为
\[
 H_\epsilon(x)=(2\pi)^{-6}\sum_R\alpha_R
 \operatorname{DD}_5\!\left(\arctan(\cdot/\epsilon);b_{R,1},\ldots,b_{R,6}\right). \tag{A}
\]
这是由
\[
 \int_\epsilon^\infty\Re\prod_{j=1}^6(t-i b_j)^{-1}dt
 =\operatorname{DD}_5\!\left(\arctan(\cdot/\epsilon);b_1,\ldots,b_6\right)
\]
得到的。若 \(g=(\epsilon-i\cdot)^{-1}\)，则
\(g^{[5]}=i^5\prod_jg(b_j)\)，故右边对 \(\epsilon\) 的导数恰为负的积分
被积函数。

重复节点由 confluent Newton 表处理：若首末节点重合，k 阶差商取
\(F_\epsilon^{(k)}(b)/k!\)。\(\epsilon>0\) 使零节点和任何重数都有限。
在全部节点不同且非零的开层，\(\epsilon\downarrow0\) 给
\[
 G(x)={1\over128\pi^5}\sum_R\alpha_R\operatorname{DD}_5
 (\operatorname{sgn};b_{R,1},\ldots,b_{R,6}). \tag{B}
\]
零/重节点不能逐 cone 使用 (B)：它们可产生发散 principal parts，必须先以
(A) 对完整 cone 和进行处理。

实际校准覆盖6的11种重数分割（含零节点），并有全零
\(1/(5\epsilon^5)\) 与六重节点2的
\(\Re[1/(5(\epsilon-2i)^5)]\) 两个解析锚点。90位直接积分的最大显示残差低于
\(5\cdot10^{-93}\)，但这是实现校准，不是 cone 证书。

## sharp body-frame 常数

对于 \(V=\{(1,s):s\in\{-1,1\}^5\}\)，\(\psi\) 是无条件、置换不变凸函数。
先取绝对值，再对坐标置换平均，得到
\[
 \psi(\xi)\ge\psi(a,\ldots,a),\qquad a=\|\xi\|_1/6.
\]
32个代表方向上的绝对符号和在全相等点为60，故
\[
 \psi(\xi)\ge10\|\xi\|_1. \tag{C}
\]
等号在六坐标相等时成立。因而
\[
p_t(0)\le(2\pi)^{-6}\prod_{j=1}^6\int_\mathbb R e^{-10t|u|}du
 ={1\over10^6\pi^6t^6}. \tag{D}
\]

在 face \(x_1=1\) 的自由变量 \(y=(x_2,\ldots,x_6)\) 中，还有锐界
\[
 \psi(\xi_1,y)\ge12\sum_{j=2}^6|y_j|. \tag{E}
\]
确实 \(\xi_1\mapsto\psi(\xi_1,y)\) 偶凸，故可降至零；对后五坐标置换平均，
而五个 Rademacher 和的绝对值总和为60。这在自由坐标相等时取等。

若 \(\max_j|\Im z_j|=s<12\epsilon\)，(E) 和 (D) 给
\[
 |H_\epsilon(1,z)|\le {1\over5\cdot10^6\pi^6(\epsilon-s/12)^5}. \tag{F}
\]
这比第001轮的 determinant-160 热核常数严格更强。

## quotient face 上界

固定一个 active body 方向 \(v_i\)，令 \(Y=\sum_{j\ne i}v_jZ_j\)，并定义
\(QY=Y_{2:6}-v_{i,2:6}Y_1\)。在坐标反射后可取 \(v_{i,2:6}=(1,\ldots,1)\)，
则31个 quotient 方向是 \(w=s-\mathbf1\), \(s\ne\mathbf1\)。它们张成
五维，且在重标号 \(r=-s\) 后
\[
 \psi_Q(\zeta)=\sum_{r\in\{-1,1\}^5}|\mathbf1\cdot\zeta+r\cdot\zeta|
 \ge\sum_r|r\cdot\zeta|\ge12\|\zeta\|_1. \tag{G}
\]
第一不等式按 \(r,-r\) 配对。故 quotient 密度 \(q_t\) 满足
\[
 \|q_t\|_\infty\le(2\pi)^{-5}(6t)^{-5}={1\over12^5\pi^5t^5}. \tag{H}
\]

对开 slab \(x_1>1-\eta\)，fiber 恒等式 \(x_1=\sum_iZ_i\) 至少给一个
\(|Z_i|\ge(1-\eta)/32\)。在卷积中以 \((Y_1,QY)\) 代换，其 Jacobian 为1；
对 \(Y_1\) 积分恰为 quotient 边际。于是 Tonelli 给 a.e.
\[
 p_t(x)\le {32^3\over12^5\pi^6(1-\eta)^2}t^{-4}.
\]
由 \(p_t\) 连续性推广至整个 slab，再令 \(\eta\downarrow0\)：
\[
 \|p_t(1,\cdot)\|_\infty\le {32\over243\pi^6}t^{-4}. \tag{I}
\]
第一坐标的尺度为 \(32t\)。其边际等式先 a.e. 成立；以区间平均并对连续的
非负密度使用 Fatou，得到固定 slice 的充分估计
\[
 \|p_t(1,\cdot)\|_1\le {32t\over\pi}. \tag{J}
\]
插值、正象限反射对称和 Minkowski 于是给
\[
 \left\|\int_0^\epsilon p_t(1,\cdot)dt\right\|_{L^{6/5}([0,1]^5)}
 \le {6\over7}\left({32\over243}\right)^{1/6}\pi^{-11/6}\epsilon^{7/6}. \tag{K}
\]

## Chebyshev 账与奇异几何

对半边长 \(h\) 的五维 box，用次数 \(d\) 的 tensor Chebyshev roots，取
\(s=h(\rho-\rho^{-1})/2<12\epsilon\)，令 \(M\) 为 (F) 的右端，
\(E_1=4M\rho^{-d}/(\rho-1)\)、\(\Lambda_d\le2d+1\)。逐坐标 telescoping
给严格解析误差
\[
 E_{an}=E_1(1+\Lambda_d+\Lambda_d^2+\Lambda_d^3+\Lambda_d^4).
\]
物理节点位移至多 \(\delta\) 时，任选 \(0<R<12\epsilon\)，
\(D=M_\epsilon(R)/R\)，安全附加误差为 \(5\Lambda_d^5D\delta\)。对于最终
\(L\le H_\epsilon\le U\)，以 Chebyshev 平均权
\(\mu_0=1\)、\(\mu_{2k+1}=0\)、\(\mu_{2k}=-1/(4k^2-1)\) 算平均，再对
\(u^{6/5}\) 用上弦 Jensen 界。短时项必须以 (K) 在范数层面加入。

body 的 \(S_2\) 还可准确描述为“坐标绝对值至多有两种取值”。因此排序正
face 中它是五条线段 \((1^k,a^{6-k})\), \(k=1,\ldots,5\), \(0\le a\le1\)。
所有0/1 binary 顶点奇异；但 \((1,a^k,0^{5-k})\) 在 \(0<a<1\) 时通常有三种
绝对值，不属于 \(S_2\)。这避免把非奇异三层锚点误删。

下一批可从低度全 tensor pilot 开始；若不通过，只按已认证误差二分 leaf。未经
mixed-derivative 证明的稀疏网格不得替代上述完整误差账。
