# 既有周期乘积来源在任意平滑宽度下的固定窗口上界

2026-10-29。本文将第28轮的固定宽度结论改进为对全部有限 \(w\ge4\) 一致的结论。证明依赖原确定来源的系数和交叉引理；来源原样保存在 `../028/source_interface_baseline.tex`，而非假设任意来源具有这些性质。下文的概率推导、局部增量和正背景比较均给出证明。

## 1. 对象与定理的精确范围

固定 \(0<p<1\)。对每个充分大的整数 \(n\)，沿用原来源的确定标签、每格一个等质量原子 \(a_j\in[j/K,(j+1)/K)\)、谱截断 \(M\) 和参数
\[
 \varepsilon=p/K,\qquad v=\frac{9p^2n}{29K^2},\qquad
 K\asymp_p\sqrt{n/\log n},\quad M=\lceil\sqrt v\rceil<K.
\]
原标签与 \(K,M,v\) 均不随下面的 \(w\) 改变。对任意有限 \(w\ge4\)，置 \(\delta=(wKv)^{-1}\)。圆周 \(\mathbb T\) 上的概率来源是
\[
 \nu_\delta(dx)=\left[\rho_{\rm bg}(x)
       +\frac pK\sum_{j=0}^{K-1}\beta_\delta(x-a_j)\right]dx,
\quad \rho_{\rm bg}=1-p+(p/K)D_M,
\quad \beta_\delta=\delta^{-1}\mathbf1_{[-\delta/2,\delta/2]},
\]
其中 \(D_M=\sum_{m\le M}\pi m c_m\cos(2\pi mx)\)，\(c_m=\sqrt{2b_m/3}\)，
\[
 b_m=4\frac{729}{29}\left[
  \left(\frac{\sin(10\pi m/27)}{\pi m}\right)^4
  +3\left(\frac{\sin(\pi m/9)}{\pi m}\right)^4\right].
\]
来源在实线上周期提升；长度超过1时按实际覆盖次数积分。对闭中心区间 \(W_L(x)=[x-L/2,x+L/2]\)，令
\[
 Y_{\delta,L}(x)=\nu_\delta(W_L(x))/L,
 \quad Z_{\delta,L}(\mathbf x)=\sum_{i=1}^n\log Y_{\delta,L}(x_i),
\]
并以 \(P\) 表示圆周的归一化Lebesgue测度。定义
\[
 \mathcal W_\delta(J)=\sup_{h\in\mathbb R}e^h
 P^{\otimes n}\{\sup_{L\in J}Z_{\delta,L}>h\}.
\]
所有坐标共用同一个 \(L\)。另设辅助原子来源
\(\nu_0=\rho_{\rm bg}dx+(p/K)\sum_j\delta_{a_j}\)，并相应定义 \(\mathcal W_0\)。

**定理。** 存在仅依赖固定 \(p\) 及本来源固定常数的 \(C_p,n_0(p)\)，使每个整数 \(n\ge n_0(p)\) 都满足
\[
 \boxed{\quad
 \sup_{4\le w<\infty}\mathcal W_{1/(wKv)}([10/9,11/9])
       \le C_p\sqrt n,
 \qquad \mathcal W_0([10/9,5/4])\le C_p\sqrt n .\quad}       \tag{T}
\]
这里 \(\sup_w\) 在概率泛函之外；不是同一事件中同时最大化不同输入。定理控制该来源在指定尺度窗口的弱尾，既不是整个尺度轴的结论，也不是任意输入的 \(C_n\) 上界。原子来源只是比较工具，不称为 \(L^1\) 测试函数。

## 2. 共同来源资格与较大窗口的粗方差

后文先在 \(J_2=[10/9,4/3]\) 工作。原来源给出
\[
 Y_{\delta,L}=1+\varepsilon F_{\delta,L}/L,
 \quad F_{\delta,L}=H_{\delta,L}+B_L,
 \quad B_L=\sum_{m\le M}c_m\sin(\pi mL)\cos(2\pi mx),
\]
其中 \(H_{\delta,L}=\mathcal P_\delta(x+L/2)-\mathcal P_\delta(x-L/2)\)，
\(\mathcal P_\delta'=\sum_j\beta_\delta(\cdot-a_j)-K\)。
\(|H|\le2\)、\(c_m\le20/(9m^2)\)，因此
\[
 |F|\le154/27,\quad |Y-1|\le(77/15)\varepsilon,
 \quad \int Y_{\delta,L}\,dP=1.
\]
共同大 \(n\) 门槛后 \(1/2\le Y\le3/2\)。本来源的全局交叉引理对全部 \(L,J\) 及允许的 \(\delta\) 给出
\[
 |\langle H_{\delta,L},B_J\rangle+
          \langle B_L,H_{\delta,J}\rangle|\le C_x,
 \quad C_x=O(K^{-1/2})+O((\log M)/K)\to0.
\]
此引理及确定标签资格是原来源输入，完整证明在第28轮归档的基线中。由于 \(\min_{L\in J_2}|\sin(\pi L)|=\sin(\pi/9)\)，正交性给
\[
 \int F_{\delta,L}^2\,dP
 \ge (b_1/3)\sin^2(\pi/9)-C_x
 \ge g_0:=\frac{b_1}{6}\sin^2(\pi/9)>0.                    \tag{2.1}
\]
最后一步只要求 \(C_x\le g_0\)，不要求把原安全尺度的平坦方差公式推广至 \(J_2\)。

设 \(\xi_L=\log Y_{\delta,L}\)、\(dQ_L=Y_{\delta,L}dP\)。独立副本公式与对数均值定理给
\[
 \operatorname{Var}_{Q_L}\xi_L
 \ge\frac19\operatorname{Var}_P Y_{\delta,L}
 \ge\frac{g_0}{16}\varepsilon^2,
 \qquad |\xi_L|\le C\varepsilon.                           \tag{2.2}
\]
其中 \(1/9\) 来自 \((\inf Y/\sup Y)^2\)，再用 \(L\le4/3\)。常数对 \(w\) 一致。

另由 \(\|D_M\|_\infty\le(20/3)(1+\log M)\)，固定 \(p<1\) 时最终
\[
 \rho_{\rm bg}\ge (1-p)/2>0.                               \tag{2.3}
\]
原子与平滑来源每周期质量均为1。既有标签、\(M<K\)、箱峰不重叠、\(K\delta\le1\)、(2.1)--(2.3)、\(n\ge2\)、\(v\to\infty\) 等资格可一次选择大 \(n\) 门槛，对全部 \(w\ge4\) 有效。没有要求增长 \(w\) 时来源密度高度仍为 \(O(\log n)\)。

## 3. 宽度1的小球引理

若 \(|\xi|\le C\varepsilon\)、\(\operatorname{Var}\xi\ge c\varepsilon^2\)，令 \(S_m\) 为独立同分布的 \(m\) 项和，则在充分小 \(\varepsilon\) 下
\[
 \sup_t\Pr\{S_m\in[t,t+1]\}\le C/(\varepsilon\sqrt m).     \tag{3.1}
\]
证明：取独立副本 \(\xi'\)。当 \(|s|\le1\) 时，最终 \(|s(\xi-\xi')|\le1\)，故
\[
 |\phi(s)|^2=E\cos(s(\xi-\xi'))
 \le1-\tfrac12s^2\operatorname{Var}\xi,
 \quad |\phi(s)|^m\le e^{-cm\varepsilon^2s^2}.
\]
使用
\[
 k(u)=\left(\frac{\sin(u/2)}{u/2}\right)^2
     =\int_{-1}^1(1-|s|)e^{isu}\,ds,
\]
并注意 \(k\) 在 \([-1/2,1/2]\) 有正下界。由Fubini及高斯积分，
\(\Pr\{S_m\in[t,t+1]\}\le C\int_{-1}^1|\phi(s)|^m ds
\le C/(\varepsilon\sqrt m)\)。不需要密度存在或局部中心极限定理。结合(2.2)，本引理对全部确定锚点 \(a\in J_2\) 下的律 \(Q_a\) 一致。

## 4. 固定锚点下的局部增量

先证不依赖 \(\delta\) 的预算。令 \(t=|L-J|<2/K\)。未平滑的 \(H_L-H_J\) 由左右两端各长 \(t/2\) 的扫入计数减去 \(Kt\) 组成。每格一个原子使每个扫入区间至多含两个原子；总计数 \(N\le4\)，且 \(EN=Kt\)。因此
\[
 \|H_L-H_J\|_2^2=E(N-Kt)^2\le EN^2\le4Kt.
\]
正概率箱卷积为 \(L^2\) 收缩，平滑后同样成立。又因 \(\sum m^2c_m^2<\infty\)，\(\|B_L-B_J\|_2\le Ct\)。计入分母 \(L\) 和对数的Lipschitz常数后，在任一固定 \(Q_a\) 下
\[
 E_{Q_a}(\xi_L-\xi_J)^2\le C\varepsilon^2Kt.               \tag{4.1}
\]
另一方面，箱峰不重叠时
\[
 \|\partial_L H_{\delta,L}\|_\infty\le1/\delta+K,
 \quad \|\partial_LB_L\|_\infty\le C(1+\log M)\le CK.
\]
前一导数是两个端点导数的**加号平均**。故
\[
 \|\xi_L-\xi_J\|_\infty
 \le C\varepsilon\min(1,t/\delta).                         \tag{4.2}
\]
这些量在整个 \(J_2\) 有效；(4.1)并不要求第29轮的特殊支持分离。

取一个左端点为 \(a\)、长至多 \(\Delta_0=(Kv)^{-1}\) 的闭cell。对坐标半块 \(A\)，始终在固定律 \(Q_a^{\otimes|A|}\) 下定义
\[
 W_A(L)=\sum_{i\in A}[\xi_L(X_i)-\xi_a(X_i)
                 -E_{Q_a}(\xi_L-\xi_a)],
 \quad M_A=\sup_{L\text{ in cell}}\sum_{i\in A}(\xi_L-\xi_a)(X_i).
\]
因为 \(E_{Q_a}e^{\xi_L-\xi_a}=\int Y_L=1\)，Jensen给锚点相对负漂移
\(E_{Q_a}(\xi_L-\xi_a)\le0\)。于是
\[
 0\le M_A\le\sup_L|W_A(L)|.                               \tag{4.3}
\]
没有对任意有向边另设负漂移。

## 5. 全部dyadic层与局部指数矩

在cell上作嵌套dyadic网格，把右端点加入第0层。第 \(j\) 层父子边长不超过 \(\Delta_0 2^{-j}\)，边数至多 \(2^{j+1}\)。由(4.1)、(4.2)及 \(n\varepsilon^2/v=29/9\)，每条边的独立中心化坐标和总方差至多 \(C2^{-j}\)，每一项的幅度至多
\(b_j=C\varepsilon\min(1,w2^{-j})\)。

这里所用Bernstein界可直接证明。独立中心变量 \(|X_i|\le b\)、\(V=\sum EX_i^2\) 满足
\[
 \log Ee^{\lambda\sum X_i}
 \le\frac{V\lambda^2}{2(1-\lambda b/3)},\qquad0\le\lambda<3/b,
\]
因为 \(E|X_i|^k\le b^{k-2}EX_i^2\)、\(k!\ge2\,3^{k-2}\)。Chernoff及负尾给
\[
 \Pr\{|\sum X_i|>\sqrt{2Vs}+(2/3)bs\}\le2e^{-s}.         \tag{5.1}
\]
\(V=0\) 时和恒为零。逐边取 \(s=u+3(j+1)\)，可数并集失败率不超过
\(He^{-u}\)，其中 \(H=4e^{-3}/(1-2e^{-3})<1\)。边之间无需独立。

固定有限参数下路径及其中心化均值连续，故dyadic链的极限覆盖所有点。记 \(A_w=1+\log w\)，累加各层得
\[
 \Pr\{M_A>C[1+\sqrt u+\varepsilon A_w^2
                            +\varepsilon A_w u]\}\le2e^{-u}. \tag{5.2}
\]
确切地，若 \(J=\lfloor\log_2w\rfloor\)，则
\[
 \sum_{j\ge0}\min(1,w2^{-j})=J+1+w2^{-J},
\]
\[
 \sum_{j\ge0}(j+1)\min(1,w2^{-j})
       =\tfrac12(J+1)(J+2)+w(J+3)2^{-J}.
\]
它们分别为 \(O(A_w)\)、\(O(A_w^2)\)；另一项
\(\sum2^{-j/2}\sqrt{u+3(j+1)}\le C(1+\sqrt u)\)。因此无限细层费用全部计入。

为说明(5.2)足以给指数矩，写 \(D=C(1+\varepsilon A_w^2)\)、\(E=C\)、\(F=C\varepsilon A_w\)。若 \(T\sim\mathrm{Exp}(1)\)，尾比较给
\(M_A\preceq D+E\sqrt{T+\log2}+F(T+\log2)\)。用
\(4E\sqrt T\le T/4+16E^2\)、\(4F\le1/4\)，再积分 \(Ee^{T/2}=2\)，得到
\[
 E_{Q_a^{\otimes|A|}}e^{4M_A}\le C_{C_*}
 \quad\text{若}\quad \varepsilon A_w^2\le C_* .            \tag{5.3}
\]
统一大 \(n\) 门槛合法，因为
\(\varepsilon A_w\le\sqrt{C_*\varepsilon}\to0\)。该结论至少同时覆盖 \(4\le w\le n\)，因为
\(\varepsilon(1+\log n)^2=O_p((\log n)^{5/2}/\sqrt n)\to0\)。

## 6. 两块小球与每cell的弱尾

分坐标为大小与 \(n\) 可比的独立块 \(A,B\)。在固定乘积律 \(\mathbf Q_a\) 下记
\(S_A=\sum_A\xi_a\)、\(S_B=\sum_B\xi_a\)、\(S=S_A+S_B\)。由(3.1)，两个块和的宽度1区间集中概率均不超过
\(Cq\)，\(q=(\varepsilon\sqrt n)^{-1}\asymp v^{-1/2}\)。

逐点 \(\sup_{\rm cell}Z\le S+M_A+M_B\)，而确定锚点的换测度密度是 \(e^S\)，故
\[
 e^hP\{\sup_{\rm cell}Z>h\}
 \le E_{\mathbf Q_a}[e^{h-S}\mathbf1_{S+M_A+M_B>h}].        \tag{6.1}
\]
对 \(S\ge h\) 用条带 \([h+k,h+k+1)\)：每层概率至多 \(Cq\)，权重至多 \(e^{-k}\)，合计 \(Cq\)。对 \(S<h\) 用 \([h-j-1,h-j)\)：权重至多 \(e^{j+1}\)。\(j=0\)仍由小球直接处理；\(j\ge1\)时必须有 \(M_A>j/2\) 或 \(M_B>j/2\)。条件于整个块 \(A\)，独立的 \(S_B\) 落入一个宽度1区间的概率至多 \(Cq\)，所以
\[
 \mathbf Q_a\{S\in[h-j-1,h-j),M_A>j/2\}
 \le Cq\,\mathbf Q_a\{M_A>j/2\}
 \le C_{C_*}q e^{-2j}.                                    \tag{6.2}
\]
另一块同理。于是 \(\sum_{j\ge1}e^{j+1}e^{-2j}<\infty\)，导出
\[
 \sup_h e^hP\{\sup_{\rm cell}Z>h\}\le C_{C_*}q.           \tag{6.3}
\]
没有把同一块的和与局部最大值当成独立量。两端点已包含在局部最大值中，也不另收每cell一个常数费用。

不超过 \(1+|J_2|Kv\) 个cell覆盖 \(J_2\)。由(6.3)及union bound，
\[
 \mathcal W_\delta(J_2)
 \le C_{C_*}\frac{1+|J_2|Kv}{\varepsilon\sqrt n}
 \le C_{C_*,p}\sqrt n.                                    \tag{6.4}
\]
这里 \(K\sqrt v=3p\sqrt{n/29}\)。特别地，对全部 \(4\le w\le n\)，(6.4)由共同常数和共同 \(n_0\) 控制。至此有限增长宽度范围的证明闭合。

## 7. 正背景区间夹逼与任意薄的平滑

对 \(L>\delta\)，正背景与逐原子比较给
\[
 \nu_0(W_{L-\delta}(x))\le\nu_\delta(W_L(x))
              \le\nu_0(W_{L+\delta}(x)).                  \tag{7.1}
\]
确切地，若原子 \(b\) 位于小闭窗口，其整个箱支撑都包含于 \(W_L\)，贡献为1；若箱与 \(W_L\) 的重叠长度为正，原子必在大闭窗口内，贡献至多1。恰好接触端点不破坏不等式。背景部分由非负性和窗口包含关系比较。有限窗口内只需加有限个周期原子。

除以窗口长度，再在左侧把 \(L\) 替换为 \(L+\delta\)，得到
\[
 Y_{\delta,L}\le(1+\delta/L)Y_{0,L+\delta},
 \qquad Y_{0,L}\le(1+\delta/L)Y_{\delta,L+\delta}.          \tag{7.2}
\]
无需错误地将整个 \(\nu_\delta\) 写成 \(\nu_0\) 的卷积；背景并未平滑。

置 \(a_0=10/9\)。若 \(J+[0,\delta]\subset J'\)，(7.2)和共同尺度平移给
\[
 \mathcal W_\delta(J)\le e^{n\delta/a_0}\mathcal W_0(J'),
 \quad \mathcal W_0(J)\le e^{n\delta/a_0}\mathcal W_\delta(J'). \tag{7.3}
\]
这里用到的简单事实是：\(A\le B+d\) 蕴含
\(\sup_he^hP(A>h)\le e^d\sup_te^tP(B>t)\)。严格阈值无需连续分布。

取 \(I=[10/9,11/9]\)、\(J_1=[10/9,5/4]\)、上述 \(J_2\)，并取
\(w_*=n\)、\(\delta_*=(nKv)^{-1}\)。最终 \(\delta_*\le1/36\)，故
\(J_1+[0,\delta_*]\subset J_2\)，而任意 \(0<\delta\le\delta_*\) 有
\(I+[0,\delta]\subset J_1\)。由(6.4)在 \(w=n\) 的特例及(7.3)，
\[
 \mathcal W_0(J_1)\le C_p e^{1/(a_0Kv)}\sqrt n.
\]
对任意 \(w\ge n\)，\(\delta\le\delta_*\)，再次用(7.3)：
\[
 \mathcal W_\delta(I)\le C_p e^{n(\delta+\delta_*)/a_0}\sqrt n
              \le C_p e^{2/(a_0Kv)}\sqrt n.                 \tag{7.4}
\]
两次倍率都已计入，且其积趋于1。\(4\le w\le n\)直接用(6.4)。这就证明(T)，常数不依赖 \(w\) 的增长速度。

## 8. 可测性、边界与适用限制

有限 \(\delta>0\) 的路径连续。原子来源采用闭窗口；固定 \(x\) 和紧尺度区间时只有有限个周期原子可能被扫入，计数路径为右连续的向上阶跃，加上连续背景并除以正的 \(L\)，仍右连续。\(Z_{0,L}\)同理。因此
\[
 \sup_{L\in[a,b]}Z_{0,L}
 =\sup_{L\in((a,b)\cap\mathbb Q)\cup\{a,b\}}Z_{0,L}
\]
可测。右端点明确加入。不存在极限与上确界的未证交换。尤其对称箱在端点处的极限可能只有半个原子，故不能用逐点平滑极限代替(7.1)。

定理是充分大 \(n\) 的陈述。\(n=1\)若另有同样幅度资格，直接 \(\sup_LY_L\le3/2\) 即可；不声称原渐近参数已覆盖所有小维数。\(p\uparrow1\)时(2.3)的门槛不一致，故保留固定 \(p<1\)。任意增大平滑峰高虽不改变质量，也不保留固定倍 \(\log n\) 的高度限制。

结论排除了在**同一确定来源、共同原子与背景、固定窗口**内单靠平滑变薄取得无界正对数增益。这不排除不同来源、非乘积结构、改变窗口或其他机制；也不证明 \(C_n\asymp\sqrt n\)。
