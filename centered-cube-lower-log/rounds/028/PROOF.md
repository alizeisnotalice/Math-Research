# 既有周期乘积来源的固定窗口平方根上界

结论：对原文件指定的确定标签来源、固定 (0<p<1)、固定有限 (w\ge4)，在原参数选择下，对充分大的每个整数 (n)，有
\[
 \sup_{h\in\mathbb R}e^hP^{\otimes n}\!\left\{
 \sup_{L\in I}\sum_{i=1}^n\log Y_L(X_i)>h\right\}
 \le e+C|I|\sqrt{K/\delta}\le C_{p,w}\sqrt n,
 \qquad I=[10/9,11/9].                                      \tag{1}
\]
这里 (C) 只依赖本来源中固定的低频系数和固定窗口；最终常数允许依赖固定的 (p,w)。这不是所有输入的 (C_n) 上界，也不是允许 (w=w_n\to\infty) 的统一结论。

核心修正是使用宽度为 1 的 occupation strip，只要求和分布的区间集中概率，不要求整个卷积分布有有界密度。因此正密度混合成分留下的残余尖峰不构成这条证明的缺口。

以下除已在原来源文件证明的交叉界与系数界外，概率论步骤全部展开证明；不依赖 Rice 公式或局部中心极限定理。源文件的原样归档为本目录 `source_interface_baseline.tex`；出处与哈希见 `provenance.json`。

## 1. 来源接口，以及从 safe 集提升到整个 I 的关键补充

令 (F_L=H_{L,\delta}+B_L)、\(\varepsilon=p/K\)，则
\[
Y_L=1+\varepsilon F_L/L,
\quad \int Y_L=1,
\quad |F_L|<154/27,
\quad |Y_L-1|\le(77/15)\varepsilon=:r_n.
\]
最终 (r_n\le1/2)，故 (1/2\le Y_L\le3/2)。原来的核引理给 safe 集上的 (G(L,L)=\int F_L^2\ge1/4)，并不直接给整个 (I) 的该下界。

但原文件的“真实交叉与平滑预算”引理是对**全部** (L,J) 成立的，且 (C_x\to0)。由于
\[
\begin{aligned}
G(L,L)
 &=\|H_{L,\delta}\|_2^2+\|B_L\|_2^2
     +2\langle H_{L,\delta},B_L\rangle\\
 &\ge \frac{c_1^2}{2}\sin^2(\pi L)-C_x
 \ge \frac{b_1}{3}\sin^2(\pi/9)-C_x,
\end{aligned}
\]
其中 (c_1^2=2b_1/3)、(b_1>0)，最终在整个 (I) 上有
\[
G(L,L)\ge g_0:=\frac{b_1}{6}\sin^2(\pi/9)>0.                 \tag{2}
\]
这一粗正下界无需安全尺度的平坦方差近似。原参数 (M=\lceil\sqrt v\rceil<K) 满足交叉引理前提，而且
\[
C_x=O(K^{-1/2})+O((\log M)/K)\longrightarrow0.
\]
因此 (1) 的全窗口结论可以闭合。若不调用这个全窗口交叉引理，则下面相同论证只能先应用于 safe 集，且必须支付每个连通分支的左端点费用；见第 6 节。

以下统一取充分大参数，明确要求 n≥2、ε=p/K≤15/308、C_x≤g_0、Kδ≤1、箱峰不重叠、1≤M<K，以及既有来源非负且周期质量1。它们在原参数选择中最终同时成立。

## 2. tilted 方差与 score 二阶矩

对每个固定 (L)，定义概率测度
\[
dQ_L(x)=Y_L(x)dP(x),\qquad
d\mathbf Q_L(\mathbf x)=\prod_iY_L(x_i)dP^{\otimes n}(\mathbf x).
\]
在 \(\mathbf Q_L\) 下，各坐标仍独立同分布。令单坐标
\(\xi_L=\log Y_L\)、\(\eta_L=\partial_L\log Y_L\)，全和
\(Z_L=\sum_i\xi_L(X_i)\)、\(T_L=\sum_i\eta_L(X_i)\)。

取两个独立的 (P) 坐标，利用对数的均值定理及 (Y\in[1-r_n,1+r_n])，得到
\[
\begin{aligned}
\operatorname{Var}_{Q_L}(\xi_L)
 &=\frac12\iint[\log Y_L(x)-\log Y_L(y)]^2
                 Y_L(x)Y_L(y)dP(x)dP(y)\\
 &\ge\left(\frac{1-r_n}{1+r_n}\right)^2
                 \operatorname{Var}_{P}(Y_L)\\
 &=\left(\frac{1-r_n}{1+r_n}\right)^2
                  \varepsilon^2G(L,L)/L^2
 \ge c_0\varepsilon^2.                                    \tag{3}
\end{aligned}
\]
此外 \(|\xi_L|\le C\varepsilon\) 一致成立。

score 均值**准确为零**：
\[
E_{Q_L}\eta_L=\int\partial_LY_L=0.                         \tag{4}
\]
注意这不是 \(\partial_LE_{Q_L}Z_L=0\)，后者一般不成立，因为测度也随 (L) 变化。

平滑原函数满足 \(\mathcal P_\delta'=K(\mu_{K,\delta}-1)\)。箱峰互不重叠，故
\[
\|\mu_{K,\delta}-1\|_2^2=1/(K\delta)-1.
\]
几乎处处
\[
\partial_LH_{L,\delta}(x)
 =\tfrac12[\mathcal P_\delta'(x+L/2)+\mathcal P_\delta'(x-L/2)],
\]
由平移保持范数和凸性，
\[
\|\partial_LH_{L,\delta}\|_2^2\le K/\delta-K^2\le K/\delta.
\]
同时
\[
\partial_LB_L=\sum_{m\le M}\pi m c_m\cos(\pi mL)\cos(2\pi mx),
\quad
\|\partial_LB_L\|_2^2
 \le\frac{\pi^2}{2}\sum_{m\ge1}m^2c_m^2\le C.
\]
由 \(\partial_LY_L=\varepsilon(F_L'/L-F_L/L^2)\)、(L\in I)，并使用 \(K\delta\le1\)，得到
\[
J_L:=E_{Q_L}\eta_L^2
 =\int\frac{(\partial_LY_L)^2}{Y_L}
 \le C\varepsilon^2K/\delta
 =Cp^2/(K\delta).                                         \tag{5}
\]
积分交换与微分合法性：固定的有限 (K,M,\delta>0) 下，\(\mathcal P_\delta\) 为 Lipschitz 函数，(Y_L(x)) 对 (L) 一致 Lipschitz，且一致正下界成立。其导数有界，Fubini 和支配收敛适用；导数恒等式只需在 (dL\,dP) 几乎处处成立。

## 3. 自含的宽度 1 集中概率估计

以下引理避免任何密度假设。

设实随机变量 \(\xi\) 满足 \(|\xi|\le C_1\varepsilon\)、\(\operatorname{Var}\xi\ge c_1\varepsilon^2\)，且 \(\varepsilon\to0\)。若 (S_m) 是 (m) 个独立同分布副本的和，则
\[
\sup_aP\{S_m\in[a,a+1]\}\le\frac{C}{\varepsilon\sqrt m}.   \tag{6}
\]
右端大于 1 时仍有效；也可取与 1 的较小值。

证明：记特征函数为 \(\phi(t)=Ee^{it\xi}\)，取独立副本 \(\xi'\)。对 \(|t|\le1\)，最终 \(|t(\xi-\xi')|\le1\)，所以
\[
|\phi(t)|^2=E\cos(t(\xi-\xi'))
 \le1-\tfrac14t^2E(\xi-\xi')^2
 =1-\tfrac12t^2\operatorname{Var}\xi.
\]
因此 \(|\phi(t)|^m\le\exp(-cm\varepsilon^2t^2)\)。使用恒等式
\[
k(u):=\left(\frac{\sin(u/2)}{u/2}\right)^2
 =\int_{-1}^1(1-|t|)e^{itu}dt,
\]
并在 (u=0) 连续延拓。(k\ge c>0) 于 \(|u|\le1/2\)，故对区间中心 (a+1/2)，
\[
\begin{aligned}
P\{S_m\in[a,a+1]\}
 &\le C E k(S_m-a-1/2)\\
 &\le C\int_{-1}^1|\phi(t)|^m dt
 \le C\int_{\mathbb R}e^{-cm\varepsilon^2t^2}dt
 =C/(\varepsilon\sqrt m).
\end{aligned}
\]
Fubini 仅在有限区间上交换有界被积函数，因而无正则性缺口。将 (3) 代入，(6) 对所有 (L\in I) 的 tilted 单坐标律一致成立。

## 4. 两块独立消除 score 与 strip 的相关性

把坐标分为 (A,B)，两块大小均与 (n) 可比。令 (S_A=\sum_A\xi_i\)、(T_A=\sum_A\eta_i\)，另一块同理。对任意宽度为 1 的区间 (J)，条件于块 (A)，由块间独立性与 (6)，
\[
E_{\mathbf Q_L}[\mathbf1_{S_A+S_B\in J}|T_A|]
 \le\frac{C}{\varepsilon\sqrt n}E_{\mathbf Q_L}|T_A|
 \le\frac{C}{\varepsilon\sqrt n}\sqrt{|A|J_L}.
\]
最后一步使用 (4)，所以块内 score 方差按坐标相加。交换两块后，由三角不等式得到
\[
E_{\mathbf Q_L}[\mathbf1_{Z_L\in J}|T_L|]
 \le C\varepsilon^{-1}\sqrt{J_L}
 \le C\sqrt{K/\delta}.                                    \tag{7}
\]
这里没有把全和的集中概率与 score 的期望错误地相乘；条件化所用的集中概率总来自另一独立块。

## 5. 固定宽度 occupation 不等式及完整结论

对任意绝对连续实路径 (z:[a,b]\to\mathbb R)，
\[
\mathbf1_{\sup z>h}
 \le\mathbf1_{z(a)>h-1}
   +\int_a^b\mathbf1_{h-1<z(L)<h}(z'(L))_+\,dL.             \tag{8}
\]
若左端事件成立而首项不成立，路径从 \(z(a)\le h-1\) 到达 (h)，截断函数
\(u\mapsto\min(1,\max(0,u-h+1))\) 至少增加 1；链式法则给出 (8)。开区间边界的导数处理合法，因为绝对连续函数在其固定水平集上导数几乎处处为零。固定宽度无需取 \(\epsilon\downarrow0\)，因此无需 Rice 公式及连续密度。

将 (8) 用于 (Z_L\)，乘 (e^h\)，端点项由 \(E_Pe^{Z_a}=1\) 得到
\(e^hP(Z_a>h-1)\le e\)。strip 内利用 \(d\mathbf Q_L=e^{Z_L}dP^{\otimes n}\)，有
\[
\begin{aligned}
e^hE_P[\mathbf1_{h-1<Z_L<h}(T_L)_+]
 &=E_{\mathbf Q_L}[e^{h-Z_L}\mathbf1_{h-1<Z_L<h}(T_L)_+]\\
 &\le eE_{\mathbf Q_L}[\mathbf1_{h-1<Z_L<h}|T_L|]
 \le C\sqrt{K/\delta}.
\end{aligned}
\]
Tonelli 后得到 (1) 的第一个不等式。原参数精确关系为
\[
\delta=1/(wKv),\qquad v=9p^2n/(29K^2),
\]
因此
\[
\sqrt{K/\delta}=\sqrt w\,K\sqrt v
 =3p\sqrt{w/29}\sqrt n.                                   \tag{9}
\]
这是完整的既有来源、固定窗口上界。它说明仅在同一来源、该固定尺度窗口内重新组织 excursion 或选取更多尺度，不能把相应弱尾系数推到 \(\sqrt n\log n\)。它没有排除修改来源结构、扩大窗口、令 (w\) 随 (n\) 增长或使用其他输入带来不同结论。

## 6. safe 接口、原密度路线，以及未使用的断言

若仅保留原 safe 集 \(\mathcal S\) 的方差引理，不采用第 1 节的全窗口粗下界，那么 \(\mathcal S\) 是有限个相对开区间的并。对每一分支应用 (8)，都须支付一个端点费用 (e\)。闭包端点由路径与方差的连续性处理。若分支数为 (N_{\rm comp}\)，则只能先得到
\[
\sup_he^hP\{\sup_{L\in\mathcal S}Z_L>h\}
 \le eN_{\rm comp}+C|I|\sqrt{K/\delta}.
\]
原禁区由 \(|\mathcal D|=O(B^2)\) 个圆周邻域给出，在 (K(L-1)\) 的长度小于 (K) 的窗口上，其补集分支数 \(N_{\rm comp}=O(B^2)\)。固定 (w\)、(B\asymp v\asymp\log n\) 时，该端点费用仍为 (o(\sqrt n)\)。不能在未处理分支数时写成单个 (e\)。

原先“单调好区间推出正质量、密度 \(O(K/p)\) 成分”的想法可以支持部分分布的正则性，但它本身不证明全部卷积的 \(L^\infty\) 密度界：所有坐标落入剩余坏成分的混合项，即使权重指数小，仍可有无界密度。若坚持宽度趋零的 Rice/occupation 路线，这就是必须补齐的第一项。仅说“指数小坏事件”不足以控制除以 strip 宽度后的极限。

本证明并不声称建立了那条完整密度界，也不声称验证了周期系数导数的 uniform monotone interval。二者均被 (6) 的固定宽度集中估计绕开。上界结论本身没有留下该密度缺口。

## 7. 边界与计算核对

n=1时，只要同样的幅度资格成立，sup_L Y_L≤3/2直接给弱尾系数≤3/2，无需分半。原渐近构造本身只对充分大n保证全部资格；不能据此称已覆盖原小维数参数。

本机 checks.py/checks.json 核验谱导数能量、条带集中证明中的特征函数包络和参数代入；Fraction给出J系数672233/1875<360。附加有理全窗口资格见主报告。这些算术和数值证据不替代第2—5节的证明。独立审查保留于independent_audit.md。
