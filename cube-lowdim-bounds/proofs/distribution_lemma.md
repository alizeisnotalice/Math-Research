# 有限折线核的分布表示引理（传统分析证明）

状态：传统分析证明，经独立审查与主任务复核；未作 Lean 编译。本引理专为
\[
 g=H_R-\sum_{i=1}^m e_i\psi_{a_i},\qquad
 H_R(z)=\frac{(1-R/r(z))_+}{R-1},\quad
 \psi_a(z)=(1-r(z)/a)_+^2,
\]
而写；这里 \(r(z)=|z_1|+|z_2|\)、\(R>1\)、\(0<a_i\le1\)、\(e_i\ge0\)，且 \(m<\infty\)。在原点规定 \(H_R(0)=0\)。

记
\[
 A h(z)=\sum_{j=1}^2\int_{\mathbb R}
 \frac{h(z+te_j)+h(z-te_j)-2h(z)}{t^2}\,dt.
\]
这是项目源码的未归一化 `jumpGen 1`。

## 引理

存在一个局部可积函数 \(q\)（几乎处处就是上式 \(Ag\) 的逐点积分），使对每个 \(C_c^2\) 测试函数 \(\varphi\)，
\[
 \int_{\mathbb R^2}g\,A\varphi
 =\int_{\mathbb R^2}q\,\varphi.\tag{1}
\]
若又已通过标量证书证明 \(Ag\ge0\) 于坐标轴和有限折线之外，则 \(q\ge0\) 几乎处处，因而
\(\beta=q\,dz\) 是正的局部有限测度。没有额外的、支撑在折线或坐标轴上的奇异测度。

## 证明

### 1. 有限分割与正则性

令奇异折线集
\[
 S=\{z_1=0\}\cup\{z_2=0\}\cup\{r=R\}
    \cup\bigcup_{i=1}^m\{r=a_i\}.
\]
它是有限条线段、半直线和直线的有限并，故 \(|S|=0\)。函数 \(H_R\) 有界、全局 Lipschitz；每个 \(\psi_a\) 也有界、全局 Lipschitz。因此 \(g\) 有界且全局 Lipschitz，记 Lipschitz 常数为 \(L\)。

在 \(\mathbb R^2\setminus S\) 的每个 cell 内，\(r\) 是仿射函数，\(H_R\) 和每个 \(\psi_{a_i}\) 都是光滑的一元函数与该仿射函数的复合。它们的二阶导数分别由
\[
 \frac{2R}{(R-1)r^3}\le\frac{2}{R^2(R-1)}\quad(r>R),
 \qquad \frac2{a_i^2}\quad(r<a_i)
\]
控制。因此存在有限常数 \(M\)，使每个 cell 上的 Hessian 算子范数不超过 \(M\)。

### 2. 小步二阶差分的局部 \(L^1\) 界

设 \(E\subset\mathbb R^2\) 为紧集，\(j\in\{1,2\}\)，并写
\[
 \delta_{j,t}^2g(z)=g(z+te_j)+g(z-te_j)-2g(z).
\]
取 \(|t|\le t_0\)。在管状部分
\(E\cap N_{|t|}(S)\)，Lipschitz 性给出
\[
 |\delta_{j,t}^2g(z)|\le2L|t|.
\]
有限折线集在固定紧集的 \(|t|\)-管面积为 \(O_E(|t|)\)，所以该部分积分是 \(O_E(t^2)\)。若
\(z\in E\setminus N_{|t|}(S)\)，坐标线段
\([z-te_j,z+te_j]\) 留在同一 cell；二阶 Taylor 公式给出
\[
 |\delta_{j,t}^2g(z)|\le M t^2.
\]
因此
\[
 \int_E|\delta_{j,t}^2g(z)|\,dz\le C_Et^2\qquad(0<|t|\le t_0).\tag{2}
\]
对 \(|t|>t_0\)，有 \(|\delta_{j,t}^2g|\le4\lVert g\rVert_\infty\)。由 (2) 与 \(\int_{|t|>t_0}t^{-2}dt<\infty\) 得
\[
 \int_E\int_{\mathbb R}\frac{|\delta_{j,t}^2g(z)|}{t^2}\,dt\,dz<\infty.\tag{3}
\]
Tonelli 表明对几乎处处 \(z\)，两个坐标方向的逐点积分存在且其和 \(q(z)\) 属于 \(L^1(E)\)。这对每个紧 \(E\) 成立，故 \(q\in L^1_{\rm loc}\)。

### 3. Fubini 与无遗漏奇异项

对 \(\varphi\in C_c^2\)，先固定 \(j,t\)。平移不变性及有界 \(g\)、紧支撑 \(\varphi\) 给出
\[
 \int g(z)\,\delta_{j,t}^2\varphi(z)\,dz
 =\int\varphi(z)\,\delta_{j,t}^2g(z)\,dz.\tag{4}
\]
左式乘 \(t^{-2}\) 后在 \((z,t)\) 上绝对可积：小 \(t\) 时 \(\delta^2\varphi=O(t^2)\)，且其支撑只落在一个固定紧集的有界膨胀；大 \(t\) 时 \(\lVert\delta^2\varphi\rVert_1\le4\lVert\varphi\rVert_1\)，而 \(t^{-2}\) 可积。右式的局部绝对可积性由 (3) 给出。于是对 (4) 用 Fubini 并对 \(j\) 求和，得到 (1)。

这一步尤其说明：折线上的不可微性已经被 (2) 吸收为 \(L^1\) 密度，故不存在另加的折线、坐标轴或角点奇异测度。此结论不是只从“轴外 \(Ag\ge0\)”形式推出的，而是来自绝对可积 Fubini 恒等式。

### 4. 应用于 \(H_R-\sum e_i\psi_{a_i}\)

此前的多凸起标量引理在 \(S\) 外给出
\[
 Ag(z)\ge\frac4{R(R-1)}Q(r(z))\ge0
 \quad(0<r<R),
\]
而 \(r\ge R\) 时完整补丁核 \(K\) 非负且在该点为零，故其二阶差分积分非负。轴外原点势项是 \(A\)-调和的，因而此处 \(Ag=AK\)，从而仍得 \(Ag\ge0\)。由上段，\(q=Ag\ge0\) 几乎处处，故 \(\beta=q\,dz\) 为正局部有限测度。

将原始核（在原点另取任意非负有限补丁值）写成
\[
 K=\frac{2\pi R}{R-1}\,G-\frac1{R-1}+g
 \quad\text{a.e.},\qquad G(z)=\frac1{2\pi r(z)}.
\]
源码 `Delta.lean` 的基本解恒等式是 \(AG=-2\pi\delta_0\)。常数的生成元为零，(1) 因而给出
\[
 \langle K,A\varphi\rangle
 =\int\varphi\,d\beta
 -\frac{4\pi^2R}{R-1}\varphi(0).\tag{5}
\]
这里所有凸起都是有界函数，故不能改变原点点质量的系数或符号。

### 5. \(\beta\) 的总质量有限

取固定非负 \(C_c^2\) 截断 \(\chi\)，满足 \(\chi=1\) 于单位球；令
\(\chi_L(z)=\chi(z/L)\)。可取 \(0\le\chi_L\le1\)、\(\chi_L\uparrow1\)，且生成元缩放给出
\[
 \lVert A\chi_L\rVert_\infty\le C/L.
\]
将 \(\varphi=\chi_L\) 代入 (5)，并设
\(a=4\pi^2R/(R-1)\)，得到
\[
 \int\chi_L\,d\beta
 =a\chi_L(0)+\int K\,A\chi_L
 \le a+\frac{C}{L}\lVert K\rVert_1.
\]
\(K\) 已由核支配、紧支撑和径向质量公式知属 \(L^1\)。令 \(L\to\infty\)，用单调收敛得到 \(\beta(\mathbb R^2)\le a<\infty\)。因此 (5) 正好满足传递定理 `convolution_le_of_eq_zero` 对有限正测度的输入条件。

## 使用范围与限制

本证明依赖凸起数有限、所有 \(a_i>0\)。若允许无限个尺度、\(a_i\downarrow0\)，管状区和 Hessian 常数的求和需要新的可求和条件，本引理不适用。它也不替代对连续区间标量条件 \(Q\ge0\) 的认证；两者合起来才给出多凸起核的一般 \(L^1\)、全尺度弱型传递。
