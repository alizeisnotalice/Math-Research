# 小缺陷率不能恢复平均半径三角能量

记 m=2^n，Q(x,r)=x+[-r,r]^n。下面否证一个候选核能量步骤；不否证真实 Tauberian 高阶上界。

## 5. 任意小 t 的平均半径三角自能反例

拟议但错误的式子是
\[
 \int_F\!\int_F T_{r(x)}^{\rm per}(x-y)dy,dx
 \stackrel{?}{\ge}q\int_F r(x)^n dx,
 \qquad q=|F|/3^n.\tag{5.1}
\]
Terra 第 003 轮已在一维 \(t=1\) 给出反例。本节给任意 \(n\ge2\)、任意小 \(t\) 的严格连续反例。

### 构造与高密度见证

在 \(\mathbb T_3^n\) 中，令 \(0<\varepsilon<1\)，
\[
 J=\{0,1,3,4,5,6,8,9\},\quad
 A_\varepsilon=\bigcup_{j\in J}[j/4,(j+\varepsilon)/4),
 \quad F=A_\varepsilon\times\mathbb T_3^{n-1}.
\]
取 \(E=\mathbb T_3^n\setminus F\)，则 \(q=2\varepsilon/3\)。当第一坐标落在 \(j=0,1,8,9\) 的条纹时取 \(r(x)=1\)；落在 \(j=3,4,5,6\) 时取 \(r(x)=1/4\)。

注意 \(F\) 包含于完整的 \(1/4\)-周期、占空比 \(\varepsilon\) 的条纹集。对两个半径，\(2r\) 都是 \(1/4\) 的整数倍，因此每个所选见证都有
\[
 |F\cap Q(x,r(x))|\le\varepsilon(2r(x))^n.\tag{5.2}
\]
在中间条纹上可取方体第一坐标投影完全位于 \([3/4,7/4]\)，于是最大缺陷率恰为 \(\varepsilon\)。取严格参数 \(t=\varepsilon+\varepsilon^2\)（小 \(\varepsilon\) 时小于 1），即得到 \(M_c\mathbf1_E>1-t\) 的严格见证，且 \(t\to0\)。

### 精确能量

一维三角核记 \(T_r(h)=(1-|h|/r)_+\)。按两个半径分组，对上述八个小区间直接积分得到
\[
 \begin{aligned}
 E_1&=\frac{19}{32}\varepsilon^2,\\
 E_{1/4}&=\frac14\varepsilon^2-\frac1{48}\varepsilon^3.
 \end{aligned}\tag{5.3}
\]
可用两次原函数逐对核验：
\[
 P_r(z)=\frac{(z+r)_+^3-2z_+^3+(z-r)_+^3}{6r},
\]
\[
 \int_a^b\int_c^dT_r^{\rm per}(x-y)dy,dx
 =\sum_{k=-1}^1[P_r(b-c-3k)-P_r(a-c-3k)
                  -P_r(b-d-3k)+P_r(a-d-3k)].
\]
八区间的有限精确求和就是 (5.3)，没有抽样积分。

其余 \(n-1\) 个坐标积分各给因子 \(3r\)，所以左侧真实能量为
\[
 \mathcal E
 =3^{n-1}\varepsilon^2
   \left[\frac{19}{32}+4^{-n}-\frac\varepsilon{12}4^{-n}\right].\tag{5.4}
\]
拟议右侧则为
\[
 q\int_F r(x)^n dx
 =3^{n-1}\varepsilon^2\frac23(1+4^{-n}).\tag{5.5}
\]
右侧减左侧恰为
\[
 \boxed{
 3^{n-1}\varepsilon^2
 \left[\frac7{96}-\frac{4-\varepsilon}{12}4^{-n}\right]>0,
 \qquad n\ge2.}\tag{5.6}
\]
方括号在 \(n\ge2\) 时至少 \(5/96\)，因此反例对任意小正 \(\varepsilon\) 都成立。

它没有满足定理 3.1 的小 Lipschitz 常数：相邻不同半径条纹的间距固定为 \(1/4\) 量级，半径跳跃为 \(3/4\)，不会随 \(t\to0\) 消失。

### 反例的结论边界

它严格否证“高密度 \(t\to0\) 自动修复平均 \(r^n\) 三角自能”的猜想，但不否证目标 Tauberian 上界。因为
\[
 \frac q{mt}\longrightarrow\frac{2}{3\cdot2^n}<1,
\]
其首项效率远低于最优；对应 \(|F|/|E|\sim2t/3\)。若要挽救某种平均尺度公式，仍须额外使用**近首项饱和**或局部尺度几何；只写“\(t\) 足够小”不足够。

全参数系数证书：`verification/verify_endpoint_coefficients.py` 用有理多项式逐项核对0<ε<1上的2304个正部表达式；其根均不在开区间，故不是只核对离散ε值。
