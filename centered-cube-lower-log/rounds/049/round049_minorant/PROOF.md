# 第49轮：从真实 anchor-CF 比较到标量区间下界的严格桥梁

本稿在 `probe.py` 用本机 NumPy/SciPy 实际运行、保存 `RESULTS.json` 后提出。结果分两层：Fourier minorant 给出严格的区间下界转移式；在指定有限参数和宽度1/2/4下，当前CF包络收取的误差很大，Gaussian校准不能给正下界。真实人工dec模型的局部质量仍需另外证明。

## 1. 使用47轮CF命题，不替换比较模型

固定真实 Haar 查询下的来源，anchor尺度 L=7/6。设单坐标正平均为 Y_L，E Y_L=1；原anchor倾斜为 `dP₀=Y_L dP`，比较模型保留同一局部相位θ、把域标签块和Fourier块的离散状态条件解耦，记其正平均为 Ỹ_L。47轮已证明 E Ỹ_L=1，故两个anchor归一化分母精确相同，均为1。

在各自anchor倾斜下取n个 iid 坐标，定义

\[
S=\sum_{i=1}^n\log Y_{L,i},\quad
\widetilde S=\sum_{i=1}^n\log\widetilde Y_{L,i},\qquad
X=(S-v/2)/\sqrt v,\quad Y=(\widetilde S-v/2)/\sqrt v.
\]

这是共同的确定性中心化，不声称两个倾斜均值恰为v/2。所有log都取正实平均的实对数；本桥梁完全不对特征函数取复对数。

记

\[
K=2^{47}-1,\ q=K+1,\ n=10^{32},\ p=99/100,
\quad a_0=p/(KL),\quad U=24,\ V=40/9,\ W=U+V,
\]
\[
v=9p^2n/(29K^2),\quad r=Wa_0<1,\quad\lambda=\sqrt q/K.
\]

V=40/9是统一解析上界，不把有限谱当成独立Gaussian变量。把47轮(A)中η=t/√v的单标量情形直接整理，得到对**所有实数t**的有限参数界

\[
|\varphi_X(t)-\varphi_Y(t)|\le\Xi(t)
:=\min\{2,e^{\beta|t|}(\kappa_1|t|+\kappa_2t^2)\}, \tag{1}
\]

其中

\[
\beta=\frac{Wa_0}{\sqrt v(1-r)},
\]
\[
\kappa_1=\frac{n\lambda UVa_0^2}{\sqrt v}
\left[\frac{1+a_0W}{(1-r)^2}+\frac2{1-r}\right],\qquad
\kappa_2=\frac{n\lambda UVa_0^2}{v}\frac{1+a_0W}{(1-r)^2}.
\tag{2}
\]

47轮的majorant对任意有限η绝对收敛，条件r<1不依赖η。因此(1)是逐个t成立的全频界；不是将有界频率域的O符号延伸到无穷。实际数值β≈4.37651489×10⁻¹⁵、κ₁≈.00250238797、κ₂≈2.12855921×10⁻⁵。

## 2. 候选minorant与严格下界式

给定中心a、半宽d>0和0<s≤d，置

\[
f(x)=\left(1-\frac{(x-a)^2}{d^2}\right)
\exp\left(-\frac{(x-a)^2}{2s^2}\right).
\]

在区间内0≤f≤1，外部f≤0，两个端点f=0，故同时有

\[
f\le\mathbf1_{[a-d,a+d]},\qquad
f\le\mathbf1_{(a-d,a+d)}.
\tag{3}
\]

采用 `hat f(t)=∫ f(x)e^{-itx}dx`。Gaussian变换和二阶微分给

\[
\widehat f(t)=s\sqrt{2\pi}\,e^{-s^2t^2/2}
\left(1-\frac{s^2}{d^2}+\frac{s^4t^2}{d^2}\right)e^{-ita}.
\tag{4}
\]

因s≤d，括号非负，所以其绝对值是非负Gaussian权

\[
\omega_{s,d}(t)=s\sqrt{2\pi}\,e^{-s^2t^2/2}
\left(1-\frac{s^2}{d^2}+\frac{s^4t^2}{d^2}\right).
\]

准确归一化为

\[
\frac1{2\pi}\int_{\mathbb R}\omega_{s,d}(t)dt=f(a)=1.
\tag{5}
\]

f与hat f均可积。由Fourier反演和Fubini（特征函数模≤1、hat f绝对可积），

\[
|E f(X)-E f(Y)|
\le\mathcal E(s,d):=
\frac1{2\pi}\int_{\mathbb R}\omega_{s,d}(t)\Xi(t)dt.
\]

因此严格桥梁是

\[
\boxed{P(X\in(a-d,a+d))\ge E f(Y)-\mathcal E(s,d).} \tag{6}
\]

闭区间也成立。这里Y必须是上节的真实人工dec模型；不能在(6)中擅自把它换成Gaussian。该式不需要Y有密度；所缺局部输入正是对该dec模型的 E f(Y) 的下界，包括区间外的负尾贡献。

## 3. 全频误差的精确表示，没有删尾

令z=s²/d²、φ和Q为标准正态密度及右尾。换元u=st，

\[
\mathcal E(s,d)=2\int_0^\infty\phi(u)(1-z+zu^2)\Xi(u/s)du.
\tag{7}
\]

κ₁、κ₂>0、β≥0，故有唯一t_*>0满足
`exp(βt_*)(κ₁t_*+κ₂t_*²)=2`。从t_*之后精确采用帽2，得

\[
\mathcal E(s,d)=2\int_0^{st_*}\phi(u)(1-z+zu^2)
 e^{\beta u/s}(\kappa_1u/s+\kappa_2u^2/s^2)du
+4[Q(st_*)+z(st_*)\phi(st_*)]. \tag{8}
\]

末项用了 `∫_a^∞u²φ(u)du=aφ(a)+Q(a)`，是完整高频尾，不是截断后忽略。脚本对有限积分使用SciPy quad，对t_*使用brentq；这些小数求值未声称向外数值证书。(6)–(8)本身是严格解析式，且0≤E≤2。

## 4. Gaussian仅作可算校准

若单独取G∼N(0,1)，完成平方得到

\[
E f(G)=\frac{s}{\sqrt{1+s^2}}
 e^{-a^2/[2(1+s^2)]}
\left[1-\frac1{d^2}\left(
\frac{s^2}{1+s^2}+\frac{a^2s^4}{(1+s^2)^2}\right)\right].
\tag{9}
\]

本轮按原log-score区间 `[h,h+w]`，取h=v/2−2、w=1/2/4，故

\[
d=\frac{w}{2\sqrt v},\qquad
a=\frac{h+w/2-v/2}{\sqrt v}
=\frac{-2+w/2}{\sqrt v}.
\]

实际有限参数 `v≈1535.6608350596683`。在0<s/d≤1的网格加有界一维优化中得到：

| 原宽w | 数值最优s/d | E f(G) | 收费E(s,d) | Gaussian参考E f(G)−E |
|---|---:|---:|---:|---:|
| 1 | 1 | .00000207208 | .700059049 | −.700056977 |
| 2 | 1 | .00001658475 | .254547220 | −.254530636 |
| 4 | .814227603 | .014040955 | .108570535 | −.094529580 |

这些不是对所有s的优化认证；下一节另给精确证书，证明任意s都不能使这个Gaussian校准桥为正。脚本还用另一数值积分核对(9)、核对Fourier权总质量1。宽度对应的Gaussian真区间概率约.01017/.02035/.04070，并非上表minorant值。

## 5. 所有0<s≤d的Gaussian参考负性证书

本节下界的是**当前选择包络所收取的误差预算E**，不是实际CF差的下界，也不否定更精细比较或其他minorant。

对Z∼N(0,1)，(7)的权相当于 `(1-z)+z Z²`。对于任意关于|Z|单调增的g，`E Z²g≥E g`；这由两个独立拷贝的乘积差协方差非负得到。Ξ关于|t|单调增且 `Ξ(t)≥min(2,κ₁|t|)`，s≤d，故

\[
\mathcal E(s,d)\ge E\min(2,(\kappa_1/d)|Z|).
\]

本轮精确核验κ₁/d≤1。限制到|Z|≤2，得到

\[
\mathcal E(s,d)\ge\frac{\kappa_1}{d}\sqrt{\frac2\pi}(1-e^{-2}). \tag{10}
\]

另一方面，由(9)，对任意中心a，在d²≤2时

\[
E f(G)\le\max_{0<s\le d}s
\left(1-\frac{s^2}{d^2(1+d^2)}\right)
=\frac{2d\sqrt{1+d^2}}{3\sqrt3}. \tag{11}
\]

若(9)括号负则上界更直接；括号非负时依次用指数≤1、`s/√(1+s²)≤s`即可。目标d²远小于2。

为避免浮点端点，脚本直接用有理式

\[
\frac{\kappa_1}d
=\frac{2n\lambda UVa_0^2}{w}
\left[\frac{1+a_0W}{(1-r)^2}+\frac2{1-r}\right].
\]

以 `floor√q/K`、`ceil√q/K`夹住λ，给该量确切有理上下界；上端≤1。取 `d_up=w/(2 floor√v)`，则d≤d_up；并用

- `sqrt(2/π)>79/100`，由π<22/7及(79/100)²<7/11；
- `1−e^{-2}>6/7`，因为e²的第0至4阶正项和为7；
- `sqrt3>17/10`，以及 `sqrt(1+x)≤1+x/2`。

(10)的有理下界和(11)的有理上界，均以Fraction完整保存于 `RESULTS.json`。它们对任意0<s≤d给

| w | Gaussian参考桥的统一有理上界（此处为近似显示） |
|---|---:|
| 1 | <−.12777637 |
| 2 | <−.05634360 |
| 4 | <−.01306405 |

严格负号由有理数比较认证，不依赖SciPy求积及优化小数。因此在这些宽度上，当前47轮全频CF包络加这个minorant，即使拿标准Gaussian作比较基准，也不能产生正下界。该事实只诊断当前工具链的收费过大；真实dec模型尚未被Gaussian替代或否定。

## 6. 若未来获得dec局部输入，桥梁怎样使用

若能证明同一个dec模型有 `E f(Y)≥m> E(s,d)`，(6)立即给原anchor倾斜下区间概率≥m−E。对于原Haar而非anchor倾斜的log-score S，因 `dP₀^{⊗n}=e^S dP^{⊗n}`，在开区间(h,h+w)上还可准确写

\[
e^h P(S>h)\ge e^{-w}P_0(h<S<h+w)
\ge e^{-w}(m-\mathcal E).
\]

端点用开区间避免严格阈值的原子问题；minorant在端点为0，(6)仍适用。这个条件推论没有提供m，也没有改变当前数值不为正的结论，更不产生维数相关的log增益。

复现：

```sh
/opt/homebrew/bin/python3.14 probe.py
```

使用已安装NumPy2.5.4、SciPy1.18.1；无新安装、无n维Monte Carlo。只输出本目录 `RESULTS.json`。证明范围为严格scalar interval bridge和当前Gaussian校准的负性证书；实际dec局部质量/反集中及进一步弱型收益仍未闭合。
