# 真实负责混合的 Lp 谱约束：展开、可得结论与 Gamma 排除

本稿不修改已有谱障碍稿，只核对一项更强的真实几何必要条件。主结论是：`rank_budget.md` 第 6 节给出的特定 Gamma 秩混合，虽然满足其中列出的抽象熵与矩性质，却违反负责混合的真实 Lp 约束，因此不能作为满足该更强约束的抽象例子。证明使用实际信息谱 `L=T+log k(T)`，不是把 T 当作 L。

这项排除不证明一般近中央谱的正态逼近，也不闭合一般端点系数。下面同时说明：单侧矩母函数上界保留线性熵亏损时，不能直接逐阶求导得到方差不超过 mu n。

## 1. 真实矩约束与零点展开

记 `w=dP/dQ`，P 为真实负责混合，Q 为完整并的均匀律。负责分割与强矩形最大算子的对偶给出，对每个 p>1，

\[
\|w\|_{L^p(Q)}\le c_{p/(p-1)}^n,
\]

其中 c_q>1 是

\[
(q-1)c_q^q-qc_q^{q-1}-1=0
\]

的唯一正根。因此对每个 t>0，

\[
\boxed{\log E_Pe^{tL}\le nK(t),\qquad
K(t)=(1+t)\log c_{(1+t)/t}.}
\tag{1}
\]

令 `mu=1+exp(-mu)`，并设 `z(t)=log(c_{(1+t)/t})/t`。根方程等价于

\[
\frac{e^{tz}-1}{t}=1+e^{-z}.
\tag{2}
\]

左边在 t=0 按解析延拓理解。t=0 时 z=mu，关于 z 的导数为 `1+e^{-mu}=mu>0`，所以解析隐函数定理适用。写 `z=mu+a t+b t^2+O(t^3)`，比较系数得到

\[
a=-\frac\mu2,\qquad b=\frac{\mu(9-\mu)}{24}.
\]

由于 `K=t(1+t)z`，最终为

\[
\boxed{K(t)=\mu t+\frac\mu2t^2
-\frac{\mu(\mu+3)}{24}t^3+O(t^4).}
\tag{3}
\]

这个展开有一个重要的负三阶系数。只保留前两项会漏掉下面排除 scale 趋于 1 的 Gamma 模型所需的信息。

## 2. 线性熵亏损不能被删除

设

\[
M=E_PL=\mu n-\frac12\log n-A,
\qquad \delta=\mu n-M=\frac12\log n+A.
\]

以下取 n 足够大使 delta>=0；若 A 采用加上绝对常数的非负预算记号，只需相应调整常数。由 (1)—(3)，

\[
\log E_Pe^{t(L-M)}
\le\delta t+\frac{\mu n}{2}t^2
-\frac{\mu(\mu+3)n}{24}t^3+O(nt^4).
\tag{4}
\]

左侧的线性系数为零，右侧却有 delta t。因此当 delta>0 时，不能将二次系数直接比较并宣称 `Var L<=mu n`。这里两个函数在零点只接触到函数值，没有接触到一阶导数；线性余量可以容纳较大的二阶项。

### 可严格推出的粗方差界

负责权保证 w>0 Q-几乎处处，故 `E_P exp(-L)=1`。结合 (1)，可以推出

\[
\boxed{\operatorname{Var}_P L\le Cn(1+\delta)
=O\bigl(n(1+\log n+A)\bigr).}
\tag{5}
\]

证明令 `X=L-mu n`、t=1/n。由 (3)，

\[
E[e^{X/n}-1-X/n]\le C/n+\delta/n.
\]

在 L>=-n 上，有 `X/n>=-(mu+1)`，而凸性积分余项给

\[
e^y-1-y\ge\tfrac12e^{-(\mu+1)}y^2,
\quad y\ge-(\mu+1).
\]

所以这一部分的 `E X^2` 至多 `Cn(1+delta)`。在 L<-n 上，`P(L<-u)<=e^{-u}` 给出剩余贡献 `O(n^2e^{-n})`。最后用 `Var L<=E X^2` 即得 (5)。当 A=O(1) 时，这里只有 O(n log n)，没有得到 O(n)。

### 中位数与上尾

对任意满足 `P(L>=ell)>=1/2` 的中位数 ell，Chernoff 给

\[
\ell\le\inf_{t>0}\frac{nK(t)+\log2}{t}
\le\mu n+\sqrt{2\mu n\log2}+O(1).
\tag{6}
\]

最后代入 `t=sqrt(2 log2/(mu n))` 并使用 (3)。局部二次上界还给

\[
P(L-\mu n\ge x)\le
\exp\{-c\min(x^2/n,x)\},\quad x\ge0.
\tag{7}
\]

因此 `E(L-mu n)_+=O(sqrt n)`，再由均值亏损得 `E(mu n-L)_+=delta+O(sqrt n)`。对中位数的下侧应用 Markov，得到 `ell>=mu n-2delta-O(sqrt n)`。这些普适后果仍在平方根尺度，并未给出所需的对数尺度中位数结算。这里控制的是 L 的中位数；秩变量 T 的中位数 V 还涉及 `log k(T)`，不能未经证明直接替换。

### 为何 O(n log n) 不能仅由这些标量条件自动改成 O(n)

以下仅是概率论上的诊断，不是真实矩形构造。首先对所有 t>0 有 `K(t)>=mu t`：在 (2) 中测试 `z_0=mu/(1+t)`，所需不等式归结为 `e^x(1-x)<=1`，其中 `x=mu t/(1+t)`；故 z(t)>=z_0。

令 a=mu n，并选两点谱 L 取 a 和 b<0，P 概率分别为 1-p 和 p。规定

\[
b=\log p-\log[1-(1-p)e^{-a}],
\]

便精确有 `E_P e^{-L}=1`。选择 p 使 `p(a-b)=delta`，则 `M=a-delta`。对 `delta=(1/2)log n+O(1)`，可由连续性选择 `p~delta/a`、`b~log p`，从而

\[
\operatorname{Var}L=p(1-p)(a-b)^2\sim a\delta.
\]

但 L<=a 保证对全部 t>0 都有 `log E e^{tL}<=mu n t<=nK(t)`。所以全部单侧真实矩上界、归一化恒等式和近饱和熵这些标量条件，单独并不强制 O(n) 方差或正态性。它不满足已认证的矩形/Hall/LF 构造，不能冒充原几何问题的反例。

## 3. Gamma 模型的实际信息矩母函数

对 rank 密度 `f=Gam(a,theta)`，写 `L_f(T)=T+log f(T)`、T~f。直接积分得到精确公式

\[
I_{a,\theta}(t):=E_fe^{tL_f}
=\frac{\Gamma(B)}{\Gamma(a)^{1+t}}\,
\theta^{-t}[1+t(1-\theta)]^{-B},
\quad B=a(1+t)-t,
\tag{8}
\]

适用条件为 `1+t(1-theta)>0` 及 B>0。这是 `integral exp(tT) f(T)^(1+t) dT`，不是 Gamma 随机变量 T 本身的矩母函数。

设 `a theta=mu n-1`，theta 固定或趋于正常数。对固定 admissible t，Stirling 给

\[
\frac1n\log I_{a,\theta}(t)
=F_\theta(t)+o(1),\qquad
F_\theta(t)=\frac\mu\theta(1+t)
\log\frac{1+t}{1+t(1-\theta)}.
\tag{9}
\]

其零点展开为

\[
F_\theta(t)=\mu t+\frac{\mu\theta}{2}t^2
+\frac{\mu\theta(2\theta-3)}6t^3+O(t^4).
\tag{10}
\]

## 4. 主结论：原 theta_n=1+n^{-1/4} 混合违反真实 Lp 约束

原模型为

`k_n=(1-epsilon_n) f_0+epsilon_n(f_-/3+2 f_+/3)`，

其中 `f_0=Gam(a_n,theta_n)`、`theta_n=1+n^{-1/4}`、`a_n theta_n=mu n-1`、`epsilon_n=O(log n/sqrt n)`。关键是逐点 `k_n>=(1-epsilon_n)f_0`，所以其**实际**信息谱 `L=T+log k_n(T)` 满足

\[
E_{k_n}e^{tL}=\int e^{tT}k_n(T)^{1+t}\,dT
\ge(1-\epsilon_n)^{1+t}I_{a_n,\theta_n}(t).
\tag{11}
\]

固定 t>0 后，theta_n 趋于 1，因而右侧对数除以 n 的极限为

\[
F_1(t)=\mu(1+t)\log(1+t).
\]

由 (3) 和 (10)，

\[
F_1(t)-K(t)=\frac{\mu(\mu-1)}{24}t^3+O(t^4)>0
\tag{12}
\]

对所有充分小的固定正 t 成立。这里不需要任何随 n 缩小的参数，也不需要对混合的四阶累积量建立一致展开。仅 (11)—(12) 就已严格排除该模型满足 (1)。

现有有理区间证书进一步选定 `t=1/10`，此时 q=11、c_11 解多项式 `10 c^11-11 c^10-1=0`，并认证

\[
\boxed{
0.00001411891677682855533678003335204521785516464876468147205
\ <\ F_1(1/10)-K(1/10)\ <
0.00001411891677682855533678003335204521785516464876468147206.}
\tag{13}
\]

证书见 `certify_lp.json`，算法见 `certify_lp.py`。它用有理二分夹住 c_11，用已存 mu 区间并检查 `(mu-1)exp(mu)=1` 的端点符号；指数函数采用 100 项正项级数及几何余项，对数采用 80 项 atanh 级数及显式余项。所有比较在有理数中进行，小数输出向外舍入。

因此对充分大 n，(11) 的对数已经至少比 `nK(1/10)` 高出固定正数乘 n，和真实负责约束矛盾。这个结论强于“方差多出 n^{3/4} 可能不允许”的启发：即使 scale 恰取 1，Gamma 领先率的三阶项仍与真实 envelope 不相容。原构造确有额外的方差膨胀，但本证明没有通过非法比较二阶系数来排除它。

这也不意味着任何截断、修改或其它 Gamma 方向都已被排除；高阶矩对小概率远尾非常敏感，每种修改都需要重新核验。

## 5. 小 scale 的领先率确有余量，但不是有限维几何结论

主代理提出的 `0<theta<=1/2` 方向在代数上成立。令 `x=theta t/(1+t)`，则

\[
-\log(1-x)\le x+\frac{x^2}{2(1-x)}
\le x+\frac{x^2}{2\theta},
\]

其中用了 `x<theta<=1/2`。代入 (9) 得对全部 t>0，

\[
F_\theta(t)\le\mu t+\frac\mu2\frac{t^2}{1+t}.
\tag{14}
\]

另一侧，设 `r=e^{-mu}=mu-1`、`s=t/(1+t)`，以及

\[
\psi(s)=-\log(1-s)+
\log\frac{1+r e^{\mu s}}{1+r}.
\]

直接最大化

`[1/(1-s)] [1+u^(1-s)]/(1+u)`，`u>=0`，

得到 c_{1/s}；在 u=r 取值因此给

\[
K(t)\ge(1+t)\psi(t/(1+t)).
\]

有 `psi'(0)=psi''(0)=mu`，而若 p_s 表示对应倾斜 Bernoulli 概率，

\[
\psi'''(s)=\frac2{(1-s)^3}
+\mu^3p_s(1-p_s)(1-2p_s)
\ge2-\frac{\mu^3}{6\sqrt3}>0.
\]

所以

\[
\boxed{F_\theta(t)\le\mu t+\frac\mu2\frac{t^2}{1+t}
\le(1+t)\psi(t/(1+t))\le K(t).}
\tag{15}
\]

这证明了小 scale Gamma 的领先指数率方向有余量。它没有证明有限 n 时对全部 t 同时满足 (1)，也没有处理所有矩形/LF兼容性。

还必须修复 rank 单调性：theta<1 时

`d log(e^T f(T))/dT=(a-1)/T+1-1/theta`

在 `T>theta(a-1)/(1-theta)` 后变负，所以未修改的 Gamma 并不是全域合法递减秩密度。虽然该区域位于远右尾，修改远尾对所有 t 的矩约束仍不能自动忽略。因此它目前只是更强抽象模型的探索方向，不是已经构造完成的谱障碍，更不是真实近中央反例。

## 核验与结论边界

本机已有 Wolfram 引擎实际完成了精确符号计算，返回 `a=-mu/2`、`b=(9mu-mu^2)/24`，并确认 (3)、(10) 和 `(F_1-K)` 的三阶系数 `mu(mu-1)/24`。没有运行 MATLAB 或低维扫描。另已读取并逐项检查固定 t 有理证书代码的余项方向；没有修改或重新生成其文件。

本轮严格推进是：真实 Lp 约束排除了已给出的 theta_n>1 Gamma 抽象构造，并明确留下线性熵亏损、远尾与全 LF 面选择等问题。它不等价于已证明 `Sigma=O(sqrt n)`、所需 Kolmogorov 精度或一般对数端点结算。

## 附录：满足全部正参数矩界的有限三点抽象谱

下面的构造比被排除的 Gamma 例更直接。它同时满足全部正参数的矩上界 (1)、归一化、锐到常数的 M，以及自然尺度的二、三阶信息矩，却仍违反目标逆检验结算。它是一个完全明确的有限概率空间；没有构造真实矩形目录，也没有认证负责运输、Hall、LF 或逐秩反馈约束。

### A. 显式参数与概率律

令

\[
m=\mu n,\qquad \delta=\tfrac12\log n+4,
\qquad h=m-\delta,\qquad p=e^{-m}.
\]

选择 a=a_n>0 满足

\[
a^2\le m/2,\qquad a\le m/4.
\tag{A1}
\]

以下只要求 n 充分大。定义

\[
Z=1-e^{-h}\left[\tfrac12e^{-a}+(\tfrac12-p)e^a\right],
\qquad b=\log p-\log Z.
\tag{A2}
\]

因 h-a>=m/2 对充分大 n 成立，`1-e^{-h+a}<=Z<1`，所以 Z>0、p<1/2、b<0 且 b<h-a。更精确地 `b=-m+O(e^{-h+a})`。在三个点上定义 P 与信息值 L：

| 点 | P 质量 | L |
|---|---:|---:|
| + | 1/2 | h+a |
| 0 | 1/2-p | h-a |
| - | p | b |

再令 `Q_i=P_i exp(-L_i)`，得到显式三概率

\[
Q_+=\tfrac12e^{-h-a},\qquad
Q_0=(\tfrac12-p)e^{-h+a},\qquad Q_-=Z.
\tag{A3}
\]

三者严格为正，且由 (A2) 精确相加为 1。因此 L 确实等于 `log(dP/dQ)`，并精确满足 `E_P exp(-L)=1`；没有把任意随机变量错误地称为信息谱。

把 Q 质量实现为 `(0,1)` 上依次排列的三个区间，在区间上令密度分别为 `exp(h+a), exp(h-a), exp(b)`，就得到一个有限、严格正、递减的三步 rank 密度 w*。其端点没有无界或零密度问题。对应 `k(T)=exp(-T) w*(exp(-T))` 为合法 rank 密度，且 P 下的 `T+log k(T)` 恰好具有上述三点信息律。

### B. 均值与自然矩

令 `d=h-a-b>0`。直接计算为

\[
M=\tfrac12(h+a)+(\tfrac12-p)(h-a)+pb
=h-pd.
\tag{A4}
\]

由于 d=O(n)，有

\[
\boxed{M=\mu n-\tfrac12\log n-4+O(ne^{-\mu n}).}
\tag{A5}
\]

精确方差为

\[
\operatorname{Var}L
=(1-p)a^2+p(h-b)^2-p^2d^2
=a^2+O(n^2e^{-\mu n}).
\tag{A6}
\]

同理由三点逐项展开可得

\[
E|L-M|^3=a^3+O(n^3e^{-\mu n}).
\tag{A7}
\]

所以 a=O(sqrt n) 时，信息谱具有 O(n) 方差和 O(n^{3/2}) 三阶绝对中心矩；并非通过增大到 O(n log n) 方差来逃避约束。

### C. 全部 t>0 的真实矩上界同时成立

这里不使用任何小 t 展开。由 b<=h-a，精确有

\[
\begin{aligned}
E e^{t(L-m)}
&=e^{-\delta t}\cosh(at)
+p\{e^{t(b-m)}-e^{-(\delta+a)t}\}\\
&\le e^{-\delta t}\cosh(at)\le\cosh(at).
\end{aligned}
\tag{A8}
\]

对 x>=0，`log cosh x<=min(x^2/2,x)`。本稿 (15) 已经严格证明

\[
nK(t)-mt\ge\frac m2\frac{t^2}{1+t}.
\]

当 0<t<=1 时，由 (A1)，

\[
\log\cosh(at)\le\frac{a^2t^2}{2}
\le\frac m4t^2\le\frac m2\frac{t^2}{1+t};
\]

当 t>=1 时，同样有

\[
\log\cosh(at)\le at\le\frac m4t
\le\frac m2\frac{t^2}{1+t}.
\]

因此对每个充分大的 n，**同时对所有 t>0** 有

\[
\boxed{\log E_Pe^{tL}\le nK(t).}
\tag{A9}
\]

等价地，全部范数约束 `||w||_r<=c_{r/(r-1)}^n`、r>1，也同时满足。这不是只核验有限个 t，也不是只核验领先指数率。

### D. 半命中检验价仍可超出目标

最高密度点的 P 质量恰为 1/2。任何半命中检验的 Q 预算至少为 `(1/2)exp(-h-a)`，而取完整的最高点恰好达到。因此抽象 NP 预算和其对数价格精确为

\[
b_{1/2}(P,Q)=\tfrac12e^{-h-a},\qquad
\boxed{V=h+a+\log2.}
\tag{A10}
\]

在上述连续 rank 实现中，top-rank 预算端点也正是 Q_+，所以同一个 V 是 rank T 的半上尾阈值；有限步信息谱没有引入阈值歧义。

有两个直接选择：

- 给定固定 v>0，取 `a=delta+v log n`。则 (A1) 对充分大 n 成立，并且 `V=mu n+v log n+log2`、`V-M=(1/2+v)log n+4+log2+O(ne^{-mu n})`。同时方差只有 O(log^2 n)，三阶绝对中心矩只有 O(log^3 n)。
- 给定固定 c>0 且 `c^2<mu/2`，取 `a=c sqrt n`。此时 `Sigma~c sqrt n`、`E|L-M|^3~c^3n^{3/2}`，但 `V-mu n=c sqrt n-(1/2)log n-4+log2`。它在保持全部 (A9) 的同时，甚至具有平方根尺度的正检验偏移。

信息谱始终有质量 1/2 的最高原子，因此标准化之后与任一连续正态 CDF 的 Kolmogorov 距离至少为 1/4。自然方差和三阶矩本身不提供独立和结构，不能据此调用 Berry–Esseen。

### E. 结论边界与核验

这个例子严格证明：全部正参数矩族、`E exp(-L)=1`、有限正 rank、锐熵均值，以及自然二、三阶信息矩，作为一组标量条件仍不足以推出目标逆检验界。它并未证明这些标量条件穷尽真实负责几何。

尤其本例没有旧 Gamma 例对平滑比较律的小反向KL保证。更强地，[rank_budget.md附录B](rank_budget.md)已证明完整逐秩预算要求k≥c√n/(A+1)，而此例k=3、A=O(1)，所以它被这一真实必要条件严格排除。它不能被表述为满足全部已有几何比较的反例，更不能据此否定一般端点系数。这里的 V 是已明确实现的抽象 NP 价格；没有对应的真实全行源价 m_*、rho 或 G。

本机 Wolfram 对三点律的均值、全方差和 `MGF-cosh` 精确恒等式进行了独立代数化简，三个残差均为零。上述全 t 不等式与有限 rank 实现依据解析证明，没有进行低维扫描。

**联合排除的最终状态。** 本稿三阶梯构造只反驳正参数矩族与所列标量条件的充分性；完整逐秩预算已经排除了该构造，见[负责密度层数定理](rank_budget.md)。一般系数仍未确定。
