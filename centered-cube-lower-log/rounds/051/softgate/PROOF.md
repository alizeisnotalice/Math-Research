# 第051轮：合法soft split、真实query通道与指数产品倾斜

先已实际运行 `probe.py`（NumPy），72个离散正密度/真实循环box平均案例保存在 `RESULTS.json`。本稿将数值、一般恒等式和Gaussian校准分开。它不封死一般conditioned input，也不给C_n上界。

## 1. 一般source的合法c区间（exact）

在概率空间dx上取非负密度ρ、∫ρ=1。固定0<π<1，定义两子密度

\[
\rho_1=\pi\rho+c(\rho-1),\qquad
\rho_0=(1-\pi)\rho-c(\rho-1).
\]

积分分别为π和1−π。非平凡情形设
`m=ess inf ρ<1<M=ess sup ρ`。若m>0，令

\[
a_s=\frac m{1-m},\qquad b_s=\frac M{M-1},
\]

M=∞时b_s按极限取1。两子密度非负a.e.的充要条件是

\[
\boxed{-\min\{(1-\pi)a_s,\pi b_s\}\le c
\le\min\{\pi a_s,(1-\pi)b_s\}.} \tag{1}
\]

证明只须将非负性除以ρ：posterior为 `π+c(1−1/ρ)`，而 `1−1/ρ` 的ess范围为 `[−1/a_s,1/b_s]`；检查两端是否落在[0,1]。此条件也必要，因为ess端点可由正测度集合任意逼近。

若m=0，合法性强迫c=0，包括ρ只趋近0而未实际取0的情形。若ρ=1 a.e.，任何c给相同平凡split，posterior恒π；(1)的m,M表达式不用于此退化情形。

特别地，对任意合法c，不论正负，乘合相应两个端点限制得

\[
\boxed{\frac{c^2}{\pi(1-\pi)}\le a_s b_s.} \tag{2}
\]

小m会限制合法c；不是可以任取固定强相关系数的自由通道。

## 2. 必须先source合法，再对query取平均（exact）

令A为正线性、保持常数和质量的真实平均算子，Y=Aρ；假定Y>0。平均后的子密度为

\[
Y_1=\pi Y+c(Y-1),\qquad Y_0=(1-\pi)Y-c(Y-1).
\]

因此在anchor查询律 `dQ=Y dx` 下，标签Z的边际仍是Bernoulli(π)，其真实posterior是

\[
\boxed{p_Q(x)=P(Z=1\mid x)=\pi+c\frac{Y(x)-1}{Y(x)}.} \tag{3}
\]

这里用Y，不能拿source点值ρ代替。source split先满足(1)，平均后(3)自动位于[0,1]。

定义二元标签到query的最大相关系数r。对任意实平方可积g(x)，
`Cov_Q(Z,g)=Cov_Q(p_Q,g)`，故Cauchy–Schwarz给最大值，且以g=p_Q达到：

\[
\boxed{r^2=\frac{\operatorname{Var}_Q(p_Q)}{\pi(1-\pi)}
=\frac{c^2}{\pi(1-\pi)}\int\frac{(Y-1)^2}{Y}\,dx.} \tag{4}
\]

因为∫Y=1，E_Q p_Q=π。等价写法是 `c²[∫Y^{-1}dx−1]/[π(1−π)]`，但Y接近1时数值上应使用(4)的非负式，避免相消。

若 `|Y−1|≤h<1`，则

\[
r^2\le\frac{c^2}{\pi(1-\pi)}\frac{\int(Y-1)^2}{1-h}
\le\frac{a_s b_s h^2}{1-h}. \tag{5}
\]

对于log-score，只要方差非零，准确有

\[
\boxed{\operatorname{Corr}_Q(Z,\log Y)
=r\operatorname{Corr}_Q(p_Q,\log Y).} \tag{6}
\]

因此 `Corr(p_Q,log Y)≈1` 可以与极小的真实标签通道r同时发生；两者不能混同。输入点律ρdx下的通道强度另是把(4)中的Y换成ρ，不是query的r。

## 3. 已运行的离散正密度检验（numerical）

脚本使用K=31/127/511个粗格，每格16个bin、一个随机peak；每格质量相同，background m=.1/.01/10⁻⁸/10⁻¹⁴，π=.5/.01。每格密度为m加一个高度 `(1−m)*16` 的peak。选择合法上端的95%作为c；三个尺度附近的实际循环box窗口保留完整周期和余项，共72个模型。

它们是真正有限离散正密度及box平均，不是假设Gaussian相关。实际窗口长度向下取整后以其真实有理L记录；不宣称等于连续目标L。数值核验ρ₀、ρ₁非负，并验证指数倾斜恒等式。

m=10⁻¹⁴、π=.5、目标L=7/6的代表结果：

| K | 真实query通道r | Corr(p_Q,log Y) |
|---:|---:|---:|
| 31 | 1.62185×10⁻¹⁶ | .99994255 |
| 127 | 3.78140×10⁻¹⁷ | .99999649 |
| 511 | 9.12780×10⁻¹⁸ | .99999978 |

这里p_Q−π按(3)的稳定表达式计算，未先把非常小的posterior偏移加回π后再相减；后者会被浮点舍入吞掉。这些小数是计算证据，不是区间端点；严格结论由(1)–(6)给出。

## 4. 真实接受source的r=O(1/K)包络（exact envelope）

46–50轮源为 `ρ=1−p+p μ_{K,δ}+(p/K)D_M`，峰宽正、每峰内部 `μ=wv`，峰间有正测度gap且μ=0。记

\[
\ell=\frac{20p}{3K}H_M,
\quad m_{\rm up}=1-p+\ell<1,
\quad M_{\rm low}=1-p+pwv-\ell>1.
\]

gap内ρ≤m_up，peak内ρ≥M_low，所以其真实ess端点满足 `m≤m_up`、`M≥M_low`。因此(2)给每一个合法π,c的统一界

\[
\frac{c^2}{\pi(1-\pi)}
\le\frac{m_{\rm up}}{1-m_{\rm up}}
\frac{M_{\rm low}}{M_{\rm low}-1}.
\]

而真实平均满足 `Y_L=1+p(H_L+B_L)/(KL)`，`|H_L|≤2`、`|B_L|≤40/9`。在L≥L_min>0时

\[
h=\frac{p(2+40/9)}{K L_{\min}}=O(1/K).
\]

由(5)，这些合法soft splits的query最大相关系数r=O(1/K)。本式对π也统一；不需要为了小π而暗中允许常数c违反source合法性。

在p=.99、n=10³²、K=2⁴⁷−1、w64、M40、L=7/6下，脚本用Fraction保存上式给出的r²有理上界；显示平方根为

\[
r_{\rm upper}\approx 3.90525\times10^{-15}
\]

（显示值不是向外舍入端点；有理r²才是认证对象）。另选π=.5、合法上端95%的一个实际安全c时，r上帽约3.72866×10⁻¹⁶。

## 5. 指数标签倾斜只给同一方向的产品重参数（exact）

在n个独立source-label对上，以 `exp(t Σ_i Z_i)` 加权。单坐标归一化常数 `Z_t=1−π+πe^t`，归一化源是

\[
\boxed{\rho_t=\frac{\rho_0+e^t\rho_1}{Z_t}
=1+\alpha_t(\rho-1),\quad
\alpha_t=1+\frac{c(e^t-1)}{1-\pi+\pi e^t}.} \tag{7}
\]

所以n维来源精确为ρ_t的n重产品，没有新增跨坐标耦合；查询输出也是 `Y_t=1+α_t(Y−1)` 的产品。端点参数为 `α_1=1+c/π`、`α_0=1−c/(1−π)`，α_t在两者之间。它由合法子密度的正组合得到，自动合法。

若原密度写作 `ρ=1+pR`，它只把该固定方向的幅度改成p′=α_t p。在α_t≥0的常用参数范围，可解释为产品spikiness重调；一般合法c可能允许负α，故更准确的陈述是“同一仿射方向的产品参数”。

固定计数门、阈值门、空间依赖标签等不是(7)的纯指数权，不能凭本恒等式称其全部等于产品。root的固定计数比较另需其Poisson-binomial论证；本稿不重复该证明。

## 6. Gaussian rare gate为何不能直接代入：量化条件与边界

取n个独立query-label对。在实际anchor律下，标准化计数
`U_n=(ΣZ_i−nπ)/sqrt(nπ(1−π))`
和标准化加性log-score
`T_n=(Σ(log Y_i−E_Q log Y))/sqrt(n Var_Q(log Y))`
满足

\[
|\operatorname{Corr}(U_n,T_n)|
=|\operatorname{Corr}_Q(Z,\log Y)|\le r=O(1/K). \tag{8}
\]

n没有把这项相关放大；因此不能任取.9之类固定Gaussian相关，称它代表该真实soft通道。

更一般，对只依赖标签向量的任意事件G，由二元Z的一阶条件均值公式和独立坐标，准确有

\[
E[T_n\mid G]
=\operatorname{Corr}_Q(Z,\log Y)E[U_n\mid G]. \tag{9}
\]

所以，若某门已被独立证明 `|E[U_n|G]|=O(v^s)`（固定s），则在基线K超过任意v多项式时，其标准化score均值移动=o(1)。这是一条带明确前提的量化条件，不把Gaussian尾近似当成真实输入证据。

仅为Gaussian校准：若两个标准Gaussian变量相关为ρ_g，则上尾门U>u给score均值移动 `ρ_g φ(u)/Q(u)`。对u>0，Mills界给 `φ(u)/Q(u)<u+1/u`。取u=√v并以r包络代入，本轮参数移动上帽约1.53137×10⁻¹³；相应平方有理式 `r²(v+2+1/v)`保存在JSON。该数值不证明实际rare-event概率迁移。

小r**不**自动禁止所有极稀conditioning：例如所有标签均为1时，`U_n`为O(√n)，(9)的右侧未必小；该门恰给端点产品ρ₁/π。其他硬门仍需独立分析。本文没有给一般conditioned input或C_n增长作否定结论。

复现：

```sh
/opt/homebrew/opt/python@3.15/bin/python3.15 probe.py
```

只用已有NumPy和标准库；输出本目录JSON，不修改旧轮。一般公式、解析r²包络与指数恒等式是严格结果；离散相关小数和Gaussian门移动是分别标注的数值/参考模型。
