# 第 28 轮独立审查：非有界角向轮廓与有限弱型分类

## 裁决

候选放宽成立。设 (a) 是零次齐次的可测非负角向函数，

\[
0<a_*:=\operatorname*{ess\,inf}_{\sigma_K}a,
\qquad 0<E_a:=\int a\,d\sigma_K<\infty,
\]

并定义

\[
A_{\rm loc}(a):=\sup_{\substack{Q=x+rK\\r<\rho(x)/2}}
 \frac{\operatorname{avg}_Q a}{\operatorname*{ess\,inf}_Q a}<\infty.
\tag{R28.1}
\]

这里不须 (a\in L^\infty(\sigma_K))。则对每个 (0<\tau<n)，
\[
\lim_{\tau\downarrow0}\tau W^{\rm out}_{K,\rho^{\tau-n}a}
=\frac{nE_a}{a_*}.
\tag{R28.2}
\]
更强地，若另有 (0<a<\infty) \(\sigma_K\)-a.e. 且 (0<E_a<\infty)，那么
\[
W^{\rm out}_{K,\rho^{\tau-n}a}<\infty\ \text{对某个 }0<\tau<n
\quad\Longleftrightarrow\quad
a_*>0\ \text{且}\ A_{\rm loc}(a)<\infty,
\tag{R28.3}
\]
并且左边一旦成立，便对所有 (0<\tau<n) 成立。

## 充分性：三尺度证明不依赖全局上界

设 (h_\tau=\rho^{\tau-n}a)。有限 (E_a) 与极坐标公式已经给出

\[
\int_{sK}h_\tau\,dx=\frac{nVE_a}{\tau}s^\tau<\infty,
\]

所以 (h_\tau\) 局部可积；(a_*>0) 仍使有限 (L^1(h_\tau)) 成本的输入自动局部
Lebesgue 可积。对局部见证 (Q=x+rK)、(r<\rho(x)/2)，有

\[
\frac{\operatorname{avg}_Qh_\tau}{\operatorname*{ess\,inf}_Qh_\tau}
\le3^{n-\tau}\frac{\operatorname{avg}_Qa}{\operatorname*{ess\,inf}_Qa}
\le3^{n-\tau}A_{\rm loc}(a).
\tag{R28.4}
\]

因此原 (SW.5) 只需替换为

\[
\mu_\tau(H_0)\le B_KA_{\rm loc}(a)3^{n-\tau}A.
\tag{R28.5}
\]

中尺度 Tonelli 只使用 (E_a<\infty)；中央尺度只使用 (a_*>0) 和 (E_a<\infty)。
所以原 (SW.6)、(SW.7) 逐字保留。乘以 \(\tau\) 后局部项消失，故上极限仍是
(nE_a/a_*)。低角锥区有正 Lebesgue 测度，且 (h_\tau\in L^1_{\rm loc})，故原子
平滑下界仍可在其中的非零 Lebesgue 点进行；得到同一个匹配下界。此论证不使用
(a\) 的连续性或有界性。

## 方体例的精算

令 (K=[-1,1]^n)、(n\ge2)、

\[
a(x)=\left(\frac{\rho_\infty(x)}{|x_1|}\right)^\gamma,
\qquad0<\gamma<1.
\]

有 (a_*=1)。cube cone measure 在 (2n) 个面上的质量均为 (1/(2n))：
两个 (x_1=\pm1) 面的角平均为 (1)，其余 (2(n-1)) 面的平均为

\[
\frac12\int_{-1}^{1}|u|^{-\gamma}\,du=\frac1{1-\gamma}.
\]

故

\[
E_a=\frac1n+\frac{n-1}{n(1-\gamma)}
=\frac{n-\gamma}{n(1-\gamma)},
\qquad
\lim_{\tau\downarrow0}\tau W_\tau=\frac{n-\gamma}{1-\gamma}.
\tag{R28.6}
\]

对于任意一维区间 (I)，(b(t)=|t|^{-\gamma}) 满足

\[
\frac{\operatorname{avg}_I b}{\operatorname*{ess\,inf}_I b}
\le\frac2{1-\gamma}.
\]

若 (I=[-u,v]) 穿过零点，令 (q=\min(u,v)/\max(u,v))，其比值为

\[
\frac{1+q^{1-\gamma}}{(1-\gamma)(1+q)}\le\frac2{1-\gamma};
\]

若不穿过零点，凹性给出更强的 (1/(1-\gamma)) 界。局部方体中
(\rho_\infty^\gamma) 的上下比至多 (3^\gamma)，于是

\[
A_{\rm loc}(a)\le\frac{2\,3^\gamma}{1-\gamma}.
\tag{R28.7}
\]

本机 Wolfram 15.0 精确给出
\(\frac12\int_{-1}^{1}|u|^{-\gamma}du=(1-\gamma)^{-1}\)，并化简
(E_a=(n-\gamma)/(n(1-\gamma))\)；Python `Fraction` 在
(2\le n\le11\) 的有理 \(\gamma\) 样本上独立通过。它们仅作代数交叉核验。

## 必要性

设 (W=W^{\rm out}_{K,h_\tau}<\infty)。对任意 (Q=x+rK) 及任意
Lebesgue 正测度 (E\subset Q)，每个 (z\in Q) 都有 (Q\subset z+2rK)，所以

\[
M_K\mathbf1_E(z)\ge\frac{|E|}{2^n|Q|}.
\]

对严格阈值从下逼近并应用完整弱型不等式，得

\[
\operatorname{avg}_Qh_\tau
\le2^nW\operatorname*{ess\,inf}_Qh_\tau.
\tag{R28.8}
\]

令 (Q=K)。由于
\(\operatorname*{ess\,inf}_Kh_\tau=a_*\)（在半径趋于 (1) 的薄环带上逼近），
且 (\operatorname{avg}_Kh_\tau=nE_a/\tau>0)，故 (a_*>0)。对于局部 (Q)，
(a=\rho^{n-\tau}h_\tau) 和径向上下比给

\[
A_{\rm loc}(a)\le2^n3^{n-\tau}W<\infty.
\tag{R28.9}
\]

这既证明必要性，也表明条件没有偷用待证外层常数；它只是由有限弱型推出的标准
中心 (A_1) 测试。

## 关于删除 (a_*>0)

在 (R28.3) 所列的 (0<a<\infty) a.e. 假设下，若 (a_*=0)，对每个
\(\epsilon>0\) 可在 (\{a<\epsilon\}\) 中取非零 (h_\tau) Lebesgue 点。原子平滑
下界的分母含 (a(z)\le\epsilon)，而固定的中央开球有正
\(\mu_\tau\)-质量；令 \(\epsilon\downarrow0\) 即 (W=\infty)。因此这不是仅有的
小 \(\tau\) 现象，而是每个 (0<\tau<n) 都发散。

若把表述再放宽为允许 (a=0) 于正角测度集，则上述“分母含 (a(z))”的直接句子不够：
需要在零成本（或趋零成本）原子块外添加一个任意小的、正而有限成本的块，以保持定义域
中的分母严格为正。结论仍可修复，但正式分类最好保留 (0<a<\infty) a.e.，或显式写出
这个辅助块步骤。

## 标准 A1 的等价简化

令 `[a]_{A_1(K)}` 为全部平移 K 球的 `avg_Q a / essinf_Q a` 上确界。对零次齐次 a，
(R28.1) 加 `a_*>0, E_a<∞` 等价于 `[a]_{A_1(K)}<∞`：若 `r≥ρ(x)/2`，则

\[
Q=x+rK\subset3rK,\qquad \operatorname{avg}_Qa\le3^nE_a,
\qquad \operatorname*{ess\,inf}_Qa\ge a_*.
\tag{R28.10}
\]

小球情形正是 (R28.1)，故这一步未使用待证的弱型结论。

对 `h_τ=ρ^{τ-n}a`，严格成立

\[
W^{\rm out}_{K,h_\tau}<\infty
\quad\Longleftrightarrow\quad h_\tau\in A_1(K)
\quad\Longleftrightarrow\quad a\in A_1(K).
\tag{R28.11}
\]

其中有限 W 推出 `avg_Q h_τ≤2^nW essinf_Qh_τ`，故局部球上

\[
A_{\rm loc}(a)\le3^{n-\tau}[h_\tau]_{A_1(K)}.
\]

另一方面，Q=K 给出 `nE_a/(τa_*)≤[h_τ]_{A_1(K)}`；合并 (R28.10) 得

\[
[a]_{A_1(K)}\le3^n[h_\tau]_{A_1(K)}.
\tag{R28.12}
\]

反过来 `a∈A_1(K)` 蕴含本记录的充分条件，三尺度上界给 W 有限；而
`h_τ∈A_1(K)` 到 W 有限则是逐球平均比较加测度一致弱型。因而 (R28.11) 没有漏洞。

在主三尺度式中，局部项可简洁写为

\[
B_K[a]_{A_1(K)}3^{n-\tau}.
\]

所以原 SW.8 的统一类别可扩大为 `[a]_{A_1(K)}≤R`，同一显式
`C_{n,R}` 不变；原来的有界对比 `a^*/a_*≤R` 自动蕴含这个条件。

## 有限角均值不足的显式反例

方体 (K=[-1,1]^2) 上，取

\[
a=1+\sum_{k\ge1}2^k\mathbf1_{\{x_2>0,\ x_1/x_2\in I_k\}},
\qquad I_k=(2^{-k},2^{-k}+2^{-3k-4}).
\tag{R28.13}
\]

这些锥都位于 `|x_1|<x_2`。直接按上半方形锥面积计算，

\[
E_a=1+\sum_{k\ge1}\frac{2^k|I_k|}{8}
=1+\sum_{k\ge1}2^{-2k-7}=\frac{385}{384},\qquad a_*=1.
\]

然而在 (z_k=(2^{-k},1)) 取足够小的中心方体，避开各条锥的上边界及其他锥，
直线 (x_1=2^{-k}x_2) 经过中心并把方体等分，故

\[
\frac{\operatorname{avg}_{Q_k}a}{\operatorname*{ess\,inf}_{Q_k}a}
\ge1+2^{k-1}\longrightarrow\infty.
\]

于是 `[a]_{A_1(K)}=∞`，由 (R28.11) 对每个 (0<τ<2) 都有 (W_τ=∞)。
这验证仅假定 (E_a<∞,a_*>0) 的推广是错误的。并且无界例不能自动继承主稿中
只在有界角向轮廓下陈述的 AG.7/Sτ→1 端点结论。本机 Wolfram 15.0 与 Python
`Fraction` 均独立求得该级数为 `1/384`，故 `E_a=385/384`；这只核验算术，
不替代上述局部方体构造。

## 以半径比作 dyadic 分层的改进速率

此改进成立。令 `t_j=2^{j-1}` (`j=0,...,N-1`) 且
`N=ceil(log_2(2/eta))`。区间

\[
t_j\rho(x)\le r<2t_j\rho(x)
\]

覆盖中尺度 `rho(x)/2≤r<rho(x)/eta`；在恰等于分层端点或分类端点时，严格平均允许
微扰半径，或采用半开层，均不改变覆盖。见证给

\[
x+rK\subset(1+2t_j)\rho(x)K,
\qquad \int_{x+rK}f>Vt_j^n\rho(x)^n.
\]

因此对该层可用

\[
T_j(x)=\frac1{Vt_j^n\rho(x)^n}
\int_{(1+2t_j)\rho(x)K}f>1.
\]

Tonelli 的内层仍是 `rho(x)^{tau-2n} a(x)`，从而

\[
\mu_\tau(H_j)\le
\frac{nE_a}{a_*(n-\tau)}
\frac{(1+2t_j)^{n-\tau}}{t_j^n}A
\le\frac{nE_a}{a_*(n-\tau)}4^nA.
\tag{R28.14}
\]

末个不等式的完全理由是 `t_j≥1/2` 时 `1+2t_j≤4t_j`，且
`(1+2t_j)^{n-tau}≤(1+2t_j)^n`。本机 Wolfram 15.0 化简为同一不等式，
Python 在 `n=1,...,30` 的广泛参数网格作了独立数值压力检查；严格依据是前一句代数。

故归一化后中尺度总项至多

\[
\frac{4^n\tau}{n-\tau}\,N.
\tag{R28.15}
\]

取 `eta=tau`、`0<tau<min(1,n/2)`，局部项为 `O_{n,[a]A1}(tau)`，中央项不超过
`(1+tau)^n=1+O_n(tau)`，而 (R28.15) 为 `O_n(tau log(1/tau))`。结合
`1-L_0≤tau[1+log(n/tau)]`，得到上下两侧同为
`O_{n,[a]A1}(tau log(1/tau))` 的一致相对速率。该分层确实删去了旧估计中
`eta^{-n}` 的损失，适合替换而非叠加旧 SW.8。

## 一般 q 壳的对数首项

审查 `round28-rate/shell_refinement.tex` 后确认：令 `s_j=q^j/2`、
`N=ceil(log_q(2/eta))`、`0≤j<N`，即使 q 非整数也有 `q^N/2≥eta^{-1}`，
所以半开壳覆盖中尺度；严格见证允许处理全部端点。每层的精确代价是

\[
\frac{(1+qs_j)^{n-\tau}}{s_j^n}
=s_j^{-\tau}(q+s_j^{-1})^{n-\tau}.
\]

对固定 `n,q,c`、`eta=c tau`，几何级数与
`sum_{j<N} log(s_j)=O_q(N^2)` 给出

\[
\sum_{j<N}\frac{(1+qs_j)^{n-\tau}}{s_j^n}
=q^nN+O_{n,q}(1+\tau N^2),\qquad
N=\frac{\log(1/\tau)}{\log q}+O_{q,c}(1).
\]

中央项为 `1-tau log(1/tau)+O_{n,c}(tau)`，故归一化上界对全体
K 及 `[a]A1≤R` 一致为

\[
\frac{\tau W_\tau}{L}
\le1+\left(\frac{q^n}{n\log q}-1\right)\tau\log(1/\tau)
+O_{n,R,q,c}(\tau).
\]

取 `q=exp(1/n)` 后与原子下界合并，确有

\[
-1\le\liminf\frac{\tau W_\tau/L-1}{\tau\log(1/\tau)}
\le\limsup\frac{\tau W_\tau/L-1}{\tau\log(1/\tau)}\le e-1.
\]

这里不是精确高阶系数。唯一文字收紧点是：从
`O_{n,R}(tau log(1/tau))` 推出 `O_{n,R}(sqrt(tau))` 时，只能说指数 `1/2`
不依赖维数；隐含常数仍依赖 `n,R`，不能称完全“维数无关”。

## 仓库整合稿终审

已逐行审读当前 `paper/small_homogeneity_weak.tex`。未发现实质数学或符号问题。

- SA.1--SA.2 的有限性判别闭合：子集测试给 `[h_tau]A1≤2^n W`，逐球比较给反向；
  原点球得到 `a_*>0`，小/大球分别接回 `a∈A1`，故“某个 tau / 每个 tau”量词正确。
- SW.6 的逐壳 Tonelli 仍使用正确的 `rho(x)^(tau-2n)`（极坐标后为
  `s^(tau-n-1) ds`）；q 层覆盖和 SW.8 的 `4^n N_2/(n-tau)` 系数均正确。
- SW.8a 的 binomial--geometric 上界、`q=exp(1/n)` 的 `e-1` 系数、中央项负
  `tau log(1/tau)` 以及 lower/upper 的量词均一致；未误说精确二阶系数。
- SA.2 的方体积分确为 `(n-gamma)/(n(1-gamma))`，局部/大球 A1 界正确；尖峰例
  的 `385/384` 与非 A1 推论也正确。无界轮廓未被错误地接入只对有界 a 陈述的
  `S_tau→1` 结论。

本机 Wolfram 15.0 重新化简 SA.2 积分为 `(n-gamma)/(n(1-gamma))`，Python 对
q 壳的 binomial 上界作了参数压力检查；这些只作代数交叉核验。
