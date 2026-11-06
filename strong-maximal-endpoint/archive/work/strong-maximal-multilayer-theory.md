# 多层双曲轮廓：固定矩形见证的严格归约

日期：2026-11-06。本文只处理附件所述的双曲轮廓和一组固定的矩形比例。它给出的是强极大算子的认证下界，不能被解读为对所有矩形或所有输入的上界。

## 1. 固定维数的截断极限

设
\[
 S(x)=\sum_{i=1}^n-\log|x_i|,\qquad
 f_L(x)=h(S(x)-L)\mathbf 1_{[-1,1]^n}(x).
\]
这里先假定 `h` 是有限个台阶组成的非减函数：它在某个有限的 \(a_0\) 以下为零，且在某个有限的 \(a_J\) 以上为常数。这个假定同时保证所有下面的积分有限；它已覆盖三至五段阶梯的第一轮搜索。

对任意非负可测 \(F\)，对数坐标的精确公式是
\[
 \int_{[-1,1]^n}F(S(x)-L)\,dx
 =\frac{2^ne^{-L}}{(n-1)!}\int_{-L}^{\infty}
 e^{-u}(L+u)^{n-1}F(u)\,du. \tag{1}
\]
令 \(F=\Phi_n\circ h\)。由于 \(F\) 在左边最终为零、在右边有界，(1) 除以
\(2^ne^{-L}L^{n-1}/(n-1)!\) 后可由支配收敛得到
\[
 \lim_{L\to\infty}\frac{(n-1)!e^L}{2^nL^{n-1}}
 \int\Phi_n(f_L)
 = I(h):=\int_{\mathbb R}e^{-u}\Phi_n(h(u))\,du. \tag{2}
\]
这一极限的顺序是**固定 \(n,h\)** 后 \(L\to\infty\)。若允许随 \(L\) 变化的轮廓，或无界轮廓，(2) 不能直接使用，必须另给一致可积估计。

## 2. 矩形见证与认证比值

固定 \(r=(r_1,\ldots,r_n)\in[0,1]^n\)，令
\[
 T_r=\sum_{i=1}^n Z_{r_i},\qquad
 Z_r=-\log |U|,\quad U\sim\mathrm{Unif}[-r,1].
\]
在 \(x_i\ne0\) 处取矩形边 \([-r_ix_i,x_i]\)。若 \(S(x)\ge L-q\)，其平均值为
\[
 \mathbb E\,h(S(x)-L+T_r)\ge \mathbb E\,h(T_r-q), \tag{3}
\]
其中不等式只用了 \(h\) 非减。因而一旦
\[
 \mathbb E\,h(T_r-q)>1, \tag{4}
\]
就有严格见证 \(E_{L-q}\subset\{M_Sf_L>1\}\)（零测集略去）。又
\[
 |E_{L-q}|=2^ne^{-(L-q)}\sum_{k=0}^{n-1}\frac{(L-q)^k}{k!},
\]
所以 (2) 给出下列可认证的固定轮廓极限下界：
\[
 \liminf_{L\to\infty}
 \frac{|\{M_Sf_L>1\}|}{\int\Phi_n(f_L)}
 \ge \frac{e^q}{I(h)}. \tag{5}
\]
若 (4) 只有等号，可将一层高度略增后取极限，或把 \(q\) 减小一个正量；记录严格余量可避免这个技术点。

## 3. 有限阶梯的精确凸规划

取分点
\[
 a_0<a_1<\cdots<a_{J-1}<a_J=\infty,
\]
令 \(h(u)=0\) 于 \(u<a_0\)，且 \(h(u)=c_j\) 于
\([a_j,a_{j+1})\)，其中 \(0\le c_0\le\cdots\le c_{J-1}\)。定义
\[
 w_j=e^{-a_j}-e^{-a_{j+1}},\qquad
 p_j(q,r)=\Pr\{T_r-q\in[a_j,a_{j+1})\}.
\]
这里 \(e^{-\infty}=0\)，且概率端点可以直接由 Gamma--二项卷积计算。则 (4) 和 (2) 分别成为
\[
 \sum_jp_jc_j>1,\qquad I(h)=\sum_jw_j\Phi_n(c_j). \tag{6}
\]
因此，对**固定** \((q,r,a_0,\ldots,a_{J-1})\)，最大化 (5) 等价于最小化
\[
 \min_{0\le c_0\le\cdots\le c_{J-1}}
 \sum_jw_j\Phi_n(c_j)
 \quad\text{s.t.}\quad \sum_jp_jc_j\ge1. \tag{P}
\]
这是有限维凸规划：对 \(t>1\)，
\[
 \Phi_n'(t)=(1+\log t)^{n-2}(n+\log t),\quad
 \Phi_n''(t)=\frac{(n-1)(n-1+\log t)(1+\log t)^{n-3}}{t}>0,
\]
而 \(\Phi_n(t)=t\) 于 \([0,1]\)，在 \(1\) 处右导数跃升为 \(n\)。

无单调性约束时，非退化 KKT 条件是
\[
 w_j\Phi_n'(c_j)=\lambda p_j. \tag{7}
\]
单调性约束只会把违反顺序的相邻坐标合并成块；在一个正且大于 1 的块 \(A\) 上，(7) 替换为
\[
 \Phi_n'(c_A)=\lambda\frac{\sum_{j\in A}p_j}{\sum_{j\in A}w_j}. \tag{8}
\]
这给出了可审计的 active-set/pooled-adjacent-violators 求解器。落在 \([0,1]\) 的坐标须单独处理，因为目标在线性段；除非 \(\lambda p_j/w_j=1\)，最优值不会停在该开区间。

## 4. 软层并不由 Jensen 排除：一个精确反例

“分层一定不优于单示性层”的 Jensen 论证不成立。即使保持固定矩形机制，软层可严格降低 (P) 的成本。

取 \(n=2\)、\(r_1=r_2=0\)，则 \(T_r\sim\Gamma(2,1)\)，其尾为
\(\overline F(t)=e^{-t}(1+t)\)。取 \(q=4\)，并只允许两个正台阶：
\[
 h(u)=\begin{cases}
 0,&u<0,\\
 c_0,&0\le u<2,\\
 c_1,&u\ge2.
 \end{cases}
\]
本例的 (6) 数据是
\[
 (p_0,p_1)=(0.07422692920700639,\,0.01735126523666451),
 \quad(w_0,w_1)=(0.8646647167633873,\,0.1353352832366127).
\]
单层基准 \(h_\mathrm{hard}=A\mathbf1_{[0,\infty)}\)、
\(A=(p_0+p_1)^{-1}\) 的成本是
\[
 I(h_\mathrm{hard})=37.02368351072274.
\]
用约束 \(p_0c_0+p_1c_1=1\) 与 KKT 条件 (7) 得到
\[
 (c_0,c_1)=(5.490850375320629,\,34.14339126422669),
 \quad I(h_\mathrm{soft})=33.76842682633889.
\]
故相同 \((q,r)\) 的认证比值严格提高了
\[
 I(h_\mathrm{hard})/I(h_\mathrm{soft})=1.096399417749743.
\]
这只是固定二维测试族的数值证据，不改变 \(n\to\infty\) 的任何结论；但它严格否定了“Jensen 自动将最优轮廓化为单层”的预期。它说明 A1 的三至五段搜索是必要的。

## 5. 可复现的 Wolfram 验算

Wolfram 15.0.0 实际运行下列高精度输入（输出为上节数字；`FindRoot` 的精度警告来自 Piecewise 的机器精度比较，可改写为纯符号分支消除）：

```wl
phi[x_?NumericQ] := If[x <= 1, x, x (1 + Log[x])];
dp[x_?NumericQ] := If[x < 1, 1, 2 + Log[x]];
tail[t_?NumericQ] := Exp[-t] (1 + t);
q = 4; a = 2;
p1 = tail[q + a]; p0 = tail[q] - p1;
w1 = Exp[-a]; w0 = 1 - w1;
sol = FindRoot[{p0 c0 + p1 c1 == 1,
  w0 dp[c0]/p0 == w1 dp[c1]/p1},
 {c0, 5.5}, {c1, 34}];
base = phi[1/(p0 + p1)];
soft = (w0 phi[c0] + w1 phi[c1]) /. sol;
N[{p0,p1,w0,w1,c0,c1,base,soft,base/soft}/.sol,25]
```

## 6. 下一步

1. 将 (P) 实现为等边长 \(r=B^{-1}\) 的 Gamma--二项概率的对数域求解器，并在
\(n=2,4,8,16,32,64,128,256\) 与同维单层 \(H_n\) 比较；每个候选以高精度分位数残差和分段截断误差复核。
2. 对可能随 \(n\) 增长的层数、断点与高度，补充统一可积条件。固定维度的 (2) 不能推出跨维渐近；这是把本节数值现象升级为次项改进前的主要理论缺口。

## 7. 指数斜坡：一个可严格认证的增长维数候选

下面的构造比有限台阶更关键。固定 \(0<a<1\)，在未平移的轮廓坐标取
\[
 h_{n,a}(u)=e^{au/n}\mathbf 1_{\{u\ge0\}}. \tag{14}
\]
它无界，但对每一个固定 \(n\) 都是允许的：
\[
 I_{n,a}=\int_0^\infty e^{-(1-a/n)v}(1+av/n)^{n-1}\,dv
 \le\int_0^\infty e^{-(1-a)v}\,dv=\frac1{1-a}. \tag{15}
\]
不等式只用 \((n-1)\log(1+av/n)\le av-a v/n\)。所以第 1 节的截断极限仍可由支配收敛使用，且给出有限常数的模数。

取等边长 \(r=B^{-1}\)，令 \(T=\sum Z_r\)，则 \(\mathbb ET=n\mu\)、
\(\operatorname{Var}T=n\mu\)。对
\[
 q=n\mu-\sigma\sqrt n\,z,\qquad \delta_n(z)=\Pr(T<q), \tag{16}
\]
边界矩形平均为
\[
 \mathbb E\left[e^{a(T-q)/n}\mathbf 1_{\{T\ge q\}}\right]
 \ge e^{a\sigma z/\sqrt n}-\delta_n(z). \tag{17}
\]
这是 Jensen 不等式应用于未截断指数后，再用
\(e^{a(T-q)/n}\le1\) 于 \(T<q\) 扣除得到的。因此，只要右端严格大于 1，\(E_{L-q}\) 是严格见证，并且
\[
 C_n\ge \frac{e^q}{I_{n,a}}\ge(1-a)e^q. \tag{18}
\]

若使用独立同分布变量 \(Z_r=Y+\mu J\) 的标准 Cramér 中偏差定理（其矩母函数在原点邻域有限），则在 \(z=o(n^{1/6})\) 时
\[
 \delta_n(z)=(1+o(1))\overline\Phi_{\rm N}(z)
 \le (1+o(1))\frac{e^{-z^2/2}}{z\sqrt{2\pi}}. \tag{19}
\]
给定固定的 \(a\in(0,1)\)，取
\[
 z_n^2=\log n-2\log\log n+C_a, \tag{20}
\]
其中 \(C_a\) 足够大，则 (19) 小于
\(e^{a\sigma z_n/\sqrt n}-1\sim a\sigma z_n/\sqrt n\)。由 (17)--(18) 得到
\[
 \log C_n\ge n\mu-\sigma\sqrt{n\log n}
 +\sigma\log\log n\sqrt{\frac n{\log n}}
 +O_a\!\left(\sqrt{\frac n{\log n}}\right). \tag{21}
\]
因此该构造严格超过附件中单层双曲示性函数的第三阶；它不改变
\(\sqrt{n\log n}\) 的首个缺口系数。若只使用无 Mills 因子的 Chernoff 界，只能安全取
\(z_n^2=\log n-\log\log n+C_a\)，相应得到 \(\sigma/2\) 的 \(\log\log n\) 系数。

严格证明 (21) 还应把 (19) 对本分布的条件完整写出或重新证明，不能把普通 CLT 代替为相对中偏差。这里变量有解析的矩母函数
\[
 \mathbb E e^{sZ_r}=\frac1{1-s}\left(\frac1{1+r}+\frac r{1+r}e^{\mu s}\right),\quad s<1,
\]
所以标准定理的 Cramér 条件成立。

Wolfram 15 的独立数值核验取 \(a=1/2\)、\(C_a=15\)，以 Gamma--二项精确和计算 (17) 左端和 (15)：

| \(n\) | \(z\) | \(\Pr(T<q)\) | 矩形平均 | \(I_{n,a}\) | \(q-\log I_{n,a}\) |
|---:|---:|---:|---:|---:|---:|
| 32 | 3.99748 | \(2.39136\times10^{-7}\) | 1.49862 | 1.94698 | 14.67603 |
| 64 | 4.03836 | \(1.46489\times10^{-6}\) | 1.33361 | 1.97136 | 44.61383 |
| 128 | 4.08574 | \(3.27674\times10^{-6}\) | 1.22803 | 1.98506 | 110.69177 |
| 256 | 4.13755 | \(4.77655\times10^{-6}\) | 1.15814 | 1.99237 | 251.74488 |

这些是数值证据；(17)--(18) 与已引用的中偏差条件才是认证链条。

## 8. (19) 的自包含 Mills 型界

本节给出足够证明 (21) 的直接尾界，因而不需要把外部中偏差定理作为黑箱。令
\(W=\mu-Z_r\)。在 \(r=e^{-\mu}\) 时
\[
 K(\lambda):=\log\mathbb Ee^{\lambda W}
 =\lambda\mu-\log(1+\lambda)
 +\log(1+r^{1+\lambda})-\log(1+r),\quad\lambda\ge0. \tag{22}
\]
若 \(p_\lambda=r^{1+\lambda}/(1+r^{1+\lambda})\)，则 \(p_\lambda<1/2\)，并且
\[
 K'''(\lambda)=-\frac2{(1+\lambda)^3}
 -\mu^3p_\lambda(1-p_\lambda)(1-2p_\lambda)<0. \tag{23}
\]
又 \(K(0)=K'(0)=0\)、\(K''(0)=\operatorname{Var}Z_r=\mu\)。故
\[
 K(\lambda)\le\frac\mu2\lambda^2\qquad(\lambda\ge0). \tag{24}
\]

在以 \(e^{\lambda\sum W_i-nK(\lambda)}\) 倾斜的律下，连续部分
\(Y_i\) 变为参数 \(1+\lambda\) 的指数变量，离散部分仍独立 Bernoulli（参数为
\(p_\lambda\)）。于是 \(\sum W_i\) 的密度是平移/反射 Gamma\((n,1+\lambda)\)
密度的凸组合。由 Gamma 密度的众数和 Stirling 下界，
\[
 \|f_\lambda\|_\infty
 \le\frac{1+\lambda}{\sqrt{2\pi(n-1)}}. \tag{25}
\]
取 \(d=\sigma\sqrt n z\)、\(\lambda=z/(\sigma\sqrt n)\)，倾斜恒等式与 (24)--(25) 给出
\[
\begin{aligned}
 \Pr\!\left\{\sum W_i\ge d\right\}
 &=e^{nK(\lambda)-\lambda d}
   \int_0^\infty e^{-\lambda y}f_\lambda(d+y)\,dy\\
 &\le\frac{(1+\lambda)\sigma\sqrt n}
 {z\sqrt{2\pi(n-1)}}e^{-z^2/2}. \tag{26}
\end{aligned}
\]
对于 (20)，右端等于
\[
 \left(\frac{\sigma e^{-C_a/2}}{\sqrt{2\pi}}+o(1)\right)
 \sqrt{\frac{\log n}{n}}. \tag{27}
\]
因此只要 \(C_a\) 使 \(e^{-C_a/2}/\sqrt{2\pi}<a\)，(17) 对充分大的 \(n\) 严格成立。这补齐了 (21) 所需的概率估计；其余步骤只是固定 \(n\) 的截断极限和矩形见证。

## 9. 全连续问题的精确对偶与当前障碍

令 \(Q(dt)=e^{-t}\mathbf1_{t\ge0}\,dt\)，\(dP=R,dQ\)。固定比例向量的全部单调轮廓问题是
\[
 J_* =\inf_{h\uparrow,\;\int h\,dP\ge1}\int\Phi_n(h)\,dQ. \tag{28}
\]
由于第 7 节的 \(R\) 单调，点态 KKT 解已经单调。Fenchel 对偶给出可直接计算的等式
\[
 J_* = \sup_{\lambda\ge0}\left\{\lambda-
 \int_0^\infty\Phi_n^*(\lambda R(t))e^{-t}\,dt\right\}. \tag{29}
\]
这里
\[
 \Phi_n^*(y)=
 \begin{cases}
 0,&y\le1,\\
 y-1,&1<y\le n,\\
 (n-1)e^{s-1}s^{n-2},&y>n,
 \end{cases}\tag{30}
\]
其中末行的 \(s\ge1\) 由 \(y=s^{n-2}(s+n-1)\) 唯一确定。它和原问题的
\(0/1/\)指数尾 KKT 轮廓完全等价；数值时可以用 (29) 与原始积分交叉核验。

对于等边长 \(r=B^{-1}\)，二项混合可精确消去指数权重：
\[
 R_n(t)=\frac{(1+r)^{-n}}{(n-1)!}
 \sum_{k=0}^{\min(n,\lfloor t/\mu\rfloor)}{n\choose k}(t-k\mu)^{n-1}. \tag{31}
\]
故每次评价只需 \(O(n)\) 个正项的对数域求和；用 \(\lambda\) 二分求
\(\int h_\lambda dP=1\)，再比较原始目标与 (29)，可形成严格的数值残差检查。

当前没有一个一般性证明表明 (28) 的最优值必保留
\(\sigma\sqrt{n\log n}\) 的系数。仅根据激活点 \(q\) 得到
\(J\ge e^{-q}\) 不够：第 7--8 节的斜坡可把大小为
\(\sqrt{\log n/n}\) 的缺失概率分散到远尾，而只付多项式因子。因而错误的下界
\(J\ge e^{-q}\exp(cn\delta)\) 被该斜坡反例否定。要证明首系数障碍，必须对 (29) 进行统一的鞍点/拉普拉斯分析，或给出所有 KKT 尾的可比较上包络；这仍是下一轮的主要理论问题。

## 10. 单一双曲接收者机制的统一首系数障碍

本节解决第 9 节所说的障碍，但结论只针对“一个双曲轮廓、一个固定比例向量、由
\(E_{L-q}\) 接收”的机制。设任意 \(r_i\in[0,1]\)，令
\(T=\sum_iZ_{r_i}\)，并在平移坐标中考虑任意非减轮廓 \(h\)。若
\[
 J=\int_0^\infty e^{-t}\Phi_n(h(t))dt\le e^{-s_0},
 \quad s_0=n\mu-d,
 \quad d=(\sigma-\varepsilon)\sqrt{n\log n}, \tag{32}
\]
则 \(\mathbb Eh(T)<1\) 对充分大的 \(n\) 成立。故这一类机制不能把
\(\sigma\sqrt{n\log n}\) 的首个缺口系数严格减小。

由单调性，\(J\ge e^{-t}\Phi_n(h(t))\)，所以
\[
 \Phi_n(h(t))\le e^{t-s_0}. \tag{33}
\]
取固定大常数 \(c\)。若 \(t\le s_0-c\log n\)，则 \(h(t)\le n^{-c}\)。若
\(0\le x=t-s_0\le n/4\)，令 \(y=\log h(t)\)（仅在 \(h\ge1\) 时需要），则
\[
 (n-1)\frac y{1+y}\le\log\Phi_n(h(t))\le x,
\]
故 \(y\le x/(n-1-x)\)，从而 \(h(t)\le1+C x/n\)。当 \(t\le s_0\) 时 (33) 直接给出 \(h(t)\le1\)。

记 \(A=\{T\le s_0-c\log n\}\)、\(D_+=\{T>s_0+n/4\}\)。上述点态界及
\(\mathbb E|T-\mathbb ET|\le C\sqrt n\) 给出
\[
 \mathbb Eh(T)\le1-(1-n^{-c})\Pr(A)
 +C\sqrt{\frac{\log n}{n}}+\mathbb E[h(T)\mathbf1_{D_+}]. \tag{34}
\]

尾项一致地指数小。先将其分为 \(h<e^{4T/n}\) 与补集。第一部分由固定的 Chernoff 参数
\(\theta>0\) 控制：小 \(\eta\) 时一致有
\(\log\mathbb Ee^{\eta Z_r}\le\eta m(r)+C\eta^2\le\eta\mu+C\eta^2\)，而
\(D_+\subset\{T\ge n\mu+n/8\}\)，故
\(\mathbb E[e^{4T/n}\mathbf1_{D_+}]le e^{-c_1n}\)。第二部分中
\(\Phi_n(h)\ge h(4t/n)^{n-1}\)，且
\[
 R(t)=e^tf_T(t)
 \le\frac{2^nt^{n-1}}{(n-1)!}.
\]
因此
\[
 \mathbb E[h(T)\mathbf1_{D_+\cap\{h\ge e^{4T/n}\}}]
 \le J\frac{2^n(n/4)^{n-1}}{(n-1)!}=e^{-c_2n}, \tag{35}
\]
因为该括号的 \(n\) 次根趋于 \(e/2<B\)。

最后令 \(M=\mathbb ET=n\mu-D\)，其中
\(D=\sum_i(\mu-m(r_i))\)。若 \(D>2d\)，Chebyshev 给出 \(\Pr(A)=1-o(1)\)。若
\(D\le2d\)，由 \(m\) 在唯一最大点的非退化性和紧性，
\(|v(r)-\mu|\le C\sqrt{\mu-m(r)}\)，所以
\(\operatorname{Var}T=n\mu+O(\sqrt{nD})=n\mu+o(n)\)。同时各坐标的三阶中心矩一致有界且
\(v(r)\ge1\)，故 Berry--Esseen 误差一致为 \(O(n^{-1/2})\)。最坏情形 \(D=0\) 时，
\((s_0-c\log n-M)/\sqrt{\operatorname{Var}T}\) 的绝对值至多
\((1-\varepsilon/(2\sigma))\sqrt{\log n}\)。于是
\[
 \Pr(A)\gg\sqrt{\frac{\log n}{n}}, \tag{36}
\]
它压倒 (34) 的正误差，矛盾于 \(\mathbb Eh(T)\ge1\)。

## 11. A2 分块张量示性函数的精确支配

本节只处理各块都是双曲**示性函数**、接收集合是各块接收集合的直积的机制；它不涉及软轮廓、max/k-of-m 或相关构造。把坐标分为固定的块 \(b\)，维数之和为 \(n\)。令
\(S_b\) 为块内对数半径和，
\[
 E^{\rm blk}=\bigcap_b\{S_b\ge L_b\},\qquad
 g=\alpha^{-1}\mathbf1_{E^{\rm blk}},
\]
其中对固定比例向量，\(T_b\) 是第 \(b\) 块的增量，
\[
 a_b=\Pr(T_b\ge q_b),\qquad
 \alpha=\prod_ba_b,\qquad q=\sum_bq_b. \tag{37}
\]
块积矩形在接收集合 \(\cap_b\{S_b\ge L_b-q_b\}\) 上的平均至少为 \(\alpha\)，并在其内部严格大于 \(\alpha\)。令每个固定维数的 \(L_b\to\infty\)，则实际测度比趋于 \(e^q\)。注意分母必须使用原问题的完整模数，而不是块模数的乘积：
\[
 \int\Phi_n(g)=|E^{\rm blk}|\,\alpha^{-1}
 (1+\log(1/\alpha))^{n-1}.
\]
故分块机制的极限认证得分恰为
\[
 \frac{\alpha e^q}{(1+\log(1/\alpha))^{n-1}}. \tag{38}

另一方面，对同一组坐标比例，令 \(T=\sum_bT_b\)。则
\[
 a:=\Pr(T\ge q)\ge\Pr(\forall b:\,T_b\ge q_b)=\alpha. \tag{39}
\]
普通单层双曲示性函数 \(\mathbf1_{\{\sum_bS_b\ge L\}}\)，以 \(a^{-1}\) 缩放后，在相同的 \(q\) 给出得分
\[
 \frac{a e^q}{(1+\log(1/a))^{n-1}}. \tag{40}
\]
函数 \(x\mapsto x/(1+\log(1/x))^{n-1}\) 在 \((0,1]\) 上严格递增，故 (40) 至少为 (38)。这证明了此 A2 子族被单层无权双曲示性函数精确支配。Gamma--Bernoulli 和连续，边界等号只落在零测层上；取开接收集合即可满足严格超水平定义。

## 12. 单接收者单调轮廓的第三项闭合

令 \(Q_n\) 是第 10 节单一双曲接收者、任意比例向量和任意非减轮廓所得的最优极限得分。记
\(L=\log n\)、\(\rho=\sqrt{L/n}\)。第 7--8 节的指数斜坡与下面的统一上界合在一起给出
\[
 \log Q_n=n\mu-\sigma\sqrt{nL}
 +\sigma\{\log L+\log(\sigma\sqrt{2\pi})\}\sqrt{\frac nL}
 +o\!\left(\sqrt{\frac nL}\right). \tag{41}
\]

**下界。** 取 \(r=B^{-1}\)、\(h(u)=e^{u/n}\mathbf1_{u\ge0}\)。令 \(z\) 由
\[
 \overline\Phi_{\rm N}(z)=(1-L^{-1/4})\sigma\rho
\]
定义。普通 Berry--Esseen 的绝对误差是 \(O(n^{-1/2})=o(L^{-1/4}\rho)\)，故第 7 节的 Jensen 截断条件严格成立。Mills 公式给出
\[
 z^2=L-2\log L-2\log(\sigma\sqrt{2\pi})+o(1).
\]
第 7 节的 \(I_n\le2n\) 只贡献 \(O(\log n)=o(\sqrt{n/L})\) 的对数损失，遂得 (41) 的下界。

**上界。** 反设对某固定 \(C>\sigma\log(\sigma\sqrt{2\pi})\)，存在可行轮廓满足
\[
 J\le e^{-s},\qquad
 s=n\mu-\sigma\sqrt{nL}+\sigma\log L\sqrt{\frac nL}
 +C\sqrt{\frac nL}. \tag{42}
\]
第 10 节的点态界改用截止 \(s+n^{3/4}\)。对 \(0\le t-s\le n^{3/4}\)，直接反演
\(\Phi_n\) 给出
\[
 h(t)-1\le(1+O(n^{-1/4}))(t-s)/n.
\]
于是中段的期望超额至多
\((1+o(1))d/n=(1+o(1))\sigma\rho\)，其中 \(d=n\mu-s\)。其后的低高度部分以
\(h<e^{4t/n}\) 分割，在参数 \(\theta_n=cn^{-1/4}\) 的 Chernoff 界下为 \(e^{-c\sqrt n}\)；高高度部分仍由 (35) 控制。故可行性强制
\[
 \Pr\{T\le s-2\log n\}\le(1+o(1))\sigma\rho. \tag{43}
\]

沿用第 10 节的 \(D=n\mu-\mathbb ET\) 分类。若 \(D>2d\)，Chebyshev 与 (43) 矛盾。若 \(D\le2d\)，则
\(V=\operatorname{Var}T=n\mu+O(\sqrt{nD})\)，并且统一 Berry--Esseen 适用。阈值的标准化偏差满足
\[
 z_*^2\le\frac{d^2}{\mu n}+o(1)
 =L-2\log L-2C/\sigma+o(1).
\]
Mills 公式连同 Berry--Esseen 因而给出
\[
 \Pr\{T\le s-2\log n\}
 \ge\left(\frac{e^{C/\sigma}}{\sqrt{2\pi}}+o(1)\right)\rho,
\]
这与 (43) 在 \(C>\sigma\log(\sigma\sqrt{2\pi})\) 时矛盾。证明了上界，并闭合 (41)。

## 13. 锚定矩形标量轮廓的重排定理

这一节只讨论从 \([-r_ix_i,x_i]\) 形成的**锚定**矩形族；它不把任意非单调标量轮廓的真实强极大算子化为该族。令
\[
 A_rh(u)=\mathbb Eh(u+T_r),\qquad x=e^{-u},\qquad g(x)=h(-\log x).
\]
在该变量下 \(e^{-u}du=dx\)。若 \(V=e^{-T_r}\)，则其密度为
\[
 k_r(v)=R_r(-\log v),\quad0<v<1,
\]
其中 \(R_r(t)=e^tf_{T_r}(t)\) 非减，所以 \(k_r\) 非增。故
\[
 A_rg(x)=\int_0^1g(xv)k_r(v)dv
 =\frac1x\int_0^xg(y)k_r(y/x)dy. \tag{44}
\]
对固定 \((r,x)\)，(44) 的权重是 \((0,\infty)\) 上的非增函数。设 \(g^*\) 为非增重排，Hardy--Littlewood 不等式给出
\[
 A_rg(x)\le A_rg^*(x)\quad\hbox{对每个 }r,x. \tag{45}
\]
因为 \(\Phi_n\ge\mathrm{id}\)，有限模数保证 \(g\in L^1\)，且重排保持 \(\int\Phi_n(g)\)。同一个 \(g^*\) 同时适用于所有 \(r\)，所以 (45) 在取 \(\sup_r\) 后仍成立。

又 \(g^*\) 非增，故 \(A_rg^*(x)\) 及其 \(r\) 上确界都非增。其超水平接收集合因而是一个初始区间 \((0,x_0)\)（端点无需处理）。对任意 \(x<x_0\)，取一个满足 \(A_rg^*(x)>1\) 的 \(r\)，再定义
\(\widetilde g(y)=g^*(xy)\)。此时 \(A_r\widetilde g(1)>1\)，而其分母为 \(J/x\)。把
\(\widetilde g\) 在 \(y>1\) 截为零不改变这个边界平均（因 \(V\le1\)），只会减小分母；在 \(u\) 变量中所得轮廓正是非减且在 \(u<0\) 为零的第 12 节可行轮廓。因此
\[
 \frac{x}{J}\le Q_n.
\]
令 \(x\uparrow x_0\) 后，任意标量轮廓的锚定弱型比值不超过 \(Q_n\)。反向不等式来自任意边界可行的单调轮廓：其单个 \(r\) 在所有 \(x<1\) 上均严格超过阈值。故锚定标量轮廓族（包括非单调输入）的模块弱常数恰为 \(Q_n\)。

## 14. 逐坐标自适应比例的 Bellman 公式

本节不声称已经控制自适应策略。令一维算子
\(\mathcal Th(u)=\sup_{0\le r\le1}\mathbb Eh(u+Z_r)\)。迭代 \(\mathcal T^nh(0)\) 允许第 \(i\) 步的比例依赖于先前增量，故不能直接引用固定 \(r\) 的第 10--12 节结果。

对尾部值函数 \(V_0(q)=\mathbf1_{q\le0}\)，令
\[
 Gf(q)=\int_0^\infty e^{-y}f(q-y)dy.
\]
由 \(Z_r=Y+(-\log r)J\) 得到精确递推
\[
 V_k(q)=\sup_{0\le r\le1}
 \frac{GV_{k-1}(q)+rGV_{k-1}(q+\log r)}{1+r}. \tag{46}
\]
若内点可微，最优 \(r\) 满足
\[
 (1+r)V_{k-1}(q+\log r)-rGV_{k-1}(q+\log r)-GV_{k-1}(q)=0, \tag{47}
\]
其中用了 \((Gf)'=f-Gf\)。在最后一步 (47) 给出
\(r=e^{-q}\)，并重现 \(\sup_r\Pr(Z_r\ge q)=2/(1+e^q)\)。因此任何声称某个固定比例随机序支配所有策略的论证都是错误的。

另一方面，指数函数给出精确的风险敏感 Hamiltonian：对 \(0<\theta<1\)，
\[
 \mathcal Te^{\theta u}=M(\theta)e^{\theta u},\qquad
 M(\theta)=\frac1{1-\theta}\max_{0\le r\le1}
 \frac{1+r^{1-\theta}}{1+r}. \tag{48}
\]
内点最大元满足
\[
 r_\theta^\theta=1-\theta(1+r_\theta),\qquad
 r_\theta\to B^{-1},\qquad
 \log M(\theta)=\mu\theta+\frac\mu2\theta^2+O(\theta^3). \tag{49}
\]
故在中偏差共轭参数 \(\theta\asymp\sqrt{\log n/n}\) 上，指数/HJB 的最优动作接近固定的 \(B^{-1}\)，其首阶均值和方差没有改进。这支持但尚未证明自适应不能改变 \(\sqrt{n\log n}\) 系数。要完成该结论，需要对 (46) 建立统一的反集中下界或构造 Bellman/HJB 超解；固定比例的 Berry--Esseen 障碍不能直接移植。

## 15. 自适应控制的熵反集中引理

第 14 节的困难可用条件熵直接克服。以下结论适用于任意可预测比例策略
\(r_i=r_i(Z_1,\ldots,Z_{i-1})\)。令 \(P\) 为该策略的路径律，
\(d=(\sigma-\varepsilon)\sqrt{n\log n}\)，\(\omega=(\log n)^{1/4}\)，并设
\[
 \delta=\frac{d+\omega\sqrt n}{n},\qquad t=\mu-\delta.
\]
在每一步，若 \(m(r)\le t\) 则保持条件律；若 \(m(r)>t\)，将条件律按
\(e^{-\lambda Z_r}\) 倾斜，唯一选择 \(\lambda\ge0\) 使其均值为 \(t\)。把得到的路径律记为 \(Q\)。

当 \(m(r)>t\) 时，\(r\) 位于 \(B^{-1}\) 的紧邻域。设 \(\Delta=m(r)-t\le\delta\)。隐函数展开和该邻域的统一累积量界给出
\[
 \lambda=\Delta/v(r)+O(\delta^2),\qquad
 D(Q_r\Vert P_r)=\frac{\Delta^2}{2v(r)}+O(\delta^3)
 \le\frac{\delta^2}{2\mu}(1+O(\sqrt\delta)). \tag{50}
\]
这里用了 \(|v(r)-\mu|\le C\sqrt{\mu-m(r)}\)。链式法则即使策略自适应仍有效，故
\[
 D(Q\Vert P)\le\frac{n\delta^2}{2\mu}(1+o(1))
 =\frac{d^2}{2\mu n}(1+o(1)). \tag{51}
\]
在 \(Q\) 下每步条件均值不超过 \(t\)，条件方差一致有界。因此对
\(A=\{S_n\le n\mu-d-2\log n\}\)，Chebyshev 给出
\(Q(A)\ge1-O(\omega^{-2})\)。二元数据处理不等式给出
\[
 \log P(A)\ge-
 \frac{D(Q\Vert P)+\log2}{Q(A)}
 =-(1+o(1))\frac{d^2}{2\mu n}. \tag{52}
\]
特别地，\(P(A)\gg\sqrt{\log n/n}\)。

这正好与第 10 节从小模数 \(J\) 推出的 \(O(\sqrt{\log n/n})\) 失败概率上界冲突。自适应版本的两个远尾估计仍成立：每步密度满足 \(f_r(z)\le2e^{-z}\)，故联合密度及和密度分别受
\(2^ne^{-\sum z_i}\) 和 \(2^ne^{-t}t^{n-1}/(n-1)!\) 控制；而正矩母函数在所有 \(r\) 上一致地满足二阶界。这证明自适应锚定/复合轮廓也不能把首个
\(\sigma\sqrt{n\log n}\) 缺口系数降低。

这节本身不证明真实强极大算子的一般输入上界。若另有从 \(M_S\) 到该自适应复合轮廓常数的全局 rising-sun/重排归约，则它可与本节组合；该归约应作为单独引理陈述和核验。

## 16. 自适应熵界的 \(O(\sqrt n)\) 定量化

第 15 节可避免把 \(Q(A)\) 送到 1，从而得到统一的次级误差。固定足够大的常数
\(\omega\)，令
\[
 \delta=\frac{d+2\log n+\omega\sqrt n}{n},\qquad t=\mu-\delta.
\]
按第 15 节构造 \(Q\)，并记每一步的条件 KL 散度为 \(k_i\)。设
\(R_N=dQ/dP\)，则准确的似然恒等式为
\[
 \log R_N=\sum_{i=1}^nk_i-N,
 \qquad N=\sum_{i=1}^n\lambda_i(Z_i-m_{Q,i}). \tag{53}
\]
其中 \(N\) 是 \(Q\) 鞅，且 \(\mathbb E_QN^2\le Cn\delta^2\)。选定 \(\omega\) 后
\(Q(A)\ge7/8\)，其中 \(A=\{S_n\le n\mu-d-2\log n\}\)。因此在一个 \(Q\)-测度至少
\(3/4\) 的子集上，\(\log R_N\le K_{\max}+C\sqrt n\delta\)，并有
\[
 P(A)\ge\frac34e^{-K_{\max}-C\sqrt n\delta}. \tag{54}
\]

在整个范围 \(0\le d\le2\sigma\sqrt{n\log n}\) 内，(50) 的逐步界一致地给出
\[
 K_{\max}\le\frac{d^2}{2\mu n}+C\left(\frac d{\sqrt n}+1\right).
\]
这里 \(n\delta^2\sqrt\delta=o(1)\)，故 (54) 化为
\[
 P(A)\ge c\exp\left\{-\frac{d^2}{2\mu n}
 -C\left(\frac d{\sqrt n}+1\right)\right\}. \tag{55}
\]
若轮廓得分的对数大于
\(n\mu-\sigma\sqrt{n\log n}+K\sqrt n\)，则
\(d=\sigma\sqrt{n\log n}-K\sqrt n\)。(55) 的对数下界为
\[
 -\frac12\log n+(K/\sigma-C\sigma)\sqrt{\log n}+O(1),
\]
而小模数轮廓的点态分割只允许 \(O(\sqrt{\log n/n})\) 的失败概率。取固定且充分大的 \(K\) 即矛盾。故自适应锚定/复合轮廓常数的对数上界为
\[
 n\mu-\sigma\sqrt{n\log n}+O(\sqrt n). \tag{56}
\]

## 17. 停止策略与相对中偏差（待引用定理核验）

本节给出把第 16 节提升到第三项所需的自适应反集中方案。令
\(D_i=\sum_{j\le i}(\mu-m(r_j))\)，并在 \(D_i\) 首次越过
\(D_*=A\sqrt{n\log n}\) 后永久改用 \(r=B^{-1}\)，得到停止策略 \(K\)。因每步缺口不超过 \(\mu-1\)，
\(D_K\le D_*+\mu-1\)。若 \(M_K\) 是中心鞅，则
\[
 |\langle M_K\rangle_n-\mu n|
 \le C\sqrt{nD_*}=O(n^{3/4}(\log n)^{1/4}). \tag{57}
\]
各条件增量具有统一 Bernstein 矩界：指数部分具有阶乘矩，而 Bernoulli 位移满足
\(r|\log r|^j\le\sup_{a\ge0}e^{-a}a^j\le j!\)。归一化后的 Bernstein 参数为
\(O(n^{-1/2})\)。

因此，若所引用的鞅 Cramér 相对中偏差定理适用于 (57) 与这个条件 Bernstein 界，则在
\(0\le d\le2\sigma\sqrt{n\log n}\) 一致有
\[
 \Pr_K(M_K\le-d)=(1+o(1))\overline\Phi_{\rm N}
 \left(\frac d{\sigma\sqrt n}\right). \tag{58}
\]
耦合的方向如下。未停止时原策略与 \(K\) 相同，且 \(K\) 成功蕴含 \(M_K>-d\)。停止后原策略成功蕴含
\(M_\pi>D_*-d\)，其概率由一致条件矩母函数 Chernoff 界压至 \(n^{-3}\)（取 \(A\) 充分大）。故
\[
 \Pr_\pi(S_n\le n\mu-d)
 \ge(1+o(1))\overline\Phi_{\rm N}
 \left(\frac d{\sigma\sqrt n}\right)-n^{-3}. \tag{59}
\]
若 (58) 的外部定理条件被正式核验，(59) 与第 12 节的精细点态轮廓分割相结合，会把自适应复合轮廓的第三项也闭合到 (41) 的常数。当前唯一未内证的环节是 (58) 所用鞅相对中偏差定理的精确假设与常数范围。

## 18. 隐式 Gaussian 根的误差控制（依赖第 17 节的相对中偏差）

令 \(z_n>0\) 是唯一满足
\[
 \overline\Phi_{\rm N}(z_n)=\frac{\sigma z_n}{\sqrt n},
 \qquad d_n=\sigma\sqrt n\,z_n \tag{60}
\]
的根。第 17 节停止策略的方差误差给出的相对中偏差误差记为
\(\kappa_n=O(n^{-1/8}(\log n)^{5/8})\)。

对任意原策略的中心鞅 \(M\)，负倾斜累积量满足
\[
 \log\mathbb E(e^{-\lambda(Z_r-m(r))}\mid\mathcal F)
 \le\frac{V_{\max}\lambda^2}{2},\qquad
 V_{\max}=1+4/e^2. \tag{61}
\]
事实上二阶导数随负倾斜参数递减，且
\(v(r)\le1+r(\log r)^2\le V_{\max}\)。故
\[
 \mathbb E(-d-M)_+
 \le\frac{nV_{\max}}d\exp\left(-\frac{d^2}{2nV_{\max}}\right). \tag{62}
\]
在 \(d\asymp\sigma\sqrt{n\log n}\) 时，该项相对于 \(d\) 的大小是
\(n^{-\mu/(2V_{\max})}\operatorname{polylog}(n)\)，而
\(\mu/(2V_{\max})>1/8\)。这小于 \(\kappa_n\)。

于是第 12 节的中段反演（截止 \(n^{3/4}\)）和远尾界共同给出：若 \(J=e^{-s}\)、
\(d=n\mu-s\) 位于临界窗口，则
\[
 \overline\Phi_{\rm N}\!\left(\frac{d+2\log n}{\sigma\sqrt n}\right)
 \le\frac dn(1+O(\kappa_n)). \tag{63}
\]
根函数 \(H(d)=\overline\Phi(d/(\sigma\sqrt n))/(d/n)\) 在 \(d_n\) 的对数导数为
\(-\asymp\sqrt{\log n/n}\)。因此 (63) 强制
\[
 d\ge d_n-O\left(n^{3/8}(\log n)^{1/8}\right), \tag{64}
\]
其中 \(2\log n\) 平移的影响更小。

反向下界取独立的 \(r=B^{-1}\) 与指数斜坡。相对中偏差误差
\(O((\log n)^{3/2}/\sqrt n)\) 可由取 \(d=d_n+C\log n\) 吸收：该位移使正态 CDF 改变量为
\(\Theta((\log n)^2/n)\)。又 \(\log I_n=O(\log n)\)。所以在第 17 节引用的相对中偏差定理正式成立时，
\[
 \log Q_n=n\mu-d_n+O\left(n^{3/8}(\log n)^{1/8}\right). \tag{65}
\]
