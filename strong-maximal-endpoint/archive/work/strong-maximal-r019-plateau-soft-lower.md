# r019：平台加软尾构造消除 \(\sqrt n\) 分母代价

日期：2026-11-06。本文给出固定 \(c=\mu\) 的一个明确非减轮廓下界；它是实际可积函数的构造，不是仅限 KKT 的启发。

令 \(\ell=\log n\)，并定义
\[
 h_{\rm pl}(u)=
 \begin{cases}
 0,&u<0,\\
 1,&0\le u\le\ell,\\
 e^{(u-\ell)/n},&u>\ell.
 \end{cases} \tag{1}
\]

## 1. 精确 Orlicz 成本

令 \(I_n=\int_0^\infty e^{-(1-1/n)v}(1+v/n)^{n-1}dv\)，即 r018 的纯指数斜坡成本。由
\(\Phi_n(1)=1\)，(1) 的平移后成本恰为
\[
 J_{\rm pl}=\int_\mathbb R e^{-u}\Phi_n(h_{\rm pl}(u))du
 =1-e^{-\ell}+e^{-\ell}I_n. \tag{2}
\]
由于 \(e^{-\ell}=n^{-1}\) 且 \(I_n\sim\sqrt{\pi n/2}\)，
\[
 J_{\rm pl}=1+O(n^{-1/2}),\qquad\log J_{\rm pl}=o(1). \tag{3}
\]
这正是相对纯斜坡 \(I_n\asymp\sqrt n\) 的改进。

## 2. 与纯斜坡软根的严格平均比较

令 \(S_n=\sum Z_\mu^{(i)}\)，并定义纯斜坡和平台轮廓在边界的平均
\[
 \begin{aligned}
 A(q)&=\mathbb E[e^{(S_n-q)/n}\mathbf1_{S_n\ge q}],\\
 B(q)&=\Pr(S_n\ge q)+
 \mathbb E[(e^{(S_n-q-\ell)/n}-1)\mathbf1_{S_n\ge q+\ell}].
 \end{aligned} \tag{4}
\]
第二式正是 \(\mathbb Eh_{\rm pl}(S_n-q)\)。逐点比较给出
\[
 0\le A(q)-B(q)
 \le(e^{\ell/n}-1)A(q). \tag{5}
\]
令 \(q_{\rm soft}\) 是 r018 中满足 \(A(q_{\rm soft})=1\) 的唯一软根。于是
\[
 B(q_{\rm soft})\ge1-(e^{\ell/n}-1)
 =1-\frac\ell n(1+o(1)). \tag{6}

取任意固定 \(C>1\)。无需对 (4) 求导，在带 \([q_{\rm soft}-C,q_{\rm soft})\) 上直接使用 (1) 即有
\[
 B(q_{\rm soft}-C)-B(q_{\rm soft})
 \ge\Pr(q_{\rm soft}-C\le S_n<q_{\rm soft}). \tag{7}
\]
r018 的固定 \(\mu\) Fourier 密度展开给出
\[
 f_n(q_{\rm soft}-s)=\frac{\ell}{n}(1+o(1))
 \quad(0\le s\le C), \tag{8}
\]
故 (7) 的右端为 \(C\ell/n(1+o(1))\)。由 (6)，充分大 \(n\) 时
\[
 B(q_{\rm soft}-C)>1. \tag{9}
\]
这提供了严格的矩形平均余量。取固定维数后双曲截断深度趋于无穷，
\(E_{L_0-(q_{\rm soft}-C)}\) 的内部便是严格超水平见证。

## 3. 得分与 \(\log n\) 系数

由 (2)--(3)、(9)，平台构造的认证得分满足
\[
 \log L_{\rm pl}\ge q_{\rm soft}-C-o(1).
\]
使用 r018 的
\[
 q_{\rm soft}=n\mu-d_n+\frac{3-\mu}{6}\log n+o(\log n),
\]
得到固定 \(c=\mu\) 的严格下界
\[
 \log L_{\rm pl}
 \ge n\mu-d_n+\frac{3-\mu}{6}\log n+o(\log n). \tag{10}
\]

若 r017 的实际斜坡参数移动定理在 \(c=\mu+O(\sqrt{\log n/n})\) 的同一平台族上成立，则其
\((\mu-1)\log n/8\) 改善可叠加。相应的候选系数为
\[
 \frac{3-\mu}{6}+\frac{\mu-1}{8}=\frac{9-\mu}{24}. \tag{11}
\]
式 (10) 是本笔记已核验的固定 \(\mu\) 结论；(11) 仍依赖 r017 的独立定理，不能在此单独宣称。

## 4. 窄平台的乘法改进

平台宽度不必取 \(\ell\)。取
\[
 w_n=\frac12\log n+\log\log n. \tag{12}
\]
把 (1) 中的 \(\ell\) 替换为 \(w_n\)，其成本仍精确为
\[
 J(w_n)=1-e^{-w_n}+e^{-w_n}I_n
 =1+O((\log n)^{-1}). \tag{13}
\]
在软根 \(q_{\rm soft}\) 处，纯斜坡与平台的平均差满足更精确的估计
\[
 A(q_{\rm soft})-B_{w_n}(q_{\rm soft})
 =\frac{w_n}{n}(1+o(1)). \tag{14}
\]
理由是 \(S\in[q,q+w_n]\) 的贡献仅为 \(O((\log n)^3/n^2)\)，而
\(S\ge q+w_n\) 的尾部承载了 \(A(q)=1-o(1)\)，其损失因子为
\(1-e^{-w_n/n}=w_n/n+o(n^{-1})\)。

平台平均对阈值的左导数包含在 \(u=0\) 的跳跃密度及软尾导数：
\[
 -B_{w_n}'(q)=f_n(q)+\frac1n
 \mathbb E[e^{(S-q-w_n)/n}\mathbf1_{S\ge q+w_n}]
 =\frac{\log n+1+o(\log n)}{n}. \tag{15}
\]
这里使用 r018 的 \(f_n(q_{\rm soft})=(\log n/n)(1+o(1))\)。由 (14)--(15) 和单调性，平台根满足
\[
 q_{\rm pl}=q_{\rm soft}-\frac{w_n}{\log n+1+o(\log n)}+o(1)
 =q_{\rm soft}-\frac12+o(1). \tag{16}
\]
因此固定 \(\mu\) 的窄平台斜坡得分与纯斜坡得分之比为
\[
 \frac{e^{q_{\rm pl}}/J(w_n)}{e^{q_{\rm soft}}/I_n}
 \sim\sqrt{\frac{\pi n}{2e}}. \tag{17}
\]
这只是固定斜率、两段平台族的乘法常数改进；它不宣称平台宽度或一般 KKT 轮廓已经全局最优。
