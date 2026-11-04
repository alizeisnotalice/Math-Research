# 第十七轮独立审查：任意有限维范数球

## 结论与范围

候选可采用，但必须写成条件定理。令 \(K\) 是 \(\mathbb R^n\) 的范数单位球，
\(\rho(x)=\|x\|_K\)，且
\[
d\mu(x)=\rho(x)^{\tau-n}\,dx,\qquad0<\tau<n.
\]
额外假设中心 \(K\)-球极大算子对所有局部有限 Radon 测度具有统一弱型
\((1,1)\) 常数 \(B<\infty\)。则
\[
C^{\rm out}_{c,\mu,K}(1-t)-1\sim t^{\tau/n}.
\tag{17.1}
\]
并且在
\[
0<\tau<\min\left\{1,\frac{n+4-\sqrt{n^2+16}}2\right\}
\]
有
\[
C^{\rm out}_{c,\mu,K}(1-t)
=1+t^{\tau/n}+t^{2\tau/n}+o(t^{2\tau/n}).
\tag{17.2}
\]
不把 BIP 误作覆盖选择定理，也不声称任意 Besicovitch 基都有这个条件。\(B\)
只进入小半径余项，故不影响主系数。

## 径向重排与长见证

由 \(|sK|=|K|s^n\)，
\[
\mu(rK)=\kappa_{\tau,K}r^\tau,\qquad
\kappa_{\tau,K}=\frac{n|K|}{\tau}.
\tag{17.3}
\]
浴缸重排给任意可测 \(A\)
\[
\mu(A)\le\kappa_{\tau,K}
\left(\frac{|A|}{|K|}\right)^{\tau/n}.
\tag{17.4}
\]
若 \(Q=x+rK\)、\(r>C\rho(x)\)，Minkowski 不等式给
\[
(r-\rho(x))K\subset Q.
\tag{17.5}
\]
故 \(\mu(Q)\ge\kappa_{\tau,K}r^\tau(1-C^{-1})^\tau\)。又
\(|Q\setminus E|<t|K|r^n\) 蕴含
\(\mu(Q\setminus E)\le\kappa_{\tau,K}r^\tau t^{\tau/n}\)。

长类的全局包络仍须令 \(R\) 是所有实际长见证半径的上确界。于是
\(H_l\subset (R/C)K\)；取 \(r_j\uparrow R\)，因 \(C>2\)，充分大的见证球
\(x_j+r_jK\) 包含 \((R/C)K\)。这与 \(H_l\cap E=\varnothing\) 和 (17.4)
一起给出第十五轮的长类上界，未使用选择器或 BIP。

## 小半径类

在环带 \(2^k\le\rho(x)<2^{k+1}\)，短见证包含于
\[
D_k=L2^kK,\qquad L=2(1+C).
\]
由假定的弱型常数 \(B\) 及其测度扰动推论，
\[
|(H_s\setminus E)\cap\{2^k\le\rho<2^{k+1}\}|
\le B\frac{t}{1-t}|E\cap D_k|.
\tag{17.6}
\]
因为 \(\rho^{\tau-n}\le2^{k(\tau-n)}\)，且固定 \(y\ne0\) 时
\[
\sum_{k:y\in D_k}2^{k(\tau-n)}
\le\frac{L^{n-\tau}\rho(y)^{\tau-n}}{1-2^{\tau-n}},
\]
Tonelli 给
\[
\frac{\mu(H_s)}{\mu(E)}=O_{\tau,n,B}(tC^{n-\tau}).
\tag{17.7}
\]

## 任意范数的二阶薄壳界

令 \(Q_0=K\)、\(Q_u=u+K\)、\(\delta=\rho(u)<1/2\)、\(Q_s=sK\)。
层析恒等式为
\[
\mu(K)-\mu(u+K)
=(n-\tau)\int_0^\infty
\bigl(|K\cap sK|-|(u+K)\cap sK|\bigr)s^{\tau-n-1}\,ds.
\tag{17.8}
\]
括号非负：\(K,sK\) 同心嵌套，故
\(|K\cap sK|=\min\{|K|,|sK|\}\)，而任意平移交集不超过该最小值。

当 \(s\le1-\delta\)，\(sK\subset K\cap(u+K)\)；当 \(s\ge1+\delta\)，
\(K\cup(u+K)\subset sK\)。差只支撑在 \([1-\delta,1+\delta]\)。该壳内
\[
(1-\delta)K\subset(u+K)\cap sK,
\]
所以
\[
0\le |K\cap sK|-|(u+K)\cap sK|
\le |K|[1-(1-\delta)^n]=O_K(\delta).
\tag{17.9}
\]
壳宽为 \(2\delta\)，幂权在此有界，故
\[
\mu(u+K)\ge\mu(K)(1-K'_{\tau,K}\delta^2).
\tag{17.10}
\]
因此取 \(C=t^{-a}\) 时，二阶条件为
\[
2a>\tau/n,\qquad
a<\frac{1-2\tau/n}{n-\tau}.
\tag{17.11}
\]
它们可同时满足当且仅当
\(\tau<(n+4-\sqrt{n^2+16})/2\)。

## 环下界与传递

取 \(E=RK\setminus hK\)，其中
\[
h/R=\frac{t^{1/n}}{1+t^{1/n}}.
\]
对 \(\rho(x)<h\)，令 \(r=R-\rho(x)\)。在 \(h<R/3\) 时，Minkowski 不等式给
\[
x+rK\subset RK,\qquad hK\subset x+rK.
\]
缺失比例严格小于 \(t\)，故
\[
C^{\rm out}_{c,\mu,K}(1-t)-1
\ge\frac{(h/R)^\tau}{1-(h/R)^\tau}.
\tag{17.12}
\]
环和孔的边界是 Lebesgue 零集，权绝对连续，故 \(E\) 的内点与孔内点足以保证
严格光晕的这个下界。二阶时需额外 \(\tau<1\)，使
\(t^{\tau/n}t^{1/n}=o(t^{2\tau/n})\)。

不等正指数外层共同中心传递也延伸到中心 \(K\)-球：前向缺陷估计只用同一见证球
上的并集不等式；反向方向只用严格 Lebesgue 见证的有限稳定性、中心平移
\(L^1\) 连续性和有限 Boolean 原子切分。凸体 \(K\) 的边界 Lebesgue 零，因而
不依赖方体积结构。这一传递结论不需要上述弱型常数 \(B\)。

## 本机交叉核验

实际运行 Wolfram Language 15，以二维 \(\ell^1\) 单位球、
\(\tau=1/2\)、平移 \(u=(\delta,0)\) 直接对 (17.8) 积分。平移损失除以
\(\delta^2\) 在 \(\delta=1/10,1/20,1/40\) 时约为
\[
2.91486,\quad2.95380,\quad2.97596,
\]
与 (17.10) 一致。这只是独立样例核验；(17.8)--(17.10) 是统一证明。

## 审查意见

可保留为条件定理。必须保留“有限维范数球”及“该中心球族的统一弱型常数
\(B\)”；不能推广成含糊的 Besicovitch 性质或任意几何基。


---

# 第十七轮仓库整合稿终审

审查对象为仓库 paper/norm_geometry.tex。结论：**N.1--N.9 通过；N.10 建议补一条主稿反向构造接口说明。**

1. **N.7 的半径顺序已正确。** 在定义 \(z,A_C\) 并假设 \(A_C>z\) 后，先由
   \[
   \mu_{K,\tau}(E)\ge\kappa r^\tau(A_C-z)
   \]
   对每个实际长见证得出半径有界，才定义其上确界 \(R\)。随后用
   \(r_j\uparrow R\) 的实际见证覆盖共同核心。这避免了先假设最大球的循环；
   空长类也已处理。

2. **N.9 的常数与量词正确。** 壳内交叠差不超过 \(n|K|\delta\)，幂权上界为
   \(2^{n+1-\tau}\)，壳宽为 \(2\delta\)，除以
   \(\kappa=n|K|/\tau\) 正好得到
   \[
   1-A_C\le\tau(n-\tau)2^{n+2-\tau}C^{-2}.
   \]
   修正后的 Wolfram 脚本已实际复跑，窗口等价式使用正确的正二次条件
   \(\tau^2-(n+4)\tau+2n>0\)，三项断言均为 True。

3. **N.10 的非阻塞补强建议。** 第 76--81 行已准确引用主稿
   \(\ref{thm:weighted-transfer}\) 的正向/稳定 Boolean 构造；但为显式保住
   \(\cap_iE_i=F\)，建议增写一句：选有界 \(D\) 包含有限参考球及其中心扰动球，
   在 \(D\setminus F\) 的 Boolean 原子上分割，再把 \(D^c\cap F^c\) 补入各
   部分并置 \(E_i=\mathbb R^n\setminus A_i\)。同时说明以
   \(\mu_{K,\tau}\) 的 Radon 内正则性穷尽光晕。现有“其余延拓及公共交集步骤
   完全相同”可理解为该主稿步骤，故不是结论漏洞；补句会使范数球改造自足且避免
   无界全空间 Boolean 原子的歧义。

## 追加审查：N.1 可由范数几何自动获得

该改进成立，故可删除 N.1 作为额外假设，并改为范数球的引理结论
\(B=3^n\)（不声称尖锐）。

令一个有限相交 Besicovitch 闭球族的公共点平移为 \(0\)。若族不止一个球，
没有中心可以等于 \(0\)：否则它位于任一其他球中，违反 Besicovitch 条件。对
\(B(x_i,r_i)\) 置 \(a_i=\|x_i\|_K>0\)、\(u_i=x_i/a_i\)。公共点给
\(a_i\le r_i\)，而闭球的 Besicovitch 条件给
\[
\|x_i-x_j\|_K>\max(r_i,r_j)\ge\max(a_i,a_j).
\]
若 \(a_i\ge a_j\) 且 \(\|u_i-u_j\|_K\le1\)，则
\[
\|x_i-x_j\|_K
\le(a_i-a_j)+a_j\|u_i-u_j\|_K\le a_i,
\]
矛盾。因此 \(\|u_i-u_j\|_K>1\)。开集
\(u_i+\frac12\operatorname{int}K\) 两两不交，且均包含于
\(\frac32K\)，由体积比较得到族的基数至多 \(3^n\)。单球情形显然；无限族若
超过 \(3^n\)，取 \(3^n+1\) 个即与有限结论矛盾。

有限维范数空间可分，局部有限 Radon 测度自动满足 Aldaz 的
\(\tau\)-additivity条件。Aldaz 的 Definition 1.1、Theorem 1.2 表明 BIP 常数
正是所有此类测度的统一中心弱型 \((1,1)\) 常数。因此上述 packing 给出
\(B\le3^n\)，且只用于 BIP \(\Rightarrow\) 弱型；没有把 BIP 当作覆盖选择定理。

本机 Wolfram Language 15 实际核对
\(|(3/2)K|/|(1/2)K|=3^n\)（\(n=1,\ldots,12\)）通过；这是体积计算交叉检查，
不是替代上述 packing 证明。
