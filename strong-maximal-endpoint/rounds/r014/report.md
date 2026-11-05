# 第14轮：隐式正态根与所有固定阶逆对数项

日期：2026-11-06（Asia/Shanghai）。状态：严格一般渐近上、下界；本机求根只作独立复核。

## 1. 主定理

令L=log n、mu=log B、sigma=sqrt(mu)。用唯一正根定义
\[
 \overline\Phi(z_n)=\sigma z_n/\sqrt n,\qquad d_n=\sigma\sqrt n\,z_n.\tag{1}
\]
左端严格递减，由1/2降到0，右端严格递增，由0升到无穷，因此根存在唯一。则存在固定常数使充分大n时
\[
\boxed{n\mu-d_n-O(L)\le\log Q_n\le\log C_n\le\log\mathcal R_n
\le n\mu-d_n+O(n^{3/8}L^{1/8}).}\tag{2}
\]
Q_n是前轮已定义的全部单调双曲轮廓机制，mathcalR_n是一般可预测反馈上界。实际下界只需Q_n中的指数斜坡；一般上界依赖第11轮已证明的全局重排桥。
(2)的误差小于sqrt(n)/L^k，对每个固定k均成立。因此正态根的所有固定阶逆对数展开同样适用于一般log C_n。
这依然不是C_n的乘法渐近等价：对数余项仍可趋于无穷。

## 2. 停止策略的定量误差

第13轮不是只证明一个无速率o(1)。其FGL误差项逐项给
\[
 \kappa_n=n^{-1/8}L^{5/8},\quad
 P_\pi(S_n\le n\mu-d)
 \ge(1-O(\kappa_n))\overline\Phi(d/(\sigma\sqrt n))-n^{-8},\tag{3}
\]
一致于所有可预测策略和0≤d≤2sigma sqrt(nL)。
相对误差来自epsilon=32/(sigma sqrt n)、delta²=O(n^-1/4 L^1/4)，在x≤2sqrt L的FGL上、下界中最大误差O((1+x)delta)=O(kappa_n)。
这里用的是原文定理的相对尾结论，且已证明其条件二次变差的几乎处处界。不能以原策略的无条件方差替代它。

## 3. 用固定方差代理锐化中段超额

令V=1+4/e²。对每个动作r及所有lambda≥0，负倾斜的方差
\[
 v_r(\lambda)=(1+\lambda)^{-2}+c^2p_\lambda(1-p_\lambda),
 \quad c=-\log r,\quad p_\lambda=(1+e^{c(1+\lambda)})^{-1}
\]
不超过V（r=0按连续极限）。故对条件中心增量X_i=Z_i−m(r_i)，从0积分二阶导数得
\[
 \log E(e^{-\lambda X_i}\mid\mathcal F_{i-1})\le V\lambda^2/2.
\]
原策略中心鞅M=ΣX_i因此满足E e^-lambda M≤exp(nV lambda²/2)，而不是假设Σv(r_i)接近mu n。
Chernoff后对尾积分，d>0时
\[
 E(-d-M)_+\le\int_d^\infty e^{-u^2/(2nV)}du
 \le\frac{nV}{d}e^{-d^2/(2nV)}.\tag{4}
\]
因为E M=0，E(d+M)_+=d+E(-d-M)_+。又S_n=nmu−D+M、D≥0，所以在d~sigma sqrt(nL)窗口
\[
 E(S_n-(n\mu-d))_+\le d(1+o(\kappa_n)).\tag{5}
\]
证明量级：右侧(4)除以d为n^(-beta+o(1))/L，其中beta=mu/(2V)>1/8。最后一个严格不等式甚至由mu>1、V<2直接推出。
这一步替代第13轮的d+O(sqrt n)粗界；该粗界的相对L^-1/2误差不足以得到(2)。

## 4. 隐式根的一般上界

第13轮已知log mathcalR_n=nmu−sigma sqrt(nL)+o(sqrt(nL))，故只需在d~sigma sqrt(nL)窗口考察临界可行轮廓。
若J≤e^-s、s=nmu−d，则单调性给Phi_n(h(t))≤e^(t−s)。于0≤t−s≤b=n^(3/4)，有
(h−1)_+≤(1+O(n^-1/4))(t−s)/n；于t≤s−2L有h≤n^-2。
第13轮统一远尾估计原样适用：t>s+b处，h<e^(4t/n)部分≤e^-c sqrt n，其补集≤J 2^n(n/4)^(n−1)/(n−1)!≤e^-c'n。
严格可行性E h(S_n)>1及(5)从而强迫
\[
 P_\pi(S_n\le s-2L)\le\frac dn(1+O(\kappa_n)).\tag{6}
\]
这里n^-1/4、n^-2和远尾除以d/n均比kappa_n小。
将(3)应用于d+2L，再吸收n^-8，得到
\[
 \overline\Phi((d+2L)/(\sigma\sqrt n))\le\frac dn(1+O(\kappa_n)).\tag{7}
\]
Mills比在该窗口给平移2L的对数相对变化O(L^(3/2)/sqrt n)=o(kappa_n)。因此定义
\[
 H(d)=\frac{\overline\Phi(d/(\sigma\sqrt n))}{d/n}
\]
后，(7)意味着H(d)≤1+O(kappa_n)。H(d_n)=1，且在临界窗口
\[
 (\log H)'(d)=-\frac{\phi(d/(\sigma\sqrt n))}{\sigma\sqrt n\,\overline\Phi(d/(\sigma\sqrt n))}-\frac1d
 =-\Theta(\sqrt{L/n}).
\]
单调反演得d≥d_n−O(kappa_n sqrt(n/L))=d_n−O(n^(3/8)L^(1/8))。
例如可先测试s=nmu−d_n+K n^(3/8)L^(1/8)，K大于上述固定常数时即产生矛盾，从而不需要假设上确界有极值轮廓。这证明(2)右端。

## 5. 实际输入的严格指数斜坡下界

固定r=r*，T为n个独立Z_r之和。对d=d_n+C_0 L、q=nmu−d，使用h(u)=exp(u/n)1_{u≥0}。边界矩形平均满足
\[
 A_n(d)=E[e^{(T-q)/n}1_{T\ge q}]\ge e^{d/n}-P(T<q).\tag{8}
\]
同一FGL定理在iid情形的delta=0，给相对尾误差O(L^(3/2)/sqrt n)。根处尾概率为d_n/n~sigma sqrt(L/n)，所以绝对误差为O(L²/n)。
将d从d_n增大C_0 L，正态尾减少
\[
 \frac{C_0L}{\sigma\sqrt n}\phi(z_n)(1+o(1))
 =\frac{C_0L^2}{n}(1+o(1)),
\]
其中phi(z_n)=sigma z_n²/sqrt n(1+O(1/L))，z_n²~L。
固定C_0充分大就使P(T<q)<d_n/n<e^{d/n}−1，因此(8)严格大于1。
第3轮的模数因子
\[
 I_n=\int_0^\infty e^{-(1-1/n)u}(1+u/n)^{n-1}du\le2n
\]
仍适用。每个固定n先对指数斜坡有界截断，保留严格余量，再取物理层深极限；这与n→∞的顺序分开。
得到真实可积输入下界log Q_n≥q−log I_n≥nmu−d_n−O(L)，完成(2)。

## 6. 新的显式项与任意固定阶展开

令a=log(sigma sqrt(2pi))、H=log L+a，w=z_n²。正态尾积分反复分部给Mills展开，余项由首个省略项控制。因此
\[
 w+2\log w=L-2a+2\log(1-w^{-1}+3w^{-2}+O(w^{-3})).
\]
逐次代入得到
\[
 z_n=\sqrt L-H/\sqrt L+(2H-1-H^2/2)L^{-3/2}
 +O(H^3L^{-5/2}).\tag{9}
\]
配合(2)，于是一般端点常数具有更精细的显式结论
\[
\boxed{\begin{aligned}
 \log C_n={}&n\mu-\sigma\sqrt{nL}+\sigma H\sqrt{n/L}\\
 &+\sigma(H^2/2-2H+1)\sqrt n\,L^{-3/2}
 +O(\sqrt n\,H^3L^{-5/2}).
\end{aligned}}\tag{10}
\]
(2)的余项被(10)所写余项吸收。
更一般地，对任意固定阶，将Mills积分分部至足够多阶，取log并逐次匹配w的逆L系数即可得到以H为变量的多项式系数。w方程左边导数1+2/w+O(w^-2)有正下界，故代入残差控制实际根误差；开平方后误差仍有可控阶数。
(2)的误差比sqrt n/L^k对每个固定k均小，故每次固定截断所得展开可转给log C_n。这不是把固定阶渐近级数断言为收敛级数。

## 7. 本机复核与复现

MATLAB R2026a以log(erfcx(z/sqrt2)/2)−z²/2计算稳定对数尾，并做70步确定性二分求根。测试维数2,4,…,256及10^4,10^8,…,10^256；另测试log n=100,1000,10000，无需实际构造超大n数组。
Wolfram15用直接Erfc公式在40位、80位分别求根，对log n=100,1000,10000均稳定至30位；其Series独立给L^-3/2系数(−2+4H−H²)/2及下一系数(9−14H+8H²−H³)/2。

|log n|正根z_n（Wolfram80位的截短显示）|(9)截断误差z_n−近似|
|--:|--:|--:|
|100|9.42969333311276|3.7274739929e-5|
|1000|31.3708609726983|−1.5665480037e-6|
|10000|99.8974458448172|−1.8571685650e-8|

三个交叉点的MATLAB/Wolfram根差至多1.43e-14；MATLAB Code Analyzer无报告。
这些是正态根和代数展开的浮点复核，不是有限维C_n的区间证书；严格保证来自第2—6节的分析估计。算法确定性，无随机种子。原始FGL存在性常数尚未数值化。
复现：MATLAB `addpath('code'); sm_root_checks;`；Wolfram `wolframscript -file code/root_checks.wl`。最小输出见同目录JSON。
本轮来源沿用sources.md登记的FGL原文，Mills余项在此由积分分部直接得到；不主张公开优先权。未发现本项目已跟踪冗余材料需要清理。

## 8. 剩余差距与下一轮两项动作

上下界仍有O(n^(3/8)L^(1/8))的对数差距，远大于O(1)，因此乘法渐近等价未解决。有限维反馈与固定核之差也未被(2)精确描述。
下一轮：① 利用均值亏损为二阶、方差变化为一阶的局部结构，估计自适应控制能否在O(L)或更小尺度产生真实差别；② 用本机Bellman递推和固定核倾斜公式校准该尺度，先找可证明的机制障碍或明确新项，不凭有限维拟合宣称极值。
