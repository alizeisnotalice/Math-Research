# 第17轮：实际指数软轮廓下界的对数级改善

日期：2026-11-06（Asia/Shanghai）。状态：实际可积输入下界族的严格改善；不声称一般C_n的对数级系数已闭合。

## 1. 结论、实际分母与严格接收

沿用Z_c=Exp(1)+c Bern((1+e^c)^-1)，S_{n,c}为n个独立和，mu=log B、sigma=sqrt(mu)、L=log n、rho=sqrt(L/n)。定义
\[
 A_c(q)=E[e^{(S_{n,c}-q)/n}1_{S_{n,c}\ge q}].
\]
对n>1，A_c连续严格递减，A_c(0)>1、A_c(q)→0，故存在唯一q_n(c)>0使A_c(q_n(c))=1。
考虑非减指数轮廓h_q(t)=exp((t−q)/n)1_{t≥q}。其真实Orlicz模数为
\[
 J(h_q)=e^{-q}I_n,\quad
 I_n=\int_0^\infty e^{-(1-1/n)u}(1+u/n)^{n-1}du
 =\frac{n e^{n-1}\Gamma(n,n-1)}{(n-1)^n}.\tag{1}
\]
I_n与c完全无关。记此固定c、仅优化阈值q的指数斜坡见证的上确界评分为
\[
 R_n^{\rm ramp}(c)=e^{q_n(c)}/I_n.
\]
这是supremal witness score，不是本项目自适应上界mathcalR_n。
实际输入合法性：取任意q<q_n(c)，就有严格A_c(q)>1；先固定n对h_q有界截断保留余量，再使用第2—3轮已证固定维数深层极限，得到真正紧支撑可积输入下界。让截断高度及物理层深分别趋极限，再令q上升至q_n(c)，得到
\[
 C_n\ge Q_n\ge R_n^{\rm ramp}(c).
\]
不需要把A=1边界本身当作严格halo见证。

本轮证明，对任意固定T一致于|t|≤T，
\[
\boxed{\frac{q_n(\mu+t\rho)-q_n(\mu)}L
 \longrightarrow-\frac{p_*t^2}{2}-\frac{p_*\sigma t}{2},
 \quad p_*=(\mu-1)/\mu.}\tag{2}
\]
因此对明确指定c_n=mu−sigma rho/2，
\[
\boxed{\log R_n^{\rm ramp}(c_n)-\log R_n^{\rm ramp}(\mu)
 =\frac{\mu-1}{8}\log n+o(\log n).}\tag{3}
\]
这比第16轮单独的分位数比较更强：分母和严格接收均已结算。但只改善这个实际输入族的基准，未确定Q_n、反馈上界或C_n的同阶最优项。

## 2. 软边界位置与临界普通尾

所有本节估计一致于c=mu+t rho、|t|≤T。记m=m(c)、v=v(c)，则m=mu+O(rho²)、v=mu+O(rho)。令d=nm−q、eta=1/n，单步正mgf展开给
\[
 U_c(q):=E e^{(S-q)/n}
 =\exp(d/n+O(n^{-1}))=1+d/n+O(L/n)
\]
于d=O(sqrt(nL))。又A=U−G，其中G=E[e^{(S-q)/n}1_{S<q}]，F(q)=P(S<q)。
对d=(1+epsilon)sigma sqrt(nL)，FGL相对尾给F=o(rho)，所以A≥U−F>1。
对d=(1−epsilon)sigma sqrt(nL)，FGL给F/rho→∞；而
0≤F−G≤E(q−S)_+/n≤sqrt(nv)/n=O(n^-1/2)，因为q<nm且中心和均值零。因此A<1。
固定任意epsilon>0夹逼并用严格单调性，得到软根处
\[
 d=nm-q_n(c)=\sigma\sqrt{nL}(1+o(1)).\tag{4}
\]

进一步，第14轮的固定V=1+4/e²负中心mgf界给
\[
 E(q-S)_+\le(nV/d)e^{-d^2/(2nV)}=o(\sqrt{nL}).
\]
故根处F−G=o(rho)。而A=1意味着G=U−1=d/n+O(L/n)，从而
\[
 F(q_n(c))\sim d/n\sim\sigma\rho.\tag{5}
\]
FGL只用于(4)的粗位置，不承担微小尾差比较。

## 3. 尾密度比与实际边界导数

置N=n−1、a=q_n(c)−c，选唯一负倾斜lambda>0使N步倾斜均值为a。
由(4)与一步均值Taylor，lambda~rho/sigma。第16轮已经直接证明的统一Fourier引理给倾斜居中密度g_N：
\[
 g_N(0)\asymp N^{-1/2},\quad
 \|g_N'\|_\infty=O(N^{-1}),\quad
 g_N'(0)/g_N(0)=O(N^{-1}),\quad
 \|g_N''\|_\infty=O(N^{-3/2}).
\]
恢复原密度f=f_{N,c}后，有
\[
 f'(a)/f(a)=\lambda+O(N^{-1}),\quad
 |f''(s)|/f(a)\le C e^{\lambda(s-a)}(\lambda^2+N^{-1}).\tag{6}
\]
此外，精确负倾斜公式和全局Lipschitz界给
\[
\begin{aligned}
 F_{N,c}(a)&=M_c(\lambda)^N e^{\lambda a}
 \int_{-\infty}^0e^{\lambda y}g_N(y)dy,\\
 \left|\int_{-\infty}^0e^{\lambda y}g_N(y)dy-g_N(0)/\lambda\right|
 &\le\|g_N'\|_\infty/\lambda^2.
\end{aligned}
\]
相对误差O(1/(lambda sqrt N))=O(L^-1/2)，所以f(a)/F_{N,c}(a)~lambda。
第16轮同一CDF及密度Taylor给
F_{n,c}(q)/F_{N,c}(a)=1+O(rho²)、f_{n,c}(q)/f(a)=1+O(rho²)；第一式在这里先用刚证明的f(a)/F_N(a)~lambda再得到，不循环假设F_n的大小。
结合(5)，得到关键尺度
\[
 f(a)\sim f_{n,c}(q_n(c))\sim\lambda\sigma\rho\sim L/n.\tag{7}
\]

## 4. 精确软平均参数导数

单步分布对c的有符号导数为p_c(其自身分布−delta_c)，这是第16轮已核验的Laplace导数恒等式。对n步卷积及本轮可积目标微分给
\[
 \partial_cA_c(q)=np_c[A_c(q)-B_{N,c}(q-c)],\qquad
 B_{N,c}(x)=E[e^{(S_{N,c}-x)/n}1_{S_{N,c}\ge x}].\tag{8}
\]
分母仍是原来的n，而非N；这一点不可替换。
积分尾的直接微分给partial_q A=−A/n−f_{n,c}(q)。因此在软根上
\[
 q_n'(c)=np_c\frac{A_c(q)-B_{N,c}(a)}{1/n+f_{n,c}(q)}.\tag{9}
\]
导数与积分交换可从有限Gamma–二项混合和一致小正指数矩验证；原目标的正指数率1/n严格小于1。

## 5. 卷积Taylor、相同首系数与收益

下列所有B及导数均指B_{N,c}，x位于a附近。精确关系为
\[
 B'=-B/n-f,\quad B''=B/n^2+f/n-f',\quad
 B'''=-B/n^3-f/n^2+f'/n-f''.
\]
在a处B(a)≤E exp((S_N−a)/n)=O(1)，由(6)—(7)，
\[
 B'(a)=-f(a)(1+O(L^{-1})),\qquad
 B''(a)=-\lambda f(a)(1+o(1)).\tag{10}
\]
令X=Z_c−c。卷积A=E B(a−X)给
\[
 A-B(a)=(c-m)B'(a)+\tfrac12 E X^2 B''(a)+\mathcal E.
\]
由(6)的全局界和上面B'''公式，Taylor积分余项满足
\[
 |\mathcal E|\le C\{f(a)(\lambda^2+n^{-1})+n^{-3}\}.\tag{11}
\]
其证明保留指数权重：在连接a与a−X的线段上，f及其导数由e^{lambda|X|}乘相应局部尺度控制，B由O(e^{|X|/n})控制；E|X|³(e^{lambda|X|}+e^{|X|/n})一致有界。
由于c−m=(c−mu)+O(rho²)、E X²=mu+O(rho)，(10)—(11)推出
\[
 A-B(a)=-f(a)\{(c-\mu)+\mu\lambda/2\}+o(f(a)\rho).
\]
(7)使(9)的分母为f(a)(1+o(1))，因为1/n比f(a)小1/L。于是对c=mu+t rho一致有
\[
 q_n'(c)/(n\rho)=-p_*\{t+\sigma/2\}+o(1).
\]
沿有界t区间积分，n rho²=L，便得(2)，再取t=−sigma/2得(3)。

## 6. 本机复核与数值范围

MATLAB R2026a对软根做55步二分，正倾斜后Gamma–二项混合以众数两侧12sd权重递推。稳定单步log mgf使用
-log1p(−1/n)+log1p(p_c expm1(c/n))，避免大数相减。

|n|指定偏移软根收益/log n|
|--:|--:|
|256|−0.00999260989|
|1024|−0.00389253839|
|4096|0.000564022649|
|1e4|0.00287604154|
|1e5|0.00753210873|
|1e6|0.01094897063|
|1e7|0.01357600032|

极限(mu−1)/8=0.03480806785。小n负收益不与渐近定理冲突。
脚本还记录q根减1/n后的数值严格平均余量；这些双精度值不是定向证书，严格性来自A的严格单调性。
I_n以(1)的稳定Stirling形式计算；只截断logGamma至−1/(360n³)时，解析余项在0与1/(1260n^5)之间，但Gamma尾的浮点误差未认证，不能把该余项界冒充总误差证书。
Wolfram15用原始未倾斜二项权重和上不完全Gamma函数全量求和，在n=256、1024以60/100位独立求软根及I_n。它与MATLAB的正倾斜实现独立，仍共享同一分布模型。
复现：MATLAB `addpath('code'); sm_soft_shift_checks;`，再运行Wolfram `wolframscript -file code/soft_shift_checks.wl`。
Wolfram两档精度的独立软根吻合至25位；MATLAB Code Analyzer无报告。算法确定性，无随机种子；双精度二项截断与Gamma计算未作定向认证。

## 7. 差距与下一轮两项动作

现已证明同一实际软轮廓族可从比例偏离中获得正的对数级收益，超越“均值最大比例永远最优”的错误推断。
一般C_n的隐式根误差仍为第14轮的O(n^(3/8)L^(1/8))，不能由本轮家族比较直接缩小。
下一轮：① 用一项Edgeworth的可控绝对误差计算固定比例软根相对nmu−d_n的首个log n项，并结算I_n；② 比较该下界与可预测反馈模型在同尺度的损失，寻求一般上界改进或机制障碍。
FGL来源沿用sources.md；密度估计来自第16轮自包含证明。无新增优先权主张，无冗余已跟踪材料需清理。
