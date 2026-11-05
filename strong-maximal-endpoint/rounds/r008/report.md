# 第 8 轮：完整单调轮廓机制的三项展开

日期：2026-11-06（Asia/Shanghai）。状态：第5轮Q_n的严格精化及一般C_n的新下界；两名独立子代理审查通过。

## 1. 定理与适用范围

仍以第5轮Q_n表示任意非减轮廓、任意比例向量、单一双曲接收见证的最优固定维数深层极限值。
写 ell=log n、rho=sqrt(ell/n)、c0=sigma log(sigma sqrt(2pi))，则
\[
 \boxed{\log Q_n=n\mu-\sigma\sqrt{n\ell}
 +\sigma\{\log\ell+\log(\sigma\sqrt{2\pi})\}\sqrt{n/\ell}
 +o(\sqrt{n/\ell}).}\tag{1}
\]
c0=1.17791928933...。因此一般输入也有相应严格下界：log C_n≥右侧的某个带o项的序列。
这里的等号仅适用于Q_n，不适用于C_n；一般上界仍为nmu+mu/2+o(1)。
相对附件单层基准表达式，本轮机制新增sigma loglog n sqrt(n/logn)，并改进固定系数；没有用该单层族全局三项上界作为前提。

## 2. 下界：正态分位数和指数轮廓

固定r=B^-1，T=ΣZ_r，均值nmu、方差nmu。
令eta=ell^-1/4，并用标准正态上尾定义正数z_n：
\[
 \overline\Phi_N(z_n)=(1-\eta)\sigma\rho.
\]
定义q=nmu−sigma sqrt(n)z_n。普通Berry–Esseen给
\[
 P(T<q)\le(1-\eta)\sigma\rho+O(n^{-1/2})
 =(1-\eta+o(\eta))\sigma\rho.
\]
因为 eta rho=ell^(1/4)/sqrt(n)，其量级严格大于绝对逼近误差。
Mills反演给
\[
 z_n^2=\ell-2\log\ell-2\log(\sigma\sqrt{2\pi})+o(1),
\]
进而z_n/sqrt(ell)=1−O(logell/ell)=1−o(eta)。
第3轮Jensen估计仍给指数轮廓平均≥exp(sigma z_n/sqrt(n))−P(T<q)>1，对充分大n具有严格余量。
其模数因子I_n≤2n，固定n先做有界截断再取层深极限，得到log Q_n≥q−log(2n)。
展开sqrt(z_n²)并用(logell)²/ell→0、log(2n)=o(sqrt(n/ell))，即得(1)的下界方向。
本段只在尾概率约rho的窗口使用绝对Berry–Esseen；不误用到附件单层的更小临界尾窗口。

## 3. 上界：把逆模数的斜率保留为1

固定常数C>c0，测试
\[
 s=n\mu-\sigma\sqrt{n\ell}+\sigma\log\ell\sqrt{n/\ell}+C\sqrt{n/\ell},
 \qquad d=n\mu-s.
\]
假设存在某个r向量及非减h，使J≤e^-s且E h(T)≥1。
第5轮点态界Phi_n(h(t))≤e^(t−s)保持。
将其中间区域上端从s+n/4改为s+b_n，b_n=n^(3/4)。对0≤v=t−s≤b_n，若h>1，则
\[
 y=\log h\le\frac{v}{n-1-v},\qquad
 h-1\le(1+O(n^{-1/4}))v/n.
\]
令M=E T≤nmu、V≤Cn，中间区域的期望超额因此≤
\[
 (1+o(1))E(T-s)_+/n\le(1+o(1))(d+\sqrt V)/n=(1+o(1))\sigma\rho.\tag{2}
\]

远尾依旧可统一忽略。高振幅h≥exp(4t/n)部分仍≤J×2^n(n/4)^(n−1)/(n−1)!=exp(−cn)。
低振幅部分用theta_n=c n^-1/4和第5轮一致二阶mgf界；远尾阈值s+b_n≥nmu+b_n/2给
\[
 E[e^{4T/n}1_{T>s+b_n}]\le\exp(-c'\sqrt n).
\]
具体指数是 −(c/2)sqrt(n)+C1 c²sqrt(n)+O(1)，固定c足够小即为负。
于T≤s−2log n有h≤n^-2，故E h≥1与(2)迫使
\[
 P(T\le s-2\log n)\le(1+o(1))\sigma\rho.\tag{3}
\]
所有o项一致于r向量与h。

## 4. 低尾的精确常数完成矛盾

写D=nmu−M。D>2d时Chebyshev使(3)左侧趋1，立即矛盾。
D≤2d时第5轮方差稳定性给
\(V=nmu+O(\sqrt{nD})\)。阈值低于M的标准化距离z最多为(d+2ell)/sqrt(V)，并且
\[
 z^2\le\frac{d^2}{\mu n}+o(1)
 =\ell-2\log\ell-2C/\sigma+o(1).
\]
这里平方标准化误差为O(ell sqrt(D/n))+O(ell^(3/2)/sqrt(n))=o(1)，不能只写相对o(1)后忽略指数中的常数。
如果阈值高于均值，左尾至少1/2−o(1)，更直接矛盾。
一致Berry–Esseen绝对误差O(n^-1/2)=o(rho)，结合正态Mills公式得
\[
 P(T\le s-2\log n)
 \ge\left(\frac{e^{C/\sigma}}{\sqrt{2\pi}}+o(1)\right)\rho.
\]
C>sigma log(sigma sqrt(2pi))时其系数严格大于sigma，与(3)矛盾。
所以对每个固定C>c0，充分大n时log Q_n≤s；下界方向已给相同常数，由此证明(1)。
这依赖第5轮已证明的均值缺口/方差稳定性、统一三阶矩和密度比界，不依赖任何待补的全空间halo比较。

## 5. 大维度本机校准

MATLAB R2026a将二项混合限制在众数两侧12个标准差，以众数递推相对权重，避免大阶乘相减。
正倾斜mgf使用log1p/expm1稳定计算；I_n用下降阶乘递推和递减几何尾上界。
五个维数从1万到1亿全部实际运行，总耗时约17秒；不使用Monte Carlo。
以下残差为
\((\log L_n-nmu+sigma\sqrt{nlogn})/\sqrt{n/logn}-sigma\log\log n\)，理论极限为c0。

|n|实际二项求和项数|数值第三常数残差|估计省略二项质量|
|--:|--:|--:|--:|
|1e+04|993|0.92999441|4.68e-32|
|1e+05|3135|1.05883161|5.89e-33|
|1e+06|9909|1.12944474|3.66e-33|
|1e+07|31329|1.16969468|3.53e-33|
|1e+08|99065|1.19417741|3.54e-33|

有限n的残差可从任一方向接近极限，不能以这五点拟合替代证明。
二项省略质量约1e-32，由binocdf估计；这些大维度值是双精度诊断，未作定向认证。
Wolfram 15在n=10000以60位对同一截断混合及完整I_n有限和独立核对：平均余量约1.0000043e-10，对数评分差2.03e-12，省略质量4.67658e-32。
第4轮的有限阶梯定向证书不自动认证本节大维度数值。

## 6. 一般上界线的定量障碍

Hagelstein–Parissis第2节Lemma 1、Lemma 2和Theorem 2给坐标复合及强矩形示性halo的显式界
\[
 T_n(\alpha)\le U_n(\alpha)^n,\qquad
 U_n(\alpha)=1+\frac{4}{\exp[-\log(1-\alpha)/n]-1}.
\]
将第5轮有界核心归约的H=exp(A rho)、tau=1−2/sqrt(n)、delta=sqrt(H/tau−1)代入，接收示性阈值alpha=1/(1+delta)。
其余量1−alpha约sqrt(A)(logn/n)^(1/4)，因此该显式根为
\(U_n=1+16n/\log n\,(1+o(1))\)，发散而非接近B。
要使这个特定显式根低于B，必须1−alpha<((B−1)/(B+3))^n，要求指数小余量。
所以现有固定维数Solyanik定理的这一显式版本不能解决当前核心；这不是说真实T_n必有大底数。
原始源码位置和实际MATLAB代入核验见sources.md。

## 7. 差距与下一轮

Q_n的展开现在闭合到sqrt(n/logn)尺度，并给一般C_n匹配该表达式的严格下界。
一般上界与该下界仍差sigma sqrt(nlogn)量级，首系数是否为一般最优仍未知。

下一轮两项动作：① 用加权测度下的重排检验非单调标量轮廓能否逃离Q_n障碍；② 转向能利用跨分量质量的新矩形与不等振幅，而不重复已排除的共同膨胀机制。
复现：MATLAB `rows=sm_large_soft;`；Wolfram `wolframscript -file code/large_soft_check.wl`。
本轮无已跟踪冗余文件需要清理，不作公开首创性声明。
