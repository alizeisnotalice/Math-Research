# 从锐强型范数导出的精确几何相位不等式

本稿先由已知一般强矩形 $L^p$ 范数推出几何极大函数的 $L^1$ 上界，再用显式正函数截断构造导出相位泛函不等式。第二步是这里的推导，并非引用论文已有的相位定理。维数 $n$ 固定，随后取截断尺度极限；所得结论对每个 $n$ 均为精确不等式。

## 1. 原始文献输入与几何极大函数

Grafakos–Montgomery-Smith 的 [作者公开 PDF](https://stephenmontgomerysmith.github.io/preprints/uncentered.pdf) 主定理及 strong maximal Corollary 给出
\[
\|M_S\|_{L^p(\mathbb R^n)\to L^p(\mathbb R^n)}=c_p^n,\qquad
(p-1)c_p^p-pc_p^{p-1}-1=0,\quad p>1. \tag{1}
\]
本次打开的六页 PDF 中，主定理在打印第 1 页，强矩形推论在打印第 4 页（零起始页号 3）；项目 `sources.md` 的 L01 是同一文献。这里只有 (1) 来自该文献。

令 $y_p=(p-1)(c_p-1)$。方程等价于
\[
y_p=1+\left(1+\frac{y_p}{p-1}\right)^{-(p-1)},\qquad 1<y_p<2.
\]
故 $y_p\to\mu$，其中 $\mu=1+e^{-\mu}$ 是唯一正根，并且
\[
p\log c_p\longrightarrow\mu. \tag{2}
\]
对非负 $f\in L^1(\mathbb R^n)$，定义
\[
Gf(x)=\sup_{R\ni x}\exp\left(\frac1{|R|}\int_R\log f\right),\qquad e^{-\infty}=0.
\]
这里 $R$ 是有界轴向矩形，可用内点包含 $x$ 的约定，使各上水平集为开集。$\log^+f$ 局部可积，故每个平均有明确的扩展实值。对任意 $p>1$，Jensen 不等式给
\[
Gf\le [M_S(f^{1/p})]^p,\qquad
\int Gf\le c_p^{np}\int f.
\]
令 $p\to\infty$，由 (2) 得
\[
\boxed{\displaystyle \int_{\mathbb R^n}Gf\le e^{\mu n}\int_{\mathbb R^n}f.} \tag{3}
\]
对含零点的 $f$，Jensen 可先将 $\log f$ 从下截断，再单调取极限，因而无需假定 $f$ 处处正。

## 2. 精确相位泛函的陈述

令 $F:\mathbb R/(2\pi\mathbb Z)\to(0,\infty)$ 连续，$\xi\in\mathbb R^n$ 任意。令 $c(\theta)=(c_1(\theta),\ldots,c_n(\theta))$ 是周期 Borel 可测场，值域包含于某个紧集；连续性不是必需的。给定 $c\in\mathbb R^n$，令各坐标独立且
\[
Z_{c,i}=W_i+c_iJ_i,\quad W_i\sim\operatorname{Exp}(1),\quad
J_i\sim\operatorname{Bern}(q(c_i)),\quad q(t)=\frac1{1+e^t}.
\]
写 $M(t)=EZ_t=1+t/(1+e^t)$、
\[
(K_c H)(\theta)=E[H(\theta+\xi\cdot Z_c)],\qquad
\langle H\rangle=\frac1{2\pi}\int_0^{2\pi}H(\theta)\,d\theta.
\]
则有
\[
\boxed{\displaystyle
\left\langle\exp\left\{\sum_{i=1}^nM(c_i(\theta))+
(K_{c(\theta)}\log F)(\theta)\right\}\right\rangle
\le e^{\mu n}\langle F\rangle.} \tag{4}
\]
期望使用该接收相位所指定的整个动作向量。定理不要求 $\sum_i\xi_i=0$。当 $\xi=0$ 时，$K_c\log F=\log F$，且一元微分给 $M(t)\le M(\mu)=\mu$，所以 (4) 逐点成立。以下证明 $\xi\ne0$。

## 3. 正的原点截断与源质量

对 $x\in[-1,1]^n$ 且各坐标非零，置 $u_i=-\log|x_i|\ge0$，定义
\[
g_T(x)=\exp\left(\sum_i\min(u_i,T)\right)F(\xi\cdot u), \tag{5}
\]
立方体外置零，立方体内的坐标超平面上任取正的有界值，例如 $e^{nT}\min F$。这些零测集不影响积分。于是 $g_T$ 在整个立方体内几乎处处满足 $\min F\le g_T\le e^{nT}\max F$；没有在原点附近挖去零区域。特别是所有位于立方体内的见证矩形上，$\log g_T$ 有界且可积。

取绝对值坐标及其 $2^n$ 个符号分支，得
\[
\int g_T=2^n\int_{[0,\infty)^n}e^{-\sum_i(u_i-T)_+}F(\xi\cdot u)\,du.
\]
立方体 $[0,T]^n$ 外的贡献至多
\[
2^n\|F\|_\infty[(T+1)^n-T^n]=O_{n,F}(T^{n-1}). \tag{6}
\]
对任意有界可测周期函数 $H$，只要 $\xi_j\ne0$，沿第 $j$ 个坐标作完整周期划分就给
\[
\frac1{L^n}\int_{[a,a+L]^n}H(\xi\cdot u)\,du
=\langle H\rangle+O\left(\frac{\|H\|_\infty}{|\xi_j|L}\right), \tag{7}
\]
误差对位置 $a$ 一致；其它坐标只造成相位平移。余下不足两个周期的区间长度有界，因此 (7) 不需 $H$ 连续。用 $H=F$ 与 (6) 得
\[
\frac{\int g_T}{2^nT^n}\longrightarrow\langle F\rangle. \tag{8}
\]

## 4. 实际矩形平均、统一截断损失与接收积分

取 $C_0$ 使所有 $|c_i(\theta)|\le C_0$，并固定 $A>C_0+1$。当 $T>2A$、接收点满足 $u\in[A,T-A]^n$ 时，令 $\theta=\xi\cdot u$，第 $i$ 个区间的两个端点为 $x_i$ 与 $-e^{-c_i(\theta)}x_i$（按大小排序）。所得真实轴向盒 $R_x$ 完全包含于 $[-1,1]^n$，因为两端绝对值均不超过 $e^{C_0-A}<1$。

盒内独立均匀取样 $y$ 给出精确表示
$-\log|y_i|=u_i+Z_{c(\theta),i}$：同号分支概率 $1-q(c_i)$，反号分支概率 $q(c_i)$，两者上的对数增量分别为 $W_i$ 与 $c_i+W_i$。于是
\[
\frac1{|R_x|}\int_{R_x}\log g_T
=\sum_i u_i+\sum_iM(c_i(\theta))+
(K_{c(\theta)}\log F)(\theta)-\mathcal L_T(u), \tag{9}
\]
\[
\mathcal L_T(u)=\sum_iE[(Z_{c_i(\theta)}-(T-u_i))_+].
\]
这里相位仍用未截断的真实 $-\log|y_i|$，因此核平均项没有近似误差。若 $a\ge\max(0,c)$，指数分布的直接积分给
\[
E[(Z_c-a)_+]=(1-q(c))e^{-a}+q(c)e^{-(a-c)}
=2(1-q(c))e^{-a}.
\]
两个积分及最后恒等式另经本机 Wolfram 核验，见 `cap_loss_check.wl`、`cap_loss_check_output.txt`。因 $T-u_i\ge A>C_0$，(9) 的总损失满足
\[
0\le\mathcal L_T(u)\le2n e^{-A}. \tag{10}
\]
接收点可能位于 $R_x$ 的端点；将盒各端稍微向外扩大仍留在立方体内，因 $\log g_T$ 在该立方体上有界，其平均趋于 (9)。故采用内点约定的 $Gg_T(x)$ 仍不小于 (9) 右端的指数。

置有界可测周期函数
\[
H(\theta)=\exp\left\{\sum_iM(c_i(\theta))+(K_{c(\theta)}\log F)(\theta)\right\}.
\]
乘以 $dx=e^{-\sum u_i}du$，(9)—(10) 中的 $\sum u_i$ 精确抵消，故
\[
\int Gg_T\ge2^ne^{-2ne^{-A}}\int_{[A,T-A]^n}H(\xi\cdot u)\,du.
\]
由 (7)，固定 $A$ 后令 $T\to\infty$，再结合 (3)、(8)，得到
\[
e^{-2ne^{-A}}\langle H\rangle\le e^{\mu n}\langle F\rangle.
\]
最后令 $A\to\infty$，即得 (4)。先固定 $n,c,F,\xi$，依次令 $T\to\infty$、$A\to\infty$；所有有界性与可测性已在各步说明。$F$ 为常数且全部 $c_i=\mu$ 时 (4) 取等。

## 5. 非有界但逐点有限动作场的延拓

式 (4) 也适用于任意逐点有限的实值 Borel 场。先将每个动作截断到 $[-C,C]$，用已证明的有界场结论；固定相位下，有限维的动作向量随 $C$ 增大最终完全相同。被积函数始终介于 $0$ 与 $e^{\mu n}\max F$ 之间，故支配收敛即可取极限。这个延拓在有界场证明之后进行，不把非有界动作直接代入固定 $A$ 的矩形构造。

还可以允许 $c_i=+\infty$，按 $Z_{+\infty}=W_i$、$M(+\infty)=1$ 定义。因为 $C\operatorname{Bern}(q(C))\to0$ 依概率，且 $\log F$ 有界连续，核期望与均值均收敛，再用同一支配函数。这里不引入 $c_i=-\infty$。

对全部有限实数 $c$、$a\ge0$，还有统一尾界 $E(Z_c-a)_+\le2e^{-a}$。当 $c\le a$ 时第4节恒等式直接给出；当 $c>a$ 时，置 $s=c-a>0$，则
$$E(Z_c-a)_+=(1-q(c))e^{-a}+q(c)(c-a+1)
\le e^{-a}+e^{-a}(s+1)e^{-s}\le2e^{-a}.$$
这条尾界并不取消在第4节中保证实际矩形位于立方体内所需的有界动作假设；非有界场仍由截断结论延拓。

## 6. 与真实指标中位数之间尚缺的连接

式 (4) 控制几何平均。取 $F(\theta)=e^{-\epsilon\cos\theta}$ 时，其指数中的相位项精确为 $-\epsilon E\cos(\theta+\xi\cdot Z_c)$，但真实指标见证的边界是 $S_c-\epsilon\cos(\theta+\xi\cdot Z_c)$ 的中位数。尚未证明在 $\epsilon=O(\log n)$、动作均值损失受控的相关全相位场上，该中位数一致不超过 $\sum_iM(c_i)-\epsilon E\cos(\theta+\xi\cdot Z_c)+O(1)$ 或足够小的 $o(\log n)$ 余项。因此本稿不据 (4) 宣称指标 halo 的非局部上界，也不宣称一般端点系数闭合。
