# 非切向单调轮廓的精确归约与深度自适应见证

本稿使用已审 r728 的切向全场评分定理。归约保留真实指标源，允许动作依赖接收深度、相位乃至完整接收点。结论针对端点为 $x_i$ 与 $-e^{-c_i}x_i$ 的跨原点角点矩形族；它没有把这族矩形等同于一般强极大算子的全部成功矩形。区间比例导数、余弦轮廓的 Jacobian 抵消以及一般 Gamma 尾比另经本机 Wolfram 核验，见 `reduction_checks.wl`、`reduction_checks_output.txt`；其余论证在下文完整给出。

## 1. 陈述

令 $h\in C^1(\mathbb R/2\pi\mathbb Z)$ 为实函数，$\xi\in\mathbb R^n$，分解
\[
\lambda=\frac1n\sum_i\xi_i,\qquad a=\xi-\lambda\boldsymbol1,\qquad \sum_i a_i=0.
\]
假设
\[
1-\lambda h'(\phi)>0\qquad\text{对所有 }\phi. \tag{1}
\]
对 $x\in[-1,1]^n$，写 $u_i=-\log|x_i|$、$\ell=\sum_i u_i$，定义
\[
E_T=\{x\in[-1,1]^n:\ell>T+h(\xi\cdot u)\}.
\]
坐标超平面上的零测点任意定义。对任意非零坐标接收点 $x\in\mathbb R^n$ 与有限实动作 $c\in\mathbb R^n$，$R(x,c)$ 为以 $x_i,-e^{-c_i}x_i$ 作两端点的轴向矩形，端点按大小排序。令
\[
\Omega_T^{\rm cross}=\{x:\exists c\in\mathbb R^n,\ |E_T\cap R(x,c)|>|R(x,c)|/2\}.
\]
这个集合允许每点独立选择动作，故包含任何依赖 $T$、相位和深度的动作场所给出的接收集。

对每个固定 $B<\infty$，若 $\inf_b\|h-b\|_\infty\le B\log n$，则当 $a\ne0$ 时
\[
\limsup_{T\to\infty}\frac{|\Omega_T^{\rm cross}|}{|E_T|}
\le\exp\{\mu n+\varepsilon_n(B)\log n\},\qquad \varepsilon_n(B)\to0. \tag{2}
\]
误差统一于满足 (1) 的轮廓与频率；但每个维数先取 $T\to\infty$，不声称相位混合对联合 $(T,n)$ 有统一速率。当 $a=0$ 时有独立的径向估计
\[
\limsup_{T\to\infty}\frac{|\Omega_T^{\rm cross}|}{|E_T|}\le e^{\mu n+C}, \tag{3}
\]
其中 $C$ 为绝对常数。给 $h$ 减去常数可相应平移深度参数 $T$，故以下可设 $\|h\|_\infty\le A\le B\log n$。

## 2. 单调隐式轮廓与源的精确变换

定义圆周映射的提升 $Q(\phi)=\phi-\lambda h(\phi)$。由 (1)，其导数在紧圆周上有正下界，且 $Q(\phi+2\pi)=Q(\phi)+2\pi$；故 $Q$ 是严格递增的 $C^1$ 微分同胚提升。写其逆为 $\Theta$，并令
\[
\widetilde h(\psi)=h(\Theta(\psi)),\qquad
\widetilde h(\psi)=h(\psi+\lambda\widetilde h(\psi)). \tag{4}
\]
$\widetilde h$ 是 $C^1$ 周期函数，与 $h$ 具有相同值域，特别 $\|\widetilde h\|_\infty\le A$。对每个 $\psi$，函数 $t\mapsto t-h(\psi+\lambda t)$ 严格递增，唯一零点为 $\widetilde h(\psi)$，因此
\[
t>h(\psi+\lambda t)\quad\Longleftrightarrow\quad t>\widetilde h(\psi). \tag{5}
\]
置 $v=T-\ell$、$\psi=a\cdot u+\lambda T$，则 $\xi\cdot u=\psi-\lambda v$。由 (5)，原指标源精确等于
\[
E_T=\{x\in[-1,1]^n:v<-\widetilde h(a\cdot u+\lambda T)\}. \tag{6}
\]
圆周换元 $\psi=Q(\phi)$ 还给出
\[
\langle e^{-\widetilde h}\rangle
=\frac1{2\pi}\int_0^{2\pi}e^{-h(\phi)}(1-\lambda h'(\phi))\,d\phi
=\langle e^{-h}\rangle, \tag{7}
\]
因为 $\int e^{-h}h'=0$。相位中的 $\lambda T$ 不被丢弃；后文使用 Haar 的平移不变性处理它。条件 (1) 不要求 $\lambda$ 或振幅单独很小；余弦轮廓 $h(\phi)=A\cos\phi$ 的充分必要条件为 $|\lambda|A<1$。

## 3. 成功概率连续性、有理动作与可测上确界

本节假设 $a\ne0$，因此 $n\ge2$。给定动作 $c$，独立增量为 $Z_{c,i}=W_i+c_iJ_i$，其中 $W_i\sim\operatorname{Exp}(1)$、$P(J_i=1)=(1+e^{c_i})^{-1}$。令
\[
Y_{c,\psi}=\sum_iZ_{c,i}-\widetilde h(\psi+a\cdot Z_c),\qquad
p(c,\psi,v)=P(Y_{c,\psi}>v),
\]
并定义左中位数 $r_c(\psi)=\inf\{v:P(Y_{c,\psi}\le v)\ge1/2\}$。此时
\[
p(c,\psi,v)>1/2\quad\Longleftrightarrow\quad v<r_c(\psi). \tag{8}
\]
证明所需的连续性如下：$Z_c$ 有 $n$ 维 Lebesgue 密度，$\boldsymbol1$ 与 $a$ 线性无关，所以 $(\sum Z_{c,i},a\cdot Z_c)$ 有二维密度。对第二变量作任意固定值的一维平移，Fubini 给 $Y_{c,\psi}$ 无原子。再用共同指数变量与共同均匀变量实现 Bernoulli 阈值，若 $(c_k,\psi_k,v_k)\to(c,\psi,v)$，则 $Y_{c_k,\psi_k}\to Y_{c,\psi}$ 几乎处处；极限无原子使尾概率 $p$ 联合连续。因此 (8) 成立，且严格成功事件有开放的动作裕量。

特别，对整数 $C\ge1$，
\[
r_C(\psi):=\sup_{c\in[-C,C]^n}r_c(\psi)
=\sup_{c\in\mathbb Q^n\cap[-C,C]^n}r_c(\psi),\qquad
r_*(\psi)=\sup_Cr_C(\psi). \tag{9}
\]
等式无需预设中位数关于动作连续：若 $v<r_c(\psi)$，严格成功概率的连续性可在同一闭动作盒内找到有理动作仍严格成功。$r_c$ 的可测性也直接由 $r_c(\psi)=\sup\{v\in\mathbb Q:p(c,\psi,v)>1/2\}$ 给出。故 $r_C,r_*$ 都是周期 Borel 函数。

所有实动作的标量和满足统一 $\operatorname{med}(\sum Z_{c,i})\le\sum M(c_i)+C_0\le\mu n+C_0$，来自 r727 的统一中心矩与标量 Berry–Esseen（有限小维数可由方差界吸收）。于是
\[
r_c(\psi)\le\mu n+A+C_0. \tag{10}
\]
另一方面 $c=0$ 给所有相位一个有限的统一下界，因此 $r_C,r_*$ 均有限且上下有界。

枚举 $\mathbb Q^n\cap[-C,C]^n$。对 $\epsilon>0$，取第一个满足 $r_c(\psi)>r_C(\psi)-\epsilon$ 的动作，便得到有界 Borel 场 $c_{C,\epsilon}(\psi)$。应用 r728 的**全场一致**切向评分上界于频率 $a$、轮廓 $\widetilde h$，再令 $\epsilon\downarrow0$，得
\[
\frac{\langle e^{r_C}\rangle}{\langle e^{-\widetilde h}\rangle}
\le e^{\mu n+\varepsilon_n(B)\log n}.
\]
其误差不依赖动作盒 $C$；令 $C\to\infty$ 并用单调收敛，就得到
\[
\frac{\langle e^{r_*}\rangle}{\langle e^{-\widetilde h}\rangle}
\le e^{\mu n+\varepsilon_n(B)\log n}. \tag{11}
\]
这证明优化所有动作的包络仍受控，包含原动作 $c(\xi\cdot u)=c(\psi-\lambda v)$ 的深度依赖性。

## 4. 任意深度自适应动作的逐点包络与体积

对任意接收点（允许在单位立方体外），取 $u_i=-\log|x_i|\in\mathbb R$，盒内均匀样本的对数坐标为 $u+Z_c$。记 $S_c=\sum_iZ_{c,i}$。由 (5)，样本进入源的条件精确为
\[
u+Z_c\in[0,\infty)^n\quad\text{且}\quad
S_c-v>\widetilde h(\psi+a\cdot Z_c).
\]
去掉第一条件只会增大命中率，所以任何实际严格成功都推出
\[
p(c,\psi,v)>1/2,\qquad v<r_c(\psi)\le r_*(\psi). \tag{12}
\]
这对每个 $x,T,c$ 成立，不要求动作只依赖相位，也不要求给定动作场位于统一紧集。实际成功集本身可测：固定 $x$ 的矩形指标平均连续依赖有限动作，故严格成功可由可数有理动作捕获。

先在单位立方体内计算包络 $\{u\ge0:T-\sum u_i<r_*(a\cdot u+\lambda T)\}$。总深度 $l$ 的体积密度为 $2^ne^{-l}l^{n-1}/(n-1)!$；给定 $l$，$u=l\mathsf D$ 为均匀单纯形。因 $a\ne0$ 且切向，$a\cdot\mathsf D$ 有 $L^1$ 密度。其 $l$ 倍绕圈律以总变差趋 Haar：用有限区间阶梯密度作 $L^1$ 逼近，每个区间放大后绕圈的 TV 误差趋零，再用推前的 TV 收缩。加上任意平移 $\lambda T$ 不改变到 Haar 的 TV 距离。

令 $b_T=2^ne^{-T}T^{n-1}/(n-1)!$。用 $v=T-l$，对每个固定 $v$ 应用上述 TV 收敛；以 $e^v(1+|v|)^{n-1}$ 控制低侧，(10) 控制高侧。源亦由 (6) 同样计算，得到
\[
|E_T|/b_T\to\langle e^{-\widetilde h}\rangle,\qquad
|\text{单位立方体内的包络}|/b_T\to\langle e^{r_*}\rangle. \tag{13}
\]
因此即使动作随 $T$ 任意改变，点态上包 (12) 与固定包络的 (13) 仍成立。

立方体外接收点也可统一处理。若某坐标 $|x_i|=s>1$，其相反端点绝对值为 $t>0$，该区间落在 $[-1,1]$ 内的比例至多
\[
\frac{1+\min(t,1)}{s+t}\le\frac2{s+1}.
\]
源包含于单位立方体，严格命中率 $>1/2$ 因此强迫每个 $|x_i|<3$。结合 (10)、(12)，所有外部成功点满足
$u_i\ge-\log3$、至少一个 $u_i<0$、$\sum_i u_i>T-(\mu n+A+C_0)$。
固定那个负坐标并积分其有限区间，其余 $n-1$ 轴的指数尾给出体积 $O_{n,A}(e^{-T}T^{n-2})=o(b_T)$。当 $n=1$ 时此区域对大 $T$ 为空。源至少包含 $\{u\ge0:\sum u_i>T+A\}$，其体积不小于 $c_{n,A}e^{-T}T^{n-1}$，所以外部成功体积与源之比为 $O_{n,A}(T^{-1})$，对频率、轮廓及动作选择统一。由 (11)—(13) 即证 (2)，而且控制的是整个指定跨原点矩形族的成功集，不仅是某个预先固定动作场。

## 5. 纯径向分支：独立结算，不使用 Haar

当 $a=0$，$\xi=\lambda\boldsymbol1$。由 (1)，$l\mapsto l-h(\lambda l)$ 严格递增且趋于两端无穷，故存在唯一 $\tau_T$ 满足
\[
\tau_T=T+h(\lambda\tau_T),\qquad |\tau_T-T|\le A.
\]
源精确为 $E_T=\{x\in[-1,1]^n:\sum_i u_i>\tau_T\}$。令 $q_n^*=\sup_{c\in\mathbb R^n}\operatorname{med}(\sum_iZ_{c,i})$，由统一标量估计 $q_n^*\le\mu n+C_0$。任意实际成功动作满足
\[
P(S_c>\tau_T-\ell)>1/2\quad\Longrightarrow\quad \ell>\tau_T-q_n^*.
\]
在单位立方体内，深度尾为 Gamma$(n,1)$ 尾（乘 $2^n$）。对固定 $n$，因 $\tau_T\to\infty$，
\[
\frac{\int_{\tau_T-q_n^*}^{\infty}e^{-l}l^{n-1}\,dl}
{\int_{\tau_T}^{\infty}e^{-l}l^{n-1}\,dl}\longrightarrow e^{q_n^*}\le e^{\mu n+C_0}.
\]
外部成功点仍有 $|x_i|<3$，同第 4 节只有低一阶的深度尾体积。故得 (3)。这里没有相位均匀化步骤，亦没有把径向源的单一相位值替换为 Haar 平均。

## 6. 范围

在条件 (1) 下，本稿把一般非切向频率的这类指标源精确归约为切向轮廓，并通过可测动作包络控制任意深度自适应选择。若 $h=A\cos$、$A=O(\log n)$、$\lambda=O(n^{-1/2})$，(1) 最终自动成立；频率的切向部分可任意大或很小。若 (1) 失效，隐式阈值可能多支，本归约不适用。即使 (1) 成立，这仍是指定跨原点角点矩形族的非局部深层上界，没有覆盖一般输入与一般成功目录，故不确定一般端点的精确对数系数。
