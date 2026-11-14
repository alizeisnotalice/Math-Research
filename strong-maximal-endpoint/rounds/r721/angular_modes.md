# 深层指标源的角向 Fourier 扰动与接收侧比响应

固定 $n\ge2$，在 $[-1,1]^n$ 内写 $u_i=-\log|x_i|$，$L=\sum u_i$。取实的非零切向频率 $\xi$，$\sum_i\xi_i=0$，以及真实指标源 $E_{T,\epsilon}=\{L>T+\epsilon\cos(\xi\cdot u)\}$。动作 $c_i\in(0,\infty)$，其穿接收点 $x$ 的区间为 $[-e^{-c_i}x_i,x_i]$（两端按大小排序）。盒内均匀取样得到源坐标 $u+Z$，$Z_i=W_i+c_iJ_i$，$W_i$ 是独立单位指数，$J_i$ 是独立成功率 $q_i=e^{-c_i}/(1+e^{-c_i})$ 的 Bernoulli，两个向量独立。实际角相位为 $\theta=\xi\cdot u\pmod{2\pi}$。设 $S=\sum Z_i$、$\varphi=\xi\cdot Z$，$r_\theta(\epsilon,c)$ 为 $S-\epsilon\cos(\theta+\varphi)$ 的中位数，$r=r_\theta(0,c)$。以下所有结论针对所选实际见证族，不假设一般极大算子的最优盒均有这种形状。

## 深层极限与二阶系数

允许 $c=c(\theta)$ 连续随接收相位改变且处于固定紧集。固定 $n,\xi$ 和充分小 $|\epsilon|$，这些动作证明的接收体积对源体积之比在 $T\to\infty$ 时趋于
\[
 \frac{\langle\exp r_\theta(\epsilon,c(\theta))\rangle}{\langle\exp[-\epsilon\cos\theta]\rangle}, \tag{1}
\]
其中 $\langle\cdot\rangle$ 对周期相位均匀平均。事实上盒内成功判据为 $T-L<r_\theta$。总深度 $l$ 的体积密度是 $2^ne^{-l}l^{n-1}/(n-1)!$，给定 $l$ 时 $u=lD$，$D$ 在标准单纯形均匀。因有两个 $\xi_i$ 不等，$\xi\cdot D$ 有绝对连续紧支撑分布（以该线性函数作一个单纯形内坐标）。Riemann–Lebesgue 使每个非零整数 $h$ 的 $\mathbb E e^{ih l\xi\cdot D}\to0$，故相位弱收敛至圆周均匀律。令 $v=T-l$，公共归一化后积分权为 $e^v(1-v/T)^{n-1}1_{v\le T}$。阈值 $r_\theta$ 是有界连续函数；除至多可数个具有正质量水平集的 $v$ 外，可对事件 $v<r_\theta$ 用弱收敛。在 $v$ 低侧以 $e^v(1+|v|)^{n-1}$ 支配，超过共同阈值上界时被积函数零；分母同理。就得到 (1)。这里先固定维数与频率、再令深度趋无穷，极限属于实际指定盒带来的下界。

写 $f(s)$ 为 $S$ 密度，$m(s)=\mathbb E[e^{i\varphi}\mid S=s]$，$m_2(s)=\mathbb E[e^{2i\varphi}\mid S=s]$。相位二阶展开在正规点进行：例如对 $n\ge2$ 假设 $r$ 不在 Bernoulli 子集和 $\sum_{i\in J}c_i$ 上即可。这个充分条件避免了在二维混合密度折点误用二阶可导性。条件于 Bernoulli 模式及指数向量的 Dirichlet 比例 $D$ 后，$S=d+W$，$\varphi=\varphi_0+W\xi\cdot D$，$W\sim\Gamma(n,1)$ 独立；$W\mapsto d+W-\epsilon\cos(\theta+\varphi_0+W\xi\cdot D)$ 在 $|\epsilon|\|\xi\|_\infty<1$ 时严格递增。中位数附近的逐条件隐式求导和紧单纯形下积分由均匀导数界给出，$f(r)>0$。由此 $r_\theta$ 在正规点连续、二次可导，以下余项对 $\theta$ 一致（正规参数的紧邻域亦可）。

记 $h(s,\theta)=\Re[e^{i\theta}m(s)]$、$j(s,\theta)=[1+\Re(e^{2i\theta}m_2(s))]/2$。移动门槛的尾概率展开为 $\bar F(s)-\epsilon f h-\epsilon^2(f j)'/2+o(\epsilon^2)$。可直接以固定模式、$D$ 的标量根 $W=(s-d)+\epsilon C+\epsilon^2 C\partial_s C+o(\epsilon^2)$（$C=\cos(\theta+\varphi_0+(s-d)\xi\cdot D)$）代入 Gamma 生存函数核实。设 $r_\theta=r+\epsilon a+\epsilon^2 b+o(\epsilon^2)$，解尾概率 $1/2$ 得
\[
 a=-h,\qquad b=-\frac{f'}{2f}a^2-a\frac{(fh)'}f-\frac{(fj)'}{2f}.
\]
平均后 $\langle a\rangle=0$、$\langle h^2\rangle=|m|^2/2$、$\langle j\rangle=1/2$。因此 (1) 的对数相对 $r$ 的二阶系数是
\[
 C_0(\xi)=\tfrac14\Big[(1+f'/f)(|m|^2-1)+(|m|^2)'\Big]_{s=r}
 =-\frac{(\partial_s+1)[f(1-|m|^2)](r)}{4 f(r)} . \tag{2}
\]

## 固定动作与一个允许适应的坐标：二阶增益均非正

单轴密度为 $k_i(z)=(1-q_i)e^{-z}(1_{z\ge0}+1_{z\ge c_i})$。消去任一个坐标 $i$ 后，约束和为 $s$ 的切片测度乘 $e^s$，在其余坐标上的表示就是 $e^{\sum_{j\ne i}z_j}\prod_{j\ne i} k_j(z_j)dz_j$ 乘以
\[
 (1-q_i)(1_{\sum_{j\ne i}z_j\le s}+1_{\sum_{j\ne i}z_j\le s-c_i}).
\]
这是随 $s$ 非减的测度。切片上相位 $e^{i\xi\cdot Z}$ 除去共同因子 $e^{i\xi_i s}$ 后等于 $\exp(i\sum_{j\ne i}(\xi_j-\xi_i)z_j)$，不依赖 $s$。于是 $e^sf(s)(1-|m(s)|^2)$ 是对所有复数 $z$ 取下确界的积分 $\int|\exp(i\sum(\xi_j-\xi_i)z_j)-z|^2d\nu_s$，这里 $\nu_s$ 为上述非减测度；该量单调非减。故 (2) 在所有正规点给 $C_0(\xi)\le0$。

更精确，记 $f_{-i}$ 为删去第 $i$ 坐标后的和密度。由上式对 $s$ 求导，仅两个边界平面贡献，积分平方距离的最优常数处其梯度为零，故只对测度求导。恢复共同相位后得到
\[
 D_\xi:=(\partial_s+1)[f(1-|m|^2)](r)
 =(1-q_i) f_{-i}(r)\,\mathbb E_{0,i}|e^{i\varphi}-m|^2
 +q_i f_{-i}(r-c_i)\,\mathbb E_{c_i,i}|e^{i\varphi}-m|^2. \tag{3}
\]
$\mathbb E_{u,i}$ 表示设置 $Z_i=u$，并在其余坐标的和等于 $r-u$ 的条件律下取期望，无密度的层贡献零。第二系数是 $(1-q_i)e^{-c_i}=q_i$；这一步不能把面上的相位当作独立的另一次抽样。

若只让第 $i$ 个 $c_i$ 随 $\theta$ 微调，且 $c_i$ 是未扰动 $r(c)$ 的内点驻点，令 $H_{ii}=-\partial_{c_i}^2r(c)$。以单轴密度导数 $\partial_{c_i}k_i=q_i k_i-q_i\delta_{c_i}$（分布意义）对命中率求导，可得 $\partial_{c_i}\bar F(r)=q_i[1/2-P(S_{-i}>r-c_i)]$。驻点要求 $P(S_{-i}>r-c_i)=1/2$，从而 $h_i:=f_{-i}(r-c_i)>0$，$H_{ii}=q_i h_i/f(r)$。同样对条件相位的密度分子与分母求导得
\[
 \partial_{c_i}m(r)=\frac{q_i h_i}{f(r)}(m-m_i),\qquad m_i=\mathbb E_{c_i,i}e^{i\varphi}.
\]
这里总导数沿 $s=r(c)$，与定 $s$ 偏导相同，因为该坐标分位梯度为零。令 $c_i(\theta)=c_i+\epsilon v(\theta)$，对任意有界实 $v$，二阶额外项为 $\langle-v\Re(e^{i\theta}\partial_{c_i}m)-H_{ii}v^2/2\rangle$。逐相位配方所得最大值为 $|\partial_{c_i}m|^2/(4H_{ii})$。而由 (3) 的第二项与 Jensen 不等式，
\[
 f(r)|\partial_{c_i}m|^2/H_{ii}=q_i h_i |m-m_i|^2\le D_\xi .
\]
故允许任意**一个**实际侧比随接收相位调整，对数评分的总二阶响应仍然非正。这个结论对任意非零切向频率与满足上述正规性、驻点条件的基准成立，不需对动作曲率数值拟合。

## 多侧比同时响应的精确有限维判别

在所有动作的内点驻点且 $H=-\nabla_c^2r$ 正定时，置 $g=\nabla_c m(r(c))$；同时选择各坐标 $c(\theta)=c+\epsilon v(\theta)$ 的最优二阶系数恰为
\[
 C_{\rm all}=\frac14\left(g^* H^{-1}g-\frac{D_\xi}{f(r)}\right). \tag{4}
\]
这是逐相位对 $v$ 取极值的 Schur 修正，实际可由连续 $v(\theta)=-H^{-1}\Re(e^{i\theta}g)$ 实现。非对角元由同一尾导数给出
\[
 H_{ij}=\frac{q_iq_j}{f(r)}\left[\tfrac12-P(S_{-i,-j}>r-c_i-c_j)\right],\qquad i\ne j.
\]
仅用 (3) 时每个坐标分别消费同一个 $D_\xi$，不能作为控制全部 $g^*H^{-1}g$ 的证明。这里既未证明 (4) 对所有频率非正，也未从指定盒的局部展开推出一般 $H_n$ 或 $C_n$ 的对数系数。

## 大维数、固定有限频率支撑，仅支撑坐标允许响应

固定 $d$，对每个 $n>d+2$ 取一个使**共同参数**的 $n$ 步中位数最大的 $c_n$（在紧化参数 $[0,\infty]$ 中取），所有坐标初始均用 $c_n$。设 $\xi=(t_1,\ldots,t_d,0,\ldots,0)$、$\sum_{j=1}^d t_j=0$，仅前 $d$ 个侧比允许随接收相位调整。以下极限定理不声称未扰动基准已在异构动作全集中全局最优。

写 $M(c)=1+c/(1+e^c)$。动作增量有一致三阶绝对矩、方差一致正且有界；中位数的 Berry–Esseen 展开统一为 $nM(c)+O(1)$。在紧化端点 $c=\infty$ 的单步增量依分布趋至单位指数，$n$ 固定时和的中位数连续，故极大值可取。与 $c=\mu$ 比较并利用 $M$ 的唯一非退化极大得 $c_n=\mu+O(n^{-1/2})$，且大 $n$ 时为内点。写 $r_n$ 为这个共同动作的中位数，$r_n=nM(c_n)+O(1)$。对称性及共同参数的驻点给每个坐标的分位梯度均为零。此处即使在子集和接缝也可二阶求导：对 $n\ge4$，分模式的 Gamma 分布函数在零端点延零后为 $C^2$，对相位扰动的条件标量根即使过零，积分用该延拓，Gamma密度及一阶导数在零均为零。需要的删除一、二轴分布密度和尾亦连续。（前面的有限维判别也可在这些条件下使用。）

设 $f_l^{(n)}$ 为 $l$ 个共同动作的和密度。对固定 $j$、$l=n-j$，任意有界 $z$，有 $f_l^{(n)}(r_n-z)/f_n^{(n)}(r_n)\to1$，对紧 $z$ 区间一致。其局部中心极限证明可直接对特征函数反演：近零时一致三阶展开给标准化局部高斯密度误差 $O(n^{-1/2})$，远离零处有公共指数卷积的衰减界 $(1+u^2)^{-l/2}$，$r_n-lM(c_n)-z=O_j(1+|z|)$。用同一密度上界（或 Gamma 卷积峰值）得该密度比对所有 $z$ 一致有界。

于是对前 $d$ 个坐标在 $S=r_n$ 条件下的联合律，与其未条件律之总变差为 $o(1)$：其密度比是 $f_{n-d}^{(n)}(r_n-\sum_{j\le d}Z_j)/f_n^{(n)}(r_n)$，分成有界深度及一致指数尾。固定某个响应坐标于 $0$ 或 $c_n$ 后，删去该坐标作同一论证亦得其余固定数目坐标趋于未条件律。故 $m\to\prod_{j\le d}\phi(t_j)$ 及 $m_i\to e^{it_i\mu}\prod_{j\ne i}\phi(t_j)$。这里 $r_*=e^{-\mu}$、$p_*=r_*/\mu$，$\phi(t)=(1+r_*e^{it\mu})/[\mu(1-it)]$。收敛对任何给定紧频率集合一致，因上述条件化界为总变差界而 $c_n\to\mu$。将 (3) 的消去坐标选择在支撑外，其 $0$ 面、$c_n$ 面条件相位极限分布相同，上述密度比趋一；所以 $D_\xi/f(r_n)\to1-|\prod\phi(t_j)|^2$。

在前 $d$ 轴曲率矩阵中，每个对角元趋于 $p_*$。需注意，非对角元出现除以 $f_n^{(n)}(r_n)\asymp n^{-1/2}$，仅证明删二轴尾趋于 $1/2$ 并不足够。利用一轴驻点可得更精确的消去。令 $F_{n-2}$ 为删二轴的分布函数、$x=r_n-2c_n$，则 $\mathbb E F_{n-2}(x+c_n-Z)=1/2$，其中 $Z$ 为单步共同动作。局部密度比和全局密度峰值给
\[
 \frac{1/2-F_{n-2}(x)}{f_n^{(n)}(r_n)}
 =\mathbb E\int_0^{c_n-Z}\frac{f_{n-2}^{(n)}(x+y)}{f_n^{(n)}(r_n)}\,dy
 =c_n-M(c_n)+o(1)=o(1).
\]
最后等式对有界 $Z$ 用一致局部极限，截尾用 $Z$ 的一致可积性。因而固定 $d$ 的非对角块元均趋零，$H_d\to p_*I_d$ 且最终正定。驻点处的 $g_i\to p_*[\phi(t_i)-e^{it_i\mu}]\prod_{j\ne i}\phi(t_j)$。这里只对 $d$ 维响应子空间作逆，不要求 $n$ 维 Hessian 的算子范数极限。

定义 $A(t)=|\phi(t)|^2$、$B(t)=p_*|\phi(t)-e^{it\mu}|^2$。支撑内各侧比同时响应后的二阶系数因此趋于
\[
 \frac14\left[\prod_j A(t_j)-1+\sum_i B(t_i)\prod_{j\ne i}A(t_j)\right]. \tag{5}
\]
用 $\mu=1+r_*$ 化简可得
\[
 1-A(t)-B(t)=\frac{\mu t^2+2r_*t\sin(\mu t)+2r_*^2(1-\cos(\mu t))/\mu}{\mu^2(1+t^2)}
 \ge\frac{1-2r_*}{\mu}\frac{t^2}{1+t^2}.
\]
这里 $r_*<1/2$。恒等式已另用本机 Mathematica 的精确 `FullSimplify` 以实未知量 $t$、$0<r_*<1/2$ 核验，见 [WL 源码](single_factor.wl)和[原调用输出](single_factor_output.txt)。因为 $A,B\ge0$，(5)方括号不超过 $\prod_j(A(t_j)+B(t_j))-1$。对每个固定非零切向模式，极限严格为负；对离开零向量的紧频率集，费用一致为负。该有限支撑结论已允许支撑内多个侧比同时调整，未覆盖支撑外 $n-d$ 个坐标同时随相位响应或频率支撑随维数增长的情形。
