# 切向全频率的条件相位总变差

本稿始终要求频率切向 $\sum_i\xi_i=0$。不使用二阶公共频移不变性来改变非线性相位，不处理带任意径向分量的频率。未找到满足下面假设的反例；本稿给出一条可独立审查的全频率证明。关键是 Gamma–Dirichlet 的精确结构、卷积将弱收敛提升为总变差，以及条件切片的 $L^1$ Lipschitz **差分**界。正文定理及全场评分结算均已通过独立严格审查。Gamma 导数极值方程及协方差比值局部常数另经本机 Wolfram 核验，见 `kernel_checks.wl`、`kernel_checks_output.txt`；计算仅核验其中代数步骤，不替代以下概率证明。

## 1. 精确陈述与全实动作协方差

置 $q(c)=(1+e^c)^{-1}$、$M(c)=1+cq(c)$、$v(c)=1+c^2q(c)(1-q(c))$，$\mu=1+e^{-\mu}$。给定任意有限实动作 $c_1,\ldots,c_n$，独立取 $W_i\sim\operatorname{Exp}(1)$、$J_i\sim\operatorname{Bern}(q(c_i))$，并写
\[
Z_i=W_i+c_iJ_i,\quad m=\sum_iM(c_i),\quad S=\sum_iZ_i,\quad
D=\sum_i[\mu-M(c_i)].
\]
**条件投影定理。** 对每个 $D_n=o(n)$ 和固定 $R<\infty$，
\[
\sup_{\substack{D\le D_n,\ \sum a_i=0,\ \sum a_i^2=1\\ |s-m|\le R\sqrt n}}
\left\|\mathcal L\left(\sum_i a_i(Z_i-M(c_i))\mid S=s\right)
-\mathcal L\left(\sum_i a_i(Z_i-M(c_i))\right)\right\|_{\rm TV}
\longrightarrow0. \tag{1}
\]
条件律取后文联合密度给出的版本，中央窗口内其分母最终一致正。通过缩放、确定性平移和模 $2\pi$ 映射的总变差收缩，(1) 对任意非零切向 $\xi$ 的相位 $\theta+\xi\cdot Z\pmod{2\pi}$ 同样成立，且统一于频率大小及支撑。$\xi=0$ 的相位为常量，比较恒等。

首先对所有实 $c$ 有
\[
(v(c)-\mu)^2\le C[\mu-M(c)]. \tag{2}
\]
证明：比值在 $c\ne\mu$ 连续。在 $c=\mu$，$M''(\mu)=-p_*$、$v'(\mu)=\mu-1$、$p_*=(\mu-1)/\mu$，比值极限为 $2\mu(\mu-1)$；在 $+\infty$ 趋于 $\mu-1$，在 $-\infty$ 趋零。分母仅在 $\mu$ 为零，故比值全局有界。因此
\[
\sum_i(v(c_i)-\mu)^2\le CD,\quad
\left|\sum_i a_i(v(c_i)-1)\right|\le C\sqrt D,\quad
\frac1n\sum_i(v(c_i)-1)\to\mu-1 \quad(D=o(n)). \tag{3}
\]
这里第二式使用切向性。r727 的中心矩估计统一于全部实动作，包括稀疏的很负动作与任意大的正动作。

## 2. 一维对数凹密度的三个初等界

若 $k$ 为均值零、方差在 $[1/2,1]$ 的一维对数凹密度，则存在通用常数使
\[
\|k\|_\infty+|Dk|(\mathbb R)\le C,\qquad
\int |y|\,d|Dk|(y)\le C. \tag{4}
\]
$Dk$ 是 BV 分布导数，可包含支撑端点的跳跃。简要证明：令峰值为 $K$、众数为 $b$。总质量一与对数凹性给 $k(y)\le eK\exp[-K|y-b|/e]$（在众数两侧各到 $K/e$ 水平，距离至多 $e/K$，其后由对数凹性控制）。因此 $E|Y-b|^2\le C/K^2$，方差下界推出 $K\le C$；均值零又给 $|b|K\le C$。单峰性给 $|Dk|(\mathbb R)=2K$，并由 Stieltjes 分部积分给 $\int|y-b|\,d|Dk|(y)=1$，从而得到 (4)。支撑有硬边界时按 BV 测度理解即可。

令 $\mathcal D_t k(y)=t^{-1}k(y/t)$。由 (4) 对尺度参数积分，得到
\[
\|\mathcal D_tk-\mathcal D_uk\|_1\le C|\log(t/u)|,\qquad
\|\mathcal D_tk-k\|_1\le C|t-1|,\quad t,u>0. \tag{5}
\]
这是一条 $L^1$ 差分估计：可先平滑，利用尺度变化的分布导数总变差至多 $1+\int|y|d|Dk|$，再取 BV 近似极限；不声称原密度的尺度映射处处存在 $L^1$ 导数。第二式在 $|t-1|\le1/2$ 用第一式，其余用概率密度间 $L^1$ 距离至多二。

此外这族密度在 $L^1$ 中相对紧：方差上界给统一尾质量，BV 界给统一平移连续性。具体可在 $[-L,L]$ 上用固定小网格的区间平均逼近，误差由网格宽乘 BV 控制；区间平均向量处于有限维紧集，随后令网格变细、$L$ 变大。

## 3. 指数部分的全频率平滑核

对任意切向单位向量 $a$，令 $R_n=\sum_iW_i$、$U_W=(R_n-n)/\sqrt n$。精确 Gamma–Dirichlet 分解给 $W=R_n\mathsf D$，其中 $R_n\sim\Gamma(n,1)$ 与均匀单纯形向量 $\mathsf D$ 独立。置
\[
V_0=n\,a\cdot\mathsf D,\qquad V_W=\sum_i a_i(W_i-1)=(R_n/n)V_0.
\]
$V_0$ 均值零，方差为 $n/(n+1)$，其密度 $k_{n,a}$ 对数凹：单纯形切片体积的 Brunn–Minkowski 不等式给出此性质，$n=2$ 时直接为区间均匀密度。故 (4)—(5) 统一适用。条件于 $R_n$ 后用 (5)，
\[
\|\mathcal L(U_W,V_W)-\mathcal L(U_W)\otimes\mathcal L(V_0)\|_{
m TV}
\le C E|R_n/n-1|\le C/\sqrt n. \tag{6}
\]
标准化 Gamma 密度由 Stirling 公式逐点趋标准正态密度；两者质量均为一，Scheffé 引理给 $L^1$ 收敛。于是 (6) 的右侧再加一个仅依赖 $n$ 的 $o(1)$，即可将 $U_W$ 替换为与 $V_0$ 独立的标准正态 $G$。这是对全部切向单位向量一致的二维总变差近似。

还需要精确核的切片控制。写 $\gamma_n(r)=e^{-r}r^{n-1}/(n-1)!$，$r=n+\sqrt n x$。$(U_W,V_W)$ 的联合密度为
\[
K_{n,a}(x,y)=\sqrt n\gamma_n(r)\,\mathcal D_{r/n}k_{n,a}(y)\quad(r>0),
\]
在 $r\le0$ 时置零。若 $r_1,r_2>0$，由 (5) 得两个切片的 $L^1$ 距离不超过
\[
\sqrt n|\gamma_n(r_1)-\gamma_n(r_2)|
+C\sqrt n\min(\gamma_n(r_1),\gamma_n(r_2))|\log(r_1/r_2)|.
\]
Gamma 密度单峰，故区间内部的密度不小于两端最小值。将上式分别以 $\int|\gamma_n'|$ 与 $\int\gamma_n(r)/r\,dr$ 控制，得到
\[
\|K_{n,a}(x,\cdot)-K_{n,a}(x',\cdot)\|_1\le C|x-x'|. \tag{7}
\]
这里使用 $n\sup_r|\gamma_n'(r)|+n\sup_r\gamma_n(r)/r\le C$（$n\ge3$）。后一项由 $\gamma_n(r)/r=\gamma_{n-1}(r)/(n-1)$ 及 Stirling 得到；前一项在 $r=n-1\pm\sqrt{n-1}$ 达到候选极值，代入 Stirling 即得。跨过 $r=0$ 的差分同样成立，因为 $\gamma_n(0)=0$。这完全是 $L^1$ 差分论证，允许 $k_{n,a}$ 带硬边界。

## 4. Bernoulli 部分的联合弱极限

置 $B_i=c_i(J_i-q(c_i))$，$U_B=n^{-1/2}\sum_iB_i$、$V_B=\sum_i a_iB_i$。中心三阶矩统一有界；$V_B$ 方差统一有界，故其分布族紧。由 (3) 和一维 Berry–Esseen，$U_B\Rightarrow N(0,\mu-1)$。

证明 $U_B$ 与 $V_B$ 渐近独立时可合法使用 Bernoulli 的连续对数特征函数。事实上 $B_{-c}$ 与 $B_c$ 同分布；对 $c\ge0$，
\[
\beta_c(t)=e^{-icq(c)t}[1-q(c)+q(c)e^{ict}]
\]
在实轴从不为零：$c>0$ 时 $q(c)<1/2$，$c=0$ 时恒等于一。在每个固定紧实 $t$ 区间，$\ell_c=\log\beta_c$ 的前面所需三阶导数对**所有**实 $c$ 一致有界。证明可分为有界 $|c|$ 的连续紧性，以及 $|c|\to\infty$ 时各高阶项被 $|c|^ke^{-|c|}$ 控制。这里不对任意离散分布假定无零点，只使用这组明确的 Bernoulli 特征函数。

令 $b_i=\operatorname{Var}B_i=v(c_i)-1$。有 $\ell_i'(0)=0$、$\ell_i''(0)=-b_i$。对固定有界 $u,t$，Taylor 展开给
\[
\begin{aligned}
\mathcal R_n&=\sum_i[\ell_i(u/\sqrt n+t a_i)-\ell_i(u/\sqrt n)-\ell_i(t a_i)]\\
&=\frac u{\sqrt n}\sum_i\ell_i'(t a_i)
+\frac{u^2}{2n}\sum_i[\ell_i''(t a_i)-\ell_i''(0)]+O(n^{-1/2}).
\end{aligned}
\]
第一项中 $\sum_i\ell_i'(t a_i)=-t\sum_i b_i a_i+O_t(\sum a_i^2)=O_t(\sqrt D+1)$；第二项绝对值至多 $C_{u,t}n^{-1}\sum_i|a_i|=O(n^{-1/2})$。故
\[
\mathcal R_n=O_{u,t}(\sqrt{D/n}+n^{-1/2})\to0. \tag{8}
\]
指数化后这意味着联合特征函数与边际特征函数乘积之差趋零。因而沿任何 $V_B\Rightarrow\nu$ 的子列，
\[
(U_B,V_B)\Rightarrow N(0,\mu-1)\otimes\nu. \tag{9}
\]
这里允许极限 $\nu$ 非 Gaussian，未假设最大投影系数趋零。

## 5. 卷积升级为联合 TV，再升级为条件切片 TV

写 $U=(S-m)/\sqrt n=U_W+U_B$、$V=\sum_i a_i(Z_i-M(c_i))=V_W+V_B$。由 (6)，其联合律在 TV 中与 $(G+U_B,V_0+V_B)$ 相差 $o(1)$，统一于 $a$；这里 $G,V_0$ 与 Bernoulli 向量独立。

考虑任意满足 $D=o(n)$ 的参数序列。抽取子列，使 $V_B\Rightarrow\nu$，并由第 2 节紧性使 $k_{n,a}\to k$ 于 $L^1$。由 (9)，$(U_B,V_B)$ 弱收敛到乘积律。一个固定的 $L^1$ 概率密度与弱收敛概率测度卷积后，密度在 $L^1$ 收敛：先用连续有界概率密度逼近该核，再用弱收敛逐点、Scheffé 引理及卷积的 $L^1$ 收缩。将此用于核 $\phi_1\otimes k$，并先吸收 $k_{n,a}-k$，得
\[
J_n(x,y)\longrightarrow\phi_\mu(x)\,(k*\nu)(y)\quad\text{于 }L^1(\mathbb R^2), \tag{10}
\]
其中 $J_n$ 是实际 $(U,V)$ 的联合密度，$\phi_\mu$ 为方差 $\mu$ 的正态密度。记实际 $V$ 的边际密度为 $p_n$，由边际化的收缩性，亦有
\[
\int\|J_n(x,\cdot)-\phi_\mu(x)p_n\|_1\,dx\to0. \tag{11}
\]
这一步才是联合总变差；(9) 本身并不足够。

实际 $J_n$ 是精确核 $K_{n,a}$ 与 Bernoulli 联合律的卷积，故 (7) 给
$\|J_n(x,\cdot)-J_n(x',\cdot)\|_1\le C|x-x'|$。因此
$A_n(x)=\|J_n(x,\cdot)-\phi_\mu(x)p_n\|_1$ 是统一 Lipschitz 的非负函数。由 (11)，$\sup_xA_n(x)\to0$：若某点值为 $b$，Lipschitz 性会迫使 $\int A_n\ge cb^2$。于是实际 $U$ 的边际密度 $f_{U,n}(x)$ 一致趋 $\phi_\mu(x)$，且
\[
\sup_{|x|\le R}\left\|\frac{J_n(x,\cdot)}{f_{U,n}(x)}-p_n\right\|_1\to0. \tag{12}
\]
每个被考虑的参数序列均有上述子列论证；若 (1) 的上确界不趋零，选择一条违反收敛的序列，再作该子列抽取，即与 (12) 矛盾。这证明定理所声明的**一致性**，不是只对预先固定频率的逐点结论。

## 6. 相位均值与真实中位数的连接

这一节始终保留切向频率。令 $h$ 为实连续周期函数，中心化后 $\|h\|_\infty\le A_n$，其中 $A_n=O(\log n)$；动作满足 $D\le D_n=O(\log n)$。由 (1)，在 $|s-m|\le\sqrt n$ 内，条件相位律与未条件律的 TV 差至多某个统一的 $\epsilon_n\to0$。记 $H=h(\theta+\xi\cdot Z)$、$T=S-H$。因投影二维联合密度存在，$T$ 有连续分布；以下以左 $1/2$ 分位数定义其中位数即可。

(7) 对边际化给原尺度的密度差分界 $|f_S(s)-f_S(t)|\le C|s-t|/n$，并且 $f_S(m)\asymp n^{-1/2}$。全实动作统一中心三阶矩又给 $F_S(m)-1/2=O(n^{-1/2})$。对 $|x|+A_n\le\sqrt n$，直接对 $s=m+x+z$ 积分并只保留 $|z|\le A_n$，得
\[
F_T(m+x)=F_S(m+x)+f_S(m)EH
+O\left(\frac{A_n\epsilon_n}{\sqrt n}+\frac{A_n(|x|+A_n)}n\right). \tag{13}
\]
符号为正，因为 $T\le m+x$ 等价于 $H\ge S-m-x$，且
$\int[1_{H\ge z}-1_{z\le0}]dz=H$。把条件相位替换为未条件相位的误差由 TV 控制，再以密度差分界冻结 $f_S$，就是 (13)。

再展开 $F_S(m+x)=1/2+f_S(m)x+O(n^{-1/2}+x^2/n)$，在 $x=-EH$ 两侧取足够大的统一误差，得到
\[
\left|\operatorname{med}(S-H)-(m-EH)\right|
\le C\left(1+A_n\epsilon_n+\frac{A_n^2}{\sqrt n}\right)=o(\log n). \tag{14}
\]
当 $A_n$ 有界时右侧可保留为 $O(1)$。所有常数与相位、动作配置、频率大小及支撑无关；$o(\log n)$ 的一致性来自 (1)，不声称已有显式速率。

## 7. 任意全场的近远相位结算

固定 $B<\infty$。对每个维数，允许任意非零切向频率 $\xi\in\mathbb R^n$、任意实连续周期轮廓 $h$，以及任意有限实值 Borel 动作场 $c(\theta)\in\mathbb R^n$。给 $h$ 减去一个常数不改变评分，故假设 $\|h\|_\infty\le A\le B\log n$。令 $r_\theta$ 为 $S_{c(\theta)}-h(\theta+\xi\cdot Z_{c(\theta)})$ 的左中位数，则
\[
\boxed{\displaystyle
\log\frac{\langle e^{r_\theta}\rangle}{\langle e^{-h}\rangle}
\le\mu n+o_B(\log n),} \tag{15}
\]
其中余项统一于频率大小及支撑、轮廓和全部动作场。这里没有要求整个动作场的均值损失处处小；条件总变差定理 (1) 的 $D=o(n)$ 假设仍然保留。

证明：全实动作的中心三阶矩统一有界，且 $n\le\operatorname{Var}S\le2n$。标量 Berry–Esseen 因而给 $\operatorname{med}S\le m+C$，常数绝对。又 $S-h\le S+A$，故对每个相位，无论损失大小，均有
\[
r_\theta\le\mu n-D(\theta)+A+C. \tag{16}
\]
固定 $K>2B+1$，分成近场 $\mathcal N=\{\theta:D(\theta)\le K\log n\}$ 与其余远场。第 6 节对近场统一给出某个 $\delta_n\ge0$、$\delta_n\to0$，使
\[
r_\theta\le\sum_iM(c_i(\theta))-K_{c(\theta)}h(\theta)+\delta_n\log n
\quad(\theta\in\mathcal N).
\]
将正被积函数从近场扩展到整个圆周，再用 r725 的精确几何相位不等式，得到
\[
\frac{\langle e^{r_\theta}1_{\mathcal N}\rangle}{\langle e^{-h}\rangle}
\le e^{\mu n+\delta_n\log n}. \tag{17}
\]
远场则用 (16) 及 $\langle e^{-h}\rangle\ge e^{-A}$，
\[
\frac{\langle e^{r_\theta}1_{\mathcal N^c}\rangle}{\langle e^{-h}\rangle}
\le e^{\mu n}\exp[C+2A-K\log n]\le e^{\mu n}n^{-1}
\]
对充分大 $n$ 成立。相加并取对数即得 (15)。这里不把远场用于条件 TV 比较。

对于有界 Borel 场，r727 第 5 节的真实指标构造与深度极限仍适用：非零切向频率必不是常向量，所以单纯形投影具有 $L^1$ 密度，其放大绕圈律以总变差趋 Haar。$T=S-h(\theta+\xi\cdot Z)$ 的连续分布由本稿二维联合密度得到，不需保留非活跃坐标。以左中位数定义成功门槛，严格命中率 $>1/2$ 恰对应深度偏移小于 $r_\theta$。因此指定实际接收集与指标源的深层体积比趋于 (15) 左侧取指数的评分。对无界但有限实值场，双边截断每个动作后，在每个相位最终恢复全部原坐标；由 (16) 的统一上界 $e^{r_\theta}\le e^{\mu n+A+C}$，支配收敛给评分泛函的极限。r725 的几何不等式对同一截断亦由 $M(c)\le\mu$ 支配通过。故 (15) 统一覆盖这些有界见证族及其评分极限。

式 (15) 是**切向相位指标源与指定穿点侧比场的深层评分上界**。它没有把该接收集等同于完整 halo，也没有移除非线性问题中的切向限制。

## 8. 范围与下一步核验

定理的输入是 $D=o(n)$ 与严格切向的归一化投影。它涵盖近均值 $O(\log n)$ 窗口、全支撑及稀疏异常动作。证明没有用‘协方差趋零即独立’或‘弱收敛即 TV’，相应升级分别在 (8)—(9)、(10)—(12) 完成。需要独审的核心技术点是统一 BV 尺度差分、Gamma 核切片差分以及 Bernoulli 对数特征函数在紧实频率区间的全动作导数界；正文已给具体证明。一般非切向频率、任意实际成功目录和一般端点结论均不在本稿定理范围内。
