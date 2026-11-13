# 第707轮证明：独立标量机制的绝对 1/log n 展开

2026-11-14。记号沿用第706轮。Q_n 优化所有确定独立比例向量与所有非负非减标量轮廓；一般强极大算子常数记 C_n。设
$$ \mu=1+e^{-\mu},\quad \sigma=\sqrt\mu,\quad L=\log n,\quad
 \overline\Phi_G(z)=\sigma z/\sqrt n,\quad d_n=\sigma\sqrt n\,z, $$
$$ c_Q=(9-\mu)/24,\quad B_Q=(13\mu-15)/12,\quad
 A_Q=(\mu+3)/8-\log(2\pi\mu). $$
此文中的 z 始终是精确定义的 z_n，不是用有限项渐近式替代的近似值。

## 定理

完整独立确定比例标量机制满足
$$ \boxed{\log Q_n=\mu n-d_n+c_Qz_n^2+B_Q+
 \frac{-4\log z_n+A_Q}{z_n^2}+o(z_n^{-2}).} \tag{1} $$
等价地，在 c_Qz_n² 项仍保持精确隐式 z_n 的前提下，末项可写成
$$ \frac{-2\log L+A_Q}{L}+o(L^{-1}). \tag{2} $$
一般 C_n 获得(1)右侧的同精度下界，未得到相应匹配上界。

## 1. 全轴二阶 Edgeworth 的统一余项

写 S=ΣZ_i，Z_i=Exp(1)+c_i Bern((1+e^(c_i))^(−1))，无穷位移坐标定义为 Exp(1)。令 M=ES、V=Var S、K_j=cum_j(S)。每个坐标任意固定阶绝对矩全域一致有界；V 与 n 同阶（上下比值有统一正界），且指数部分给特征函数的平滑控制。第22轮第5节的 Fourier 证明向后保留一阶，得到
$$ F_S(M+\sqrt V x)=\Phi_G(x)-\phi(x)\left[
 \frac{K_3}{6V^{3/2}}H_2(x)+\frac{K_4}{24V^2}H_3(x)
 +\frac{K_3^2}{72V^3}H_5(x)\right]+O(n^{-3/2}), \tag{3} $$
全轴误差一致于所有独立比例向量。H_j 为概率论 Hermite 多项式。

说明统一余项：五阶绝对矩的一致有界性，在小频率区给 log 特征函数展开至四阶、余项 O(|t|⁵/n^(3/2))。指数化后保留三阶项、四阶项和三阶项的平方；三阶项纯虚，用其三阶 Taylor 余项控制。四阶实部及更高余项的指数增长由充分小固定频率窗中的负二次指数吸收。误差被 Cn^(−3/2)(|t|⁵+|t|⁷+|t|⁸+|t|⁹)exp(−ct²) 控制。窗外抽出指数因子的半数幂，除以 |t| 后的积分指数小。全线 CDF Fourier 反演即给(3)。该论证不要求同分布或单坐标位移有界。

## 2. 保留精确 Mills 比的软根

第22、705轮允许将 Q_n 及最佳临界斜坡的近极值限制到 D=nμ−M≤KL，其中 K 足够大且固定。此类有
$$ |V-\mu n|+|K_3-n\mu(3-\mu)|\le C_K\sqrt{nL}. \tag{4} $$
令 a_n=K_3/(6V)，对当前 V 定义
$$ \overline\Phi_G(z_V)=\sqrt V\,z_V/n,\quad d(V)=\sqrt V\,z_V,\quad
 R_V=\overline\Phi_G(z_V)/\phi(z_V),\quad H_V=(1+R_V/z_V)^{-1}. $$
设 f_G=φ(z_V)/√V。纯单位斜率软根是 A(q_s)=1，其中
$$ A(q)=E[e^{(S-q)/n}\mathbf1_{S\ge q}]. $$
考虑 q=M−d(V)+Δ，|Δ|≤C L。(3)中的两个二阶 CDF 项在此点均为 o(1/n)，因为 φ(z_V)=O(L/√n)，其多项式次数固定。Gaussian 位移的二阶误差与第一 Edgeworth 项的位置变化也都是 o(1/n)。因而
$$ F(q)=\frac{\sqrt V z_V}{n}+
 [\Delta+a_n(1-z_V^2)]f_G+o(1/n). \tag{5} $$

写 U(q)=E exp((S−q)/n)、G(q)=E[e^((S−q)/n)1_{S<q}]。累积量与指数展开保留绝对 1/n 阶：
$$ U(q)-1=\frac{d(V)}n-\frac\Delta n+
 \frac V{2n^2}(1+z_V^2)+o(1/n). \tag{6} $$
这里与第706轮相比，同时保留 V/(2n²)、−Δ/n，不能只保留其中的 log n/n 量级部分。

下侧超额的一阶矩必须计算而非只作 O(1) 估计：
$$ E(q-S)_+=\sqrt V\,\phi(z_V)[1-z_VR_V]+o(1). \tag{7} $$
证明：S≥0，q>0 且 q=O(n)，所以该期望是 ∫_0^q F(s)ds。逐项积分(3)，其统一余项积为 O(n^(−1/2))。Gaussian 和 Hermite 项向负无穷延伸的误差指数小；第一修正的原函数由 H_1φ 给出，在测试点为 O(z_Vφ(z_V))=o(1)，其余显式项更小。把终点 q 换回 M−√V z_V 的超额差至多 O(L)barΦ_G(z_V−O(L/√n))=O(L^(3/2)/√n)=o(1)。Gaussian 超额的显式值即为(7)。

第22轮倾斜密度还给 E(q−S)_+²≤C f(q)/λ³=O(√n/√L)=o(n)。因此指数加权低尾可以线性化到所需精度：
$$ F(q)-G(q)=E(q-S)_+/n+o(1/n). \tag{8} $$
把(5)—(8)代入 A=U−G，解掉 Δ 的系数 f_G+1/n。由于 1/(nf_G)=R_V/z_V，得到
$$ T(V,a_n)=H_V\left[a_n(z_V^2-1)+\frac Vn
 \left(1+\frac{(1-z_V^2)R_V}{2z_V}\right)\right]. \tag{9} $$
在 Δ=T±ε/L 处，平均减一的主项具有相反符号，量级为 ε/n；所有误差是 o(1/n)。单调性夹逼严格给
$$ q_s=M-d(V)+T(V,a_n)+o(L^{-1}), \tag{10} $$
一致于整个定位亏损类。

## 3. 软根与比例优化到 1/z² 阶

记 a_0=(3−μ)/6、b_0=(5μ−6)/6。Mills 比的分部积分给
$$ R(z)/z=z^{-2}-z^{-4}+3z^{-6}+O(z^{-8}),\qquad
 H=1-z^{-2}+2z^{-4}+O(z^{-6}). $$
令 V/n=μ、a_n=a_0 代入(9)，得到
$$ T=a_0z^2+b_0+\frac{3a_0+\mu/2}{z^2}+O(z^{-4})
 =a_0z^2+b_0+\frac{3}{2z^2}+O(z^{-4}). \tag{11} $$
(4)造成的替换误差 O(L^(3/2)/√n)=o(1/L)，所以(11)在(10)中可统一使用。

第706轮隐式导数的精确表达式为
$$ \eta_n=d'(\mu n)=z_nH_n/(2\sigma\sqrt n). $$
该轮二阶 Taylor 误差、每坐标全局分离优化的三阶误差都为 O(L^(3/2)/√n)=o(1/L)。因此对 D≤KL 取上确界的纯软根，比例优化增益仍为
$$ \frac{\mu-1}{8}z_n^2H_n^2+o(1/L). $$
所需展开现在是
$$ z_n^2H_n^2=z_n^2-2+5z_n^{-2}+O(z_n^{-4}). $$
共同 c_eta=μ−μη_n+O(η_n²) 实现该上确界到 o(1/L)；甚至共同 μ−μη_n 已足够，因为由二阶参数误差造成的评分损失 O(nη_n⁴)=O(L²/n)=o(1/L)。无须假设一般向量预先等比例。

遂有
$$ \sup_{D\le KL}q_s=n\mu-d_n+c_Qz_n^2+
 b_0-\frac{\mu-1}4+\frac{7+5\mu}{8z_n^2}+o(z_n^{-2}). \tag{12} $$

## 4. 最佳临界斜坡：保留平均导数中的 1

固定 A>√μ、C>(1+A²)/2，临界斜坡参数为 β=1+a√(L/n)、0≤a≤A、0≤w≤CL。第701、705轮已给唯一最佳宽度、最佳 a→√μ 和 q−q_s=O(1)，全部一致于亏损类。故优化时可把 a 限制到 √μ 的任意固定正邻域。

定义
$$ R=\overline\Phi_G(z_n)/\phi(z_n),\quad
 N=z_n/R=n\phi(z_n)/\sqrt{\mu n},\quad b=(\beta-1)\sqrt n. $$
需将旧估计 nf(q_s+u)/L→1 锐化为
$$ nf(q_s+u)=N+o(1) \tag{13} $$
对有界 u 一致。证明：在负倾斜中心作局部 CLT，第22轮的中心相对误差为 O(n^(−1/2))；勒让德指数与原方差的二次近似差、倾斜方差比及 O(L) 测试位移共同给
log[f(q_s+u)/f_G]=O_K(L^(3/2)/√n)。因此绝对误差 nf−N_V=O(L^(5/2)/√n)=o(1)。第706轮已有 z_V²−z_n²=O(√(L/n))，而 dN/d(z²)=O(1)，故 N_V=N+o(1)。

纯软平均的精确导数是
$$ nA'(q)=-A(q)-nf(q). $$
于是对 q−q_s 有界，n[A(q)−1]=−(1+N)(q−q_s)+o(1)。斜坡平均相对纯软平均的差商公式保留绝对 o(1)，得到
$$ (q-q_s)(1+N)=b\sigma z_n-w+o(1). \tag{14} $$
误差说明：乘 n 后的指数二阶差商余项由
(β−1)E[u_+²exp(2u_+/n)]/n+wE[u_+exp(2u_+/n)]/n 控制，为 O(L^(3/2)/√n)=o(1)。把截断 min(u_+,w) 换成 w 的误差至多 wP(S<q+w)=O(L^(3/2)/√n)。又有 E(S−q)_+=M−q+O(1)=σ√n z_n+O_K(L)，乘 β−1 后仍有 o(1) 精度。

最佳宽度的精确驻点关系是 w=log(I(β)−1)+log[nf(q)/(βA_s)]。此时 βA_s=1+o(1)，而最佳 a 附近 b→∞，Gamma 截尾概率趋一。Stirling 给
$$ \log I(\beta)=\tfrac12(L+b^2+\log(2\pi))+o(1), $$
故
$$ w=\tfrac12(L+b^2+\log(2\pi))+\log N+o(1),\qquad
 \log K=N^{-1}+o(L^{-1}). \tag{15} $$
代回(14)，对 b 完成平方 bσz_n−b²/2=μz_n²/2−(b−σz_n)²/2。上界直接由此得到；下界取窗口内可行的 b=σz_n 及其唯一最佳宽度。因此
$$ \log R_n^{crit}(\mathbf c)-q_s(\mathbf c)=
 \frac{\mu z_n^2/2-L/2-\log N-\log(2\pi)/2}{1+N}
 -\frac1N+o(L^{-1}). \tag{16} $$
利用 N=z_n²+1+O(z_n^(−2)) 与 L−z_n²=4log z_n+log(2πμ)+o(1)，得到一致展开
$$ \log R_n^{crit}(\mathbf c)-q_s(\mathbf c)
 =\frac{\mu-1}2+
 \frac{-4\log z_n-\log(2\pi\sqrt\mu)-\mu}{z_n^2}+o(z_n^{-2}). \tag{17} $$
分母 1+N 中的 1 对最终常数有贡献，不能删除。该段在 L 记号下贡献 −2log L/L，没有 (log L)²/L。

## 5. 完整轮廓与系数合并

第705轮已在所有独立比例向量上证明
$$ \log Q_n-\log Q_n^{crit}=D_0/L+o(1/L),\qquad
 D_0=(\mu-1-\log\mu)/2. $$
由定位，(12)、(17)可以分别优化比例并保留 o(1/L) 误差。因 1/L−1/z_n²=o(1/L)，把(12)、(17)及 D_0/z_n² 相加，常数项为 B_Q，1/z_n² 系数为
$$ \frac{7+5\mu}{8}-\mu-\log(2\pi\sqrt\mu)+D_0
 =\frac{\mu+3}{8}-\log(2\pi\mu)=A_Q. $$
这证明(1)。式(2)只转换最后一项，并未替换前面的 c_Qz_n²。

若把前面 z_n² 也展开，需要额外保留
$$ z_n^2=L-2\log L-\log(2\pi\mu)+
 \frac{4\log L+2\log(2\pi\mu)-2}{L}+o(1/L). $$
此时完整公式变成
$$ \begin{aligned}
 \log Q_n={}&\mu n-d_n+c_QL-2c_Q\log L+B_Q-c_Q\log(2\pi\mu)\\
 &-\frac{(\mu+3)\log L}{6L}
 +\frac{(5\mu-9)/24-(\mu+3)\log(2\pi\mu)/12}{L}+o(1/L).
 \end{aligned} \tag{18} $$
d_n 仍保持精确隐式值。一般 C_n 只获得同式下界；本轮没有进行一般上界闭合。
