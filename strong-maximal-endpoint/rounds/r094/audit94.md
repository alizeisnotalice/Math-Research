# 第94轮独立审查

**结论。** 逐式核对 core94 (1)–(15) 后，完整无条件对数凸下集族的上系数 $c_Q$ 通过审查。匹配下界沿用第93轮等权产品子水平集。单坐标 Steiner 对称化在高密度区可严格减少整个真实 halo；任意指标没有被归约到这个已解的族。下述各渐近均先固定 $B$ 和量词中出现的其他常数，再令整数维数趋于无穷，未给有限维认证起点。

## 1. 原始物理和受限面

- 斜带左侧接收点 $(y-\delta,y)$ 的指定边长 $1/4$ 正方形只缺面积 $\delta^2/2$ 的三角形。水平跨度不足以碰另一斜边，垂直范围位于 $(0,1)$。源内部加上该外片的体积至少 $1+\sqrt{2\eta}/8$。真单位矩形的全 halo 为 $1+4A\log A$，其中 $A=(1-\eta)^{-1}$。$0<\eta\le10^{-4}$ 时 $(\sqrt2-1)/8>4\sqrt\eta/(1-\eta)^2$，故 (1) 的方向和常数均正确。
- 乘后缀 $D$ 时两组坐标的极大函数精确相乘。对于 $b=M1_D>p$ 的每个后缀截面，均有阈值 $p/b\in[p,1)$，可应用整段严格差，得到 (2)。在 $b>p+L/N$ 子集保留 $\sqrt{1-p/b}\ge\sqrt{(L/N)/(p+L/N)}$。后缀维数 $N-2$ 的临界门槛与 $N$ 的门槛仅差 $o(L/N)$，隐式根导数和固定系数两侧夹逼合法；深度 $T$ 可依维数选定。全体积差下界为 $Q_N N^{-1/2+o(1)}$ 的对数意义，各源本身达到 $q_0+c_QL+o(L)$。没有据此声称二者相差已占正比例。
- 对数图像的凸下集具有有限的 Lipschitz 凹屋面。正体积给非退化原点盒，沿固定对角方向无上界会包含任意大原点盒，与总体积有限矛盾。切向导数给 $w=N^{-1}\mathbf1-\nabla_V u$，$\sum w_i=1$；向下减小各坐标的点仍在下集，迫使各 $w_i\ge0$。几乎处处同一射线的支持平面限制所有真正成功矩形采样点，侧边在单矩形内固定、可随接收点改变。变元 Jacobian 为常数乘 $e^{Ns}$，证明 (5) 的体积上界；可积 Lipschitz 权 $e^{Nu}$ 的弱切向导数可积，截断分部积分给所列均值。此处只对几乎处处方向取法向，无穷个面不需并算。
- 有界系数族中 $\sum\lambda_i=N$ 且 $\lambda_i\le A$ 给正比例数目 $\lambda_i\ge1/2$。加权密度和 CDF 展开可用 r093 的 Fourier 界；$V_2=\mu(N+H)+O_A(\sqrt{ND})$，竞争向量有 $D=O(L)$、$H=O(\sqrt{NL})$。r022 逐坐标最小化使用 $t=\ell\lambda_i$ 再乘 $\lambda_i$，因此负二次项含 $\sum\lambda_i^3=N+O_A(H)$，正项 $\mu\ell H/2$ 足以吸收该误差。有限稀疏面有固定正概率全部取小于 $\mu/2$ 的未平移指数增量，给所述分位上界。普通无条件凸体的支撑泛函和坐标绝对值混合满足 $EV_i\ge\sqrt2-1$、$\operatorname{Var}V_i\ge1/24$，四阶中心矩和插值给 $P(Y>0)\ge1/(4C)$，由此得 (9)。凸性、无条件性在该结论里同时使用。

## 2. 全权重尾比较

**优化 Laplace 曲率。** 写 $u=s\lambda$，单坐标对数 Laplace 变换 $F$ 如正文，极小动作满足 $e^{uc}=1+u(1+e^{-c})$。在此处 $p_{(1+u)c}=p_c/(1+u)$ 使 $F_{cc}=up$，交叉项和 $F_{uu}$ 代回 Schur 补得到 §6.1 的表达式。关键为 $p(1-y)^2\le (1/2)(v/2)\le u/2$，其中 $(1-y)^2\le\min(v/2,1)^2\le v/2$。因此 $\ell''\ge1/[2(1+u)^2]$，且 $\ell+(1/2)\log(1+u)$ 凸。原始的 Jensen 差 $P$ 和动作差 $A$ 都非负，这是**率**的比较；仅凭此还不能比较概率。

**任意权重下的倾斜前因子。** 倾斜后两独立部件确为 $a_i Y_i$ 和 $a_i d_i B_{d_i}$，其中 $a_i=\lambda_i/(1+u_i)$、$d_i=(1+u_i)c_i$，零权及无平移动作均合法。对纯指数和以 $v^2=\sum a_i^2$ 归一、对权重降序，任意序列有坐标收敛子序列。去掉头 $m$ 项的尾在每个固定 Fourier 频率上与同方差高斯对数特征函数相差 $O_t(m^{-1/2})$，因为最大剩余归一权不超过 $1/\sqrt m$，平方和不超过一。极限为可能退化的高斯与独立中心化指数级数之和。无高斯时 $\sum b_j\ge1$，任一固定开子区间紧含于 $(-1,\infty)$ 可选有限部分使其支撑左端严格低于目标区间，再用独立剩余尾的方差趋零；有限部分在其支撑内部有正密度。故任意极限对目标区间有正质量。分别固定一个区间在 $(-1,-1/2)$，另一个在 $(C,C+1)$；否则取质量趋零子序列，即与 Portmanteau 矛盾。纯指数和的密度对数凹，这两个质量下界给出 intervening $[-1/2,C]$ 的统一密度正下界。零权可直接丢掉。此论证没有假设完整 Bernoulli 混合对数凹。

Bernoulli 部件的方差至多 $C_0 v^2$。Cantelli 先给 $D-ED\le v/4$ 的固定正概率，再以 Chebyshev 丢掉 $D-ED<-C_1v$ 的小概率。因独立，在 $0\le h\le v/8$，$q-h-D$ 落入纯指数和已控制的密度区间，给完整倾斜密度在 $q-h$ 至少 $c/v$。积分 $e^{-sh}$ 得 (11)。假如 $q\ge ES$，同样的统一四阶矩与插值用于 $ES-S$ 给固定左尾下界，与竞争假设矛盾，故正鞍点确实存在。

**由粗到细的局部化。** 在真实鞍点 $\Psi(s)+sq\le0$，以 $N\ell(s)+sq$ 的曲率先得 $s=O(\rho)$。在 $s_*$ 附近的基准最小值是 $I_*=-z^2/2+O(\rho L)$。$s^2V_s\lesssim Ns^2+P$，由 $d(u,s)=(u-s)/(1+s)-\log[(1+u)/(1+s)]$ 分别控制 $|u-s|\le1$ 和其补集。将 (11) 与严格竞争尾预算合并后，
\[
 cN(s-s_*)^2+A+P/2\le C+\tfrac12\log(1+P/L).
\]
所以这三项均有统一上界，$s\asymp\rho$，每个 $u_i$ 都在一个固定紧区间，$\sum(u_i-s)^2\lesssim P$。特别有 $\sum(\lambda_i-1)^2=O(N/L)$，足量坐标的 $\lambda_i$、$a_i$ 同时有正的上下常界；重权坐标仍可随维数增长，不能在此直接套用有界权 Edgeworth。

**微分估计和精确 Mills。** $D(u,c)=(F(u,c)-\ell(u))/u$ 在零处光滑延拓，极小曲线 $D=D_c=D_u=0$，而 $D_{cc}=p>0$，紧区间内一致；$c\to\infty$ 端的导数由 $c^j e^{-c}$ 控制且极小点在内部。因此 (12) 的二次动作亏损及相应的一阶、二阶 $u$ 导数界均成立。由 $\sum u_i^3\lesssim Ns^3+P^{3/2}$ 和 Cauchy–Schwarz，
\[
 \sum u_i^2\min(1,|c_i-c(u_i)|)\lesssim A+P+\sqrt{A Ns^3}=R .
\]
在 $u\ell'(u)$、$u^2 F_{uu}(u,c(u))$ 处作二阶 Taylor，线性项因 $\sum(u_i-s)=0$ 而消失。利用 $s$ 是**真实**驻点，进一步得到 $Ns|s-s_*|\lesssim R$，故 $|b^2-b_*^2|\lesssim R$。若只用前一步 $N(s-s_*)^2=O(1)$，精度不足。

倾斜中心和与等均值、等方差高斯的逐因子 Fourier 望远镜中，即使一个坏的大权重被删掉，每项乘积仍有 $cN-1$ 个正则权因子，且已替换的高斯因子也提供同样的衰减。三阶绝对矩 Taylor 界在全实轴有效，积分 $|t|^3(1+c't^2/N)^{-(cN-1)/2}$ 一致有界，故得 (14) 的**密度**绝对误差 $O(N^{-1/2}+P^{3/2}/L^{3/2})$。直接积分指数权所付绝对误差为该量除以 $b$，而精确 Mills 因子 $M(b)\asymp1/b$，故相对误差同阶。其对数导数为 $O(1/b)$，两精确因子之比只损失 $O(R/L)$，未使用可能遗留 $O(1/L)$ 的各自渐近截断。指数率保留 $A+P/2$ 后，吸收 $\sqrt{ANs^3}/L$ 余 $O(1/\sqrt{NL})$，且 $P=O(1)$ 可吸收 $P^{3/2}/L^{3/2}$。因此同点真实 iid CDF 对照的对数比至少 $-C/\sqrt N$，即 (15)。

最后 $c(s_*)=\mu+O(\rho)$，基准 iid 均值 $\mu N+O(L)$（利用最优平均一阶导数为零），方差 $\mu N+O(\sqrt{NL})$。在长度 $\epsilon L/2$ 的目标区间，其原始密度至少 $cL/N$；这里 r093 的密度绝对误差 $O(1/N)$ 小于主密度。iid 命中分位已经由 r092 固定向量结论控制，所以 $F_*(q)\ge\tau+BL/N+c_\epsilon L^2/N$。乘 $e^{-C/\sqrt N}$ 在这个**下界**上仅损失 $O(\sqrt L/N)$，确实严格小于固定 $\epsilon$ 的余量。由此对每个固定 $\epsilon>0$ 排除全部竞争向量；等权向量给反向，(10) 成立。其应用精确止于真实无条件对数凸下集的全 halo 上系数；一般 $H_N$ 与 $C_N$ 仍开放。

## 3. 独立轻量本机核验

为检查倾斜后三个尺度（指数权、平移长度和 Bernoulli 概率）是否同时正确变换，我另外取 $\lambda=(1/8,3/8,1,5/2)$、侧比 $(.2,.7,.4,1)$、$s=.12,.5,1.1$。按 $2^4$ 个平移组合枚举，剩余不同尺度的指数和用精确部分分式 $d_i=\prod_{j\ne i}a_i/(a_i-a_j)$，分别计算原分布 CDF 和倾斜分布指数加权的截断积分。三个 $q$ 分别为 $3.6270838950,2.4799515720,1.7131336521$；原 CDF 为 $0.4857630783,0.2650019891,0.1211809773$；两种求值绝对差为 $6.11\times10^{-16},1.67\times10^{-16},2.78\times10^{-17}$。这是有限数据的浮点恒等式检查，未代替上面任何统一估计。零平移情形以 $r=1$ 而非 $r=0$ 测试，避免程序中使用无穷动作。

```python
from itertools import product
from math import exp, expm1, fsum, log, log1p, prod

weights = (.125, .375, 1., 2.5)
r = (.2, .7, .4, 1.)
c = tuple(-log(x) for x in r)
p = tuple(x/(1+x) for x in r)

def mixture(scale, prob, cut, integral):
    coeff = [prod(a/(a-v) for j, v in enumerate(scale) if j != i)
             for i, a in enumerate(scale)]
    terms = []
    for mask in product((0, 1), repeat=4):
        t = cut - sum(i*w*z for i, w, z in zip(mask, weights, c))
        if t <= 0:
            continue
        mass = prod(p if i else 1-p for i, p in zip(mask, prob))
        terms.append(mass*sum(d*integral(t,a) for a,d in zip(scale,coeff)))
    return fsum(terms)

for s in (.12, .5, 1.1):
    a = tuple(w/(1+s*w) for w in weights)
    u = tuple((1+s*w)*z for w,z in zip(weights,c))
    pt = tuple(1/(1+exp(v)) for v in u)
    q = sum(v+w*z*p for v,w,z,p in zip(a,weights,c,pt))
    psi = sum(log1p(exp(-d))-log1p(r0)-log1p(s*w)
              for d,r0,w in zip(u,r,weights))
    orig = mixture(weights,p,q,lambda t,a: -expm1(-t/a))
    tilt = mixture(a,pt,q,lambda t,a:
                   exp(-s*t)*(-expm1(-(1/a-s)*t))/(1-s*a))
    print(s, q, orig, abs(orig-exp(psi+s*q)*tilt))
```

## 4. 定稿增补 (16)–(18)

对非凹的 Lipschitz 屋面，原采样条件是 $s\le\bar Z+u(\theta-Z_V)=u(\theta)+w\cdot Z+D_\theta(Z)$。删掉题设的 $D_\theta>b_N$ 事件，对每个严格见证仍剩下严格大于 $p_N-AL/N$ 的概率。$w\cdot Z$ 无原子，可在全权重结论中用固定 $B+A$。所得 (16) 只在给定角向概率假设及有限矩费用下有内容，未为一般非凹源自动提供矩上界。

(17) 的比较可比 (10) 的系数精度更强。对 $q=s_p^{\rm det}+J/\sqrt L$，在 $J$ 选择以前先固定一个有余量的 $q_0+O_B(L)$ 窗口；每一固定 $J$ 最终都在其中。匹配行动虽随 $q$ 变，但它在 $s_p^{\rm det}$ 处的 CDF 至少为 $1-p_N$。相同窗口内稳健密度下界使 $F_*(q)$ 至少增加 $cJ\sqrt L/N$。对这一显式下界使用 (15) 后损失至多 $C'\sqrt L/N$，取固定 $J$ 足够大产生严格矛盾。等权可选、r093 下界及 (5) 合起来给真实最优体积比相对 $\exp(s_p^{\rm det})$ 为 $1+O_B(1/\sqrt L)$。这里不要求两者逐维精确相等。

对达到 $q=q_0-CL$ 的选择，匹配 iid 的一项 CDF 展开给 $F_*(q)=\tau+O_{B,C}(L^2/N)$，故原假设给两 CDF 对数比的**上界** $O(\rho L)$。用 (15) 第一行的**下界**，保留正的 $A+P/2$ 并吸收低阶量，可得 $A+P=O(\rho L)$。于是 $H\lesssim P/s^2=O(\sqrt{NL})$。$μ-m(c_i)\lesssim(d_i')^2+u_i^2$ 跟随最优动作处的一阶导数为零及紧区间光滑性，在远离最优处使用有限矩。由 (12)，$sD\lesssim A+\sum u_i^3=O(\rho L)$，所以 $D=O(L)$。法向满足 $H=N^2|\nabla_Vu|^2$；超出相应界的方向既无行动达到测试秩，其径向 halo 体积根据切平面上界直接至多是自身角向源体积的 $e^{q_0-CL}$ 倍。这是对方向的上界，没有声称所有好方向都实现高接收体积。
