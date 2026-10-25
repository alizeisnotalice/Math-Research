# alpha1-029：五维 Cauchy 模型完整上界资格独立重建

**理论资格通过，Astra014的完整 (I) 预算通过；同模型新 (g) 尚未锁定，联合 (C_5) 结论待核。** 本报告从五维Cauchy半群及原 (alpha=1) 原件重建一般 (L^1)、全部正尺度、全空间有限正源与cube domination，不以替换 (alpha=3/2) 指数作为证明。

已核完整面预算

\[
 I_1=\int_{[0,1]^4}G_1(1,y)^{5/4}dy
 <I_*:=\frac{36112162239}{500000000000000}.
\]

它与下述传递定理兼容，包括 (S_2) 上 (G_1=\infty)。若后续**同一个模型**的essential cube下界 (g_*>0) 被认证，则可以直接使用

\[
 C_5\le5^{5/4}I_*/g_*^{5/4}.
\]

使用该有理 (I_*) 达到成本32的准确门槛在

\[
 [0.000151983629808272418348117992983612438676962688957798593869,
  0.000151983629808272418348117992983612438676962688957798593870].
\]

本批未读取或联合未锁新 (g)，未预宣 (C_5<32)。主/他人文件只读，仅写独立clone的 `upper_alpha1_029/`，未提交、未推送、未用MATLAB。

## 1. 模型与锁定来源

固定Lebesgue、中心轴平行共同半径方体 (Q=[-1,1]^5)，

\[
 V=\{(1,\epsilon_2,\ldots,\epsilon_5):\epsilon_j=\pm1\},\quad
 \widehat p_t(\xi)=e^{-t\psi(\xi)},\quad
 \psi(\xi)=\sum_{v\in V}|v\cdot\xi|.
\]

16条方向不除以 (\sqrt5)，时间权均为1；Fourier逆因子 ((2\pi)^{-5})。审阅原件 `primary/finite_source_transfer.md`、`upper/round001_singular_strata.md`，并对照 `alternative_alpha/transfer.md` 的域/尺度差异。原 (alpha=1) 原件的主链正确；本报告把奇异点的正部步骤写成不含 (infty-infty) 的形式。

Astra014 manifest精确为 `7a65aa85088c1beaaaa8a6ac60974449d9175f5f0d37ef9c0c6327207dd8bea3`。其报告、概率类型、预算、脚本均锁入 `registration.json`；manifest七项原始字节只作hash检查，不执行诊断脚本。数值/QMC点不参与本次数学验收；未重跑Kanter或 (eta=3/4) 控制。

## 2. Cauchy半群、域及cutoff

取16个独立标准Cauchy变量，密度

\[
 s_t(z)=\frac{t}{\pi(t^2+z^2)},\qquad
 X_t=\sum_vvZ_t^v.
\]

选任一det-48五方向基，其五项和有正乘积密度，余十一项为概率卷积。所以 (p_t>0)、质量1、自伴，且

\[
 p_t(x)=t^{-5}p_1(x/t),\qquad
 \|p_t\|_\infty\le C t^{-5},\quad C=1/(48\pi^5).
\]

frame (\sum_vvv^T=16I) 给

\[
 4|\xi|\le\psi(\xi)\le16|\xi|.
\]

因此 (e^{-t\psi}) 及其任意多项式频率权都可积；正时间 (p_t) 实际为有限的 (C^\infty\cap C_0) 函数，((t,x)) 联合连续。不存在正时间无限密度。质量、缩放及近似恒等性给 (P_t\to I) 于 (L^1,L^2,C_0)。Fourier乘子直接给

\[
 \operatorname{Dom}_{L^2}L=H^1(\mathbb R^5),\qquad
 \operatorname{Dom}((-L)^{1/2})=H^{1/2},\quad L=-\psi(D).
\]

不是 (alpha=3/2) 的 (H^{3/2}/H^{3/4}) 域。对 (phi\in C_c^\infty)，标准Cauchy归一化给

\[
 L\phi(x)=\frac1{2\pi}\sum_v\int_\mathbb R
 \frac{\phi(x+sv)+\phi(x-sv)-2\phi(x)}{s^2}ds.
\]

二差分的 (L^p) 范数被 (min(s^2|v|^2\|D^2\phi\|_p,4\|\phi\|_p)) 支配，(p=1,\infty)，故 (L\phi\in L^1\cap C_0\cap L^\infty)。积分无须在零处作未控主值。其符号常数由 (int_0^\infty(1-\cos s)s^{-2}ds=\pi/2) 核对。

取光滑径向**非增** (0\le\chi\le1)，原点邻域为1，紧支撑，(chi_R(x)=\chi(x/R)\uparrow1)。生成元的一阶齐次性给精确式

\[
 \|L\chi_R\|_\infty=R^{-1}\|L\chi\|_\infty,\qquad
 \|L\chi_\varepsilon\|_1=\varepsilon^4\|L\chi\|_1.
\]

这两个幂分别用于无穷远质量与原点原子，不能换成 (R^{-3/2},\varepsilon^{7/2})。

## 3. 五维奇异集：S1 power、S2 log

令 (S_k) 为至多 (k) 条方向span的有限并。离 (S_{k-1}) 正距离的紧集要求非小Cauchy参数的方向秩至少 (k)：否则大参数和落在 (S_{k-1})，小参数和的范数不足以跨过固定距离。选 (k) 个独立非小方向补成五基；其受限密度各 (O(t))，其余 (5-k) 个峰值各 (O(t^{-1}))，未选参数为概率卷积。因此

\[
 p_t(x)\le C_Et^{2k-5}\quad(E\Subset\mathbb R^5\setminus S_{k-1},0<t\le1).
\]

先在开邻域作受限law的a.e.密度界，再由正时间连续性延至每一点。无需假定任意四方向独立。

若固定 (x\in S_k\setminus S_{k-1})，最小表示的 (k) 个系数非零。含它们的基中，限制余十一变量到 (O(t)) 邻域的事件具有固定正概率。前 (k) 个基密度至少 (ct)，其余至少 (c/t)，故同阶下界成立。尤其

\[
 p_t(x)\asymp_xt^{-3}\ (x\in S_1\setminus\{0\}),\qquad
 p_t(x)\asymp_xt^{-1}\ (x\in S_2\setminus S_1).
\]

截断 (int_\varepsilon^1p_tdt) 因而分别按 (arepsilon^{-2}) 和 (log(1/\varepsilon)) 发散；原点按 (t^{-5}) 缩放发散。离 (S_2) 的紧集有 (p_t\le C_Et)，与 (t\ge1) 的 (Ct^{-5}) 合成局部一致可积支配，所以 (G=\int_0^\infty p_tdt) 在 (mathbb R^5\setminus S_2) 有限且连续，真实无限集为 (S_2)。这不是四维“只在线上无限”的结论。

确切地，

\[
 (1,1,1,1,0)=\tfrac12(1,1,1,1,1)+\tfrac12(1,1,1,1,-1)\in S_2\setminus S_1,
\]

所以其Green为无穷，不能被有限样本均值当成有限点Green用于筛选。各pair plane在cube面的截面至多一维，在五维体积和四维面测度下均为零集。

## 4. 障碍分解：有限质量及H1极限

给有界紧支撑 (f\ge0)、(kappa>0)，零输入单独取零解。设 (F=\|f\|_\infty,m=\|f\|_1)。定义

\[
 V_t=t\sum_{j\ge0}P_{jt}f,\quad
 u_{t,0}=0,\quad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+.
\]

由 (|P_sf|_\infty\le\min(F,Cms^{-5}))，右Riemann和给

\[
 \|V_t\|_\infty\le tF+\frac54(Cm)^{1/5}F^{4/5}\le A
 \quad(0<t\le1).
\]

(V_t-tf\in C_0)，(P_tV_t=V_t-tf)。单调迭代的极限满足固定点式及 (u_t\le(V_t-\kappa t)_+)。固定t时这个上帽有本质紧支撑，故 (u_t\in L^1\cap L^2)；没有把可能不可积的 (V_t) 当作L1控制，也不需要随t的统一支撑。

令 (L_t=(P_t-I)/t)、(sigma_t=f+L_tu_t)。两分支给 (0\le\sigma_t\le\kappa)，在 (u_t>0) 处等于(kappa)。现在因 (u_t\in L^1)，保质量才合法给 (int\sigma_t=m)。能量身份

\[
 \langle u_t,-L_tu_t\rangle+\kappa\|u_t\|_1=\int fu_t
\]

给统一界 (|u_t|_1\le mA/\kappa)、(|u_t|_2^2\le mA^2/\kappa)、(|\sigma_t|_2^2\le\kappa m)、(|L_tu_t|_2\le\|f|_2+\sqrt{\kappa m})。

设 (psi_t=(1-e^{-t\psi})/t)。固定Q后 (t\le Q^{-1},\psi\ge Q) 时 (psi_t\ge(1-e^{-1})Q)，因此高频 (u_t) 的L2尾一致为 (O(Q^{-1}))。固定低频部分全局H1有界，局部Rellich配合高频尾及对角抽取给局部强L2、a.e.和全局弱L2子列。Fatou保留 (u\ge0,u\in L^1\cap L^2)。

对Schwartz测试 (L_t\phi\to L\phi) 于L2，故 (Lu=\sigma-f\in L^2)，椭圆符号给 (u\in H^1)。a.e.收敛使 (u>0) 处最终 (u_t>0)，用有界支撑测试及DCT得 (sigma=\kappa) 于 (D=\{u>0\})，且 (0\le\sigma\le\kappa)。先从弱极限及大cutoff得 (int\sigma\le m)，再对 (Lu=\sigma-f) 测试 (chi_R)：

\[
 |\langle u,L\chi_R\rangle|\le\|u\|_1R^{-1}\|L\chi\|_\infty\to0.
\]

所以回收完整质量 (int\sigma=m)，(kappa|D|\le m)。未以局部紧性替代质量守恒，也未在非局部接触集上假设 (Lu=0)。

## 5. 全空间正部源，含S2上的无限值

写 (G=U_T+W_T)，则 (|U_T|_1=T)、(|W_T|_\infty\le C/(4T^4))。取 (T=(C/(2\lambda))^{1/4}) 得

\[
 |\{G>\lambda\}|\le2(C/2)^{1/4}\lambda^{-5/4}.
\]

另对有限体积E，(int_EG\le\int_0^\infty\min(1,C|E|t^{-5})dt<\infty)。因此 (G\in L^1_{loc})、a.e.有限，而每个 (h>0) 的 (H=(G-h)_+\in L^1)。这些为Lebesgue类；canonical Green仍允许在S2逐点为无穷。

对 (phi\in C_c^\infty)，(int G|L\phi|\le\int_0^\infty\min(\|L\phi\|_\infty,Ct^{-5}\|L\phi\|_1)dt<\infty)。Fourier支配收敛给 ((P_t\phi-\phi)/t\to L\phi) 于L∞；积分半群导数、使用大时间 (P_T\phi\to0)，得到

\[
 \langle G,L\phi\rangle=-\phi(0),\qquad LG=-\delta_0
\]

于整个空间分布。

Tonelli和半群律给 (P_tG=\int_t^\infty p_sds<\infty)（处处有界），以及extended身份 (G=P_tG+U_t)。为避免S2上的 (infty-infty)，直接写正数形式

\[
 H=(G-h)_+\le(P_tG-h)_++U_t\le P_tH+U_t.
\]

它对extended非负代表有意义；与紧支撑测试配对的三项积分均有限。于是对所有非负 (phi\in C_c^\infty(\mathbb R^5))，取生成元极限得

\[
 \mu:=LH+\delta_0\ge0.
\]

这是**全空间**正分布，故为Radon测度。未删去S2、level界面或无穷尾的测试函数。大cutoff给

\[
 \int\chi_Rd\mu=1+\int HL\chi_R\to1,
\]

且 (chi_R\uparrow1)，故 (mu(\mathbb R^5)=1)。令 (a=\min(G,h))，由 (H=G-a) a.e.得 (mu=-La)。因 (a) 有界、(L\chi_\varepsilon\in L^1)，

\[
 0\le\int\chi_\varepsilon d\mu
 \le h\|L\chi_\varepsilon\|_1
 =h\varepsilon^4\|L\chi\|_1\to0,
\]

所以 (mu(\{0\})=0)。奇异Green值只是零Lebesgue集上的extended代表；不妨碍这个L1正部、分布方程或有限源测度论证。

## 6. Cube domination、同一接收集、全部尺度和一般L1

暂假设一个同模型认证essential cube下界 (G\ge g>0)。对 (0<h<g)，置

\[
 N=32(g-h),\quad K=H/N\in L^1,\quad M=\|K\|_1.
\]

则 (K\ge1_Q/32) a.e.，故其共同半径dilation逐尺度支配相同半径的中心cube平均。

全空间Fourier方程为 (-\psi\widehat K=(\widehat\mu-1)/N)。先在Schwartz函数上配对，再用L1 kernel/有限测度的L2 Young和H1稠密性，得到

\[
 K*Lu=(\mu*u-u)/N\quad(u\in H^1).
\]

因为 (psi(r\xi)=r\psi(\xi))，(K_r(x)=r^{-5}K(x/r))、(mu_r=(x\mapsto rx)_\#\mu) 满足

\[
 \boxed{K_r*Lu=(\mu_r*u-u)/(Nr).}
\]

此处是一阶尺度 (r^{-1})。对第4节障碍取 (kappa=\lambda/M)，同一个 (u,\sigma,D) 对每个固定r都有效。由 (f=\sigma-Lu)、(u=0) 于D外、(mu_r*u\ge0) a.e.，得

\[
 K_r*f\le K_r*\sigma\le\kappa M=\lambda\quad\text{a.e. on }D^c.
\]

先合并正有理r的可数零集，再用 (|K_r*f-K_s*f|_\infty\le\|f|_\infty\|K_r-K_s|_1\to0) 延至**全部**r。因此

\[
 \lambda|\{\sup_{r>0}K_r*f>\lambda\}|\le\lambda|D|\le M\|f\|_1.
\]

一般 (f\in L^1) 用 (f_j=\min(|f|,j)1_{[-j,j]^5}\uparrow|f|)。非负卷积MCT及 (sup_r\sup_j=\sup_j\sup_r) 同时给极大函数可测、严格超水平集的递增极限和质量极限。故所有一般L1输入及全部正尺度的中心cube maximal常数不超过M。没有改成非恒定背景、非中心cube或独立坐标尺度。

## 7. 面积分、层蛋糕及资格成本

Canonical Green满足 (G(rx)=r^{-4}G(x))，含无穷值。令 (V(\lambda)=|\{G>\lambda\}|)，弱L5/4先给 (V(1)<\infty)，齐次性给 (V(\lambda)=\lambda^{-5/4}V(1))。用10个cube面的锥参数及Jacobian (r^4)，非负Tonelli给

\[
 \sum_F\int_FG^{5/4}=5V(1)=160I_1,\qquad V(1)=32I_1.
\]

这本身即证明a.e.面L5/4有限，未假定在S2连续或处处有限。essential cube下界与essential boundary下界由齐次射线的正体积扫掠相匹配；不升级成任意逐点边界陈述。

层蛋糕给 (|H|_1=128I_1h^{-1/4})，从而

\[
 M=\frac{4I_1h^{-1/4}}{g-h}.
\]

对数导数为 ((5h-g)/(4h(g-h)))，唯一极小在 (h=g/5)，因此

\[
 \boxed{C_5\le5^{5/4}I_1/g^{5/4},\qquad
       \text{cost}^4=3125I_1^4/g^5.}
\]

这条传递只需同模型essential下界与完整I上界，**不要求G在S2有限**。

## 8. 新Gaussian shape预算怎样覆盖全部退化邻域

取 (eta=1/2) 的独立正stable (S_i)，(A=\sum S_iv_iv_i^T)。条件Gaussian covariance是 (2t^2A)，所以密度为 ((4\pi)^{-5/2}t^{-5}\mathbb E[(\det A)^{-1/2}e^{-x^TA^{-1}x/(4t^2)}])。直接非负时间积分

\[
 \int_0^\infty t^{-5}e^{-a/(4t^2)}dt=8/a^2
\]

给 (kappa_1=1/(4\pi^{5/2}))。A=TM带 (T^{-1/2})，且 (mathbb ET^{-1/2}=1/(8\sqrt\pi))，因此

\[
 G_1=c_1\mathbb E_{R_{1/2}}H_1,\quad
 c_1=1/(32\pi^3),\quad H_1=d^{-1/2}Q^{-2}.
\]

正时间mixture与第2节有限连续热核是同一代表：在零点的Gaussian期望由非负Fourier积分有限，提供 (F_x\le F_0) 的可积支配；时间Tonelli允许Green为无穷。

(H_1^{5/4}=d^{-5/8}Q^{-5/2})。固定正定M的球面积分与160面片轨道给 (L^{5/4}) 范数界 ((\pi^2/60)^{4/5}d^{-1/10})。R在全部坐标符号/置换共轭下不变，所以在形状期望用Minkowski：

\[
 I_1\le c_1^{5/4}\frac{\pi^2}{60}
      [\mathbb E_Rd^{-1/10}]^{5/4}.
\]

这里不在单个M的一个面片擅用全球面系数。旗分解是纯整数矩阵身份，不依赖beta；但原概率和的scale必须在beta1/2重算。五类型准确给

\[
 R_f^2=(81/8,27,36,36,64),\quad Z_f=2661/64,
\]
\[
 q_f^{\rm type}=(27,108,144,96,512)/887.
\]

由 (5(1/10)=\beta) 的精确总尺度抵消和全权熵最优，

\[
 \mathbb E_Rd^{-1/10}\le4\Gamma(3/2)
 [\Gamma(6/5)/\Gamma(11/10)]^5(2661/64)^{-1/10}<\infty.
\]

这个有限负矩保证形状Lq范数在R下可积；积分Minkowski可由非负截断及极限完成。因此它直接控制所有近退化M的尾，不删角点管或S2邻域。它的Lq结论是a.e.陈述；canonical Green仍可在S2取无穷，二者完全兼容。

## 9. 本机必要准确归一化批

本批只执行alpha1的有限检查：SymPy核时间积分、Cauchy符号常数、五类型平方根的准确有理平方、概率和、最优权和、成本导数及11110半和。未执行014的diagnostic、Kanter点或beta3/4控制。

独立常数路径取未移位Mellin比 (m=2\Gamma(1/5)/\Gamma(1/10))，并用

\[
 J^{10}=2^{10}\pi^5m^{50}/(2661/64),\qquad
 U_1^8=\frac{m^{50}}{2^{34}\,2661\,60^8\pi^9}.
\]

512-bit Arb包Gamma/π，150位Fraction整数幂夹8次/10次根。独立m、J、c1、I及精确公式的g门槛五个区间全部包含于014原384位锁定球，确认

\[
 U_1\in[
 0.000072224324476276829882877307945281140961817175378511064079,
 0.000072224324476276829882877307945281140961817175378511064080]
 <I_*.
\]

复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_alpha1_029/check_normalization.py
```

实际退出码0，约0.444秒；python-flint0.9.0、SymPy1.14.0，具体记录及精确端点见 `check_results.json`。没有几何搜索、概率/空间采样、旧β3/4Gamma基线重算或MATLAB。

## 10. 等待锁定 g 后的最终资格

目前完整有限源/全尺度/一般输入传递已通过，新完整I也通过。最后只可联合一个已经锁定并验收的 **alpha1、beta1/2、16未单位化权1、同一R/G/c1** 的essential单位cube下界 (g_0)。旧alpha3/2的g、有限点或必要筛选、未完成Bernstein/Hermite残差均不能替代该证书。

若采用有理 (I_*)，最终无须再跑Gamma；用原g证书的严格有理下端做准确比较即可：

\[
 \boxed{3125I_*^4<32^4g_0^5.}
\]

这里需要的安全g门槛是本报告开头的 (0.0001519836298082724183\ldots)，而不是精确Gamma上界表达式对应的略低 (0.0001519836298053715270\ldots)。在g原件锁定前，本批不执行该联合比较，不宣称C5突破。

必要白名单只有本报告、登记、归一化脚本、结果、manifest五项；不复制源文件或诊断缓存。后续g联合若获新指令另记源hash/比较结果，不改写已锁理论证明。
