# R6-006-C：α=3/2 body32 stable frame 独立模型首批

本报告只研究一个替代势模型。原问题仍为 Lebesgue 测度、中心轴平行、六坐标共同半径、一般 L1 输入及全部正尺度的 cube maximal operator。冻结的适用上界仍为 C6≤64。本报告没有认证新的 C6 上界，也不把本模型的势、极小值或成本与 α=1 Cauchy body32 的 g/I 混用。

已读主 clone 的 AGENTS、coordination/phase_20261025/protocol.md、frozen_thresholds.json、roles.json 及 STATE.md。只写本目录；未改旧 audit、未 stage/commit/push、未调用 MATLAB、未派生代理。以下解析论证供独立审计；数值批没有连续误差证书。

## 1. 固定模型、生成元与齐次性

记 a=3/2，d=6，h=d−a=9/2，β=a/2=3/4。取

\[
 V=\{(1,s_2,\ldots,s_6):s_i\in\{-1,1\}\},\qquad
 \psi_a(\xi)=\sum_{v\in V}|v\cdot\xi|^a.
\]

所有 32 个权重严格为 1。在完整坐标反射和置换群下，body 方向线是一个轨道，所以完整对称性要求这些权重相同；共同乘权只重标时间，最终成本不变，本批没有伪造权重优化自由度。

独立一维对称 stable 变量的特征函数为 exp(−t|η|^a)，故 X_t=ΣvZ_{v,t} 给出真实保正、保质量、自伴卷积半群，符号为 exp(−tψ_a)。也可用正 β-stable subordinator 对 Gaussian 作独立混合直接构造这些概率分布。

生成元在光滑紧支撑 φ 上为

\[
 L_a\phi(x)={c_a\over2}\sum_{v\in V}\int_{\mathbb R}
 {\phi(x+sv)+\phi(x-sv)-2\phi(x)\over |s|^{1+a}}\,ds,
 \qquad c_a={\Gamma(1+a)\sin(\pi a/2)\over\pi}.
\]

这里 1/2 与全实线积分同时保留。常数可直接核：将 s^(−1−a) 写为 Gamma 的 Laplace 积分，Tonelli 给

\[
 \int_0^\infty(1-\cos s)s^{-1-a}ds
 ={1\over\Gamma(1+a)}\int_0^\infty{r^{a-1}\over1+r^2}dr
 ={\pi\over2\Gamma(1+a)\sin(\pi a/2)}.
\]

故 L_a 的 Fourier 符号确为 −ψ_a。二差分在 0 附近为 O(s²)，在无穷远由 4||φ|| 控制；因此 L_aφ∈L1∩L∞∩C0，且

\[
 ||L_a\chi_R||_\infty=R^{-a}||L_a\chi||_\infty,
 \qquad ||L_a\chi_\epsilon||_1=\epsilon^{d-a}||L_a\chi||_1.
\]

Σvvᵀ=32I，ℓa/ℓ2 范数比较给

\[
 32^{a/2}|\xi|^a\le\psi_a(\xi)\le32|\xi|^a.
\]

于是生成元 L2 域为 H^(3/2)，形式域为 H^(3/4)。代码用 SymPy 准确核了 Gram、496 对均独立、编号 [0,3,5,14,22,25] 六行行列式为 +160。记

\[
 M={\Gamma(1/a)\over a\pi},\qquad C_0={M^6\over160}.
\]

用该六方向密度上界、其余 26 项概率卷积得 p_t≤C0t^(−d/a)=C0t^(−4)。缩放直接给

\[
 p_t(x)=t^{-4}p_1(t^{-2/3}x),\qquad
 G_a(rx)=r^{-9/2}G_a(x),\quad G_a=\int_0^\infty p_tdt.
\]

## 2. 奇异层：rank2 已可积，rank1 仍发散

定义 S1=∪_{v∈V}span(v)。本模型的结论是

\[
 G_a(x)<\infty\text{ 且连续于 }\mathbb R^6\setminus S_1,
 \qquad G_a(x)=+\infty\text{ 于 }S_1.
\]

先给不依赖多维渐近定理的上界论证。一维 stable 密度 ℓ 满足 ℓ≤M，并有显式尾包络

\[
 \ell(s)\le A|s|^{-5/2},\qquad
 A={\Gamma(5/2)\over\pi\sin(\pi/3)^{5/2}}.
\]

其推导：对 x>0 的 Fourier 积分旋转到角 θ∈(0,π/(2a))；旋转后的常数项积分为 i/x，实部为 0；|exp(−z)−1|≤|z| 当 Re z≥0，故 ℓ(x)≤Γ(1+a)/(π(x sin θ)^(1+a))。令 θ↑π/(2a) 得上式。轮廓弧在任意严格小于该角的扇区因 exp(−R^a cos(aθ)) 衰减，最后由上述可积包络取极限。于是 ℓ_t(s)≤At|s|^(−5/2)。下述线发散只需一维双侧尾可比性，不需该上包络的锐常数。

若 K 与 S1 距离正，取足够小 δ。任何落到 K 邻域的参数 fiber 至少有两个不同参数绝对值≥δ；否则去掉最多一个方向后，其余 31 项的向量和离该方向线过近。任意两条不同 body 方向独立，可补成六方向基。两个受限密度各 O(t)，其余四个基密度各 O(t^(−1/a))，对余参数只作概率积分。有限事件求和给

\[
 \sup_{x\in K}p_t(x)\le C_Kt^{2-4/a}=C_Kt^{-2/3}\quad(0<t\le1).
\]

事件可各用不同的补基；对开放邻域先以 Tonelli 得到 a.e. 密度界，再由真实 p_t 的连续性升级逐点。短时 t^(−2/3) 与长时 C0t^(−4) 均可积，DCT 给线外连续性。

对 x=bv≠0，选含 v 的基。令剩余 26 个缩放参数的基坐标和位于一个固定小盒；该概率严格正且不随 t 改变。指定的大坐标密度在 b+O(t^(1/a)) 上下有界于正倍 t，其余五个基坐标在 O(t^(1/a)) 内密度为正倍 t^(−1/a)。故 p_t(bv)≥ct^(1−5/a)=ct^(−7/3)，Green 发散；原点由 p_t(0)=t^(−4)p1(0)>0 同样发散。

这也核对了题设猜算的适用限定：k 是最少需要的独立 signed support ray 数，不是任意选出的 k 个表示方向。Watanabe 的有限谱定理对相应最小 stratum 给 p1(rξ)≍r^(−(1+a)k)，缩放后恰为 t^[k−(6−k)/a]。本批准确指数表为 −7/3、−2/3、1、8/3、13/3、6。尤其 rank2 短时可积；不能因某点同时属于更大 span 而给它错误的更快指数。

实际 SciPy 在 determinant-160 六方向子模型上核了 t=10^(−2),…,10^(−5) 的正乘积密度。rank1 最后对数斜率 −2.33320858，rank2 为 −0.66641716，分别趋近 −7/3、−2/3。它们仅校验子模型与缩放实现；完整 frame 的解析结论由前面的事件证明支持，不由这几个浮点点值证明。

## 3. 正 face 只有一个奇异角，且有严格幂次帽预算

在 x=(1,y)，y∈[0,1]^5 上，S1 只在 y=(1,…,1) 与 face 相交。记 ρ=dist(x,S1)>0。任取充分小邻域 U，使 x1>1/2 且 dist(U,S1)>ρ/2。每个落到 U 的 fiber 有某参数≥1/64；删除这项后，又必须有另一参数≥ρ/(62√6)。任意非零 6×6 body 行列式是 32 的整数倍（减去首行后五行均为偶数），所以任何补基的绝对行列式≥32。

单大参数事件和双大参数事件分别给

\[
 p_t(x)\le C_1t^{-7/3},\qquad
 p_t(x)\le C_2\rho^{-5/2}t^{-2/3},
\]

\[
 C_1=AM^5 64^{5/2},\qquad
 C_2=31A^2M^4 64^{5/2}(62\sqrt6)^{5/2}.
\]

这两个界均对全部 t>0 成立，源于同一真实密度。按 t=ρ^(3/2) 拆分得

\[
 G_a(x)\le D\rho^{-2},\qquad D=3C_2+\tfrac34C_1.
\]

对非负 x，最近的 body 线是全正方向线。令 r=||1−y||2，则

\[
 \rho^2=\sum_{i=1}^5(1-y_i)^2-{(\sum_{i=1}^5(1-y_i))^2\over6}
 \ge r^2/6.
\]

因此 G_a(1,y)≤6D r^(−2)。取 q=4/3，δ≤1，五维正正交象限球积分准确给

\[
 \int_{||1-y||_2<\delta}G_a(1,y)^{4/3}dy
 \le {(6D)^{4/3}\pi^2\over28}\delta^{7/3}.
\]

这是明确覆盖 singular corner 的解析预算，不是丢弃零测奇异点。其常数非常粗：本机浮点 D≈2.1706×10^9，帽系数≈1.0801×10^13；这些小数未作区间认证，实际定理使用上面的准确 Gamma 表达式。此预算当前不足以形成竞争性积分证书，但确认 rank2 poles 消失后，待支付的 face 奇异结构只剩一个角点。

## 4. 正部源、保守卷积、一般 L1 与全部尺度重推

下面给替代模型自己的条件性接口，不继承 α=1 接受标签。

由 p_t≤C0t^(−4)，对任意有限测集 E，∫E G_a≤∫0∞min(1,C0|E|t^(−4))dt<∞。令 U_T=∫0^T p_tdt、W_T=∫T∞p_tdt，则 ||U_T||1=T、||W_T||∞≤(C0/3)T^(−3)。取 T=(2C0/(3λ))^(1/3) 得

\[
 |\{G_a>\lambda\}|\le2(2C_0/3)^{1/3}\lambda^{-4/3}.
\]

故 H_c=(G_a−c)+∈L1，G_a∈L1loc。对 φ∈C∞c，L_aφ∈L1∩L∞ 使 ∫G_a|L_aφ| 可积，半群积分给 L_aG_a=−δ0。Tonelli 及半群律给 G_a=P_tG_a+U_t a.e.；正部凸性和 1-Lipschitz 性给 P_tH_c≥H_c−U_t。因此

\[
 L_aH_c=\mu_c-\delta_0,\qquad\mu_c\ge0.
\]

此分布正性处理了全部方向线及交叠，不依赖线外点值延拓。用大 cutoff 和 ||L_aχR||∞=O(R^(−3/2))，得 μc(R6)=1。又令 h=min(G_a,c)，μc=−L_ah；小 cutoff 的 L1 缩放为 ε^(9/2)，故 μc({0})=0。

对非负、有界、紧支撑 f，令 F=||f||∞、M_f=||f||1，并取 κ>0。对 t≤1 定义

\[
 V_t=t\sum_{j\ge0}P_{jt}f,\quad
 u_{t,0}=0,\quad u_{t,k+1}=(P_tu_{t,k}+t(f-\kappa))_+.
\]

右 Riemann 和给 ||Vt||∞≤F+(4/3)(C0 M_f)^(1/4)F^(3/4)=:A_f。Vt−tf∈C0，所以递增迭代的最小固定点 ut≤(Vt−κt)+ 对固定 t 有本质紧支撑。令 Lt=(Pt−I)/t、σt=f+Ltut，则 0≤σt≤κ，σt=κ 于 ut>0，且质量守恒给 ∫σt=M_f。能量等式

\[
 \langle u_t,-L_tu_t\rangle+\kappa||u_t||_1=\int fu_t
\]

给 ||ut||1≤M_f A_f/κ，||ut||2²≤M_f A_f²/κ，以及 ||Ltut||2≤||f||2+√(κM_f)。符号 ψt=(1−exp(−tψa))/t 在 t≤Q^(−1)、ψa≥Q 时至少为 (1−e^(−1))Q，给统一高频 L2 尾 O(Q^(−2))。频率截断的局部紧性使某 t↓0 子列局部强 L2/a.e. 收敛到 u，σt 弱收敛到 σ。测试传极限、Fatou 及大 cutoff 质量回收得到

\[
 f=\sigma-L_au,\quad u\in L^1\cap L^\infty\cap H^{3/2},
 \quad0\le\sigma\le\kappa,\quad\sigma=\kappa\text{ a.e. on }D=\{u>0\},
 \quad\int\sigma=\int f.
\]

H^(3/2) 来自 ψa û∈L2 与 ψa≍|ξ|^(3/2)，不是照搬 H1。特别是 κ|D|≤||f||1。

如果未来独立认证 g_a≤ess inf_{[-1,1]^6}G_a 且 g_a>0，取 0<c<g_a，N=64(g_a−c)、K=H_c/N，则 K≥1/64 于该 cube。令 Kr(x)=r^(−6)K(x/r)，准确的 Fourier 缩放是

\[
 K_r*L_au={1\over Nr^{3/2}}(\mu_{c,r}*u-u).
\]

因此在同一个 D 的补集，Kr*Lau≥0，继而 Kr*f≤κ||K||1。取 κ=λ/||K||1；正有理尺度的可数零集先合并，再利用 bounded f 下 ||Kr−Ks||1→0 所给的卷积连续性覆盖全部 r>0。cube 平均被 Kr*f 控制，故坏集测度≤||K||1||f||1/λ。最后以 fj=min(|f|,j)1_{[-j,j]^6}↑|f|，单调收敛扩展到一般 L1。这一推导保持中心 cube 的共同半径，未把各坐标尺度分离。

## 5. α=3/2 自己的成本式

令 I_a=∫[0,1]^5 G_a(1,y)^(4/3)dy。以上 weak Lq 与齐次性已保证其有限；仍没有认证数值上界。12 faces、各面 32 象限、坐标置换对称给 J=6·64 I_a。直接沿 cube 的径向坐标积分（Jacobian r5）得

\[
 \int(G_a-c)_+={h\over a d}Jc^{-a/h}.
\]

所以条件性模型成本为

\[
 \mathcal C_a(c)={h\over a}{I_a c^{-a/h}\over g_a-c}.
\]

令 c=a g_a/d=g_a/4；导数为零且端点均发散，所以这是唯一最小值：

\[
 \boxed{\mathcal C_{3/2}=4^{4/3}I_a/g_a^{4/3}}.
\]

SymPy 准确验证了导数与化简。只有同一 G_a 的连续 g_a 下界、I_a 上界、严格成本<64 及独立审计都闭合后，才能用上一节传给原 C6。当前未做到。

## 6. 有限数值批、正混合公式与停止结果

β=3/4，令 Sv 独立、E exp(−λSv)=exp(−λ^β)，A=ΣSv vvᵀ。条件 Gaussian 密度为 (4π)^(−3)det(A)^(−1/2)exp(−xᵀA^(−1)x/4)。Tonelli 后把时间解析积分，得到

\[
 G_a(x)={a\,\Gamma((6-a)/2)\over2^{a+1}\pi^3}
 \mathbb E\left[\det(A)^{-1/2}(x^TA^{-1}x)^{-(6-a)/2}\right].
\]

该公式允许值为 +∞，尤其不可在 body 线上把有限样本均值当有限 G。它不含交错 Fourier 锥求和或时间求积。Sv 的 Kanter 采样式为

\[
 S={\sin(\beta U)\over\sin(U)^{1/\beta}}
 \left({\sin((1-\beta)U)\over E}\right)^{(1-\beta)/\beta},
 \quad U\sim\mathrm{Unif}(0,\pi),\ E\sim\mathrm{Exp}(1).
\]

登记的单一模型、16 空间参数包括 8 个固定 face Sobol 点、7 个固定边界探针、1 个解析 singular corner。只对 15 个有限点运行，两份独立 scramble，各 4096 个协方差样本。64 维 QMC 是随机变量积分维度，不是 64 个拟合参数。未优化 α、权重或点位，未增加样本规模。补入解析常数后原输入原 seed 原样本作一次确定性复现，没有扩大有限批。

两副本形式成本为 16.4804206208、16.2734107534；经验最小点均为 (1,.99,.9,.8,.45,.1)，势均值 2.24216e−5、2.27127e−5。8 个 face 点的相对副本差不超过 1.74%；近 body 线点差为 39.0% 与 37.9%，单个最大样本可占均值分子的 47.3%。最大协方差条件数约 1.50e5。Laplace 值与负矩检查约 1e−4 相对精度支持归一化实现，但不给 QMC 误差保证。代码字段明确将传统 sample standard error 标为不是 QMC error bound。

因此该模型值得后续解析/区间机制研究，但当前 MC 已触发停止，不扩大重尾采样。上述低于 64 的数值不可登记为突破：有限点 min 不是 cube g 下界，8 点平均不是 I 上界，角点帽粗常数极大，且随机混合积分没有误差证书。

下一有限批的合理验收对象是：先缩小显式角点帽常数；再在离角点正距离的固定空间盒，对正 Gaussian 混合的参数尾部给确定误差预算。若该尾预算无法做到竞争性，应停止该数值表示，不能只扩大样本来获得更漂亮的形式成本。

## 7. 来源、hash、真实运行句柄

原始 stable 核估计实际核读 [Watanabe 2007 原论文](https://www.ams.org/journals/tran/2007-359-06/S0002-9947-07-04152-9/)：一维尾可比性、strict stable 缩放、有限谱支持分层分别见 Theorem1.2、式(2.3)、Theorem3.10/Remark3.11。官方 PDF 返回403，实际读到的是原论文公开镜像，不以搜索摘要代替定理核验。[可复取镜像](https://artefacts-discovery.researcher.life/full_text/DA-2/03/03ea54ee0d103fab984b9a92829cb8af/full_text/04d1990cdd74cf68186f726d0f836c18.pdf) 的426814字节 SHA256 为 29369f085260e598d859f1533485f996c4d5ecf770535ffb7949d0590c30e926。未留存 PDF 副本。

采样式实际核读 [Simon 原论文 arXiv:1310.1888v2](https://arxiv.org/pdf/1310.1888v2) 的(1.1)、第3节(3.1)，由 Mβ=Sβ^(−β) 代回核到上述式子；332399字节 SHA256 为 1c47ee810f3f963ac894a22ecac0f7e26f2765c4dc84bb1085e106dd870334c3。它明确溯源 Kanter1975 Corollary4.1；Kanter原文未获取，不能称其原文已经核读。源明细含在 sources.json。

复现命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
  cube-dim-6-bounds/upper/round006/alternative/stable_frame_batch.py
```

工作目录为 /Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit。最终实际调用退出0，程序内部 wall_seconds≈0.1461s，工具输出 chunk85a9c1；无存活计算句柄。Python3.15.0、NumPy2.5.4、SciPy1.18.1、SymPy1.14.0。主输入 hash：eae97103aca20226fb87a7736eb923afd45dab024c4708634b019123aa839271。源码和注册 hash 由脚本写入 stable_frame_results.json，文件清单 hash 另见 manifest.json。
