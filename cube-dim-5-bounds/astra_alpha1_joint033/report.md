# alpha1-joint033：Levy–Gamma–Dirichlet 降维与默认预算排除

本轮把 \(\beta=1/2\) 的三个 rank‑2 联合负矩降为 2、3、4 维 Dirichlet 积分，其中三线类进一步成为 Gamma 闭式。只计算三个已有代表，没有重扫 560 三子集，没有重跑旧 Green 点值。

有限诊断给默认 joint-Hölder 预算
\(\mathcal I_{\rm joint}\approx7.24619395\times10^{-5}\)，比现有 014/029 预算 \(7.22243245\times10^{-5}\) 更弱。随后用 Dirichlet 多项式准确均值及 Jensen，严格证明该**默认精确 joint-Hölder 预算**至少

\[
\boxed{\mathcal I_{\rm joint}\ge
6.9820795642834831209183330858867\ldots\times10^{-5}
>6.89\times10^{-5}>6.8\times10^{-5}.}
\]

因此这套固定三子集权重和外层 Hölder，即使把全部联合矩精确算出，也不能达到本轮目标。**这不是实际 \(I\) 的下界，也不排除新权重、更大联合块或其它估计。** 未认证 m=4、5 的联合矩上界，没有新可采用的完整 \(I\) 上界或 \(C_5\) 界。

## 1. 输入、来源和完整接口

固定 \(\alpha=1,\beta=1/2,b=1/10\)，原独立正 stable 满足 \(\mathbb Ee^{-zS_i}=e^{-\sqrt z}\)。16 条方向仍为未单位化 \((1,\pm1,\pm1,\pm1,\pm1)\)。本轮只读取 `upper_rank2_interface` 的锁定整数商几何，其 manifest 四项 hash 全部通过；**没有沿用该来源中的 \(\beta=3/4\) Laplace 尺度**。

源归一化为：三个模型分别有 \(m=3,4,5\) 个独立共线类，

\[
Q_t(X)=\sum_{i<j}c_{ij}X_iX_j,\qquad
K_t(s)=\Lambda_t^s\mathbb E Q_t(X)^{-s}.
\]

\(c_{ij}>0\) 是 primitive 整数商方向的平方面积，\(\Lambda=16^2\det(NN^T)\)；不能省略非正交 annihilator 的 Gram 因子。读取的 \(X_i=\sum_{k\in\ell_i}a_kS_k\) 原系数是精确整数。对新 \(\beta=1/2\) 必须重算

\[
\kappa_i=\sum_{k\in\ell_i}\sqrt{a_k},\qquad
\mathbb Ee^{-zX_i}=e^{-\kappa_i\sqrt z}.
\]

实际精确结果为

|Gram 轨道|代表|\(m\)|\(\Lambda\)|新 \(\kappa\)|\(h_i=\kappa_i^2/4\)|\(\theta_t\)|
|---|---|---:|---:|---|---|---:|
|48|(0,1,2)|3|768|(8,8,8)|(16,16,16)|3/16|
|64|(0,1,6)|4|1024|(8,4,8,4)|(16,4,16,4)|3/8|
|112|(0,3,5)|5|1792|(12,2,6,6,2)|(36,1,9,9,1)|7/16|

程序另精确检查“所有原 \(S_i=1\)”时 \(Q=\Lambda\)。这不同于错误地令合并后的所有 \(X_i=1\)。全部正系数、补集单项式和来源 hash 保存在 `results.json`。

三子集真实 scope 的占用率仍为 \(\rho=247/256\)。置
\(s=b\rho=247/2560\)、\(m(s)=\Gamma(1+2s)/\Gamma(1+s)\)，已有 AM–GM/分数 Hölder 接口在新 beta 上为

\[
J\le\mathcal J_{\rm joint}:=
4\Gamma(3/2)m(s)^{3/\rho}
\prod_{t=1}^3K_t(s)^{\theta_t/\rho}.
\tag{1}
\]

前因子来自 \(16\Gamma(3/2)(16^5)^{-1/10}=4\Gamma(3/2)\)。对应完整面预算为

\[
I\le\mathcal I_{\rm joint}:=
c_1^{5/4}\frac{\pi^2}{60}\mathcal J_{\rm joint}^{5/4},
\qquad c_1=\frac1{32\pi^3}.
\tag{2}
\]

014/029 的全时间、全单位面、一般 \(L^1\) 接口是本轮输入，不重算它们的旧数值批，也不把旧 beta 的指数换名充当新证明。025、026、030 均不修改。

## 2. Reciprocal-Gamma 与全部前因子

对当前 Levy 归一化，独立
\(Y_i\sim\mathrm{Gamma}(1/2,1)\) 给 \(X_i=h_i/Y_i\)。定义次数 \(m-2\) 的正平方自由多项式

\[
\Psi(y)=\frac1\Lambda\sum_{i<j}c_{ij}h_ih_j
\prod_{k\notin\{i,j\}}y_k.
\]

于是**精确**有

\[
\frac{Q(X)}\Lambda=\frac{\Psi(Y)}{\prod_iY_i},\qquad
K(s)=\mathbb E\left[\Big(\prod_iY_i^s\Big)\Psi(Y)^{-s}\right].
\tag{3}
\]

若采用父任务的 \(Q^\star(y)=\sum c_{ij}\prod_{k\ne i,j}(y_k/h_k)\)，则
\(\Psi=(\prod h_i)Q^\star/\Lambda\)。这里把 \((\prod h_i)/\Lambda\) 完整吸收到 \(\Psi\) 中，因此下式不再另乘 \(\prod h_i^{-s}\) 或 \(\Lambda^s\)，也没有漏掉它们。

乘上 \(\prod Y_i^s\) 将各 Gamma shape 从 \(1/2\) 改为

\[
a=\frac12+s=\frac{1527}{2560},
\]

并产生因子 \([\Gamma(a)/\Gamma(1/2)]^m\)。在这个新乘积 Gamma 概率中，令 \(R=\sum Y_i\)、\(Z_i=Y_i/R\)。径向换元的 Jacobian 为 \(R^{m-1}\)，所以 \(R\sim\Gamma(ma,1)\)、\(Z\sim\mathrm{Dirichlet}(a,\ldots,a)\)，二者独立。由于 \(\Psi(RZ)=R^{m-2}\Psi(Z)\)，

\[
\mathbb E R^{-(m-2)s}
=\frac{\Gamma(ma-(m-2)s)}{\Gamma(ma)}
=\frac{\Gamma(m/2+2s)}{\Gamma(m/2+ms)}.
\]

因此完整公式为

\[
\boxed{K(s)=G_m(s)\,mathbb E_{\mathrm{Dir}(a^m)}\Psi(Z)^{-s},\quad
G_m(s)=\left[\frac{\Gamma(1/2+s)}{\Gamma(1/2)}\right]^m
\frac{\Gamma(m/2+2s)}{\Gamma(m/2+ms)}.}
\tag{4}
\]

Dirichlet 概率密度包含 \(\Gamma(ma)/\Gamma(a)^m\)；用概率期望写 (4) 时不能再额外乘或除这些因子。单纯形规范化也与 [NIST DLMF 5.14.2](https://dlmf.nist.gov/5.14.E2) 的积分恒等式一致。

全部 \(s>0\) 时径向 shape \(m/2+2s>0\)。单纯形边界可能使 \(\Psi\to0\)，不能当成处处有正下界；但原负矩有限，因为 \(Q\ge c_{ij}X_iX_j\)，独立 Levy 变量的负矩全部有限。正积分和换元均由 Tonelli 合法化，不需要把边界截掉。

## 3. 三线类精确闭式及另外两型的低维诊断

对 m=3，\(c_{ij}=1,h_i=16,\Lambda=768\)，故

\[
\Psi(z)=\frac{256}{768}(z_0+z_1+z_2)=\frac13
\quad\text{在单纯形上}.
\]

因此无需数值一维积分，直接有

\[
\boxed{K_{48}(s)=3^s
\left[\frac{\Gamma(1/2+s)}{\Gamma(1/2)}\right]^3
\frac{\Gamma(3/2+2s)}{\Gamma(3/2+3s)}.}
\tag{5}
\]

256-bit Arb 对精确有理 \(s\) 给

\[
K_{48}=0.6548927205346375423170642162655410573160\ldots.
\]

这是已认证的精确矩包围。m=4、5 只执行预登记的每型 65536 点、seed 33033 的一个 scrambled Sobol Dirichlet 批，未追加种子。用独立 stick-breaking
\(B_j\sim\mathrm{Beta}(a,(m-j-1)a)\)、\(Z_j=B_j\prod_{k<j}(1-B_k)\) 将积分维数降到 m−1，即 3、4；程序通过逆不完全 Beta 函数产生节点。

m=4 的归一多项式可简写为

\[
\Psi_{64}(z)=\frac1{16}\bigl(e_2(z_0,z_1,z_2,z_3)+3z_1z_3\bigr).
\]

m=5 的十个三次正单项式完整保留在 JSON，系数和为 \(225/448\)，没有用 m=4 的尺度代替。实际结果：

|型|\(K_t\)|证据|
|---|---:|---|
|48|0.6548927205346375|Gamma 闭式，Arb 认证|
|64|0.6577835822205504|65536 点单批 Dirichlet QMC，非认证|
|112|0.6519902417518968|65536 点单批 Dirichlet QMC，非认证|

m=4 的八块均值范围约 [0.65777348,0.65780198]；m=5 约 [0.65195316,0.65202863]。这些块仅是描述性检查，不是置信区间或严格误差界。有限样本所有 \(\Psi\) 正，无裁剪；最小样本单纯形坐标分别约 \(1.79\times10^{-10}\)、\(3.24\times10^{-13}\)，归一和最大浮点误差分别为 \(2.22\times10^{-16}\)、\(4.44\times10^{-16}\)。

代入 (1)–(2) 的非严格预算诊断为

\[
\mathcal J_{\rm joint}\approx2.0500703459633103,
\qquad
\mathcal I_{\rm joint}\approx7.246193949909558\times10^{-5}.
\]

这是当前默认 rank‑2 分组，和已优化的 014 flag 分组并非同一个预算表达式；保留块内联合信息并不自动保证它优于另一个已优化分组。本次没有竞争力，未启动高维积分认证。

## 4. Jensen 给默认预算的严格下界

为避免仅凭 QMC 停留在推测，本轮另登记一个有限低成本认证批。\(x\mapsto x^{-s}\) 在正半轴凸，所以

\[
\mathbb E\Psi(Z)^{-s}\ge(\mathbb E\Psi(Z))^{-s}.
\tag{6}
\]

\(\Psi\) 的每个单项式都为平方自由、次数 \(k=m-2\)。Dirichlet 的准确单项式矩为

\[
\mathbb E\prod_{j\in A}Z_j=\frac{a^k}{(ma)_k}\quad(|A|=k).
\]

于是只需把所有正系数求和，不需任何空间或概率样本。实际精确有理结果为

|型|系数和|\(\mathbb E\Psi\)|\(K_t\) 的 Jensen 下界（显示近似）|
|---|---:|---:|---:|
|48|1|1/3|0.6548927205346375|
|64|9/16|13743/554752|0.6488425629888607|
|112|225/448|20985561/11651335360|0.6178546279299093|

把 \(G_m(s)(\mathbb E\Psi)^{-s}\) 代入 (1) 的正幂右端，得到**方法预算**的下界，而非 \(J\) 本身的下界。256-bit Arb 实际给

\[
\mathcal J_{\rm joint}\ge
1.9900712939259878804209653886551359791689\ldots,
\]
\[
\mathcal I_{\rm joint}\ge
6.9820795642834831209183330858867283154751\ldots\times10^{-5}.
\]

代码用完整球分别与精确有理 \(68/10^6\)、\(689/10^7\) 比较，均严格为大于；没有从展示小数生成端点。`jensen_floor.json` 保留全部球、源 hash、运行 PID 和真值断言。

必须保持方向：已知 \(I\le\mathcal I_{\rm joint}\)，而本节证明 \(\mathcal I_{\rm joint}>6.89\times10^{-5}\)。这仅排除通过这个默认预算证明 \(I<6.89\times10^{-5}\)，**不推出真实 \(I>6.89\times10^{-5}\)**。它也未优化新 type 权重，未排除任何不同外层 Hölder 或更大联合块。

## 5. 复现、锁与未做事项

首批实际约 0.231 秒（PID 41082），Jensen 的固定 Arb 运算约 0.00050 秒（PID 41127），两脚本均退出 0。120 秒 checkpoint 按代表保存，完成后无活动句柄。所有计算使用本机 Python/SymPy/SciPy/python‑flint，不使用 MATLAB；没有增加随机点、全空间网格或旧 beta=.75 试验。

复现（在本目录）：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B levy_dirichlet.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B jensen_floor.py
```

只写 `astra_alpha1_joint033/`。旧 025、026、030 的白名单 hash 在最终交付前核验不变；未读取或重跑父任务提到的 `green_generic.json` 点值。014/029 报告只用于接口与基线比较，不重复其积分或 Gamma 批。

必要归档为报告、两个登记、两个最终脚本、三型结果、Jensen 证书和 manifest。没有全量 QMC 节点、重复输入快照或 checkpoint 副本。未提交、未推送。可复用的新成果是 (4) 的完整 Levy–Gamma–Dirichlet 公式、三线类闭式 (5)，以及默认机制的严格预算排除；不提供新的完整 \(I\) 上界。
