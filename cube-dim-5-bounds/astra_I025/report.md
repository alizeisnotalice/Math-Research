# I025：rank‑2 联合负矩、正确 trace 倾斜与默认预算排除

本轮结论有三种不同证据等级，不能混用。

1. **严格接口与公式**：三子集 rank‑2 联合矩的 AM–GM / 分数 Hölder 接口正确；真实占用率仍是 \(\rho=247/256\)。建立了适用于 \(0<s<1/2\) 的正级数，使用正确的 \(T^{-2s}\) trace 倾斜，没有调用此范围不成立的正 Wishart 积分。
2. **已实际认证的排除**：仅正级数的第 0、1 项就证明，默认 \(\theta=(3/16,3/8,7/16)\) 的**精确 rank‑2 Hölder 预算**大于 \(1.18\times10^{-5}\)。第 1 项预算常数约 \(1.26188590217706\times10^{-5}\)。这是该固定方法预算的下界，**不是真实 \(J\) 或 \(I\) 的下界**，也不排除改变分组、权重或 Hölder 机制后的方法。
3. **非严格诊断**：\(-1/2\) moment lift 得 \(I\approx1.38084415\times10^{-5}\)，劣于 I024；直接 joint 诊断约 \(1.29555211\times10^{-5}\)，只比父任务提供的 I024 \(1.29676880\times10^{-5}\) 略低。未认证 joint 级数的奇异尾，所以没有新可采用的辅助 \(I\) 上界，更没有新 \(C_5<32\)。

因此本轮不继续高精度认证被目标预算排除的 moment lift，也不为这个默认 rank‑2 预算作大规模尾积分。完整 \(C_5\) 上界仍为 32；本轮没有修改下界 \(g\)。

## 1. 固定输入、所有权与来源

\(\beta=3/4\)，\(b=3/20\)，16 个未单位化方向按
\(v=(1,\pm1,\pm1,\pm1,\pm1)\)、后四坐标以 \(-1\) 先的字典序编号 0–15。原变量 \(S_i\) 独立，\(\mathbb Ee^{-\lambda S_i}=e^{-\lambda^\beta}\)。设

\[
A=\sum_iS_iv_iv_i^T,\qquad
m(t)=\frac{\Gamma(1+t/\beta)}{\Gamma(1+t)}.
\]

沿用父任务和 I021/I022 已审的接口

\[
J=16\Gamma(7/4)\mathbb E(\det A)^{-b},\qquad
I\le c^{10/7}\frac{\pi^2}{60}J^{10/7},\quad
c=\frac3{2^{15/2}\pi^{5/2}}.
\]

此处原 \(S\) 概率和原问题的 shape 倾斜概率不混用。本轮后面引入的 \(T^{-2s}\) 是各二维商的 **trace 倾斜**，与原问题 \(T^{-\beta}\) 也不同。

独立 clone：`/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I`；clone 起始 HEAD 为 `9bea6c1e502ef8d9b9529f960cd620236a11978c`。已读实际根 `AGENTS.md`，`cube-dim-5-bounds/` 无额外 AGENTS。本代理只写 `astra_I025/`，不提交、不推送、不改其他项目；遵从父任务的专项限制。原 3/4 维项目与 ZIP 未访问或修改。I021、I022 的实际来源路径和 SHA‑256 固定在 `diagnostic.json`；I024 数字仅来自父任务消息作比较，未读取或重复优化它的权重。

## 2. 三子集接口与有限几何

对三子集 \(B\) 令 \(G_B=V_B^TV_B\)、\(d_B=\det G_B\)，并把所有方向正交投影到其二维正交补，记投影 \(w_i\)。设

\[
A_B=\sum_iS_iw_iw_i^T,\quad Q_B=A_B/16,\quad D_B=\det Q_B.
\]

原 frame 满足 \(\sum v_iv_i^T=16I_5\)，故 \(\sum w_iw_i^T=16I_2\)。Cauchy–Binet 中每个五方向基包含十个三子集，因此

\[
\det A=\frac1{10}\sum_Bd_B\Big(\prod_{j\in B}S_j\Big)\det A_B,
\qquad
\frac{\det A}{16^5}=\sum_B\theta_BF_B,
\]
\[
\theta_B=\frac{d_B}{10\,16^3},\qquad
F_B=\Big(\prod_{j\in B}S_j\Big)D_B,\qquad \sum_B\theta_B=1.
\]

三种 Gram 轨道确实覆盖所有三子集，不以 Gram 相同自动假定等价。证明如下：把一条 projective cube 方向通过坐标变号送到全 1；其余方向取距它至多 2 的反极代表，记变号集合 \(U,V\)。可能的 \((|U|,|V|,|U\cap V|)\) 为 \((1,1,0),(1,2,0),(1,2,1),(2,2,0),(2,2,1)\)。其中第一、第三种通过换锚点等价，Gram 行列式 48；第二、第四种通过换锚点和反极代表等价，行列式 64；最后一种行列式 112。坐标排列使每种标准化类型内部等价。固定锚点的有序另外两方向计数分别为 60、90、60，除以每个三子集的 6 个有序锚定表示，得 \(16(60,90,60)/6=(160,240,160)\)。

本轮仅核三个代表的精确投影，没有扫五方向基，也没有空间网格。

| \(d_B\) | 本轮代表 | 数量 | 总 \(\theta\) | 非零商方向 | 非零商线类 |
|---|---|---:|---:|---:|---:|
|48|(0,1,2)|160|3/16|12|3|
|64|(0,1,6)|240|3/8|12|4|
|112|(0,3,5)|160|7/16|13|5|

三被选变量的商投影为零；48、64 型还有一个额外零投影。故 \(F_B\) 的真实 scope 大小为 15、15、16。cube 对称性使各原坐标占用相同，

\[
\rho=\frac1{16}\sum_B\theta_B|\mathrm{scope}(F_B)|
=\frac{15(3/16+3/8)+16(7/16)}{16}=\frac{247}{256}.
\]

AM–GM 后用产品概率的 fractional Hölder（权 \(\theta_B/\rho\) 每坐标占用 1）得到，令 \(s=b\rho=741/5120\)，

\[
J\le\mathcal J_{\rm joint}:=
2\Gamma(7/4)m(s)^{3/\rho}
\prod_{d\in\{48,64,112\}}
\left(\mathbb ED_d^{-s}\right)^{\theta_d/\rho}.
\tag{1}
\]

选中的三变量与 \(D_B\) 独立，才可抽出 \(m(s)^3\)。这里没有把 shape 变量误当独立。普通 Hölder 的版本取 \(\rho=1\)、\(s=b\)。负矩有限：任取两条非平行商方向，\(D_B\ge aS_iS_j\)、\(a>0\)，所以所有所用负矩有限；截断和单调收敛足以合法化各正积分及 Hölder。

## 3. 一维正角积分与失败诊断

任意二维 SPD 矩阵 \(Q\) 对角化后，由 \(\tan\theta\) 换元得

\[
\det Q^{-1/2}=\frac1{2\pi}\int_0^{2\pi}(e_\theta^TQe_\theta)^{-1}\,d\theta.
\]

对 \(Q_B\)，Tonelli 与一维 stable 负矩给

\[
K_B:=\mathbb ED_B^{-1/2}
=\frac{m(1)}{2\pi}\int_0^{2\pi}
\left[\sum_i\left(\frac{(w_i\cdot e_\theta)^2}{16}\right)^\beta\right]^{-1/\beta}d\theta.
\tag{2}
\]

所有角上 \(\sum_i(w_i\cdot e)^2/16=1\)。因此方括号至少 1，无分母退化；各个投影零点只造成低阶光滑性损失，不造成积分发散。对 \(0<s<1/2\)，凹幂 Jensen 给 \(\mathbb ED_B^{-s}\le K_B^{2s}\)。代入 (1) 得

\[
\mathcal J_{\rm lift}=
2\Gamma(7/4)m(b\rho)^{3/\rho}
\prod_d K_d^{2b\theta_d}.
\]

登记后实际 SciPy 求积，按每个投影零角分段，仅三个代表：

|型|\(K_d\)，非认证诊断|
|---|---:|
|48|0.5684940575684908|
|64|0.5766494221814679|
|112|0.5709533990035084|

普通 Hölder：\(\mathcal J\approx1.48183631198009\)、\(\mathcal I\approx1.38322865863015\times10^{-5}\)；分数 Hölder：\(\mathcal J\approx1.48004770600045\)、\(\mathcal I\approx1.38084415485443\times10^{-5}\)。SciPy 误差估计不是向外舍入证书。这条路线没有竞争力，按父任务指示未作 Arb 角积分加密。

## 4. 保留真实 joint 矩的正级数

合并二维商内每条无向共线类。写

\[
Q=\sum_{\ell=1}^L X_\ell u_\ell u_\ell^T,
\quad \|u_\ell\|=1,\quad
X_\ell=\sum_{i\in\ell}\frac{\|w_i\|^2}{16}S_i.
\]

\(X_\ell\) 独立，Laplace 指数为 \(c_\ell\lambda^\beta\)，其中
\(c_\ell=\sum_{i\in\ell}(\|w_i\|^2/16)^\beta\)。令

\[
C=\sum_\ell c_\ell,\quad p_\ell=c_\ell/C,\quad
T=\operatorname{tr}Q=\sum_\ell X_\ell,\quad
z_\ell=(u_{\ell,1}+iu_{\ell,2})^2,
\quad Z=\sum_\ell(X_\ell/T)z_\ell.
\]

直接写出矩阵的两个对角元及非对角元可得

\[
\det Q=\frac{T^2}{4}(1-|Z|^2),\qquad |Z|<1\quad\text{a.s.}
\]

对 \(r>0\) 定义 **正确的 trace 倾斜**

\[
\frac{d\mathbb P_r}{d\mathbb P}
=\frac{T^{-r}}{m(r)C^{-r/\beta}}.
\]

取 \(r=2s\)、\(\mu_k(r)=\mathbb E_r|Z|^{2k}\)，二项式正展开和单调收敛给精确恒等式

\[
\boxed{\mathbb E(\det Q)^{-s}
=2^{2s}m(2s)C^{-2s/\beta}
\sum_{k=0}^{\infty}\frac{(s)_k}{k!}\mu_k(2s).}
\tag{3}
\]

无需 \(s\ge1/2\)。\(\mu_0=1\)，所有项非负；有限和只提供负矩的**下界**，不能误称上界。

有限矩可以通过稳定线性和计算，不必积分 stable 密度。对任意角 \(\phi\)，记
\(a_\ell=\Re(e^{-i\phi}z_\ell)\)、\(Y_\phi=\Re(e^{-i\phi}Z)\)。对实数 \(|t|<1\)，线性系数全部为正，精确有

\[
\mathbb E_r(1-tY_\phi)^{-r}
=\left[\sum_\ell p_\ell(1-ta_\ell)^\beta\right]^{-r/\beta}
=:g_\phi(t).
\tag{4}
\]

因为左边乘上 \(\mathbb ET^{-r}\) 后就是
\(\mathbb E[\sum X_\ell(1-ta_\ell)]^{-r}\)，其 stable 尺度已知。(4) 也可在小复邻域写成 \(uZ+v\overline Z\) 的双变量生成式；这里只需实变量版，避免额外分支约定。

对 (4) 在零处求导，并用圆周上 \(\cos^{2k}\) 的平均，得

\[
\mu_k(r)=\frac{4^k(k!)^2}{(r)_{2k}}
\frac1{2\pi}\int_0^{2\pi}[t^{2k}]g_\phi(t)\,d\phi.
\tag{5}
\]

如果 \(h_\phi(t)=\sum p_\ell(1-ta_\ell)^\beta=\sum h_nt^n\)、\(g=h^{-q}\)、\(q=r/\beta\)，则

\[
h_n=(-1)^n\binom\beta n\sum_\ell p_\ell a_\ell^n,\quad
g_0=1,\qquad
ng_n=-\sum_{j=1}^n(n-j+qj)h_jg_{n-j}.
\]

\([t^{2k}]g_\phi\) 是次数不超过 \(2k\) 的三角多项式。因此 \(2N+1\) 等距角节点在**精确算术**下足以一次计算所有 \(k\le N\) 的平均。当前 N=64、128 程序使用浮点/long double，故只作诊断，不宣称这些高阶系数已认证。

## 5. 实际联合矩诊断与严格低阶预算排除

独立正 stable 使用 Kanter 表示；每个代表仅合并后的 3/4/5 个变量，用两个预登记 scramble（1025、2025），各 65536 点。为避免接近 rank‑1 时的相减消失，程序按正项
\(\det Q=\sum_{\ell<j}X_\ell X_j(1-\Re(z_\ell\overline z_j))/2\)
求行列式。

|型|\(\mathbb ED^{-b\rho}\) 两次 QMC 均值，非严格|N=128 正部分和，浮点诊断|
|---|---:|---:|
|48|0.8134558807553008|0.8132125996413652|
|64|0.8167429694171863|0.8165061991181519|
|112|0.8144039447622610|0.8141649510167158|

QMC 的 \(\mathbb ED^{-1/2}\) 与 (2) 的独立角积分吻合到约 \(2\times10^{-5}\)；两次 QMC 的 \(-b\rho\) 矩相差约 \(4\times10^{-6}\)。这些是诊断交叉检查，不能代替严格误差界。代入 (1) 的诊断为

\[
\mathcal J_{\rm joint}\approx1.41544457190783,
\qquad \mathcal I_{\rm joint}\approx1.29555210621909\times10^{-5}.
\]

严格排除不需要尾估计，也不依赖 QMC。令 \(h=\sum p_\ell z_\ell\)，从 (4) 的二阶系数直接得到

\[
\mu_1(r)=\frac{1-\beta+(r+\beta)|h|^2}{1+r}.
\tag{6}
\]

本轮三个代表的 \(h\) 虚部精确为零。48、64 型 \(h=0\)；112 型 \(|h|^2\approx2.07670569316795\times10^{-6}\)。令

\[
L_{d,0}=2^{2s}m(2s)C_d^{-2s/\beta},\qquad
L_{d,1}=L_{d,0}(1+s\mu_{d,1}(2s)).
\]

由 (3)，\(\mathbb ED_d^{-s}\ge L_{d,1}\ge L_{d,0}\)。所以 (1) 的**右端预算常数**满足

\[
\mathcal J_{\rm joint}\ge
2\Gamma(7/4)m(s)^{3/\rho}\prod_dL_{d,n}^{\theta_d/\rho}
=:\mathcal J_{{\rm floor},n},\quad n=0,1.
\]

使用同一单调换算定义 \(\mathcal I_{{\rm floor},n}\)。256-bit Arb 在精确有理输入、平方根和有理幂上实际核得：

|正级数截断|预算 \(\mathcal J\) 下界常数（显示近似）|预算 \(\mathcal I\) 下界常数（显示近似）|
|---|---:|---:|
|N=0|1.35030607818643369596|1.21122657869464907942e−5|
|N=1|1.38959585468501465413|1.26188590217706455899e−5|

`budget_floor.py` 用完整 Arb 球分别与精确有理数 \(118/10^7\) 比较，两个比较均严格成立；不是从上表展示小数构造证书。完整球及输入/代码 hash 位于 `budget_floor.json`。所用 `arb` 的包围语义与比较算子语义已对照 [python-flint 官方说明](https://python-flint.readthedocs.io/en/latest/general.html)：比较返回真要求对整个输入球成立，输入必须避免非预期的二进制 float 转换。本证书全部数学输入来自整数、有理数及精确符号投影。这已经排除默认 (1) 到达 \(1.18\times10^{-5}\)，从而也排除它到达 \(1.1\times10^{-5}\)。**由于 \(J\le\mathcal J_{\rm joint}\)，对右端的下界绝不推出真实 \(J\) 或 \(I\) 的下界。**

## 6. 未闭合项和可恢复方向

未认证 (3) 的无限尾，故不能把约 \(1.29555\times10^{-5}\) 当成新上界。一个方向明确但可能偏弱的尾包络可保留供后续使用：设 \(x=|Z|^2\)、\(0<s<1/2\)。因为 \((s)_{N+1+j}/(N+1+j)!\le(s)_j/j!\)，

\[
\sum_{k>N}\frac{(s)_k}{k!}x^k\le x^{N+1}(1-x)^{-s}.
\]

所以 (3) 在原概率下的负矩尾至多

\[
\mathbb E D^{-s}x^{N+1}
\le K^{2s}\left(\mathbb E x^{(N+1)/(1-2s)}\right)^{1-2s}.
\]

右边最后的期望在**原概率**下，不能替换成 \(T^{-2s}\) 倾斜矩。它的整数矩可由 (4) 的 \(r\downarrow0\) 极限计算：
\(\mathbb EY_\phi^n=-(n/\beta)[t^n]\log\sum p_\ell(1-ta_\ell)^\beta\)。此尾式本轮未用于数值认证，亦不声称它能给出足够精度。

默认 rank‑2 精确矩的目标预算已经被低阶严格排除；继续突破需要保留更大联合块，或直接降低外层 AM–GM / Hölder 损失。此次排除不覆盖任意新权重、新分组或不经过 (1) 的方法。

## 7. 实际运行、复现与精简归档

所有三批在执行前分别登记输入、有限范围、120 秒 checkpoint、资源和句柄。本机 Python 使用 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python`，SciPy 1.18.1、SymPy 1.14.0、python‑flint 0.9.0。没有调用 MATLAB，没有使用全域空间网格。首次角积分运行约 0.109 秒；联合有限批约 0.280 秒；Arb 低阶预算认证约 0.118 秒，均退出 0。clone session 55393、joint session 11159 均完成；无遗留运行句柄。

在本目录依次运行：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B rank2_angle.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B rank2_joint_probe.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B budget_floor.py
```

必要归档为本报告、三个最终脚本、三个预登记、三份结果及 manifest。checkpoint 仅用于运行恢复，完成后无需归档；没有全量样本、重复输入快照或临时日志。未提交、未推送，交父任务独立审核及整合。
