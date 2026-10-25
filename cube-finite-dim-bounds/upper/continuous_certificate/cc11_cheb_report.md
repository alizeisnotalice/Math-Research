# 11J：固定三维 Green 的连续全域下界证书

2026-10-25。固定 (\tau=1/3)，真实生成元为

\[
L=-\sum_{v=(1,\pm1,\pm1)}|v\cdot D|
-\tfrac13\sum_{i<j,\ \sigma=\pm1}|(e_i+\sigma e_j)\cdot D|.
\]

**本地完整证书 PASS，已送原 independent_audit，独审结论尚待回传。** 136 个闭 patch 全部在根节点通过 Bernstein 下界；没有使用子域分支、没有变更初始精度或表示、没有扩大网格。依赖锁定真实热核的 96 锥表示，证书证明

\[
G(x)\ge13/2500=0.0052\qquad(x\in[-1,1]^3).
\]

本批只认证 minimum 接口，不认证 J、最终 maximal 常数、一般 L1 或共同全部尺度迁移；当前有效 (U_3=8,U_4=16) 不因本报告单独改变。未提交或推送。

## 1. 认证结果与运行

执行环境为已有 venv，Python 3.15.0、python-flint 0.9.0、Arb 384 bits。四代表 patch 先端到端通过，耗时 2.02165 秒；随后固定完整 136 patch 批耗时 70.22595 秒，复用其中 4 个已计算节点缓存。其余 132 patch 均实际新算，未作失败后参数调整。

全域最小认证下界来自 patch 

\[
y\in[15/16,1],\qquad z\in[7/16,1/2],
\]

其数值为

\[
g_{\rm cert}=0.0052426479787179679751287605981257529687\ldots.
\]

这是 Bernstein 多项式下界扣除统一插值余项后的严格 dyadic 下界，不是 sampled minimum。全部 patch 的准确 dyadic 数在 cc11_cheb_all_result.json 的 `H_lower` 四元组中，显示小数仅供阅读。对目标的余量大于 (4.26479\times10^{-5})。同一证书也严格高于有理数 (131/25000=0.00524)；是否采用更强 glow 由主在独审后选择。

最大节点球半径小于 (2.296921\times10^{-38})，最大 Bernstein 系数球半径小于 (2.042669\times10^{-26})。节点误差已通过整个矩阵转换传播；下界中仅额外扣一次解析插值余项。

## 2. 正则化表示及所有退化投影

取

\[
H_\epsilon(x)=\int_\epsilon^\infty p_t(x)\,dt,\quad\epsilon=1/50.
\]

真实热核非负，故 (G\ge H_\epsilon) 于所有点；没有短时间误差需要支付。对锁源 48 个正负配对锥，记 (b_j=r_j\cdot x/k_j)、(k_j=\psi(r_j))、(W=2/(k_0k_1k_2))。则

\[
H_\epsilon(x)=\frac1{8\pi^3}\sum_C W f_\epsilon[b_0,b_1,b_2],
\qquad f_\epsilon(b)=\tfrac12\log(\epsilon^2+b^2).
\]

准确相等节点由 Fraction 判断，分别使用 Hermite 规则

\[
f[a,a,b]=\frac{f(b)-f(a)-f'(a)(b-a)}{(b-a)^2},\quad
f[a,a,a]=f''(a)/2,
\]

其中 (f'(a)=a/(\epsilon^2+a^2))、(f''(a)=(\epsilon^2-a^2)/(\epsilon^2+a^2)^2)。零节点完全正常；不会用宽球相交来猜测准确重节点。所有 signed 锥项均保留，任何近重节点相消均由 Arb 球计算包围。原 (G) 的 body/face 方向奇线在截断核中也被直接覆盖。

公式和 rational-point evaluator 来自首批已锁 cc11_regularized.py，本批不覆盖其文件。该 evaluator 已在首批对 5 个含退化情形的点与独立形式的 Arb 正时间有理积分重叠校准；本批的共享理论依赖仍是锁源 96 锥表示对应真实指定生成元。

## 3. 全纯界与有理节点替代

程序以准确 Fraction 逐 ray 验证

\[
k_r\ge3(|r_y|+|r_z|).
\]

由锥上线性及三角不等式，得到所有频率上的 

\[
\psi(\xi)\ge3(|\xi_y|+|\xi_z|).
\]

因此 (F(y,z)=H_\epsilon(1,y,z)) 在 

\[
\max(|\Im y|,|\Im z|)<3\epsilon
\]

全纯，且虚部上界为 (s<3\epsilon) 时

\[
|F|\le\frac{15}{1088\pi^3(\epsilon-s/3)^2}.
\]

此式直接由 Fourier 绝对值与 (p_1(0)=15/(544\pi^3)) 积分得出，实部任意。

每个 patch 半宽 (h=1/32)。21 个 Chebyshev roots 为 (t_j=\cos((2j+1)\pi/42))。正 roots 用 Arb 确认其 (2^{128}t_j) 的唯一 floor，再取精确 dyadic (q_j)；负 roots 使用准确反对称值，中央 root 为准确零。由显式检查，

\[
|q_j-t_j|<2^{-128},\qquad|x_j^{\rm rational}-x_j^{\rm exact}|<h2^{-128}.
\]

在任意实点，取 Cauchy 半径 (R=1/50)，每个坐标导数满足

\[
|\partial_yF|,|\partial_zF|\le
D=\frac{15}{1088\pi^3(\epsilon-R/3)^2R}
<125.056215.
\]

所以每个准确 rational-point enclosure 加半径 (2Dh2^{-128}<2.296921\times10^{-38})，就包围真实 Chebyshev 节点的函数值。此步骤避免了无理节点之间可能的代数相关性丢失，也没有假定节点误差仅经一个不适用的低精度 Bernstein 放大界。

## 4. 插值余项的严格支付

固定每维 degree (n=20)，Bernstein ellipse 参数 (\rho=3)。每个 ellipse 的虚部高度

\[
s=h(\rho-\rho^{-1})/2=1/24,\quad
\epsilon-s/3=11/1800,
\]

故全体 patch 统一可用 (M<11.906178395825)。一维 Chebyshev 系数满足 

\[
|a_k|\le2M\rho^{-k}\quad(k\ge1).
\]

roots 插值把每个 (T_k) alias 到 (\pm T_l) 或零，区间绝对值不超过 1，因此一维插值余项不超过 (4M\rho^{-n}/(\rho-1))。从 cardinal 的 DCT 展开

\[
\ell_j(t)=\frac1{n+1}\left[1+2\sum_{k=1}^n\cos(k\theta_j)T_k(t)\right]
\]

直接得到保守 Lebesgue 界 (\Lambda_n\le2n+1)。用 (I-I_yI_z=(I-I_y)+I_y(I-I_z))，二元插值误差最多

\[
E=(2n+2)\frac{4M\rho^{-n}}{\rho-1}
<2.868313237155\times10^{-7}.
\]

程序计算的是此量的 Arb enclosure，每个 patch 下界减去整个 enclosure 后再取向下端点；不使用四舍五入小数作为误差常数。

## 5. 连续多项式下界和边界覆盖

Arb DCT 从包围准确节点值的矩阵 V 产生 Chebyshev 系数。对 (T_k(2u-1)=\sum_l a_{kl}u^l)，degree-20 Bernstein 转换矩阵准确取

\[
B_{jk}=\sum_{l\le\min(j,k)}a_{kl}\frac{\binom jl}{\binom{20}l}.
\]

系数由整数递推构造，使用 Fraction 保存，直到转 Arb；全部四则和 DCT 三角函数均包围计算。二元 Bernstein 系数矩阵为 (BAVA^TB^T)。在闭单位方形上，每个 Bernstein 基函数非负且总和为一，所以多项式的整个闭域下界为系数的最小下端点。不是节点值比较。

本地登记允许至多 binary depth 8 的 de Casteljau 分支；实际 136 patch 均在根节点通过，完全未使用该余量。

使用准确闭矩形

\[
[i/16,(i+1)/16]\times[j/16,(j+1)/16],\quad0\le j\le i\le15.
\]

它们覆盖 (0\le z\le y\le1)，含所有边界与 diagonal。diagonal patch 验证整个方形，故有额外覆盖而无缺口。对非零 cube 点，cube 符号/置换对称和 (G(rx)=r^{-2}G(x)) 将其送到此三角面；缩小半径只增加 G。原点有 (G(0)=\infty)。因此下界是全 cube 的逐点结论，不是 a.e. 忽略奇线，也不假设 edge 达到最小值。

## 6. 自检、重现、锁定与审核边界

新增 cc11_cheb_selfcheck.py 实际完成：

- 21 个准确有理点上逐项验证 21 个 Bernstein 基底转换，共 441 个准确等式；两侧均 degree≤20，因而也证明了转换的多项式恒等式。
- 441 个 Arb DCT 逆关系 entry 包含对应单位阵值，半径小于 (2^{-300})。
- 四代表 patch 在非插值 rational 点比较完整 Bernstein 多项式和直接 H evaluator，差值上界均远小于统一 E。
- 精确核对所有 136 patch index，无遗漏或重复；全部根证书高于 13/2500。

这仍是作者自检，不替代独立审核。已经把 checker、结果和关键证明接口交原 independent_audit。

最小完整重算：

```text
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cc11_cheb_check.py --stage all --fresh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cc11_cheb_selfcheck.py
```

在本目录执行；`--fresh` 强制重算全部 59976 节点，避免信任外部缓存。节点缓存位于仓库外 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/cc11_cheb_cache/`，保存精确 midpoint/radius dyadic 四元组，且配置/输入 hash 全检查。缓存不是独审必需材料。锁证据为脚本、报告、完整 patch/leaf 结果与 cc11_cheb_manifest.json。所有新文件名均为 cc11_cheb_*，首批已锁材料未变。

后续尚需独审验收，以及主协调的 J 上界和最终成本/迁移接口。本批没有自行宣称 (U_3<8)。
