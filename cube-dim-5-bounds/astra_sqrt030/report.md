# sqrt030：倒数线性形状的闭式负项与两个参数的严格排除

本轮固定 \(\alpha=1,\beta=1/2\)，只检查预登记的 \(\rho\in\{7/8,31/32\}\)，形状为
\(Z=a_0(1-\rho u)^{-3/2}\)、\(a_0\ge0\)。得到新的闭式 Cauchy 上界，已用 256-bit Arb 严格排除这两个形状达到所需角点预算：

|\(\rho\)|真实 \(U\) 的严格上界（显示近似）|达到目标所需 \(U_{\rm crit}\)|所有 \(a_0\ge0\) 的角点 Young 值严格上界（显示近似）|
|---|---:|---:|---:|
|7/8|0.618595049860920|0.900316410087420|4.92981515746951e−5|
|31/32|0.829345772938888|0.975051128551463|9.35236723064249e−5|

两个上界均由完整 Arb 球严格比较为小于较弱的有理阈值
\(15198/10^8=0.00015198\)，从而也小于父任务给定的 \(0.00015198363\)。不是从展示小数制造证书。这里上界作用于 **这两个固定形状的 Young 证书值**，不是真实 Green 函数的上界，也不排除其他 \(\rho\)、其他 \(\lambda\) 或其他机制。

首批一个新 seed、32768 样本的真实正矩诊断也未显示竞争力。随后闭式 Cauchy 认证将筛选提升为严格排除，所以没有启动二维积分 pilot、空间网格或追加随机样本。无新完整 \(g\)，无新 \(C_5\) 突破；025、026 全部锁定文件保持不变。

## 1. 尺度消去与两个正矩

16 条方向 \(v_i=(1,\pm1,\pm1,\pm1,\pm1)\) 不单位化，后四坐标按 \(-1\) 先的字典序排列。固定角点 \(x=(1,1,1,1,1)\)。原变量 \(S_i\) 独立，\(\mathbb Ee^{-tS_i}=e^{-\sqrt t}\)。令

\[
A=\sum_iS_iv_iv_i^T,\quad T=\sum_iS_i,\quad P_i=S_i/T,\quad
M=A/T,\quad d=\det M,\quad b_A=x^T\operatorname{adj}(A)x.
\]

原 shape 概率为 \(dR=Z_{\rm tilt}^{-1}T^{-1/2}d\mathbb P\)，其中
\(Z_{\rm tilt}=1/(16\Gamma(3/2))=1/(8\sqrt\pi)\)。角点系数为

\[
r_i=\frac{(v_i\cdot x)^2-1}{24}\in\{1,1/3,0\},\qquad
u=\sum_i r_iP_i,
\]

三类数量依次为 1、5、10。置

\[
\lambda_i=1-\rho r_i>0,\quad H=\sum_i\lambda_iS_i,\quad
L=\sum_i\lambda_iP_i=H/T=1-\rho u.
\]

对 \(Z=a_0L^{-3/2}\)，Young 两项的矩为

\[
U=\mathbb E_R[\sqrt d\,L^{-3}],\qquad
V=\mathbb E_R[B L^{-9/2}],\quad B=x^T\operatorname{adj}(M)x.
\]

因为 \(\sqrt d=\sqrt{\det A}\,T^{-5/2}\)、\(B=b_AT^{-4}\)，恰有

\[
\boxed{U=Z_{\rm tilt}^{-1}\mathbb E[\sqrt{\det A}\,H^{-3}],\qquad
V=Z_{\rm tilt}^{-1}\mathbb E[b_AH^{-9/2}].}
\tag{1}
\]

第一式的 \(T\) 指数为 \(-1/2-5/2+3=0\)，第二式为 \(-1/2-4+9/2=0\)。没有假定 \(T\) 与 \(P\) 独立。

角点最优 Young 值为

\[
\sup_{a_0\ge0}c_1(3a_0^2U-2a_0^3V)
=c_1\frac{U^3}{V^2},\quad
c_1=\frac1{32\pi^3},\quad a_0^{\rm opt}=U/V.
\tag{2}
\]

本轮门槛按父任务提供的精确十进制有理数
\(g_*=15198363/10^{11}\) 定义，\(U_{\rm crit}=(g_*V^2/c_1)^{1/3}\)。严格排除实际比较更弱阈值 \(15198/10^8\)，不依赖最后几位舍入。

## 2. \(V\) 的完整闭式与角点化简

定义

\[
C=\sum_i\sqrt{\lambda_i},\qquad
D=\sum_i\lambda_i^{-1/2}v_iv_i^T.
\]

Laplace 一阶矩是

\[
\mathbb E[S_ie^{-t\lambda_iS_i}]
=\frac{e^{-\sqrt{t\lambda_i}}}{2\sqrt{t\lambda_i}}.
\]

\(b_A\) 是 \(S\) 中四次、每变量次数至多 1 的 Cauchy–Binet 多项式；也可由 \(\det(A+qxx^T)\) 对 \(q\) 的一次系数取得。因此独立性给

\[
\mathbb E[b_Ae^{-tH}]
=\frac{e^{-C\sqrt t}}{16t^2}x^T\operatorname{adj}(D)x.
\]

用 \(H^{-9/2}\) 的正 Mellin 积分、\(r=\sqrt t\) 换元，

\[
\boxed{V=\frac{128}{35}\frac{x^T\operatorname{adj}(D)x}{C^5}.}
\tag{3}
\]

常数实际由 SymPy 精确核验：
\(Z_{\rm tilt}^{-1}\,2\Gamma(5)/(16\Gamma(9/2))=128/35\)。所有积分为非负，Tonelli 可用。

角点三类方向的外积和分别为

\[
\mathbf1\mathbf1^T,\qquad 4I+\mathbf1\mathbf1^T,\qquad
12I-2\mathbf1\mathbf1^T.
\]

程序以精确整数逐项核对这三个矩阵。令 \(a=1-\rho\)、\(b=1-\rho/3\)，则

\[
C=\sqrt a+5\sqrt b+10,\quad
D_\perp=12+4b^{-1/2},\quad
D_\parallel=2+9b^{-1/2}+5a^{-1/2},
\]
\[
\det D=D_\perp^4D_\parallel,\qquad
x^T\operatorname{adj}(D)x=5D_\perp^4,
\qquad V=\frac{128}{7}\frac{D_\perp^4}{C^5}.
\tag{4}
\]

\(\rho=0\) 时 \(C=16,D=16I\)，(3) 精确给
\(V=(8/35)|x|^2=8/7\)，校准通过。两个所需值经 256-bit Arb 实际计算为

\[
V_{7/8}=2.19985747450644752694029533462356128\ldots,
\]
\[
V_{31/32}=2.47937877597652496722867533702009053\ldots.
\]

完整球在 `exact_V.json`。

## 3. 新闭式 Cauchy 上界及严格排除

本轮新增的关键不是高阶矩，而是对 (1) 在原概率下选取恰好可积分的两个 Cauchy 因子：

\[
\mathbb E[\sqrt{\det A}\,H^{-3}]
=\mathbb E\left[\sqrt{\det A\,H^{-5}}\sqrt{H^{-1}}\right]
\le\sqrt{\mathbb E[\det A\,H^{-5}]\,\mathbb E[H^{-1}]}.
\]

五次多线性行列式的 Laplace 式为

\[
\mathbb E[(\det A)e^{-tH}]
=\frac{\det D}{32t^{5/2}}e^{-C\sqrt t}.
\]

于是正 Mellin 积分直接给

\[
\mathbb E[(\det A)H^{-5}]=\frac{\det D}{16C^5},\qquad
\mathbb E[H^{-1}]=\frac2{C^2}.
\]

由于 \(Z_{\rm tilt}^{-2}=64\pi\)，得到对任意正 \(\lambda_i\) 都成立的新上界

\[
\boxed{U^2\le\frac{8\pi\det D}{C^7}.}
\tag{5}
\]

右端有限，亦验证 Cauchy 两因子的可积性。式 (5) 的尺度正确：同时把所有 \(\lambda_i\) 乘 \(k>0\)，右端乘 \(k^{-6}\)，与 \(U^2\) 一致。

将 (4)、(5) 代入 (2)，新的有限认证批只计算两点，没有参数优化或扫描。256-bit Arb 结果为

\[
U_{7/8}\le0.61859504986091963173225495225059017\ldots
<U_{{\rm crit},7/8},
\]
\[
U_{31/32}\le0.82934577293888825227593912937062036\ldots
<U_{{\rm crit},31/32}.
\]

两族最佳角点 Young 值依次不超过

\[
4.92981515746951149624481300024540066\ldots\times10^{-5},
\quad
9.35236723064249191534440859071175459\ldots\times10^{-5}.
\]

`cauchy_upper.py` 对完整球断言其严格小于 \(15198/10^8<g_*\)，全部通过。证书输入、源码 hash、完整球及 PID 在 `cauchy_upper.json`。这是对两个固定形状、全部 \(a_0\ge0\) 的排除；即使允许 \(a_0\) 随空间自适应，也无法绕开该角点的所有 \(a_0\) 上限。它不声称真实 Green 值小于这些数。

## 4. 已完成的一次正矩数值筛选

新 seed 30030，32768 个 IID 样本，两个 \(\rho\) 共用同一批数据。正 stable 通过 \(Y_i\sim\Gamma(1/2,1)\)、\(S_i=1/(4Y_i)\) 生成，统计量为

\[
\frac{T^{-1/2}}{Z_{\rm tilt}}\frac{\sqrt{\det M}}{(1-\rho u)^3}.
\]

这是有限正矩，不是点值 Green 的奇异统计量。理论 tilt 常数固定，没有样本权重自归一化，也没有裁剪。每个 4096 样本 chunk 保存 RNG 与数组以供 120 秒 checkpoint 恢复。两个结果为

|\(\rho\)|\(U\) 样本均值|样本标准误，非置信区间|均值/所需阈值|换算的 \(g_{\rm opt}\) 诊断，非界|
|---|---:|---:|---:|---:|
|7/8|0.519647254611742|0.00516698448611|0.5771829201|2.92238362844e−5|
|31/32|0.636277271626265|0.00717605742625|0.6525578536|4.22331889619e−5|

tilt 权重均值 0.9900472468，样本标准误 0.0041259848，有效样本数 20883.4571。最小特征值 \(1.63305\times10^{-10}\)，最大条件数 \(3.06175\times10^{10}\)；条件数中位数 61.33、99% 分位约 200312.8。未出现非正浮点行列式，未触发条件数 \(10^{12}\) 的预登记高精度回退，无裁剪。完整记录在 `mc_result.json`。

这些数字独立于后续严格排除，只是初筛证据。没有把一次样本标准误当严格误差界，没有追加 seed 或重采样。

## 5. 四子集正二维积分常数核验

为保留可复用的新正接口，本节核对父任务提出的下界；本轮没有执行其求积。

记 \(W_0=16^5\)。对秩 4 的四子集 \(J\)，令

\[
L_J=\sum_i\det(v_j:j\in J;v_i)^2S_i,\quad
C_J=L_J(\mathbf1),\quad p_{Ji}=\det(v_j:j\in J;v_i)^2/C_J,
\quad \theta_J=\frac{C_J}{5W_0}.
\]

\(p_{Ji}=0\) 对 \(i\in J\)，\(\sum_J\theta_J=1\)。Cauchy–Binet 和平方根凹性给

\[
\sqrt{\det A}\ge\sqrt{W_0}\sum_J\theta_J
\Big(\prod_{j\in J}S_j^{1/2}\Big)\sqrt{\sum_i p_{Ji}S_i}.
\]

使用 \(H^{-3}=\Gamma(3)^{-1}\int t^2e^{-tH}dt\)、平方根的正 Frullani 积分，及
\(\mathbb E[S^{1/2}e^{-qS}]=K_0(\sqrt q)/\sqrt\pi\)，得

\[
\begin{aligned}
U\ge\frac{2048}{\pi^2}\sum_J\theta_J
\int_0^\infty t^2\prod_{j\in J}K_0(\sqrt{t\lambda_j})
\int_0^\infty s^{-3/2}\bigg[&
e^{-\sum_{i\notin J}\sqrt{t\lambda_i}}\\
&-e^{-\sum_{i\notin J}\sqrt{t\lambda_i+sp_{Ji}}}
\bigg]ds\,dt.
\end{aligned}
\tag{6}
\]

因子经 SymPy 精确核为

\[
\frac{Z_{\rm tilt}^{-1}\sqrt{W_0}}{\Gamma(3)\,2\sqrt\pi\,\pi^2}
=\frac{2048}{\pi^2}.
\]

括号对全部 \(s,t>0\) 非负，因为 \(p_{Ji}\ge0\)；各 \(K_0\) 为正。原选中四变量与余下变量确实独立，且 \(p_{Jj}=0\)，所以乘积化合法。所有 \(\lambda_i>0\)，\(L\ge\min\lambda_i\)，原 \(U\) 有限；下界项非负且受原期望控制，Tonelli 及整体有限性无缺口。实际求积认证仍需防小 \(s\) 相减消失和无限尾界，不能只依赖形式正性。

\(K_0\) 的正积分由 [NIST DLMF 10.32.10](https://dlmf.nist.gov/10.32.E10) 经正实参数换元得到，本地 026 已核此表示。本轮没有额外导数或 \(U_{\rm raw}^6\) 展开。若未来应用 (6)，必须按选中方向的 \(r_i\) 类别和余下联合 \((r_i,p_{Ji})\) profile 聚合，不能仅按法向类型忽略角点相关性。

但本轮两形状已由 (5) 严格排除，即使把 (6) 的下界算得任意精确也不可能达目标，故不执行二维 pilot。

## 6. 实际运行、锁和复现

两批均在运行前登记：`registration.json` 固定两 \(\rho\) 和唯一 MC seed；`cauchy_registration.json` 固定新上界、两点、精度和阈值。实际脚本均退出 0，PID 与耗时在 `run_manifest.json`、`cauchy_upper.json`。未调用 MATLAB，无高阶图矩计算、无参数优化、无全域空间网格，无活动句柄。

复现命令（在本目录）：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B reciprocal_linear_screen.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cauchy_upper.py
```

只写独立 clone 的 `cube-dim-5-bounds/astra_sqrt030/`。交付前按原 manifest 核验 025 和 026 白名单 hash 全部不变。完成的 checkpoint 仅供恢复，不归档；无全量样本副本。必要文件是本报告、两个登记、两个脚本、闭式结果、MC 结果、Cauchy 证书、运行来源 manifest 和最终 manifest。未提交、未推送。
