# sqrt026：真实角点正矩筛选与 Gauss 上界

本轮完成用户指定的一个 seed、32768 样本的有限正矩筛选，并独立核验锁定 016 矩后计算了真正的 \(\sqrt d\) 两节点 Gauss 上界。

- **非严格诊断**：真实角点 \(U\approx0.0004409408677683\)，样本标准误约 \(1.4843\times10^{-5}\)；所需 \(U_{\rm crit}\approx0.0005242970629252\)。诊断低约 16%，约 5.62 个样本标准误，固定 \(u^3\) 形状没有显示竞争力。没有将这句话当严格排除。
- **严格上界**：锁定 \(\mu_0,\ldots,\mu_3\) 给 \(U\le0.0006987893728633\ldots\)，仍高于阈值；相应全 \(a\ge0\) 的角点 Young 值上界为 \(0.0003598350082697\ldots\)，不能严格排除当前目标。
- **理论核验**：四子集 Cauchy–Binet 后用平方根凹性，方向确实给正下界；它可以化为两维正 Mellin–Frullani 积分。完整式见下，不假定 \(T\) 与归一形状独立。由于首批真实 \(U\) 诊断未通过，不启动该积分 pilot 或加密。

没有新完整 \(g\)，没有新 \(C_5\) 突破；025 锁定文件未修改。本轮只写 `astra_sqrt026/`，无 MATLAB、无空间网格、未提交推送。

## 1. 模型与角点权重

固定 \(\alpha=1\)、\(\beta=1/2\)，16 条未单位化方向
\(v_i=(1,\pm1,\pm1,\pm1,\pm1)\)，后四个坐标按 \(-1\) 先的字典序编号。角点 \(x=(1,1,1,1,1)\)。独立正 stable \(S_i\) 满足

\[
\mathbb Ee^{-tS_i}=e^{-\sqrt t},\qquad
A=\sum_iS_iv_iv_i^T,\quad T=\sum_iS_i,\quad M=A/T,\quad d=\det M.
\]

令 \(W_c=\sum_{a<b}M_{ab}x_ax_b\)（下文用 \(W_0=16^5\) 表示 frame 常数，避免混淆）。因为 \(M_{aa}=1\)，

\[
u=\frac{W_c+2}{12}
=\sum_i r_i\frac{S_i}{T},\qquad
r_i=\frac{(v_i\cdot x)^2-1}{24}.
\]

正好有一条方向 \(r_i=1\)、五条 \(r_i=1/3\)、十条 \(r_i=0\)，故 \(0\le u\le1\)。记 \(U_{\rm raw}=\sum r_iS_i\)。

shape 概率为

\[
\frac{dR}{d\mathbb P}=\frac{T^{-1/2}}{Z_{\rm tilt}},\qquad
Z_{\rm tilt}=\mathbb ET^{-1/2}=\frac1{16\Gamma(3/2)}.
\]

所需真实正矩及准确另一矩为

\[
U=\mathbb E_R[\sqrt d\,u^6]
=Z_{\rm tilt}^{-1}\mathbb E[\sqrt{\det A}\,U_{\rm raw}^{6}T^{-9}],
\qquad
V=\mathbb E_R[Bu^9]=\frac{965344619}{31225862828160}.
\tag{1}
\]

\(T^{-9}\) 的幂来自 \(-1/2-5/2-6=-9\)，没有漏掉原 shape 倾斜。

Young 形式为 \(H=d^{3/2}B^{-2}\ge3\sqrt d\,Z^2-2BZ^3\)。取 \(Z=au^3\)、\(a\ge0\)，角点的表达式是 \(c_1(3a^2U-2a^3V)\)，其中 \(c_1=1/(32\pi^3)\)。其最大值及必要阈值为

\[
g_{\rm corner}^{\rm best}=c_1\frac{U^3}{V^2},\qquad
U_{\rm crit}=\left(\frac{g_{\rm needed}V^2}{c_1}\right)^{1/3},
\quad g_{\rm needed}=\frac{15198363}{10^{11}}.
\]

这里的门槛有理数取自父任务给定的 \(0.00015198363\)，不是声称重建完整 \(I\) 预算中的未舍入门槛。本轮没有发现靠最后几位舍入改变判断的情形。

## 2. 一个 seed 的真实正矩诊断

登记后按 \(Y_i\sim\mathrm{Gamma}(1/2,1)\)、\(S_i=1/(4Y_i)\) 采样；对应密度

\[
f_S(s)=\frac1{2\sqrt\pi}s^{-3/2}e^{-1/(4s)},\quad s>0,
\]

其 Laplace 积分为 \(e^{-\sqrt t}\)。seed 为 26026，32768 个 IID 样本，8 个 4096 样本 chunk，每个 chunk 保存 RNG 和累加状态，可在 120 秒处恢复。没有增加第二 seed 或追加样本。

样本量是

\[
X=\frac{T^{-1/2}}{Z_{\rm tilt}}\sqrt{\det M}\,u^6,
\qquad \widehat U=\frac1N\sum X_j.
\]

使用准确的理论 \(Z_{\rm tilt}\)，没有以样本权重总和作自归一化。计算全部特征值、条件数和行列式；预先设置条件数大于 \(10^{12}\) 或非正/非有限浮点行列式时，以原始 \(S\) 重建 60 位矩阵复核，不裁剪特征值或行列式。本次无需触发回退。

实际运行结果：

|量|实际值|
|---|---:|
|\(\widehat U\)|0.0004409408677683214|
|样本标准误（非严格区间）|0.0000148425562192457|
|\(U_{\rm crit}\)|0.0005242970629252041…|
|\(\widehat U/U_{\rm crit}\)|0.8410134234|
|用样本 \(U\) 换算的最佳角点值（不是界）|0.0000904077364532705|
|tilt 权重样本均值|1.000543275309977|
|tilt 有效样本数|20906.8036729|
|最小特征值|1.8924012412e−9|
|最大条件数|2.6421457885e9|
|行列式范围|[3.6266059175e−33, 0.9200650714]|
|非正/非有限行列式、裁剪数、回退数|全部 0|

条件数的 50%、90%、99%、99.9% 分位约为 59.3、2055.5、193242.8、17446378.4。所有明细在 `mc_result.json`。\(X\le T^{-1/2}/Z_{\rm tilt}\)，右侧具有全部正矩；这不是带极点的点值 Green 采样。尽管如此，单 seed 样本标准误和诊断差距都不能作为严格误差认证。

真实 \(U\) 诊断已低于门槛约 16%，因此按父任务指示，不再为必然低于真实 \(U\) 的四子集正下界启动二维积分。

## 3. 正确测度的两节点 Gauss 上界

首先逐项核验 `astra_upper016/manifest.json` 的全部 11 个白名单 hash，全部匹配。016 锁状态为 `locked_two_necessary_screens_and_selected_m3_cubic_minorant_ceiling`；本轮只读这组锁定来源，不运行其写回脚本、不使用未锁矩。

本次取正测度

\[
d\mu(t)=u^6\,dR\quad\text{经 }t=d\text{ 推前},\qquad
\mu_j=\mathbb E_R[d^ju^6],\quad j=0,1,2,3.
\]

注意它与 016 的 \(d\,u^6dR\) 测度不同。锁定原始 \(d^jW_c^n\) 矩经精确有理二项式展开给

\[
\begin{aligned}
\mu_0&=\frac{112733}{14370048},&
\mu_1&=\frac{3231145}{33941155248},\\
\mu_2&=\frac{896669090318}{69083562226503825},&
\mu_3&=\frac{256667069791797248}{74931831186788506305375}.
\end{aligned}
\]

\(d\in[0,1]\)：\(M\) 正半定且所有对角元为 1，Hadamard 不等式给 \(\det M\le1\)。令 \(\Delta=\mu_0\mu_2-\mu_1^2>0\)，Gauss 节点为

\[
z_{1,2}=\frac{a\mp\sqrt{a^2-4b_0}}2,\quad
a=\frac{\mu_0\mu_3-\mu_1\mu_2}{\Delta},\quad
b_0=\frac{\mu_1\mu_3-\mu_2^2}{\Delta}.
\]

权为 \(\lambda_1=(\mu_0z_2-\mu_1)/(z_2-z_1)\)、\(\lambda_2=(\mu_1-\mu_0z_1)/(z_2-z_1)\)。256-bit Arb 实际核验 \(0<z_1<z_2<1\)、两权正，且对 0–3 次矩均重现其准确有理值的包围。

对 \(f(t)=\sqrt t\)，\(f^{(4)}(t)=-15/(16t^{7/2})<0\)。两节点双 Hermite 三次多项式 \(P\) 满足

\[
f(t)-P(t)=\frac{f^{(4)}(\xi)}{4!}(t-z_1)^2(t-z_2)^2\le0
\]

对 \(t>0\) 成立，\(t=0\) 由连续性延拓。因此 \(P\) 是全 \([0,1]\) 上界，Gauss 对三次多项式准确，得到真实矩的严格上界

\[
\boxed{U\le\lambda_1\sqrt{z_1}+\lambda_2\sqrt{z_2}.}
\tag{2}
\]

节点约为 \((0.00629157742768,0.27007548111674)\)，权约为 \((0.00767121565345,0.000173782498257)\)。实际完整 Arb 球给

\[
U\le0.0006987893728632779408994991678327874307930875\ldots,
\]
\[
c_1\frac{U^3}{V^2}\le
0.0003598350082697121377665386634275804668438460\ldots.
\]

完整球、准确 Gauss 二次多项式和比较结果在 `gauss_upper.json`；`strictly_excludes_required_g=false`。这确实是严格上界，但**不能**据它排除当前 \(u^3\) 形状。诊断低于阈值和严格上界未闭合，分别保留。

## 4. 四子集平方根正下界及两维正积分核验

以下为理论核验，未运行该积分的数值 pilot。对秩 4 的四子集 \(J\)，设

\[
a_{Ji}=\det(v_j:j\in J;v_i)^2,\quad
L_J=\sum_i a_{Ji}S_i,\quad C_J=\sum_i a_{Ji},\quad
W_0=16^5.
\]

\(a_{Ji}=0\) 对所有 \(i\in J\)，且 \(C_J=16\det(V_J^TV_J)>0\)。每五方向基在四子集展开中计五次，所以

\[
\frac{\det A}{W_0}=\sum_J\theta_JF_J,
\quad \theta_J=\frac{C_J}{5W_0},\quad
F_J=\left(\prod_{j\in J}S_j\right)\frac{L_J}{C_J},\quad
\sum_J\theta_J=1.
\]

秩不足 4 的子集贡献为零，可直接省略。由平方根凹性，正确方向是

\[
\boxed{\sqrt{\det A}\ge\sqrt{W_0}\sum_J\theta_J
\left(\prod_{j\in J}\sqrt{S_j}\right)\sqrt{L_J/C_J}.}
\tag{3}
\]

不能把 (3) 换成 AM–GM 导致的错误方向。选中的四个 \(S_j\) 与 \(L_J\) 独立，因为其系数确实为零；但 \(U_{\rm raw}\)、\(T\) 与这些量不独立，必须保留在后续共同期望里。

令 \(p_{Ji}=a_{Ji}/C_J\)。使用正公式

\[
T^{-9}=\frac1{\Gamma(9)}\int_0^\infty t^8e^{-tT}\,dt,
\quad
\sqrt{\ell}=\frac1{2\sqrt\pi}\int_0^\infty(1-e^{-s\ell})s^{-3/2}\,ds.
\]

将 (3) 代入 (1)，Tonelli 给

\[
U\ge\frac{32}{315}\sum_J\theta_J
\int_0^\infty\!\int_0^\infty
t^8s^{-3/2}\,mathcal K_J(t,s)\,ds\,dt,
\tag{4}
\]

其中常数精确化简为
\(Z_{\rm tilt}^{-1}\sqrt{W_0}/(2\sqrt\pi\Gamma(9))=32/315\)。核为

\[
\mathcal K_J(t,s)=
\mathbb E\left[\left(\prod_{j\in J}S_j^{1/2}\right)
U_{\rm raw}^6e^{-tT}\left(1-e^{-s\sum p_{Ji}S_i}\right)\right]\ge0.
\]

展开 \(U_{\rm raw}^6\) 时仅六个正权方向参与。对 \(\sum n_i=6\) 的正系数多项式项，定义

\[
M_a(q)=\mathbb E[S^ae^{-qS}]
=\frac1{\sqrt\pi}(4q)^{-(a-1/2)/2}K_{a-1/2}(\sqrt q),\qquad q>0.
\tag{5}
\]

特别地 \(M_0(q)=e^{-\sqrt q}\)，\(M_{1/2}(q)=K_0(\sqrt q)/\sqrt\pi\)，一般选中方向
\(M_{n+1/2}(q)=(2\sqrt q)^{-n}K_n(\sqrt q)/\sqrt\pi\)。公式 (5) 直接来自上述正 stable 密度及 [NIST DLMF 10.32.10 的 Bessel 正积分](https://dlmf.nist.gov/10.32.E10)，所用参数 \(q>0\) 满足该表示的条件。

于是

\[
\begin{aligned}
\mathcal K_J(t,s)
={}&\sum_{\sum n_i=6}\frac{6!}{\prod_i n_i!}\prod_i r_i^{n_i}
\left(\prod_{j\in J}M_{n_j+1/2}(t)\right)\\
&\qquad\times\left[
\prod_{i\notin J}M_{n_i}(t)
-\prod_{i\notin J}M_{n_i}(t+sp_{Ji})\right].
\end{aligned}
\tag{6}
\]

这里每个括号都非负：\(M_n(q)\) 是非负函数的 Laplace 变换，随 \(q\) 递减，且所有 \(p_{Ji}\ge0\)。因此不是把任意 Laplace 差当作正数；逐个单项式都有正性依据。小 \(s\) 的数值实现仍需避免相减损失，认证时要用正积分或可靠差值包围。

聚合四子集时必须同时保留所选四方向的 \(r_i\) 类别和剩余方向的联合 \((r_i,p_{Ji})\) profile；仅按四维法向量类型聚合，会丢失与角点 \(U_{\rm raw}\) 的相关性。本轮未枚举、聚合或积分该核。

积分有限性由 (3) 的每项非负且总和不超过原 \(\sqrt{\det A}\)，再结合 \(0\le\sqrt d\,u^6\le1\) 得到。实际无限域认证仍需要显式尾界；本轮没有声称完成它们。

## 5. 复核、停止条件与交付

Python 脚本 `corner_screen.py` 实际完成一次登记批；MC 和 Gauss 阶段合计约 0.083 秒，PID 记于 `run_manifest.json`，退出码 0，无遗留句柄。程序不写任何 016 或 025 文件。025 的全部 manifest hash 在交付前再次逐项确认未变。

复现：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B corner_screen.py
```

证据边界：真实 \(U\) 的单批诊断没有竞争力，所以停止二维下界积分；这属于预登记筛选决策。现有四矩的严格上界还不能排除固定 \(u^3\)，故不能把“停止加密”写成该形状的数学不可能性。父任务若改形状、改 \(\alpha\)、获得更高精度正矩证据或更强解析上界，应另开新批，不重写本轮结果。

白名单：报告、登记、最终脚本、MC 结果、Gauss 结果、运行来源 manifest、最终 manifest。完成的 checkpoint 数组和 RNG 状态不归档；没有全量样本、额外随机种子、空间产物或重复源文件。未提交、未推送。
