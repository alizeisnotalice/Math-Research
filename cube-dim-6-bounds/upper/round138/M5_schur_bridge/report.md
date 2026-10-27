# C138：五维具体矩证书经Schur--Holder桥的预算资格

## 结论

Schur桥及同law来源核对均成立，但当前具体五维矩上包配此Holder桥的上界预算，严格劣于已保存round128完整矩上包，也不能达到
\(g_0=1/65000\) 或 \(g_*=1/63500\) 的对应门槛。唯一新三比较Fraction批及独立准确算术复核均通过；没有运行Bessel/Gamma、积分、优化或旧基线。

这里的严格下界属于一个指定的**上包预算表达式**，不属于实际
\(M_6\)、最优可达矩预算或 \(C_6\)。更强的具体 \(M_5\) 证书或利用相关性的方法仍未被排除；没有新的global g或 \(C_6\) 上界。

## 1. 具体概率模型与来源

令
\(V_5=\{(1,\pm1,\pm1,\pm1,\pm1)\}\)，
\(V_6=\{(u,\pm1):u\in V_5\}\)，独立
\(S_v=1/(4Y_v)\)，\(Y_v\sim\Gamma(1/2,\text{rate }1)\)。其Laplace变换为 \(e^{-\sqrt\lambda}\)。定义
\[
\Sigma_d=\sum_{v\in V_d}S_vvv^T,
\qquad M_d=\mathbb E(\det\Sigma_d)^{-1/(2d)}.
\]
五维只读源 `astra_upper014/report.md` 第1节正是此16方向、此独立stable归一化和此底层矩阵 \(A\)；源中的条件Gaussian协方差
\(2t^2A\) 没有被误当成 \(\Sigma_5\)。其第3节式(4)证明
\[
J_{1/2}=16\Gamma(3/2)M_5=8\sqrt\pi M_5.
\]
`astra_upper014/budget.json` 的alpha1/beta1/2/b1/10条目保存准确向上有理端
\[
J_U=2044690556653767/10^{15},\qquad M_5\le J_U/(8\sqrt\pi).
\]
`upper_alpha1_029/check_results.json` 独立包围同一表达式，并记录包含于锁定来源球。独审还核查了原writer的有向Arb有理端比较与hash资格。输入是具体determinant moment证书，并非从五维最大算子常数反推；没有读取或修改五维live计算或状态。

## 2. Schur与Holder推导

将 \(\Sigma_6\) 写成块矩阵 \(\begin{psmallmatrix}A&b\\b^T&c\end{psmallmatrix}\)。按末坐标正负配对，
\[
A=\sum_{u\in V_5}(S_{u+}+S_{u-})uu^T
\overset d=4\Sigma_5.
\]
对任意 \(z\in\mathbb R^5\)，逐项完成平方给
\[
z^TAz-2b^Tz+c
=\sum_u[S_{u+}(u\cdot z-1)^2+S_{u-}(u\cdot z+1)^2]
\ge\sum_u\frac{4S_{u+}S_{u-}}{S_{u+}+S_{u-}}=H.
\]
其中 \(H=\sum_{j=1}^{16}E_j^{-1}\)，独立
\(E_j=Y_{j+}+Y_{j-}\sim\mathrm{Exp}(1)\)。最小化左侧得到
\(c-b^TA^{-1}b\ge H\)，故 \(\det\Sigma_6\ge\det A\,H\)。Holder指数6/5与6给
\[
M_6\le2^{-5/6}M_5^{5/6}B^{1/6},
\qquad B=\mathbb EH^{-1/2}.
\]
此步骤不假设 \(A,H\) 独立。具体 \(M_5\) 上包有限，而
\(H^{-1/2}\le\sqrt{\sum_jE_j}/16\) 保证 \(B<\infty\)。

可选标量恒等式为
\[
B=\pi^{-1/2}\int_0^\infty[uK_1(u)]^{16}\,du.
\]
它由独立Laplace积分及Mellin公式直接推出；换元
\(u=2\sqrt t\) 满足 \(t^{-1/2}dt=du\)，无额外因子2。该积分及附带尾界本轮没有数值执行，因为下面的廉价严格诊断已分离预算。

## 3. 唯一新准确三比较

令 \(E_{\min}=\min_jE_j\sim\mathrm{Exp}(16)\)。点态
\(H\le16/E_{\min}\)，于是
\(B\ge\mathbb E\sqrt{E_{\min}}/4=\sqrt\pi/32\)。指定输入
\(M_{5,U}=J_U/(8\sqrt\pi)\) 产生的桥预算六次方为
\[
Q_U=2^{-5}M_{5,U}^5B
=\frac{J_U^5B}{2^{20}\pi^{5/2}}
\ge\frac{J_U^5}{2^{25}\pi^2}
>L:=\frac{J_U^5}{2^{25}(22/7)^2}.
\]
标准库Fraction实际证明 \(L>10^{-7}\)，并证明以下三个准确目标均小于 \(10^{-7}\)：

| 比较目标 | 仅展示近似值 | \(L/\text{目标}\)，仅展示 |
|---|---:|---:|
| \(\Theta(1/65000)\) | 5.863551139121912107973e-8 | 1.838984486755113363687 |
| \(\Theta(1/63500)\) | 6.745257062289188138183e-8 | 1.598601726007602999221 |
| round128完整 \(M_U^6\) | 6.711146303358518834960e-8 | 1.606726942719792467237 |

这里 \(\Theta(g)=(8/9)^6\,24576^5g^6\)，第一个目标与round128保存的完整准确threshold相等。准确 \(L\)、源哈希、三项真值及展示值分离存于 `candidate_certificate.json`，没有用展示小数作决定。

作者新批为1进程、nice10、无时间上限，内部wall约0.01455秒。未运行特殊函数、积分、优化、native或旧基线。来源及推导详见 `proof.md` 与 `source_manifest.json`；准确比较入口为 `check_budget.py`。

## 4. 独审与证明状态

独审原件：
`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-cover090-audit/cube-dim-6-bounds/audit/upper138_review/`。
其 `manifest.json` 状态为
`SEALED_PASS_C138_SOURCE_QUALIFIED_SCHUR_CURRENT_BUDGET_REVIEW`，SHA256为
`738510ec37a7827830d554cc8ff7d6140e79ae6cd89364216dfcfe404da42116`。
`proof_and_source_review.md` 接纳同law识别、Schur系数、Holder指数、Mellin常数及当前Jcap来源资格。独立三比较进程exit0、stderr空并已reap，完整wall0.039811666996683925秒、CPU0.032219秒；新增特殊函数、数值积分和旧基线调用为0。

因此本轮严格确认的是当前证书与该桥组合的预算不足。这个结论既不排除更小的真实五维矩，也不排除更好的五维矩上包、相关性估计或其他六维方法。本报告补齐归档说明；未修改锁定证书、代码或来源，也未启动后续科学计算。
