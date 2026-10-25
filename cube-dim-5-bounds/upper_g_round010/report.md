# 第010批：alpha=3/2 全域正 Gaussian 混合下界

本批取得新的、覆盖整个闭立方体的严格下界

\[
 \boxed{G(x)\ge \frac1{65000}\quad(x\in[-1,1]^5).}
\]

固定16条未归一化 body 方向、权1共同stable时间，Fourier归一化与 `alternative_alpha/transfer.md` 相同。该结论只是一条较弱的新全域基准：没有达到条件目标 `g=1/10000`，没有认证面中段积分，也没有改进当前 `C5<=32` 上界。原点或其他奇异处若Green核取无穷大，以下不等式仍成立。

方法与已排除的 cosine 尾积分不同：本批从**正**的条件 Gaussian 密度出发，通过指数倾斜、熵形式 Jensen 和 tight-frame 矩阵界获得每个时间的全域下界，最后对所有正时间积分。没有丢弃短时段，也不涉及有限点值推广。

## 1. 正混合及矩阵归一化

令 `beta=3/4`，各 `S_v(t)>0` 独立且

\[
 \mathbb E e^{-\lambda S_v(t)}=e^{-t\lambda^\beta}.
\]

正stable subordinator 的存在和此 Laplace 归一化可由 Lévy 测度 `beta/Gamma(1-beta) s^(-1-beta)ds` 构造；参见 [Lalley 讲义 Proposition 5.1、Corollary 5.2，PDF第11–12页](https://galton.uchicago.edu/~lalley/Courses/385/Old/LevyProcesses.pdf)。本处只使用该存在性与 Laplace 变换；下面的混合公式与所有不等式均直接推导。

给定 `S_v`，取独立中心 Gaussian，其各方向增量方差为 `2*S_v`。其无条件方向特征函数是 `exp(-t*|xi|^(2*beta))`，正好给所需 `alpha=3/2`。令

\[
 V=[v]_{v\in\mathcal V},\qquad
 A=\sum_v S_v vv^T=V\operatorname{diag}(S_v)V^T.
\]

条件密度及无条件密度为

\[
 F_x(A)=(4\pi)^{-5/2}(\det A)^{-1/2}
      e^{-x^TA^{-1}x/4},\qquad p_t(x)=\mathbb E_P F_x(A).
\]

16方向有限精确枚举验证 `VV^T=16I`。对任何 `x`，加权最小二乘恒等式和可行选择 `z=V^Tx/16` 给

\[
 x^TA^{-1}x
 =\min_{Vz=x}\sum_v\frac{z_v^2}{S_v}
 \le\frac1{256}\sum_v\frac{(v\cdot x)^2}{S_v}. \tag{1}
\]

这是矩阵逆的上界，方向与后续下界所需一致。该不等式可逐个正 `S_v` 使用，没有概率高置信事件或例外方向。

## 2. 固定参数指数倾斜

对每个时间 `t>0`，固定 `u=1/10`、`lambda=u*t^(-4/3)`，定义乘积概率测度

\[
 \frac{dQ}{dP}=\exp\{-\lambda\sum_vS_v+16t\lambda^\beta\}.
\]

由 Laplace 变换求导及 Tonelli 的 `1/S=int_0^infty exp(-rS)dr`，令

\[
 m=\beta u^{\beta-1},\qquad
 n=\frac{e^{u^\beta}}\beta\Gamma(1/\beta,u^\beta),
 \qquad D=16(1-\beta)u^\beta,
\]

则

\[
 \mathbb E_QS_v=m t^{4/3},\qquad
 \mathbb E_QS_v^{-1}=n t^{-4/3},\qquad
 D(Q\Vert P)=D. \tag{2}
\]

例如逆矩公式的中间式为

\[
 \mathbb E_QS_v^{-1}
 =e^{t\lambda^\beta}\int_\lambda^\infty e^{-tr^\beta}dr,
\]

换元后恰得 (2) 的上不完全Gamma式。这里使用有限一阶矩只限于倾斜测度；没有误用原stable变量的无穷均值。

## 3. 熵 Jensen 给每个时间的全域下界

`log det` 凹性和 (2) 给

\[
 \mathbb E_Q\log\det A\le\log\det\mathbb E_QA
 =5\log(16m t^{4/3}). \tag{3}
\]

由 (1)、(2) 及 `sum_v(v dot x)^2=16|x|^2`，

\[
 \mathbb E_Q x^TA^{-1}x
 \le \frac n{16}t^{-4/3}|x|^2
 \le\frac{5n}{16}t^{-4/3}\quad(x\in[-1,1]^5). \tag{4}
\]

再对 `F_x*dP/dQ` 使用 Jensen：

\[
 \log p_t(x)\ge\mathbb E_Q\log F_x-D(Q\Vert P).
\]

其可积性可由 (2) 检查：`log S` 的正负部分分别受 `S`、`1/S` 控制；任选五个独立方向的 Cauchy–Binet 项给 `det A` 的正下界，故 `log det A` 可积，(4) 控制逆二次型。因此该 Jensen 使用没有未检查的无穷差。

合并 (3)–(4)，定义

\[
 C=(4\pi)^{-5/2}(16m)^{-5/2}e^{-D},\qquad a=5n/64,
\]

得到

\[
 p_t(x)\ge C t^{-10/3}e^{-a t^{-4/3}}
 \quad\text{对每个 }t>0,\ x\in[-1,1]^5. \tag{5}
\]

这个证明唯一的空间放宽是 `|x|^2<=5`，所以同时覆盖所有面、棱、角、方向奇异集合及内部。无需排除任何边界邻域。若以a.e.密度代表元定义 `p_t`，正Gaussian混合本身提供与该分布一致的逐点代表；essential下界当然也成立。

## 4. 完整时间积分和实际计算

直接换元给

\[
 \int_0^\infty t^{-10/3}e^{-a t^{-4/3}}dt
 =\frac34\Gamma(7/4)a^{-7/4}.
\]

于是 (5) 对所有时间积分即得

\[
 G(x)\ge C\frac34\Gamma(7/4)a^{-7/4}.
\]

现有本机 Python 环境的256位 Arb 实际运行给常数包围

\[
 C=[4.16968695366106699065726719765273\ldots\times10^{-7}],
 \quad a=[0.102680957659234242602334203657114\ldots],
\]

\[
 C\frac34\Gamma(7/4)a^{-7/4}
 =[1.54312969885831833948392973526252\ldots\times10^{-5}
 \ \mathbin{+/-}\ 5.89\times10^{-80}]
 >\frac1{65000}.
\]

Arb比较同时验证该公式值小于 `1/10000`；因此已知的固定参数解析包围本身不足以达到条件目标。不能由此推断真实 `inf G` 小于目标。

复现命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/upper_g_round010/tilted_gaussian_bound.py
```

实际成功运行约0.0066秒，未调用MATLAB。代码SHA-256：`f62764cf49abcc10f66efecc926c51e52be1680bd9b67dd823468e781310550e`；完整Arb球及常数见 `tilted_gaussian_bound.json`。

## 5. 有限候选账本与停止范围

候选1是原测度下的log-Jensen基准：`log det A<=5 log(sum S_v)` 与 (1) 可得浮点诊断约 `4.70705e-6`，本批没有将该数值另外整理为严格证书。候选2是上述固定倾斜参数方案，已完成解析证明与本机区间认证。登记中的候选3（方向依赖倾斜或正积分覆盖）未执行，后续突破由另一个研究者接续，避免重复探索。

本批至此停止。已完成的是一个全域 `g>=1/65000` 基准；尚未完成 `g>=1/10000`、面中段 `I` 认证及 `C5<32`，也没有声称满足用户要求的新上下界改善门槛。
