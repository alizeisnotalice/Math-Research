# upper_expect058：径向消去后的全域 Jensen 预算及准确障碍

## 结论与范围

本批只检验一个新机制：先准确消去 determinant 和 stable 总尺度倾斜，
再消去 Gamma 径向，最后对**整个 Dirichlet 单纯形这一个完整概率 cell**
使用五次根 Jensen。得到有限、覆盖全部形状的上界，但证明这个上界表达式
本身无法达到目标。没有增加 MC、固定形状、求积点、图矩或参数网格。

令
\[
 J_{\mathrm L}=\mathbb E_R[d^{-1/10}r_{\mathrm L,8}(M)^{1/5}].
\]
本批推导出一个有限预算 \(B_J\)，满足
\[
 J_{\mathrm L}\le B_J,\qquad
 B_0\le B_J\le\sqrt3 B_0,
\]
其中
\[
 B_0=1.711781284229948446351219887123095300798661169\ldots,
 \quad \sqrt3B_0=2.964892155731772084154696048441153981308965462\ldots.
\]
055 的准确 \(g_0=1/10000\) 门为
\[
 J_{\rm gate}=\frac1{10000}
 \left(\frac{60^4 32^9\pi^7}{3125}\right)^{1/5}
 =1.345336046567759917886552188938232715163990807\ldots.
\]
整数算术证书证明
\[
 \boxed{B_0-J_{\rm gate}>0.366445237662188528464667698184862585634670361.}
\]
所以本轮这个**全域单 cell Jensen 预算**必然不过门。
这不是 \(J_{\mathrm L}\ge B_0\)，不是实际 \(I\) 或 \(C_5\) 的排除结论，
也不排除条件分区、较精细混合矩或直接认证形状期望。
新粗上包围比已锁定的 \(J_{\rm flag}<2.044691\) 更弱，不能替换现有最优预算。

另外得到两个可复用的严格接口：

* 在原倾斜概率 \(R\) 下，直接用于五次根 Jensen 的变量 \(Q_M\) 有
  \(\mathbb E_R Q_M^p=\infty\) 对所有 \(p\ge1/2\)。有限全域切线与仅依赖
  \(Q_M\) 的有理上函数不能避开这个问题。
* 正确消去径向后，对应变量 \(Q_Z\) 满足
  \(0<Q_Z\le2^{-20}\)。所有正矩有限，并有无奇异端点的条件分区公式。

## 1. 固定模型与精确 determinant 消去

沿029、038、049、051的同一对象：16个未单位化方向
\(v=(1,\varepsilon_2,\ldots,\varepsilon_5)\)，独立 positive stable
\(S_v\) 的 Laplace transform 为 \(e^{-\sqrt u}\)。
\[
 A_S=\sum_v S_vvv^T,\quad T=\sum_vS_v,\quad M=A_S/T,
 \quad d=\det M,\quad Z_{\rm tilt}=\mathbb ET^{-1/2}=\frac1{8\sqrt\pi}.
\]
\(R\) 是由 \(T^{-1/2}/Z_{\rm tilt}\) 倾斜后推出的形状概率。
不使用 \(T\) 和形状独立的假设。

记 \(\mathcal A\) 为完整1920元素轨道与051正 Lobatto 权重的归一化平均。
对 \(N=UMU^T\)，有
\[
 M[tM^{-1}+(1-t)N^{-1}]N=(1-t)M+tN.
\]
因此
\[
 [d\det(tM^{-1}+(1-t)UM^{-1}U^T)]^{-1/2}
 =\sqrt d\,[\det((1-t)M+tUMU^T)]^{-1/2}.
\]
定义对任意正定 \(X\)
\[
 Q_X=\mathcal A_{U,t}[\det((1-t)X+tUXU^T)]^{-1/2}.
\]
于是精确逐点恒等式是
\[
 r_{\mathrm L,8}(M)=\sqrt d\,Q_M,
 \qquad d^{-1/10}r_{\mathrm L,8}(M)^{1/5}=Q_M^{1/5}. \tag{1}
\]
同时 \(Q_M=T^{5/2}Q_{A_S}\)，从而 **倾斜先准确取消**：
\[
 \boxed{J_{\mathrm L}=8\sqrt\pi\,\mathbb E[Q_{A_S}^{1/5}].} \tag{2}
\]
非负性允许 Tonelli；这里尚未做任何 Jensen。

## 2. 直接在 \(R\) 下取矩的真实发散障碍

051端点的总质量是 \(a=1024/51051>0\)，两端均给
\(\det M^{-1/2}\)，因此 \(Q_M\ge a d^{-1/2}\)。
写 \(S_i=(4Y_i)^{-1}\)，其中独立 \(Y_i\sim\Gamma(1/2,1)\)。
固定正概率事件
\[
 0<Y_1=y\le1,\qquad 1\le Y_j\le2\quad(2\le j\le16).
\]
在此事件上，\(P_j\le y/Y_j\le y\) 对 \(j\ne1\)，所以
\(M=P_1v_1v_1^T+N\)，且 \(\operatorname{tr}N\le75y\)。
由 rank-one min–max 以及 \(\operatorname{tr}M=5\)，
\[
 d\le5(75y)^4.
\]
又 \(T\le4/y\)，故 \(T^{-1/2}\ge\sqrt y/2\)。
\(Y_1\) 的密度在 \(0<y\le1\) 下界为
\(e^{-1}/\sqrt{\pi y}\)。其余15变量落在 \([1,2]\) 的概率
\(p_*\) 严格为正。因此对每个 \(p>0\)
\[
 \mathbb E_RQ_M^p\ge
 4e^{-1}p_*a^p5^{-p/2}75^{-2p}
 \int_0^1y^{-2p}\,dy.
\]
右端在 \(p\ge1/2\) 时发散，包括临界对数发散。
这是一个正概率邻域，不是只使用概率零的退化矩阵。

特别地，直接 \(\mathbb E_RQ_M^{1/5}\le(\mathbb E_RQ_M)^{1/5}\)
只给无穷。任意有限个正节点的五次根切线的最小值，在无穷处仍至少是
正斜率线性增长，因此其 \(R\) 期望也无穷。
若实有理函数 \(h(z)\) 在 \([1,\infty)\) 无极点且
\(h(z)\ge z^{1/5}\)，其分子减分母次数至少为1，故最终
\(h(z)\ge cz\) 对某 \(c>0\) 成立，也不能产生有限期望上界。
这只排除**单变量 \(Q_M\) 的这类全域 majorant**；不排除带额外形状变量
的有理函数、带非线性尾项的分段估计或先做径向消去。

## 3. 径向 Gamma 消去得到有界形状变量

令 \(R_0=\sum_iY_i\)、\(Z_i=Y_i/R_0\)。049已经准确证明
\(R_0\sim\Gamma(8,1)\)，\(Z\sim\operatorname{Dir}(1/2^{16})\)，二者独立。
定义
\[
 A_Z=\sum_i Z_i^{-1}v_iv_i^T,
 \qquad A_S=(4R_0)^{-1}A_Z.
\]
于是 (2) 给
\[
 \boxed{J_{\mathrm L}=C\,\mathbb E_Z Q_{A_Z}^{1/5},\qquad
 C=16\sqrt\pi\frac{\Gamma(17/2)}{\Gamma(8)}=\frac{6435\pi}{256}.} \tag{3}
\]
以下简写 \(Q_Z=Q_{A_Z}\)。原049表达式中的 \(H^{-1/2}\) 权重已与
\(M=A_Z/H\) 的五次根齐次性准确抵消，不能再保留一次或重复归一化。

这个变量在整个开单纯形上有显式界。令
\(B_Z=\sum_iZ_iv_iv_i^T\)，则 \(\operatorname{diag}B_Z=1\)，而
\[
 \begin{pmatrix}A_Z&16I\\16I&B_Z\end{pmatrix}
 =\sum_i
 \begin{pmatrix}Z_i^{-1/2}v_i\\Z_i^{1/2}v_i\end{pmatrix}
 \begin{pmatrix}Z_i^{-1/2}v_i\\Z_i^{1/2}v_i\end{pmatrix}^{T}\succeq0.
\]
Schur补给 \(A_Z\succeq256B_Z^{-1}\)。Hadamard给 \(\det B_Z\le1\)，
故 \(\det A_Z\ge256^5\)。由 log-determinant 凹性，对所有 \(U,t\)
\[
 \det((1-t)A_Z+tUA_ZU^T)\ge\det A_Z\ge256^5.
\]
所以
\[
 \boxed{0<Q_Z\le256^{-5/2}=2^{-20}.} \tag{4}
\]
这证明所有正矩有限，靠近任意 Dirichlet 边界也不产生积分尾部发散。
没有把这一有界性误写成原 \(Q_M\) 的有界性。

## 4. 唯一登记机制：完整概率单 cell Jensen

只在 (3) 后使用 Jensen：
\[
 J_{\mathrm L}\le B_J:=C(\mathbb E_ZQ_Z)^{1/5}. \tag{5}
\]
设 \(L=\mathbb E Q_{A_S}\)。齐次性与径向独立给
\[
 L=32\frac{\Gamma(21/2)}{\Gamma(8)}\mathbb E_ZQ_Z,
 \quad
 B_J=C\left[\frac{\Gamma(8)}{32\Gamma(21/2)}L\right]^{1/5}. \tag{6}
\]
这里的 \(L\) 有完整正积分表示。对每个 \(U,t\) 定义
\[
 \Psi_{U,t}(x)=\sum_v
 \sqrt{(1-t)(v\cdot x)^2+t(Uv\cdot x)^2}.
\]
Gaussian determinant 积分、独立 stable Laplace transform 与 Tonelli 给
\[
 L=\pi^{-5/2}\mathcal A_{U,t}
 \int_{\mathbb R^5}e^{-\Psi_{U,t}(x)}dx. \tag{7}
\]
由于 \(\sum_vvv^T=16I\)，Cauchy–Schwarz 给
\(\Psi_{U,t}(x)\le16|x|\)。
另一方面 \(\sqrt{(1-t)a^2+tb^2}\ge(1-t)|a|+t|b|\)。
对五个独立符号的 Rademacher 和 \(X=\sum_i\varepsilon_ix_i\)，
\[
 \mathbb EX^2=|x|^2,\quad
 \mathbb EX^4=3|x|^4-2\sum_i x_i^4\le3|x|^4.
\]
Holder给
\(\mathbb E|X|\ge(\mathbb EX^2)^{3/2}/(\mathbb EX^4)^{1/2}
\ge|x|/\sqrt3\)。16方向绝对值平均与全部32符号平均一致。
同时对 \(x\) 与 \(U^Tx\) 使用此式，得到
\[
 \frac{16}{\sqrt3}|x|\le\Psi_{U,t}(x)\le16|x|. \tag{8}
\]
这个粗夹逼无需枚举群或求积点，对每项均成立，因而完整正平均也成立。
在五维中
\(\int e^{-a|x|}dx=|S^4|\Gamma(5)/a^5\)，且
\(|S^4|=8\pi^2/3\)。于是
\[
 L_{\rm iso}:=\frac1{16384\sqrt\pi}\le L
 \le3^{5/2}L_{\rm iso}. \tag{9}
\]
将 (9) 代入 (6)，定义
\[
 B_0=C\left[\frac{\Gamma(8)}{32\Gamma(21/2)}L_{\rm iso}\right]^{1/5},
\]
便得 \(B_0\le B_J\le\sqrt3B_0\)。这是一项覆盖完整概率空间的有限上界，
也给出这个机制无法达到验收门的严格原因。
下障碍来自 (8) 的上边，故即使改善粗 Khinchine 下边、准确计算 (7)，
或把 Lobatto 换成准确 Beta 积分，**全域单 cell Jensen** 仍不可能过门。

## 5. 准确有限证书与剩余接口

半整数 Gamma 全部约为有理数乘 \(\sqrt\pi\)：
\[
 \frac{\Gamma(17/2)}{\sqrt\pi}=\frac{2027025}{256},\quad
 \frac{\Gamma(21/2)}{\sqrt\pi}=\frac{654729075}{1024},\quad\Gamma(8)=5040.
\]
因此
\[
 B_0^5=
 \frac{2(2027025/256)^5}{5040^4(654729075/1024)}\pi^4,
 \qquad
 J_{\rm gate}^5=\frac{60^4 32^9}{3125\,10000^5}\pi^7. \tag{10}
\]
程序只读取主010已锁定的完整 Machin80 有理 \(\pi\) 区间，
以 Fraction/int 对6个5次或10次根做80位十进制有理外包；每个端点都由
准确幂比较核验。另直接核验 (10) 的严格五次幂比较，不依赖展示小数。
独立回读程序再次验证全部6个幂夹逼、Gamma有理系数、门槛差和源hash。

后续若采用完整有限概率分区 \(E_1,\ldots,E_m\)（零概率 cell 按0计），
令 \(p_j=\mathbb P_Z(E_j)\)、\(m_j=\mathbb E_Z[Q_Z1_{E_j}]\)，则条件 Jensen 给
\[
 \boxed{J_{\mathrm L}\le C\sum_{j=1}^m p_j^{4/5}m_j^{1/5}.} \tag{11}
\]
(4) 保证所有 \(m_j\) 有限；(11) 无需处理 \(R\) 下的发散高矩。
**本批只认证 \(m=1\) 的全域夹逼，没有认证任何更细 cell 的 \(p_j,m_j\)。**
要推进该接口，需要共同认证完整分区的概率与截断混合行列式矩，不能只从
050样本估计这些量，也不能漏掉边界 cell。这是剩余具体缺口。

## 6. 真实运行、一次明确修复与锁定范围

首次运行 PID55221、exit1，在第5个根包围转字符串时触发 Python 默认
4300位整数序列化限制；此前没有数学断言失败。原失败receipt保留。
单独登记仅将环境 `PYTHONINTMAXSTRDIGITS=20000` 后，用同一代码、同一6个
输入、同一精度和阈值重跑；未改变数学实例。两个尝试的二分次数总上界仍
低于原登记4096。没有对失败数值裁剪、替换或继续增加新实例。

成功证书 PID55325、exit0，内部0.00393秒、1604次整数二分，外层0.0275秒。
未触及120秒checkpoint。回读运行信息另存 `verify_receipt.json` 和
`verification.json`，以实际文件为准。

```sh
PYTHONINTMAXSTRDIGITS=20000 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B certify.py
PYTHONINTMAXSTRDIGITS=20000 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B verify.py
```

049/050/051、055及主010/029均只读，完整路径与SHA256在登记/results中。
未使用050样本均值作为证书，未重跑051固定形状，未触052、共享MATLAB或
其他工作者文件，未提交/推送。没有新的 \(I\) 或 \(C_5\) 上界；正式 \(C_5\le32\) 不变。
