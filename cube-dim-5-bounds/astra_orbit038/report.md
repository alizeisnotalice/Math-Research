# astra_orbit038：先作有限轨道平均的严格 Lq 损失

## 结论

对五维、\(q=5/4\) 的 profile
\[
 H_M(\theta)=d^{-1/2}(\theta^TM^{-1}\theta)^{-2},
 \qquad d=\det M,
\]
这里 \(d\) 是行列式，维数始终是5。已得到：

1. 双 profile 的归一化交叉能量具有准确的一维 Beta–determinant 表示；所有归一化与不等式方向已核对。
2. 一个严格二元平均损失，以及更强的完整有限群 Hölder 能量界。后者保留“先平均轨道、再取范数”，能直接接回029的同一形状概率。
3. 单坐标符号翻转的交叉能量只依赖 \(\kappa=(M^{-1})_{ii}-1\)，并有有理下界；一般 sign bipartition 可用主子式和非负谱参数控制，rank 1/2 的常数明确区分。
4. 一个很粗但已经闭合的解析 \(R\) 矩修正，给原 flag 公式一个严格改进；该修正幅度不足以独自解决当前联合成本。较强的完整轨道加权期望仍待控制。

数值只验证了一个正定有理矩阵、一个固定 sign flip 和一条正项级数。没有运行 \(R\) 分布、QMC、参数网格、完整轨道积分、旧 \(I\) 重算或 \(C_5\) 联合计算；**没有得到新的 \(C_5<32\)**。

## 1. 固定模型与球面归一化

沿用冻结的 `upper_alpha1_029/report.md`，SHA-256
`298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81`。
其模型是16条未单位化 body 方向、权重1、\(\alpha=1\)、独立正 stable 指数 \(b=1/2\)，并且
\[
 G=c_1\mathbb E_R H_M,\quad c_1=(32\pi^3)^{-1},\quad
 M=\sum_vP_vvv^T,\quad P_v=S_v/T,
\]
其中 \(R\) 是原独立 stable law 经 \(T^{-1/2}\) 倾斜后推送的概率。
不假设 \(T\) 与 \(P\) 独立。半群、body、全空间正源、一般 \(L^1\) 与全部正尺度的资格均沿029，不更换算子。

记 \(d\sigma\) 为单位球面 \(S^4\) 的面积测度，\(\omega_4=8\pi^2/3\)。Gaussian 极坐标积分给
\[
 \int_{S^4}(\theta^TC\theta)^{-5/2}\,d\sigma
 =\omega_4(\det C)^{-1/2},\qquad C>0.
\]
因此
\[
 A_M:=\int_{S^4}H_M^q\,d\sigma
 =d^{-5/8}\omega_4d^{1/2}
 =\boxed{\omega_4d^{-1/8}}. \tag{1.1}
\]
\(H_M^{q/2}\in L^2(S^4)\) 且平方范数为 \(A_M\)。所有本节积分对每个正定 \(M\) 都有限。

## 2. 双 profile 交叉能量：准确的一维公式

对正交矩阵 \(U\)，设 \(A=M^{-1}\)，并定义
\[
 B_U(M)=\int_{S^4}(H_MH_{UMU^T})^{q/2}\,d\sigma,
 \qquad \rho_U(M)=B_U(M)/A_M.
\]
两个 profile 的行列式相同。Gamma 积分经 \((s_1,s_2)=(rt,r(1-t))\) 换元给
\[
 a^{-q}b^{-q}=rac1{B(q,q)}
 \int_0^1 t^{q-1}(1-t)^{q-1}[ta+(1-t)b]^{-2q}\,dt.
\]
因为 \(2q=5/2\)，再用(1.1)的球面积分：
\[
 \boxed{
 \rho_U(M)=\mathbb E_{T_0\sim\mathrm{Beta}(q,q)}
 \left[d\det\bigl(T_0A+(1-T_0)UAU^T\bigr)\right]^{-1/2}.}
 \tag{2.1}
\]
这里 Beta 变量 \(T_0\) 与原形状总尺度 \(T\) 无关，只是积分参数。
不能遗漏式中的 \(d\)；\(U=I\) 时被积函数准确等于1。

由 Cauchy–Schwarz，\(0<\rho_U\le1\)。也可由 \(\log\det\) 凹性直接核方向：两个端点的 determinant 均为 \(d^{-1}\)，混合 determinant 不小于该值。
若 \(UMU^T\ne M\)，严格凹性使 \(\rho_U<1\)；等号当且仅当两个 profile 相同。

## 3. 已证的严格二元平均不等式

对 \(1<q<2\) 及 \(u,v\ge0\)，有
\[
 (u+v)^q\le u^q+v^q+(2^q-2)(uv)^{q/2}. \tag{3.1}
\]
一个直接证明如下。对正 \(u,v\) 写 \(u=me^t,v=me^{-t}\)，可取 \(t\ge0\)。
函数 \(q\mapsto\log\sinh(qt)\) 是凹函数；比较端点 \(q=1,2\) 得
\[
 \sinh(qt)\ge2^{q-1}\cosh^{q-1}(t)\sinh(t).
\]
故 \(F_q(t)=\cosh(qt)-1-2^{q-1}(\cosh^q(t)-1)\) 的导数非负，且 \(F_q(0)=0\)，这正是(3.1)。零值情形由连续性处理。

令 \(b_q=1-2^{1-q}=1-2^{-1/4}>0\)。积分(3.1)给
\[
 \boxed{
 \int_{S^4}\left(\frac{H_M+H_{UMU^T}}2\right)^q
 \le A_M\left[1-b_q(1-\rho_U(M))\right].} \tag{3.2}
\]
右端出现明确的交叉能量差。正定 profile 处处为正；若 \(UMU^T\ne M\)，能量差严格为正。
例如范数平方也满足
\[
 \left\|\frac{H_M+H_{UMU^T}}2\right\|_q^2
 \le A_M^{2/q}[1-b_q(1-\rho_U)]^{2/q},\qquad 2/q=8/5.
\]

## 4. 完整有限轨道的 Hölder 能量界

令 \(\Gamma\) 为任意有限正交群，\(P_\Gamma\) 为轨道平均，\(h=H_M^{q/2}\)。
对有限群上的正随机变量 \(X\)，Hölder 插值给
\[
 (\mathbb EX)^q\le(\mathbb EX^{q/2})^{2q-2}(\mathbb EX^q)^{2-q}.
\]
逐点用于轨道值，再在球面上用 Hölder，得到
\[
 \int(P_\Gamma H_M)^q
 \le\left(\int P_\Gamma H_M^q\right)^{2-q}
     \left(\int(P_\Gamma h)^2\right)^{q-1}.
\]
第一积分为 \(A_M\)，第二积分恰为 \(A_Mr_\Gamma(M)\)，其中
\[
 r_\Gamma(M):=\frac1{|\Gamma|}\sum_{U\in\Gamma}\rho_U(M)
 =\frac{\|P_\Gamma h\|_2^2}{\|h\|_2^2}.
\]
因此得到
\[
 \boxed{\int_{S^4}(P_\Gamma H_M)^q\le A_Mr_\Gamma(M)^{q-1}
 =A_Mr_\Gamma(M)^{1/4}.} \tag{4.1}
\]
\(P_\Gamma\) 是 \(L^2\) 正交投影，所以
\[
 1-r_\Gamma(M)=\frac{\|h-P_\Gamma h\|_2^2}{A_M}
 =\frac1{2|\Gamma|A_M}\sum_U\|h-h\circ U^T\|_2^2.
\]
这就是可量化的完整轨道损失；所有交叉项都可由(2.1)准确计算。

现在取五维全部 signed permutations，\(|\Gamma|=3840\)。完整轨道平均不变，且 \((P_\Gamma H_M)^q\) 是 \(-5\) 次齐次函数。球面与 cube 边界的径向 Jacobian 抵消；边界分为160个同面积轨道面片。故对 \(F_+=\{x_1=1,0\le x_j\le1\}\)，
\[
 \boxed{\|P_\Gamma H_M\|_{L^q(F_+)}
 \le\left(\frac{\pi^2}{60}\right)^{4/5}
 d^{-1/10}r_\Gamma(M)^{1/5}.} \tag{4.2}
\]
其范数平方右端为 \((\pi^2/60)^{8/5}d^{-1/5}r_\Gamma^{2/5}\)。

由029的群不变性，\(G/c_1=\mathbb E_R P_\Gamma H_M\)。因此严格得到同一模型的上界
\[
 \boxed{I_1\le c_1^{5/4}\frac{\pi^2}{60}
 \left[\mathbb E_R\bigl(d^{-1/10}r_\Gamma(M)^{1/5}\bigr)\right]^{5/4}.} \tag{4.3}
\]
形状可积性没有新增要求，因为 \(0<r_\Gamma\le1\)，029已证 \(\mathbb E_Rd^{-1/10}<\infty\)。
因此(4.3)涵盖所有近退化形状与029允许的 Green 奇异面，未删去任何尺度、时间或空间管。

若只控制 sign 子群，可先得到 \(r_{\rm sign}\)，再利用完整群投影收缩得到同样的有效上界；完整群的 \(r_\Gamma\) 可能更小。不能假设 determinant 权与轨道能量独立。

## 5. 单坐标 flip：准确函数与有理损失

令 \(S\) 只翻转第 \(i\) 坐标。按该坐标分块
\[
 M=\begin{pmatrix}1&b^T\\b&C\end{pmatrix},\qquad
 s_i=b^TC^{-1}b\in[0,1),\qquad
 \kappa_i=(M^{-1})_{ii}-1=\frac{s_i}{1-s_i}.
\]
精度矩阵的 Schur 补等于 \(1/M_{ii}=1\)，故准确地
\[
 d\det(tM^{-1}+(1-t)SM^{-1}S)
 =1+4t(1-t)\kappa_i. \tag{5.1}
\]
令 \(Y=4T_0(1-T_0)\)。由 Beta 密度换元，
\(Y\sim\mathrm{Beta}(q,1/2)\)，所以 \(\mathbb EY=5/7\)，并且
\[
 \rho_i=\mathbb E(1+\kappa_iY)^{-1/2}. \tag{5.2}
\]
函数 \(z\mapsto1-(1+\kappa_i z)^{-1/2}\) 在 \([0,1]\) 上凹且零点值为0，故
\[
 \boxed{1-\rho_i\ge\frac57(1-\sqrt{1-s_i})\ge\frac5{14}s_i.} \tag{5.3}
\]
\(s_i\) 是 \(M\) 各元素的明确有理函数，非对角行不为零即给严格正损失。
由 \(\lambda_{\max}(C)\le\operatorname{tr}C=4\)，进一步有
\[
 s_i\ge\frac14\sum_{j\ne i}M_{ij}^2.
\]

另令 \(Z=1-Y\sim\mathrm{Beta}(1/2,q)\)，得到适合严格正级数的式子
\[
 \boxed{\rho_i=sqrt{1-s_i}\ {}_2F_1\left(\frac12,\frac12;\frac74;s_i\right)
 =\sqrt{1-s_i}\sum_{n\ge0}\frac{(1/2)_n^2}{(7/4)_n n!}s_i^n.} \tag{5.4}
\]
这里使用普通 \({}_2F_1\) 归一化，见 [DLMF 普通 Gauss 级数](https://dlmf.nist.gov/15.2.E1)；[Euler 积分页面](https://dlmf.nist.gov/15.6.E1) 的粗体 Olver 函数另含 \(1/\Gamma(c)\)，不能混用。
本报告的归一化由 Beta 期望直接推导，因此没有遗漏该因子。

## 6. 一般 sign bipartition 与 rank 常数

令 \(S\) 按非空划分 \(I\cup J=\{1,\ldots,5\}\) 改号，令 \(A=M^{-1}\)。
把精度矩阵的交叉块乘以 \(z=2t-1\)，Schur 分解给非负参数 \(\kappa_1,\ldots,\kappa_k\)，其中 \(k\le\min(|I|,|J|)\le2\)，使
\[
 d\det(tA+(1-t)SAS)=\prod_{\ell=1}^k(1+u\kappa_\ell),
 \qquad u=4t(1-t). \tag{6.1}
\]
具体地，它们是 \(T(I-T)^{-1}\) 的非零特征值，其中
\(T=A_{II}^{-1/2}A_{IJ}A_{JJ}^{-1}A_{JI}A_{II}^{-1/2}\)，谱属于 \([0,1)\)。

两个无需对角化的有理量是
\[
 D_{I|J}=\prod_\ell(1+\kappa_\ell)
 =\frac{\det M_{II}\det M_{JJ}}{d},
 \quad K_{I|J}=\sum_\ell\kappa_\ell
 =\operatorname{tr}(M_{II}A_{II})-|I|.
\]
因此
\[
 \rho_S\le\mathbb E(1+uK_{I|J})^{-1/2}
 \le1-\frac5{14}\frac{K_{I|J}}{1+K_{I|J}}. \tag{6.2}
\]
这是 determinant/非对角块的明确有理上包围。

对强相关区还有主子式界。因为 \(1+u\kappa\ge u(1+\kappa)\)，
\[
 \boxed{\rho_S\le\min(1,C_kD_{I|J}^{-1/2}),\quad
 C_k=\mathbb E u^{-k/2}.} \tag{6.3}
\]
必须保留 rank：
\[
 C_1=\frac{B(3/4,3/4)}{2B(5/4,5/4)},
 \qquad C_2=\frac{B(1/4,1/4)}{4B(5/4,5/4)}=3.
\]
rank 2 不能继续使用 \(C_1\)。不知道准确 rank 时，使用 \(k=\min(|I|,|J|)\) 仍有效，可能较弱。
这些公式均为理论推导，本批没有新增 split 或矩阵数值实例。

## 7. 一个闭合但幅度有限的解析 R 修正

本节是解析计算，与下一节唯一数值实例独立。记
\(w(M)=\sum_{i<j}M_{ij}^2\)。
对五个单坐标 flip 的(3.2)分别使用完整群平均的 \(L^q\) 收缩，再平均五个合法上界，结合(5.3)得
\[
 \int(P_\Gamma H_M)^q\le A_M\left(1-\frac{b_q}{28}w(M)\right).
\]
对 \(x\mapsto(1-x)^{4/5}\) 在0处用凹性切线，故
\[
 \boxed{\|P_\Gamma H_M\|_{L^q(F_+)}
 \le\left(\frac{\pi^2}{60}\right)^{4/5}d^{-1/10}
 \left(1-\frac{b_q}{35}w(M)\right).} \tag{7.1}
\]
\(0\le w\le10\)，所以括号始终为正。

为准确计算 \(R\) 的二次矩，暂写一般 stable 指数 \(0<b<1\)、方向数 \(N\)。
原变量独立且 \(\mathbb E e^{-uS_v}=e^{-u^b}\)，倾斜归一化
\(Z=\mathbb ET^{-b}=1/[N\Gamma(1+b)]\)。Tonelli 与 Laplace 导数给
\[
 \mathbb E_R P_uP_v=\frac{2b}{N^2(1+b)}\quad(u\ne v),
\]
\[
 \mathbb E_R P_u^2=\frac{1-b}{N(1+b)}+\frac{2b}{N^2(1+b)}. \tag{7.2}
\]
例如交叉项由
\[
 \frac1{Z\Gamma(b+2)}\int_0^\infty
 u^{b+1}\,b^2u^{2b-2}e^{-Nu^b}\,du
\]
直接得到；对角项再加 \(b(1-b)u^{b-2}\) 对应项。
无须正 stable 变量具有未倾斜的二阶矩；所有积分均带正 Laplace 阻尼。

在同一029模型 \(b=1/2,N=16\) 下，任意坐标对的 body 符号 \(v_iv_j\) 恰有一半正一半负。因此
\[
 \mathbb E_R M_{ij}^2
 =N\left(\mathbb E_R P_u^2-\mathbb E_R P_uP_v\right)=\frac13,
 \qquad \boxed{\mathbb E_R w=\frac{10}3}. \tag{7.3}
\]
Hadamard 给 \(d\le1\)，所以 \(\mathbb E_R(d^{-1/10}w)\ge10/3\)。令029的原 flag 上界为
\[
 J_{\rm flag}=4\Gamma(3/2)
 \left[\frac{\Gamma(6/5)}{\Gamma(11/10)}\right]^5(2661/64)^{-1/10}.
\]
从(7.1)得到严格、完整的辅助积分解析上界
\[
 \boxed{I_1\le c_1^{5/4}\frac{\pi^2}{60}
 \left[J_{\rm flag}-\frac{2(1-2^{-1/4})}{21}\right]^{5/4}.} \tag{7.4}
\]
它严格小于原 flag+Minkowski 的同一 Gamma 表达式。本批未对(7.4)运行数值预算或联合 \(g\)。

这一修正本身不足以解决目前约四分之一的积分预算差距：已知 \(J_{\rm flag}\ge1\)，且 \(b_q<1/4\)，新旧(7.4)公式之比大于
\((41/42)^{5/4}>(41/42)^2>19/20\)。所以该粗修正连5%的下降都无法给出。
真正可能显著改进的对象仍是(4.3)的完整轨道耦合期望；(6.2)–(6.3)提供下一步按主子式分层的明确接口。

## 8. 唯一预登记 SPD 与固定 flip 的准确验证

仅使用
\[
 M=\frac34I+\frac14\mathbf1\mathbf1^T,
 \quad S=\operatorname{diag}(-1,1,1,1,1).
\]
它属于实际 body-shape 内部：全正方向权 \(19/64\)，其他15个方向各 \(3/64\)。
准确矩阵数据为
\[
 d=81/128,\quad M^{-1}=\frac43I-\frac16\mathbf1\mathbf1^T,
 \quad\kappa=1/6,\quad s=1/7,
\]
\[
 d\det(tM^{-1}+(1-t)SM^{-1}S)=1+\frac23t(1-t).
\]
SymPy 只检查上述固定矩阵身份和正定特征值 \(3/4,2\)。
正级数(5.4)计算 \(n=0,\ldots,48\) 的有理项。相邻项比
\[
 \frac{a_{n+1}}{a_n}=\frac17\frac{(n+1/2)^2}{(n+7/4)(n+1)}<\frac17
\]
给遗漏尾部上界 \(a_{49}/(1-1/7)\)，再以384-bit Arb包围平方根和球面常数。
得到
\[
 \rho\in[
 0.9459152164542715526092280565685906589795630557727557293,
 0.9459152164542715526092280565685906589795630571033272456].
\]
相对 \(q\) 次积分上界的严格损失超过
\[
 b_q(1-\rho)>0.0086050829423523112741138848406740389042584363608904206.
\]
对应的 \(L^q\) 范数上界因子小于
\[
 0.9931099893629836089449832777752431980242297570954825983.
\]
这些是已证不等式的**上界表达式**，不是实际二元平均范数的数值积分。
同时验证 \(1-\rho>5/98\)，以及此例中(3.2)严格优于只用二元群的通用 Hölder 版本。

## 9. 边界、证据与交付

若 \(M=I\)，所有轨道 profile 相同，全部损失为零，故不能声称对所有形状存在统一固定下降比例。
单坐标 flip 若该行非对角元全零，也没有损失，即使其他子块仍非对角。
对强相关或退化极限，(6.3)必须按真实 rank 使用常数；本报告未把单矩阵改善外推为 \(R\) 的典型值。

完整轨道期望中 determinant 与能量的相关性、近退化主子式类型的概率成本，尚未得到足够强的量化。
(7.4)只闭合一个粗解析矩修正，不关闭所需的大幅积分改进，也没有改变029的 \(g\) 资格要求。

实际唯一运行 PID **43771**，0.0691秒，退出码0；Python-flint0.9.0、SymPy1.14.0；每阶段真实 checkpoint，120秒为可恢复单次时限。
复现：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/astra_orbit038/verify.py
```

脚本 `--output-dir` 可指向独立输出目录；结果 JSON 与运行收据分开。必要文件为本报告、登记、验证器、结果、真实收据与 checkpoint。

| 文件 | SHA-256 |
|---|---|
| `registration.json` | `b87b1e025482e1880587e0f053c172b05d4554fa937e466892ab6b52a8232e40` |
| `verify.py` | `b3a81d5a66587407991b85415760d9679dfa39d8b5d098b36dfc1dbb20ab6d19` |
| `results.json` | `d9da18cebeef923401a509c424f90a217ef0a7ccc4adcccdf09690d6888f072c` |

仅写独立 clone 的 `astra_orbit038/`，旧 `astra_cover015/` 锁保持只读；无提交、推送、共享 MATLAB 调用或他人文件修改。
