# alpha1-Lp045：一般 L^r 理论与三个指数的完整管预算

本批对全部
\[
\frac54<r<\frac53,\qquad
\gamma=\frac{r-1}{2r},\qquad k=5\gamma-\frac12
\]
证明完整正单位面上的积分上界
\[
\boxed{\int_{[0,1]^4}G_1(1,y)^r\,dy\le K_r},
\quad
K_r=c_1^r\frac{\pi^2}{60}5^{2r-5/2}B_r^r,
\quad c_1=\frac1{32\pi^3},                                  \tag{1}
\]
其中
\[
\boxed{B_r=8\sqrt\pi\,16^{5\gamma-1}
 \frac{\Gamma(1-2k)}{\Gamma(1-k)}
 \left[\frac{\Gamma(1+2\gamma)}{\Gamma(1+\gamma)}\right]^5
 \left(\frac{2661}{64}\right)^{-\gamma}}
\quad\text{满足 }\mathbb E_R d^{-\gamma}\le B_r.               \tag{2}
\]

只计算三个预登记的新指数。用 034 已证准确管体积
\(V(\delta)=108\delta^3-243\delta^4\)、目标管积分误差 \(10^{-6}\)，结果为：

| r | K_r 显示近似 | 最大通过的 dyadic 半径 | 此半径的严格误差上端，显示近似 | 对 037 已锁半径的倍率 |
|---|---:|---:|---:|---:|
| 8/5 | .00127816644883543746 | 2^-22 | 6.8898893904692881e-7 | 2 |
| 13/8 | .00291877103302939924 | 2^-22 | 8.5997743107160160e-7 | 2 |
| 33/20 | .0142834295738370512 | 2^-24 | 6.9303519207978405e-7 | 1/2 |

前两个指数将 037/034 已锁的 \(2^{-23}\) 半径提高到 \(2^{-22}\)；其中 \(r=8/5\) 在该共同半径上有更小的预算。第三个指数没有改进。每个下一更大 dyadic 半径均准确不满足本包络，并且以 \(K_r\) 的严格有理**下端**也已证明失败，故差别不是向上舍入造成的。这些是充分误差预算的阈值，不是真实管积分的必要阈值，更不是全部 \(r\) 的优化结论。

全部证明覆盖完整时间、整个正单位面和全部退化 shape 尾。没有新的真实 I、全域 g 或 C5 上界。本批没有 Green 点、空间网格、采样或连续 r 网格；037 的旧 Gamma 数值只读不重算。

## 1. 模型与完整时间归一化

固定五维、16 条未单位化方向 \(v_i=(1,\pm1,\pm1,\pm1,\pm1)\)，各方向共同时间 Cauchy 特征指数 \(t|\xi|\)。独立正变量 \(S_i\) 满足
\(\mathbb E e^{-\lambda S_i}=e^{-\sqrt\lambda}\)。写
\[
A=\sum_{i=1}^{16}S_iv_iv_i^T,\quad T=\sum_iS_i,
\quad M=A/T,\quad d=\det M.
\]
所有 \(S_i>0\) 几乎处处，方向张成全空间，故 \(M\) 几乎处处正定。其对角元均为 1，\(\operatorname{tr}M=5\)。在原产品概率上
\[
\mathbb E e^{-\lambda T}=e^{-16\sqrt\lambda},\qquad
Z_T:=\mathbb E T^{-1/2}=\frac1{8\sqrt\pi}
=\frac1{16\Gamma(3/2)}.
\]
定义 shape 概率为 \(dR/d\mathbb P=T^{-1/2}/Z_T\) 的 M 推前；不假定 T 与 M 独立。

Gaussian subordination 的正核给
\[
p_t(x)=(4\pi)^{-5/2}t^{-5}
\mathbb E[(\det A)^{-1/2}e^{-x^TA^{-1}x/(4t^2)}].
\]
由 Tonelli 与
\(\int_0^\infty t^{-5}e^{-q/(4t^2)}dt=8q^{-2}\)，得到 canonical 完整 Green 表示
\[
G_1(x)=\frac1{4\pi^{5/2}}\mathbb E[(\det A)^{-1/2}
 (x^TA^{-1}x)^{-2}]
=c_1\mathbb E_R H(M,x),                                    \tag{3}
\]
\[
H(M,x)=d^{-1/2}(x^TM^{-1}x)^{-2},\qquad c_1=1/(32\pi^3).
\]
所有核非负，式 (3) 可以取无穷值；未在此预先删掉任何奇异点，也没有时间截断。

## 2. 球面、边界与轨道平均的顺序

对固定正定 M 和 \(\omega\in S^4\)，记 \(Q=\omega^TM^{-1}\omega\)。由 \(\lambda_{\max}(M)\le5\)，有 \(Q\ge1/5\)。当 \(r>5/4\) 时，
\[
Q^{-2r}\le5^{2r-5/2}Q^{-5/2},\qquad
\int_{S^4}Q^{-5/2}d\sigma=\frac{8\pi^2}{3}\sqrt d.
\]
后一个恒等式由 Gaussian 积分的椭球变量变换及径向积分得到。因此
\[
\int_{S^4}H(M,\omega)^r\,d\sigma
\le\frac{8\pi^2}{3}5^{2r-5/2}d^{(1-r)/2}.                 \tag{4}
\]

将球面径向投影到 \(\partial[-1,1]^5\)，令
\(z=\omega/\|\omega\|_\infty\)。除零表面积的面缘方向外，
\(d\sigma=|z|^{-5}dS\)。H^r 为 \(-4r\) 次齐次，且 \(|z|\ge1\)，故
\[
\int_{\partial[-1,1]^5}H(M,z)^r\,dS
=\int_{S^4}|z(\omega)|^{5-4r}H(M,\omega)^r\,d\sigma
\le\int_{S^4}H(M,\omega)^r\,d\sigma.                      \tag{5}
\]
此处方向依赖 \(r\ge5/4\)，不能反用径向权重。

令 \(\mathcal O\) 为全部 signed permutation 群，先在固定 M 上定义
\[
\bar H_M(x)=|\mathcal O|^{-1}\sum_{U\in\mathcal O}H(M,Ux).
\]
正面片 \(F_+=\{(1,y):y\in[0,1]^4\}\) 的轨道等权覆盖全部 160 个面片。先对有限平均用凸性，再用 (4)–(5)：
\[
\|\bar H_M\|_{L^r(F_+)}
\le\left(\frac{\pi^2}{60}5^{2r-5/2}\right)^{1/r}
d^{-(r-1)/(2r)}.                                         \tag{6}
\]
方向集合在群作用下仅被排列或整体取负，独立同分布 S 与 T 倾斜均保持不变，所以 R 的 M 分布在共轭下不变，并有
\(\mathbb E_R\bar H_M=\mathbb E_RH(M,\cdot)\)。接着才用 Minkowski：
\[
\|G_1\|_{L^r(F_+)}\le
c_1\left(\frac{\pi^2}{60}5^{2r-5/2}\right)^{1/r}
\mathbb E_Rd^{-\gamma}.                                  \tag{7}
\]
没有对单个固定 M 的单独面片擅自乘全球面的 1/160。

为避免循环假设可积性，可先限制于 invariant 递增事件
\(E_n=\{\lambda_{\min}(M)\ge1/n\}\)。每个事件上核在面片上一致有界，可以直接应用以上步骤；再用下一节证明的有限 \(\mathbb E_Rd^{-\gamma}\) 和单调收敛恢复全部 M。

## 3. 原产品概率上的 Harris 分离与完整 shape 尾

准确尺度恒等式为
\[
\mathbb E_R d^{-\gamma}
=8\sqrt\pi\;\mathbb E[T^k(\det A)^{-\gamma}],\qquad
k=5\gamma-\frac12.                                       \tag{8}
\]
本轮范围 \(5/4<r<5/3\) 恰给 \(0<k<1/2\)。在原 S 产品概率中，\(T^k\) 对每个坐标递增，\((\det A)^{-\gamma}\) 对每个坐标递减。后者因为增加一个 S_i 只给 A 增加正半定秩一矩阵，正定行列式随 Loewner 序递增。

有界异向单调函数的协方差非正：一变量时，用独立副本表达协方差为
\(\frac12\mathbb E[(f(X)-f(X'))(g(X)-g(X'))]\le0\)；逐坐标条件化再归纳得到产品概率版本。对
\(f_N=\min(T^k,N)\)、\(g_N=\min((\det A)^{-\gamma},N)\) 使用该结论，再以单调收敛取极限，得
\[
\mathbb E[T^k(\det A)^{-\gamma}]
\le\mathbb ET^k\;\mathbb E(\det A)^{-\gamma}.              \tag{9}
\]
独立性只在原概率中使用，没有对倾斜的 R 误用 Harris。下列有限右端同时证明左端可积，不先假设这一点。

对于 \(0<k<1/2\)，正 Mellin 公式及 \(u=16\sqrt\lambda\) 代换给
\[
\mathbb ET^k=\frac{k}{\Gamma(1-k)}
\int_0^\infty(1-e^{-16\sqrt\lambda})\lambda^{-k-1}d\lambda
=16^{2k}\frac{\Gamma(1-2k)}{\Gamma(1-k)}.                  \tag{10}
\]
端点可积性分别由 \(k<1/2\) 和 \(k>0\) 保证。

消费 014/022 锁定的纯几何 flag 身份
\[
\det A=16^5\sum_f\theta_f^0 F_f,\qquad\sum_f\theta_f^0=1.
\]
每个 F_f 是五个正线性 stable 因子的乘积，五块支持互不相交并覆盖全部 16 个变量，其中两个是 singleton。几何身份不在本批重新生产。对一个概率块 p，令 \(s_p=\sum_i\sqrt{p_i}\)，则
\[
\mathbb E\left(\sum_i p_iS_i\right)^{-\gamma}
=\frac{\Gamma(1+2\gamma)}{\Gamma(1+\gamma)}s_p^{-2\gamma}.
\]
记 \(R_f=\prod_{\mathrm{blocks}}s_p\)。五块独立，所以
\[
\mathbb EF_f^{-\gamma}
=\left[\frac{\Gamma(1+2\gamma)}{\Gamma(1+\gamma)}\right]^5R_f^{-2\gamma}.
\]
任意正权 w 的加权 AM–GM，随后原概率上的普通 Hölder，给
\[
\mathbb E(\det A)^{-\gamma}
\le16^{-5\gamma}e^{\gamma D(w\Vert\theta^0)}
\prod_f(\mathbb EF_f^{-\gamma})^{w_f}.
\]
取 \(Z_f=\sum_f\theta_f^0R_f^2\)、\(w_f^*=\theta_f^0R_f^2/Z_f\)，则
\[
D(w\Vert\theta^0)-2\sum_fw_f\log R_f
=D(w\Vert w^*)-\log Z_f.
\]
因此对于每个 \(\gamma>0\)，相同最优权给
\[
\mathbb E(\det A)^{-\gamma}
\le16^{-5\gamma}
\left[\frac{\Gamma(1+2\gamma)}{\Gamma(1+\gamma)}\right]^5
Z_f^{-\gamma}.                                          \tag{11}
\]
014/037 已锁的五种 \(R_f^2\) 是 \((81/8,27,36,36,64)\)，对应类型权为 \((1/8,3/16,3/16,1/8,3/8)\)。其加权和准确为 \(2661/64\)，脚本只做这个有理和的检查，不重算旧 Gamma。

将 (10)–(11) 代入 (8)–(9)，注意
\(2k-5\gamma=5\gamma-1\)，恰得到 (2)。所有 Gamma 参数均正、全部 block 负矩有限，且 trace 正矩满足 \(k<1/2\)，故控制了所有 \(d\downarrow0\) 的尾部。结合 (7) 即得 (1)。当 \(r\uparrow5/3\) 时，\(k\uparrow1/2\)，本公式中的 \(\Gamma(1-2k)\) 发散；不能据“更高 r”预设预算一定更好。本批未声称真实 L^r 积分在此端点发散。

## 4. 奇异集合和完整管积分

五维 canonical Green 在 S1/S2 上的无穷值继续保留。它们是有限个线/平面与四维面片的交集，零四维测度，不影响 L^r 的 a.e. 积分结论。上面推导没有删掉任意正厚度邻域；故对任何可测 \(E\subset F_+\)，\(|E|=v\)，用 \(q=5/4<r\) 的 Hölder 得
\[
\boxed{\int_E G_1^{5/4}\le K_r^{q/r}v^{1-q/r}.}           \tag{12}
\]

034 已锁并独立提供 15 条 face 对角线的 \(\ell^\infty\) 管并集体积
\[
V(\delta)=108\delta^3-243\delta^4,\qquad 0\le\delta\le1/3.
\]
该管覆盖 S2 与正面片的相关集合。本批把这个精确几何结论作为外部输入，不重新分类 120 对或重新构管。整个管（含任意退化邻域）的贡献受 (12) 控制。

设 \(q/r=A/D\)、\(1-q/r=B/D\)。用严格有理上端 \(\widehat K_r>K_r\) 后，误差不超过 \(10^{-6}\) 的充分条件恰为
\[
\widehat K_r^A V(\delta)^B\le10^{-6D}.                    \tag{13}
\]
三个指数具体是
\[
\begin{array}{c|c}
r&\text{准确有理比较}\\\hline
8/5&\widehat K^{25}V^7\le10^{-192}\\
13/8&\widehat K^{10}V^3\le10^{-78}\\
33/20&\widehat K^{25}V^8\le10^{-198}.
\end{array}
\]
准确 \(\gamma,k\) 分别为 \((3/16,7/16),(5/26,6/13),(13/66,16/33)\)。此外
\(V'(\delta)=324\delta^2(1-3\delta)\ge0\)，所以通过一个 dyadic 半径并证明其立即上一个半径失败，即证明本包络中的最大 dyadic 半径；不需要空间网格或连续参数搜索。

## 5. 实际完整球、有理化和精确阈值

每个新 r 只运行一次 384-bit Arb 表达式 (2)、(1)。完整球与准确 dyadic 外端保存在 results.json。预先固定舍入分母为 \(10^{30}\)：若 u/l 是 K 球的外上/下端，取
\[
\widehat K=(\lfloor10^{30}u\rfloor+1)/10^{30},\qquad
\underline K=\lfloor10^{30}l\rfloor/10^{30}.
\]
实际严格核得 \(0<\underline K<K_r<\widehat K\)。三个上端按最简分数记录为
\[
\begin{array}{c|c}
8/5&639083224417718729137497879/500000000000000000000000000000\\
13/8&1459385516514699619910426879/500000000000000000000000000000\\
33/20&14283429573837051222702234887/1000000000000000000000000000000.
\end{array}
\]

半径搜索只对整数指数 n 使用阈值：先以 2,4,8,16,32 找到两端，再二分，测试 \(\delta=2^{-n}\)。所有比较都是 Fraction 的整数幂；未进行超越函数阈值搜索。results.json 保留每次二进制比较及准确 lhs/rhs 分数。

独立 verify_rational.py 不导入 Arb，不调用 Gamma：它重读保存的球外端与有理 cap，用另一形式 \(V=27\delta^3(4-9\delta)\) 复核全部准确幂不等式，并以整数开根的有理外端围住误差。结果保存在 rational_verification.json。对于每个候选，严格成立
\[
\widehat K^A V(\delta)^B\le10^{-6D}
<\underline K^A V(2\delta)^B.
\]
对应较大半径的包络误差下端分别约为
\(1.0858332430\times10^{-6}\)、\(1.3896082435\times10^{-6}\)、\(1.1473242693\times10^{-6}\)。故三个准确阈值不会被上端舍入误差左右。

同半径 \(2^{-22}\) 上，r=8/5 与 13/8 的两个有理包络通过交叉整数幂准确比较，前者更小；不依靠显示小数排序。旧 037 的 \(r=3/2\)、\(K<1/3840\) 与 034 已锁 \(\delta=2^{-23}\) 只读取及核 hash，没有任何旧 Gamma 重跑。

## 6. 过程、依赖与剩余边界

运行前已写 registration.json 并向主报告源码 hash。PID45191、session95083 先处于 armed 状态；发出句柄后才输入换行启动三个标量。各耗时约 .009785、.000421、.000633 秒，均正常完成，总进程退出 0。每标量设置 120 秒超时检查点并立即保存；无超时、失败或重启。最终检查点位于 `/tmp/astra_alpha1_Lp045/checkpoint.json`，没有活动任务。

消费的锁定输入及 hash 已完整列于 registration/results：037 report/results，034 cover_results/tube_link_results，014 原概率类型。几何分解、管几何与旧 r 结果只读使用；一般 r 的归一化、球面到 face 顺序、Harris、矩有限性及 Hölder 在本报告重新证明。

复现只需本机 Python 标准库和 python-flint：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_alpha1_Lp045/certify.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B cube-dim-5-bounds/astra_alpha1_Lp045/verify_rational.py
```

第一条是三个新指数的重现入口，第二条只复核保存的有理证书。主归档前无需再次调用第一条。必要白名单为登记、两脚本、完整结果、纯有理复核结果、报告和 manifest；未提交、未更动其他目录证据。

本批放宽了一个可用充分管半径两倍，但半径仍只有约 2.38e-7。它没有证明这种极细管能完成常规区连续积分认证，没有估计真实 I 或 C5，也没有排除未测试的 r 或更强空间预算机制。现有 C5 上界 32 不变。
