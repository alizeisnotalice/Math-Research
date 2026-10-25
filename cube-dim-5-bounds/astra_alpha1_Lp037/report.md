# alpha1-Lp037：完整 L³ᐟ² 预算与奇异补集体积界

只取一个新指数 \(r=3/2\)，本轮证明并实际认证

\[
\boxed{\int_{[0,1]^4}G_1(1,y)^{3/2}\,dy
\le K_0:=\frac{\sqrt{10}\,\Gamma(2/3)^{9/2}}
{90\pi^4\,887^{1/4}}
<\frac1{3840}.}
\tag{1}
\]

\(K_0\approx0.00025860165943690446222510672064\)。对任意可测补集 \(E\subset[0,1]^4\)、\(|E|=v\)，有

\[
\boxed{\int_EG_1(1,y)^{5/4}\,dy
\le K_0^{5/6}v^{1/6}
<3840^{-5/6}v^{1/6}.}
\tag{2}
\]

因此一个完全有理的安全条件是

\[
\boxed{v\le\frac{248832}{298023223876953125}
=\frac{3840^5}{10^{30}}\quad\Longrightarrow\quad
\int_EG_1^{5/4}<10^{-5}.}
\tag{3}
\]

该体积约 \(8.35\times10^{-13}\)。若直接使用精确 Gamma 表达式 \(K_0\)，本估计允许 \(v\le10^{-30}/K_0^5\approx8.64656206\times10^{-13}\)。这只是该误差包络允许的体积，并非真实管积分的必要体积限制。\(v=10^{-12}\) 时本界约为 \(1.02453333\times10^{-5}\)，已略高于目标。预算可用于补奇异管，但需要很小的总体积；本轮不生产管几何或空间积分证书。

结论覆盖全部正时间、全部正单位面和退化 shape 尾，包括 \(S_1,S_2\) 附近的积分。它不是 \(G\) 的逐点上界；canonical \(G=\infty\) 在这些零测集上仍然保留。没有新全域 \(g\) 或 \(C_5\) 结论。

## 1. 固定模型和球面估计

固定 \(\alpha=1,\beta=1/2\)，16 个独立正 stable 变量满足 \(\mathbb Ee^{-\lambda S_i}=e^{-\sqrt\lambda}\)，方向 \(v_i=(1,\pm1,\pm1,\pm1,\pm1)\) 未单位化。写

\[
A=\sum_iS_iv_iv_i^T,\quad T=\sum_iS_i,\quad M=A/T,
\quad d=\det M,
\]
\[
\frac{dR}{d\mathbb P}=\frac{T^{-1/2}}{Z_T},\qquad
Z_T=\frac1{16\Gamma(3/2)}.
\]

采用已审 alpha1 的完整时间接口

\[
G_1(x)=c_1\mathbb E_R H(M,x),\qquad
c_1=\frac1{32\pi^3},\quad
H(M,x)=d^{-1/2}(x^TM^{-1}x)^{-2}.
\tag{4}
\]

\(M\) 几乎处处正定，且对角元全为 1，\(\operatorname{tr}M=5\)。对 \(\omega\in S^4\)，置 \(Q=\omega^TM^{-1}\omega\)，则
\(Q\ge1/\lambda_{\max}(M)\ge1/5\)。固定 \(r=3/2\) 时

\[
Q^{-2r}=Q^{-3}\le\sqrt5\,Q^{-5/2}.
\]

Gaussian 径向积分给

\[
\int_{S^4}Q^{-5/2}d\sigma=\frac{8\pi^2}{3}\sqrt d,
\]

所以逐个正定 \(M\) 有

\[
\int_{S^4}H(M,\omega)^r\,d\sigma
\le\sqrt5\frac{8\pi^2}{3}\,d^{(1-r)/2}
=\sqrt5\frac{8\pi^2}{3}\,d^{-1/4}.
\tag{5}
\]

常数一般形式为 \(5^{2r-5/2}\)，此处准确为 \(\sqrt5\)。没有把一个 cube 面上的积分直接当成全球面比例。

## 2. 从球面到面片：轨道顺序与 Minkowski

令 \(\partial Q_c\) 是 \([-1,1]^5\) 的十个单位面组成的边界。对每条非边缘射线，令 \(z=\omega/\|\omega\|_\infty\in\partial Q_c\)，则 \(|z|\ge1\)，且
\(d\sigma=|z|^{-5}dS(z)\)。由于 \(H^r\) 是 \(-4r=-6\) 次齐次，

\[
\int_{\partial Q_c}H(M,z)^r,dS
=\int_{S^4}|z(\omega)|^{5-4r}H(M,\omega)^r,d\sigma
\le\int_{S^4}H(M,\omega)^r,d\sigma,
\tag{6}
\]

其中 \(5-4r=-1\)，所以方向确实是上式的“\(\le\)”。边和角的径向图重叠为零表面积集，不产生额外系数。

令 \(\mathcal O\) 为全部 signed permutation 群，并先对固定 \(M\) 定义

\[
\overline H_M(x)=\frac1{|\mathcal O|}\sum_{U\in\mathcal O}H(M,Ux).
\]

正单位面片 \(F_+=\{(1,y):y\in[0,1]^4\}\) 的轨道等权覆盖 \(10\cdot16=160\) 个面片。对 \(r\ge1\) 先用有限平均的凸性，再用 (5)–(6)，得到

\[
\begin{aligned}
\int_{F_+}\overline H_M(x)^r,dx
&\le\frac1{160}\int_{\partial Q_c}H(M,z)^r,dS\\
&\le\frac{\pi^2}{60}\sqrt5\,d^{-1/4}.
\end{aligned}
\tag{7}
\]

原方向集在 \(\mathcal O\) 下仅发生排列和整体取负，独立同分布 \(S_i\) 及 \(T^{-1/2}\) 倾斜保持不变，所以 \(R\) 在 \(M\mapsto U^TMU\) 下不变。因此 (4) 也等于 \(c_1\mathbb E_R\overline H_M\)。**先做固定 shape 的轨道平均，再用 Minkowski**，得到

\[
\|G_1\|_{L^r(F_+)}
\le c_1\left(\frac{\pi^2}{60}\sqrt5\right)^{1/r}
\mathbb E_Rd^{-\gamma},\qquad
\gamma=\frac{r-1}{2r}=\frac16.
\tag{8}
\]

下一节证明右端有限。可先限制到递增的 invariant 事件 \(\lambda_{\min}(M)\ge1/n\)，在那里 shape 核一致有界，应用 Minkowski 后用单调收敛恢复全部正定 \(M\)。这不把 \(G\) 在奇异集上的无穷值视为有限。固定 shape 的 \(d^{-1/4}\) 在取 \(r\) 次根后成为 \(d^{-1/6}\)；两者不可混淆。

## 3. Harris 与新 trace 正矩保证完整 shape 尾有限

尺度分离现在不能沿用 \(5b=\beta\) 的抵消。准确有

\[
\mathbb E_Rd^{-\gamma}
=16\Gamma(3/2)\mathbb E[T^{5\gamma-\beta}(\det A)^{-\gamma}],
\quad 5\gamma-\beta=\frac13=:k<\beta.
\tag{9}
\]

在原独立 \(S\) 概率中，\(T^k\) 对每个坐标递增，而 \((\det A)^{-\gamma}\) 对每个坐标递减：增加一个 \(S_i\) 只给 \(A\) 增加正半定秩一矩阵。对有界异向单调函数，Harris 的方向为

\[
\mathbb E[fg]\le\mathbb Ef\,\mathbb Eg.
\]

这里可直接给出产品概率证明：一坐标时用独立副本 \(X'\)，协方差等于
\(\frac12\mathbb E[(f(X)-f(X'))(g(X)-g(X'))]\le0\)；逐坐标条件化，再对剩余坐标归纳。取
\(f_n=\min(T^k,n)\)、\(g_n=\min((\det A)^{-\gamma},n)\)，分别仍递增、递减。先用有界结论，再以单调收敛令 \(n\to\infty\)，得到

\[
\mathbb E[T^k(\det A)^{-\gamma}]
\le\mathbb ET^k\,\mathbb E(\det A)^{-\gamma}.
\tag{10}
\]

此证明不预设左端可积；下列有限右端同时证明其可积。也没有在倾斜 \(R\) 下误用独立性。

\(T\) 的 Laplace 变换为 \(e^{-16\sqrt\lambda}\)。对 \(0<k<1/2\) 的正 Mellin 积分给

\[
\mathbb ET^k=16^{2k}\frac{\Gamma(1-2k)}{\Gamma(1-k)},\quad
\mathbb ET^{1/3}=16^{2/3}\frac{\Gamma(1/3)}{\Gamma(2/3)}.
\tag{11}
\]

这是本轮新算的 trace **正矩**，不能替换成负矩，也不能沿用 \(\beta=3/4\) 数值。

另一方面，读取 014 锁定的纯几何 flag 概率。每个 flag 的五个线性因子在原 \(S\) 概率中独立、其支持覆盖全部 16 个变量。对任意正指数 \(\gamma\)，同一个加权 AM–GM / 普通 Hölder 推导给

\[
\mathbb E(\det A)^{-\gamma}
\le16^{-5\gamma}
\left[\frac{\Gamma(1+2\gamma)}{\Gamma(1+\gamma)}\right]^5
Z_f^{-\gamma}.
\tag{12}
\]

这里 \(Z_f=\sum_f\theta_fR_f(1/2)^2\)，最优 flag 权是
\(\theta_fR_f(1/2)^2/Z_f\)，不依赖负矩指数 \(\gamma\)：对数目标只有一个共同正因子 \(\gamma\)。本轮从原概率直方图重新取平方根，精确核得五个 \(R_f^2=(81/8,27,36,36,64)\)，以及 \(Z_f=2661/64\)，再代入**新**指数 \(\gamma=1/6\)：

\[
\mathbb E(\det A)^{-1/6}
\le16^{-5/6}
\left[\frac{\Gamma(4/3)}{\Gamma(7/6)}\right]^5
\left(\frac{2661}{64}\right)^{-1/6}.
\tag{13}
\]

合并 (9)–(13)，令

\[
B_0=16^{5/6}\Gamma(3/2)
\frac{\Gamma(1/3)}{\Gamma(2/3)}
\left[\frac{\Gamma(4/3)}{\Gamma(7/6)}\right]^5
\left(\frac{2661}{64}\right)^{-1/6},
\]

得到 \(\mathbb E_Rd^{-1/6}\le B_0<\infty\)。Gamma 递推、倍乘和反射恒等式

\[
\Gamma(1/6)\Gamma(2/3)=2^{2/3}\sqrt\pi\Gamma(1/3),\quad
\Gamma(1/3)\Gamma(2/3)=\frac{2\pi}{\sqrt3}
\]

上述恒等式可直接由 [NIST DLMF 的反射与倍乘公式](https://dlmf.nist.gov/5.5) 代入正参数得到，进一步化简为

\[
\boxed{B_0=\frac{64\Gamma(2/3)^3}{\sqrt3\pi\,2661^{1/6}}.}
\tag{14}
\]

原独立 stable 的全部所需负矩有限，trace 正矩因 \(1/3<1/2\) 有限；这控制了所有 \(d\downarrow0\) 的 shape 尾，没有施加人为统一条件数。

## 4. 全域常数、奇异集与补集误差

由 (8)、(14)，

\[
\int_{F_+}G_1^{3/2}
\le c_1^{3/2}\frac{\pi^2}{60}\sqrt5\,B_0^{3/2}=K_0.
\]

使用 \(2661=3\cdot887\) 恰得到报告开头的简式 (1)。这涵盖 canonical \(G\) 在 \(S_1,S_2\) 的真实无穷值：这些有限个方向线/方向平面与四维正单位面的交集均为零四维测度。\(L^{3/2}\) 界是面上的 a.e. 积分结论，不宣称点值有限，也没有删除其任意厚度邻域。任意覆盖这些集合的可测管 \(E\) 的整个贡献都受 (2) 控制。

对 \(q=5/4\)，Hölder 指数为 \(r/q=6/5\) 与其共轭 6，故

\[
\int_EG_1^{5/4}\le
\left(\int_EG_1^{3/2}\right)^{5/6}v^{1/6}
\le K_0^{5/6}v^{1/6}.
\]

因此误差预算 \(\varepsilon=10^{-5}\) 的充分条件为 \(K_0^5v\le\varepsilon^6\)。为了得到仅含有理数的可用条件，实际 Arb 严格核出 \(K_0<1/3840\)，于是 (3) 完全由有理幂比较成立。无需把显示小数当作向上端点。

## 5. 实际有限批与复现

运行前已登记唯一 \(r=3/2\)、256-bit Arb、一个完整标量批和 120 秒 checkpoint。实际只运行 `certify_Lp.py`，PID 42670，约 0.0121 秒，退出 0。对 014 manifest 的全部 7 个白名单 hash 已核验；未运行其旧 baseline、033/034、任何 QMC 或空间网格。

实际完整球给：

|量|显示近似|
|---|---:|
|\(\mathbb ET^{1/3}\)|12.5618300277584|
|\(\mathbb E(\det A)^{-1/6}\) 的显式上界|0.0440433455246175|
|\(B_0\)|7.84509372689966|
|\(K_0\)|0.000258601659436904|
|\(K_0^{5/6}\)|0.00102453332556167|
|\(10^{-30}/K_0^5\)|8.64656206081621e−13|

完整球、三个预登记有理上端的严格比较、源码/输入 SHA‑256 和精确有理安全体积在 `results.json`。原 Gamma 乘积和化简表达式通过独立 Arb 路径交叉包围；化简的严格依据是上文精确恒等式，不是浮点接近。

复现（在本目录）：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B certify_Lp.py
```

仅写 `astra_alpha1_Lp037/`，旧 025、026、030、033 锁定文件不改；没有 MATLAB、采样、大缓存或临时空间产物。必要文件为本报告、登记、最终脚本、结果和 manifest。无活动句柄，未提交、未推送。
