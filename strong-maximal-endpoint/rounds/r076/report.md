# 第76轮：近角点阈值的精确优化与源预算复用障碍

日期：2026-11-08。状态：闭合cap2固定角点模块的放大阈值优化；证明局部网格可行，并区分全局角点数与源预算复用。一般锐上界仍开放，没有新的一般上界或认证下界。

## 1. 近角点阈值已精确优化

统一J(F)=∫Φ_N(F)，Φ_N(t)=t(1+log⁺t)^(N−1)，自然对数、成功阈值1、N≥2、0≤F≤2。若实际成功矩形包含于一个固定角点矩形，体积放大≤e^τ，则后者平均严格>e^(−τ)。设

\[
 \zeta=\frac{\log(1+\log2)}{\log2},\quad p_N=1+\zeta(N-1),\qquad
 \log B_N(\tau)=\inf_{1<p\le p_N}\{Np\log[p/(p-1)]+p\tau\}.
\]

一个固定角点a的完整signed-corner halo满足|{H_aF>e^(−τ)}|≤B_N(τ)J(F)。同一a的全部2^N方向使用互不相交的源象限，没有额外2^N因子。B_N是这个Lp论证所得最优上界，未声称为算子最优弱常数。

令q_min=1/[ζ(N−1)]、f(q)=q−log(1+q)、τcrit=Nf(q_min)。τ≤τcrit时取p=p_N；τ>τcrit时有唯一q>q_min满足f(q)=τ/N，而且

\[
 \boxed{\log B_N(\tau)=N(1+q)},\qquad
 \boxed{\log B_N(\tau)=N+\sqrt{2N\tau}+O(\tau+1)\quad(\tau=o(N)).}
\]

封顶分支必须保留p_Nτ；它不是O(τ)。该分支τ=O(1/N)，故p_Nτ=O(1)。见[证明](nearcorner76.md)、[独立审查](review76.md)。

μ=1+e^(−μ)。τ=βN、角点熵为o(N)时，领先指数低于μ的阈值为β<β*=μ−1−log μ≈0.032804760864618226。

## 2. 单个几何格可共享一个角点

令ρ=√(log N/N)、c±=c_r±Kρ、r+=e^(−c−)、r−=e^(−c+)。在格∏[a^m_i,a^(m_i+1))内，所有临界向量动作矩形均包含于同一角点A_i=−r+ a^(m_i+1)的矩形；体积放大至多

\[
 \left(\frac{1+a r_+}{1+r_-}\right)^N.
\]

领先指数可支付的格比满足

\[
 1<a<a_{\max,N}=\frac{e^{\beta_*}(1+r_-)-1}{r_+},
 \qquad a_{\max,N}\longrightarrow
 \frac{e^{\mu-1}-1}{\mu-1}\simeq1.1531082515.
\]

小N可能a_max,N≤1，不能声称每个已算维数都有非平凡格。N=1048576、K=1时a_max≈1.1449412361；a=1.05给出的Hardy领先率≈1.1663181156<μ。这是固定格几何模块，全空间源成本需另外支付。[网格及紧集覆盖](cover76.md)。

## 3. 全局缺口是带成功条件的源预算分配

- **裸源重复求和发散。** F=(1/2)1_[0,1]^N在K个祖先角点盒中支付KJ(F)=K/2；全集发散。这些大尺度接收不成功，因此只否定未经成功筛选的预算求和。
- **完整高效见证需指数支付。** 已知紧支cap2见证E_N的评分≥exp(q0+αlog N+o(log N))。若以M_N个角点覆盖，或以各配一个固定角点、在分配接收集上成功的源片支付，且∑J(F_j)≤A_NJ(F)，则分别有

\[
 \log M_N,\ \log A_N
 \ge q_0+\alpha\log N-\log B_N(\tau)+o(\log N).
\]

τ=o(N)时两者的下界为(μ−1)N−o(N)。这是对完整见证的条件反推。[假设与证明](entropy76.md)。个别F_j≤F并不推出其模数和受控。

- **空间分离使角点数任意大，却不增加归一化复用。** 将一个紧支严格见证作K个远离的副本，评分不变。投影估计H_aF_K(x)≤2KL0/|x1−a1|迫使一个角点至多服务一个副本，故角点数≥K，K任意。每个副本的完整源在其M0个固定角点重复使用，复用倍数为原包的M0，与K无关；原包需M0个角点时不可误写为1。[分离副本证明](copies76.md)。低评分显式toy另有最小角点数恰为K、复用恰为1的有理证书。

有界源的任意紧严格接收集有有限精确角点覆盖：成功矩形切成以接收点为顶点的子矩形，至少一个仍成功，再用连续性和紧性提取有限子覆盖。它没有维数一致的计数结论。本轮排除了全局统一角点数这一目标，留下带成功筛选的尺度和源预算复用问题。

## 4. 有限参考前沿单独登记

令L=log N、barΦ(z_N)=√(μ/N)z_N、d_N=√(μN)z_N、q0=μN−d_N、α=(9−μ)/24。对定义的参考目标T_N=q0+αL，大N的单角点方法前沿为

\[
 \tau_*=Nf(T_N/N-1)
 =N\beta_*-p_\mu d_N+
   \left(p_\mu\alpha+\frac1{2\mu}\right)L+o(L),
 \quad p_\mu=\frac{\mu-1}{\mu}.
\]

这是参考模型，不是实际C_N有限维证书。小N若T_N<log B_N(0)，前沿不存在；若处于封顶分支则使用该分支。N=1048576的参考τ*≈33710.07404541032。保留精确d_N可把Taylor余项改进为O(L^(3/2)/√N)，但仍不改变实际常数的证明状态。[参考前沿](frontier76.md)。

## 5. 实际本机验证

确定性计算，无随机采样。Python3.12.14、SciPy1.18.1；最多12进程，BLAS每进程1线程。Wolfram15主要批次用12子内核、60/120位独立求值，结束CloseKernels[]；小型有理toy串行。本轮未启动MATLAB。主要内核计算批次不足1秒（不含启动），没有把配置12核称为持续满CPU。

|检验|范围与结果|
|---|---|
|阈值、封顶转折|N=2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576，每维10阈值；Python140、Wolfram280项|
|独立最小化|全部成功；相对闭式目标的最大绝对差2.045×10⁻⁷；验收门槛10⁻⁶|
|精度加倍|60/120位输出至35位的540字段一致；不把截取一致称为真实误差零|
|跨软件|540字段最大绝对差2.329×10⁻¹⁰，q最大相对差1.299×10⁻¹¹|
|网格|7维×3宽度：Python21、Wolfram42；跨软件42字段最大差1.999×10⁻¹⁵|
|精确预算toy|6个祖先恒等式；120个分离副本有理证书，正分离余量、角点上下数相等、复用1|

根残差、p约束、求解状态、符号凸性、内点最优值、边界分支及精度加倍通过。无截断积分；精度和容差在脚本中。浮点表不承担定向区间认证。[验证摘要](validation76.json)。

## 6. 上下界未变及下一轮两项动作

仍有Q_N≤C_N≤R_N，log Q_N=q0+αL+o(L)，log R_N=q0+(3+μ)L/12+o(L)。一般系数差(μ−1)/8未缩小，既有log C64>71.79认证未变。r073一般振幅归约已闭合。固定角点族领先指数1低于一般参考μ；高效见证排除完整低熵覆盖捷径，低评分toy只检验几何预算。

下一轮：

1. 构造有成功筛选的角点/尺度源片分配，控制∑B_N(τ_j)J(F_j)，要求分离副本不增加归一化复用；先在有限网格实现可审计优化。
2. 独立测试单包多尺度、分离高效包和临界向量动作，区分必要的指数内在复用与额外不可支付损失。

## 复现与归档

从项目根运行，Python需SciPy，Wolfram需本机许可证：

~~~bash
PY=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r076
"$PY" code/optimize76.py "$OUT/python"
"$PY" code/cover76.py "$OUT/python"
"$PY" code/copies76.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/optimize76.wl "$OUT/wolfram"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/cover76.wl "$OUT/wolfram"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/copies76.wl "$OUT/wolfram"
"$PY" code/check76.py "$OUT" "$OUT/validation76.json"
~~~

仅归档最终证明、必要代码、六份最小结果和摘要；代理草稿、缓存、日志在仓库外保留。[来源与依赖](sources76.md)。
