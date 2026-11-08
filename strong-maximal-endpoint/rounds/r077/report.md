# 第77轮：源分配的精确有限化与分数松弛障碍

日期：2026-11-08。状态：闭合固定有限角点/接收小格模型的连续源分配归约及认证计算；证明直接分数松弛存在维数损失，并验证单小格perspective修复。一般锐上界未闭合。

## 1. 本轮推进的上界接口

统一Φ_N(t)=t(1+log⁺t)^(N−1)，J(F)=∫Φ_N(F)，0≤F≤2。r076剩余问题是把整个实际halo按角点与源成本支付。本轮直接处理这个分配问题，没有用源片逐点上界代替模数和。

固定有限角点a、有限开放接收小格P和阈值σ>1。源G在有限轴平行网格上分片常数，P不跨源网格面，a_i不位于P的相应坐标开区间内部，则

\[
 Q(x)=\int_{[a,x]}G-\sigma\prod_i|x_i-a_i|
\]

在P内多仿射。检查全部2^N个闭格顶点Q(v)≥0，即认证整个P内部平均≥σ>1。退化顶点允许体积和积分均为0；内部矩形非退化。

更强地，若原F已经在包含全部角点/接收顶点坐标的公共细分上常值，那么任意可测G≤F的格平均仍≤F，保留全部顶点矩形积分，并由Jensen不增加J(G)。有界源矩形积分连续，因此原先整个开放P的约束也能延拓至顶点。故**固定分组下，任意可测源片的最小实际Orlicz成本，精确等于有限凸规划的最小值**。这是整个接收小格的精确等价，不是点云近似。[完整证明](allocation77.md)。

同一角点的源片可用pointwise max合并，J(max G_j)≤∑J(G_j)。因而有限角点覆盖可写成每角点一个源的混合整数凸规划。每接收小格至少选一个角点，顶点约束为

\[
 \sum_\ell A_{aPv,\ell}g_{a\ell}
 \ge \sigma V_{aPv}z_{aP},\qquad z_{aP}\in\{0,1\}.
\]

未选时右边为0，无需人工big-M。这个模型没有把连续全空间halo或无限角点集合有限化。

## 2. 有理原始/对偶证书

计算固定分组时，浮点线性规划仅提供候选。对Φ_2，使用t≤1时的线性支，以及t_k∈[1,2]上的下切线

\[
 L_k(t)=(2+\underline{\log t_k})t-t_k.
\]

上界用有理节点的向上Orlicz弦线。对数以20项正atanh级数和显式正尾界包围，再向外舍入至10⁻¹²。下切线确实≤Φ_2，上包络确实≥Φ_2。

浮点源舍入后与已知严格可行的F作有理凸混合，精确修复全部顶点约束。对偶乘子向下量化，检查λ_r,θ_lk≥0、∑kθ_lk≤w_l；不强迫舍入后保持驻点等式，而显式支付盒残差

\[
 d_\ell=\sum_k\theta_{\ell k}m_k-\sum_r\lambda_rA_{r\ell},
 \qquad
 \mathrm{LB}=\lambda\cdot b+\sum_{\ell k}\theta_{\ell k}b_k+
                 \sum_\ell f_\ell\min(0,d_\ell).
\]

该式是Lagrangian在0≤g≤f、非负epigraph变量上的精确下确界。因此LB≤真实有限凸最优≤UB均为有理证书；不依赖优化器宣称最优。[独立审查与基准](review77.md)。

## 3. 实际有限体积比较

源均在[0,4]²：notch在左下单位格取1、其余取1+ε，ε=1/100、1/10、1/2、1；checker在交替单位格取2和1/2；cross在中央两行或两列取2、四角余块取0。候选角点为{0,1,2,3,4}²。接收网格边界是

\[
 [-2,-1,0,1/2,1,3/2,2,5/2,3,7/2,4,9/2,5,6,8].
\]

原源的所有顶点平均至少σ+ε/10000才列入可覆盖集合，σ=1+min(ε/20,1/1000)。比较两个覆盖同一接收格集合的固定分组：greedy先选能覆盖剩余最大面积的角点；margin逐格选最小顶点平均最大的角点。没有声称这两个分组之一为全局最优。

下表用已对齐半格、64段Orlicz包络。A是该固定分组最小源成本除以原源J(F)。上/下数值均保守向外舍入：

|源|greedy可行A上界|margin最优A下界|margin/greedy成本比下界|
|---|---:|---:|---:|
|notch ε=1/2|2.785808|5.791548|2.078947|
|notch ε=1|4.624250|5.756571|1.244866|
|checker|2.155034|3.130708|1.452742|
|cross|2.978798|7.558608|2.537469|

所以“局部平均最大”与“全局源成本最小”方向可能相反，即使每组都已对源振幅做凸优化。对checker/margin，细分源格也认证了严格降成本。单位源格未包含全部半整数接收面，仅对应受限源片模型；半格满足第1节精确归约。

这些小格是实际halo内的认证接收集合，不是整个真实halo的上界；这些数值不改进C_2或一般C_N。

## 4. 两角点例排除裸分数松弛，并给出修复

令N≥2、1<σ≤2、η>0，P=(1/2,1)^N，两个真正不同的角点

\[
 a_0=0,\quad a_1=(-\eta,0,\ldots,0),\qquad
 F=2\,1_{[-\eta,1]\times[0,1]^{N-1}}.
\]

对单小格覆盖，上述整数分配的最小成本**恰为Φ_N(σ)**：至少一个角点在顶点1被选，Jensen给出此下界；G0=σ1_[0,1]^N达到。

若直接允许z0=z1=1/2，并继续收J(G0)+J(G1)，取G_j=(σ/2)1_[a_j,1]即满足缩放后的约束，成本只有σ(1+η/2)。这些源的实际平均是σ/2，不能当作阈值1的严格成功源。因此裸松弛的整数差距至少

\[
 \boxed{\frac{\mathcal I_{\rm int}}{\mathcal I_{\rm rel}}
 \ge\frac{(1+\log\sigma)^{N-1}}{1+\eta/2}.}
\]

设δ=σ−1，log(1+log(1+δ))=δ−δ²+O(δ³)。三个尺度下，对数差距下界分别为

\[
\begin{array}{c|c}
\delta&\log\text{差距下界}\\ \hline
1/N&1-\log(1+\eta/2)+o(1)\\
1/\sqrt{N\log N}&\sqrt{N/\log N}+o(\sqrt{N/\log N})\\
\sqrt{\log N/N}&\sqrt{N\log N}-\log N+o(\log N).
\end{array}
\]

固定σ>1时则是N量级损失。η=1/100、N=1048576，前三种尺度实际算得约0.99501、274.94767、3798.84354。

**修复该例：**采用perspective成本∑z_jJ(G_j/z_j)、容量G_j≤z_jF及z=0的闭延拓，Jensen给下界Φ_N(σ)，整数解达到，故该单小格的松弛恰被修复。合法的fractional source-column cover也仍可使用：每列保持真实成功阈值，以覆盖权重线性付费；不能用阈值缩放和振幅缩放冒充真实成功列。[严格证明及范围](integrality77.md)。

## 5. 本机计算与验证

Python3.12.14、SciPy1.18.1，最多12进程，每个HiGHS/BLAS计算1线程；确定性，无随机采样。Mathematica/Wolfram15用12子内核，全部批次结束CloseKernels[]。本轮未启动MATLAB。

- 24个基础分配实例、8个包络分辨率32→64的加密实例，共205组固定角点凸程序及有理证书。
- Mathematica独立重建所有矩形几何、源容量、顶点不等式和原始/对偶成本；12个LP还独立作精确有理单纯形求解，解逐项检查为Integer/Rational并精确可行，目标与SciPy最大差7.11×10⁻¹⁵。
- 原始实际Orlicz成本60/120位独立重算；精度加倍一致。基础最宽相对成本区间约1.039×10⁻⁴；8个加密实例的区间宽度比不超过0.251311。
- 高维反例56个Python值、112个Wolfram 60/120位值：维数2,4,8,16,32,64,128,256,1024,4096,16384,65536,262144,1048576；三尺度及固定δ=1/2。跨软件对数值最大差1.46×10⁻¹¹。
- 下界公式与perspective结论由解析证明支持；有限精度符合不用于证明普遍命题。本轮无截断积分，只有明确记录的有限接收家族和Orlicz包络。

计算批次本身约0.03—1.3秒，不含启动。实际使用12核并发，不将短任务描述为持续满CPU。[验证摘要和代码SHA256](validation77.json)。

## 6. 未闭合部分与下一轮两项动作

一般Q_N≤C_N≤R_N的对数系数仍为

\[
 \log Q_N=q_0+\frac{9-\mu}{24}\log N+o(\log N),\qquad
 \log R_N=q_0+\frac{3+\mu}{12}\log N+o(\log N).
\]

差距(μ−1)/8未缩小；既有log C64>71.79认证未变。r073一般振幅归约仍已闭合。本轮闭合的是有限分配的数学接口和证书，不是维数一致的全局分配预算。

下一轮两项动作：

1. 在真实成功源列上构建带线性收费的fractional cover/列生成，保留源共享并排除本轮裸松弛漏洞；比较其预算与整数分组。
2. 检验perspective在多接收小格与临界高维源中的适用性，分别审查单角点合并、跨角点相关覆盖及分离副本；只有得到维数一致预算后才尝试接入Q_N同对数系数。

## 复现

从项目根运行，Python需SciPy，Wolfram需本机许可证：

~~~bash
PY=/Users/zzh/.local/share/math-research-tools/venv/bin/python
OUT=/tmp/strong-maximal-r077
"$PY" code/allocate77.py "$OUT/python"
"$PY" code/allocate77.py "$OUT/python" refine
"$PY" code/gap77.py "$OUT/python"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/verify77.wl "$OUT/python/allocation77.json" "$OUT/wolfram"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/verify77.wl "$OUT/python/allocation77_refine.json" "$OUT/wolfram"
WolframKernel=/Applications/Wolfram.app/Contents/MacOS/WolframKernel /usr/local/bin/wolframscript -file code/gap77.wl "$OUT/wolfram"
"$PY" code/check77.py "$OUT" "$OUT/validation77.json"
~~~

仅归档最终证明、代码、有理证书和验证结果；中间失败日志、重复summary表及缓存在仓库外保留。[来源索引](sources77.md)。
