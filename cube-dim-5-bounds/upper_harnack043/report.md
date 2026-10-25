# upper_harnack043：正 Cauchy 混合的全空间 log-Hessian 界与一个局部认证

## 结果与范围

本轮固定 \(\alpha=1\)、五维16条独立、未单位化 body Cauchy 方向，Fourier 逆因子为
\((2\pi)^{-5}\)。令 \(p_t\) 为该过程的密度，
\[
 H_\epsilon(x)=\int_\epsilon^\infty p_t(x)\,dt,\qquad \epsilon>0.
\]
下文严格证明：任取由五条 body 方向组成的可逆矩阵 \(N\)，对所有 \(x,u\in\mathbb R^5\)，
\[
 -\frac{2\|N^{-1}u\|_2^2}{\epsilon^2}
 \le u^T\nabla^2\log H_\epsilon(x)u
 \le\frac{\|N^{-1}u\|_2^2/4+\|N^{-1}u\|_1^2}{\epsilon^2}
 \le\frac{21\|N^{-1}u\|_2^2}{4\epsilon^2}. \tag{1}
\]
因此中心的准确值及准确梯度确定全空间上下 log 抛物面，见 (8)。
证明不要求剩余 Cauchy 变量具有一阶或二阶矩。

唯一登记的新数值对象为 \(H_{1/64}\) 在018中心的值和四个尾坐标梯度；没有计算第五个梯度。
因此数值区域证书在仿射平面 \(x_1=1\) 上；(1) 的解析定理则在整个五维空间成立。
全部3504锥、时间 \([1/64,\infty)\)、384-bit Arb 给出
\[
 H_{1/64}(x_0)=0.00012008333756253879717053112826917\ldots,
\]
其中
\[
 x_0=\left(1,\frac{104115}{104729},\frac{103354}{104729},
                 \frac{71094}{104729},\frac{6327}{104729}\right).
\]
对原登记闭盒 \(|x_j-x_{0,j}|\le1/16384\;(j=2,3,4,5)\)，严格统一下界为
\[
 H_{1/64}(x)\ge
 [0.000120060178462742395304488412493151969916902087363374\ldots
 \ \pm2.40\cdot10^{-84}]_{\rm lower}. \tag{2}
\]
右侧含义是结果 JSON 中该下界表达式球的准确 dyadic 下端点；不是说盒内所有函数值落在该窄球。
它严格大于 \(1/10000\)，也大于比较门槛 \(11/100000\)。统一上界表达式约为
\(0.00012010914836290658210381576096286\)。

同一组中心值和梯度给出更大、明确的安全区域，未增加函数采样点：

* 以 \(x_0\) 为中心、在整个 \(x_1=1\) 平面内的闭立方体，半宽
  \(h=53571/8388608\approx0.006386160850524902\)。整个立方体上 \(H_{1/64}>1/10000\)。
* 若要求立方体包含于正单位面 \(\{1\}\times[0,1]^4\)，采用
  \(h_{\rm face}=614/104729\approx0.005862750527552063\)。
* 令 \(v=x_{2:5}-x_{0,2:5}\)、\(Q=I_4/4-J_4/36\)。闭椭球
  \[
    v^TQv\le\left(\frac{51857}{8388608}\right)^2 \tag{3}
  \]
  也满足 \(H_{1/64}>1/10000\)；如需要可与正单位面相交。

因为 \(p_t\ge0\)，扩展值 Green 势 \(G=\int_0^\infty p_tdt\) 满足 \(G\ge H_{1/64}\)，
所以以上局部下界同时适用于 \(G\)。这里没有对遗漏的 \((0,1/64)\) 时间作误差估计，
也不需要该时间段的有限上界。没有证明全域 \(g\)、全局空间覆盖或新的 \(C_5\) 上界。

## 1. 正混合、可微性及有限性

单变量 Cauchy 密度为
\[
 s_t(z)=\frac{t}{\pi(t^2+z^2)},\qquad
 \widehat s_t(\theta)=e^{-t|\theta|}.
\]
将16个独立变量中的五个放入固定可逆 body 基 \(N\)，其密度为
\[
 q_t(z)=\frac1{|\det N|}\prod_{j=1}^5s_t((N^{-1}z)_j).
\]
其余11条方向在 \(t=1\) 的和记为随机向量 \(Y\)；由 Cauchy 的尺度性质，
\[
 p_t(x)=\mathbb E q_t(x-tY),\qquad
 H_\epsilon(x)=\int_\epsilon^\infty\mathbb E q_t(x-tY)\,dt. \tag{4}
\]
可以在与 \(t\) 无关的同一个概率空间取 \(Y\)。这只用独立性、Cauchy 尺度及卷积，
不要求 \(Y\) 有密度或有限矩。

对实方向 \(u\)，写 \(a=N^{-1}u\)、\(w=N^{-1}(x-tY)\)，并令
\(\sigma_u=\partial_u\log q_t(x-tY)\)、\(\tau_u=\partial_u^2\log q_t(x-tY)\)。准确有
\[
 \sigma_u=\sum_j a_j\frac{-2w_j}{t^2+w_j^2},\qquad
 \tau_u=\sum_j a_j^2\frac{2(w_j^2-t^2)}{(t^2+w_j^2)^2}. \tag{5}
\]
标量函数 \(2(r-1)/(1+r)^2\;(r\ge0)\) 的最小值为 \(-2\)（\(r=0\)），
最大值为 \(1/4\)（\(r=3\)）；且 \(2|w|/(t^2+w^2)\le1/t\)。故
\[
 |\sigma_u|\le\|a\|_1/t,\quad
 -2\|a\|_2^2/t^2\le\tau_u\le\|a\|_2^2/(4t^2).
\]
同时对所有空间平移，
\[
 \begin{split}
 0<q_t&\le\frac1{|\det N|\pi^5t^5},\\
 |\partial_u q_t|&\le\frac{\|a\|_1}{t}q_t,\\
 |\partial_u^2q_t|&\le
 \frac{\|a\|_1^2+2\|a\|_2^2}{t^2}q_t. \tag{6}
 \end{split}
\]
混合二阶导数由极化也得到与空间位置无关的 \(O(t^{-7})\) 上界。
因此在 \([\epsilon,\infty)\times\Omega\) 上，密度、任意固定一阶及二阶空间导数分别由
常数倍 \(t^{-5},t^{-6},t^{-7}\) 可积支配。逐点导数连续，支配收敛给出
\(H_\epsilon\in C^2(\mathbb R^5)\) 及两次积分下求导；实际上该方法还给更高阶光滑性，
但本轮只用二阶。密度处处严格正，故 \(0<H_\epsilon(x)<\infty\) 对所有 \(x\) 成立。
特别地，此结论不受完整 Green 在特殊子空间可能发散的影响。

在正的、原本可能为无限总测度的 \(dt\,d\mathbb P(Y)\) 上，定义概率测度
\[
 d\nu_x(t,Y)=\frac{q_t(x-tY)}{H_\epsilon(x)}\,
             \mathbf1_{t\ge\epsilon}\,dt\,d\mathbb P(Y).
\]
上面的支配保证所有下列期望有限，并直接相减得到
\[
 u^T\nabla^2\log H_\epsilon(x)u
   =\mathbb E_{\nu_x}\tau_u+
      \operatorname{Var}_{\nu_x}(\sigma_u). \tag{7}
\]
这里是有界 score 的方差，不是重尾变量 \(Y\) 的方差。
下界丢弃非负方差；上界用
\(\operatorname{Var}\sigma_u\le\mathbb E\sigma_u^2\le\|a\|_1^2/\epsilon^2\)。
再用 \(\|a\|_1^2\le5\|a\|_2^2\)，即得 (1)。所有常数及方向的符号均已核对。

## 2. 全球 log 抛物面及固定 body 基

沿从 \(x_0\) 到 \(x_0+u\) 的线段积分 (1)，
\(\int_0^1(1-s)ds=1/2\)，得
\[
 \begin{split}
 \log H_\epsilon(x_0)+\nabla\log H_\epsilon(x_0)\cdot u
     -\frac{\|N^{-1}u\|_2^2}{\epsilon^2}
 &\le\log H_\epsilon(x_0+u)\\
 &\le\log H_\epsilon(x_0)+\nabla\log H_\epsilon(x_0)\cdot u
     +\frac{21\|N^{-1}u\|_2^2}{8\epsilon^2}. \tag{8}
 \end{split}
\]
这对任意长度的 \(u\) 有效，没有邻域半径、小时间或投影间隔假设。

本轮唯一固定矩阵的列为
\[
 (1,1,1,1,1),\ (1,-1,-1,-1,1),\ (1,-1,-1,1,-1),\
 (1,-1,1,-1,-1),\ (1,1,-1,-1,-1).
\]
准确符号计算给 \(\det N=-48\)，
\[
 N^{-1}=\frac16
 \begin{pmatrix}
 2&1&1&1&1\\1&-1&-1&-1&2\\1&-1&-1&2&-1\\
 1&-1&2&-1&-1\\1&2&-1&-1&-1
 \end{pmatrix}.
\]
对 \(u=(0,v)\)，准确有
\[
 \|N^{-1}u\|_2^2=v^TQv,\quad
 Q=I_4/4-J_4/36,\quad Q^{-1}=4I_4+(4/5)J_4. \tag{9}
\]
全部16个符号角点已准确枚举，\(\max_{\sigma\in\{-1,1\}^4}\sigma^TQ\sigma=1\)。
没有搜索或优化基 \(N\) 或 \(\epsilon\)。

数值计算得到的 \(g_j=\partial_{j+1}\log H_{1/64}(x_0)\) 为
\[
 \begin{split}
 g_1&=-1.28222064418038341573090972410139\ldots,\\
 g_2&=-1.24625287059888278355626342682936\ldots,\\
 g_3&=-0.21893608941524175236226549974291\ldots,\\
 g_4&=\phantom{-}0.19046758847326542724734917272469\ldots.
 \end{split}
\]
这四个数有384-bit Arb 球和准确有向 dyadic 端点；展示小数不是证书。
若对任意 \(v\) 使用区间梯度，线性项下界必须为
\(\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)\)，上界取相应最大值。
不能把梯度球的中点当作准确梯度，从而声称单个中点抛物面在全球成立。

## 3. 唯一中心计算：准确的一阶 log stencil

只读来源为017已独审的280真顶点/3504单纯锥表及018冻结源码、结果 JSON。
没有导入或执行它们的生产器，也未重新三角剖分。
表给 \(z_j=144r_j\in\mathbb Z^5\) 和 \(d_C=|\det(z_1,\ldots,z_5)|\)。
令 \(D=104729\)、\(k_j=Dx_0\cdot z_j\)、\(l_j=Dz_{j,2:5}\)；时间变量
\(s=144Dt\) 给出
\[
 H_\epsilon(x_0+(0,v))=
 F\operatorname{Re}\sum_Cd_C\int_e^\infty
    \prod_{j=1}^5(s-i k_j-i l_j\cdot v)^{-1}ds,
 \quad F=\frac{D^4}{144(2\pi)^5},\quad e=144D\epsilon=\frac{942561}{4}. \tag{10}
\]
积分绝对收敛：有限个锥在 \(s\ge e>0\) 上无零分母，尾部为 \(O(s^{-5})\)。
一阶空间导数也可直接用同类支配交换。

中心280个整数投影准确两两不同，因此本轮中心的一阶系数可用普通 DD 递推；
逐锥时间积分为
\[
 -\operatorname{DD}_4 f_e(k_1,\ldots,k_5),\qquad
 f_e(k)=\tfrac12\log(e^2+k^2). \tag{11}
\]
写 \(e=A/B\)，则
\(f_e(k)=\frac12\log(B^2k^2+A^2)-\log B\)。共同常数 \(-\log B\)
的四阶差商及所有一阶空间导数恒为零，故可准确删去。
每个初始节点的一阶系数准确为
\[
 \partial_{v_j}f_e(k+l\cdot v)|_{v=0}
      =\frac{l_j B^2k}{B^2k^2+A^2}. \tag{12}
\]
若差商的一步为 \(R=P/(d+\delta\cdot v)\)，只提取零阶和一阶，准确递推为
\[
 R_0=P_0/d,\qquad R_j=(P_j-\delta_jR_0)/d. \tag{13}
\]
全部计算用 `fmpq` 有理数；先把3504锥的有符号系数合并成五个完整 log stencil，
才以384-bit Arb 评估对数及 \(F\)。最终704个非零有理系数保存在外置 checkpoint。
此处“准确有理”指 stencil 系数；包含对数和 \(\pi^{-5}\) 的最终值及梯度用球包围。

由于 \(e>0\)，\(f_e\) 在整个实轴光滑，时间 resolvent 也没有实空间奇点。
若别处投影为零或重合，(4)--(8) 仍然适用；(11) 需使用连续/合流差商，不能直接除以零。
本程序没有在这些别处执行普通 DD，而是只在已验证互异的中心求一阶系数，
再用正混合定理控制全部空间。没有把“节点重合”误认为本方法的障碍或已被普通 DD 计算。

## 4. 闭盒、最大保守安全立方体及椭球的准确验收

令 \([H_L,H_U]\) 为中心球的有向 dyadic 包围。
对16个符号向量 \(\sigma\)，从梯度球得到
\(b_{\sigma,L}\le g\cdot\sigma\le b_{\sigma,U}\)，并准确取
\(c_\sigma=\sigma^TQ\sigma>0\)。下 log 抛物面是凹函数，闭盒上的最小值可在角点取到；
上抛物面是凸函数，最大值亦可在角点取到。因此半宽 \(h_0=1/16384\) 的盒上可用
\[
 H_L\exp\min_\sigma
   \{b_{\sigma,L}h_0-c_\sigma h_0^2/\epsilon^2\}
 \le H_\epsilon\le
 H_U\exp\max_\sigma
   \{b_{\sigma,U}h_0+(21/8)c_\sigma h_0^2/\epsilon^2\}. \tag{14}
\]
凹函数的角点结论先对每个真实梯度使用，随后按角点作有向区间下界；
没有把区间中点或仅16个函数样本当作整盒证明。
实际 Arb 比较严格通过两个登记门槛。

对主门槛 \(g_0=1/10000\)，准确取
\[
 \Delta_L\le\log(H_L/g_0),\qquad \Delta_L>0,
\]
其中 \(\Delta_L\) 为 Arb 下端点，约 \(0.18301579544235108\)。
每个角点的安全半宽条件为准确有理二次不等式
\[
 \Delta_L+b_{\sigma,L}h-c_\sigma h^2/\epsilon^2\ge0. \tag{15}
\]
取其16个正根的最小值，384-bit 包围给
\(h_*\approx0.0063862153390165525\)，最小根对应 \(\sigma=(1,1,-1,-1)\)。
这只是**固定保守抛物面测试**的最大半宽，不是 \(H\) 真实超水平集的最大立方体。
将根的下端点向下取至 \(2^{-24}\) 格点，得到前述
\(h=53571/8388608\)。然后对16个 (15) 准确有理代回，最小剩余量严格为正，
约 \(2.98682251704\cdot10^{-6}\)。正单位面裁剪用准确中心边距
\(\min_j\min(x_{0,j},1-x_{0,j})=614/104729\)，得到 \(h_{\rm face}\)。

对椭球，Cauchy--Schwarz 与 (9) 给
\[
 g\cdot v\ge-B\sqrt{v^TQv},\qquad
 B\ge\sqrt{g^TQ^{-1}g}.
\]
用梯度球及 \(Q^{-1}\) 计算的有向上端点 \(B_U\) 约为 \(4.284407733023408\)；
计算的是
\(4\sum_jg_j^2+(4/5)(\sum_jg_j)^2\) 的包围上端点再开方，不用中点。
令 \(\rho\) 为下式正根的向下 \(2^{-24}\) 格点：
\[
 \Delta_L-B_U\rho-\rho^2/\epsilon^2=0.
\]
得 \(\rho=51857/8388608\approx0.006181836128234863\)。准确有理代回的剩余量
严格为正，约 \(1.24786488934\cdot10^{-6}\)，从而得到闭椭球 (3) 上的严格门槛。
所有取整、二次型、面裁剪及剩余量均有准确分数记录在结果 JSON。

## 5. 实际运行、冻结来源与复现

登记后实际只运行一次新批：PID45264，全部3504锥完成，程序内墙钟约0.677745秒，
峰值 RSS76824576 bytes；python-flint0.9.0、SymPy1.14.0、384 bits。
没有触发120秒停止、没有重跑中心、没有 full \(G\) 或三阶 Taylor 基线重跑，
也没有 MATLAB、空间网格、第二个 \(\epsilon\) 或基选择优化。
018存储的中心 \(G\) 球仅作为只读数值一致性比较，未重新计算。
不能由这次局部批耗时推断全局覆盖的总成本。

冻结只读输入的 SHA-256：

| 来源 | SHA-256 |
|---|---|
| `astra_vertex017/boundary_triangulation.json.gz` | `8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9` |
| `astra_green018/kernel.py` | `5045f262ef0e00eb38746f1cf2293ee6cde138a6874e8f6850dc90d28b33ea08` |
| `astra_green018/batch_result.json` | `4eab4509995282e65c643c7423a523a9e4dbf02397060d185a0967ee5b6bd294` |

来源共同项目根是
`/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds`。
本轮代码运行前核对全部三份字节 hash。
外置完整 checkpoint：
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_harnack043/checkpoint.json`，
195466 bytes，SHA-256
`1da37f20a866551ac568558f57dffac2fe320e283aca17c7a624c3ffbe3904b6`。
checkpoint 锁定登记、当前源码及上游输入 fingerprint，每256个完整 simplex 保存；
若达到120秒则保存完整前缀、退出75，同源 `--resume` 从下个 simplex 续算。
本轮未实际触发续算，不声称已经压力测试中断恢复。

复现命令（新的仓库外工作目录；不要覆盖本轮已冻结 checkpoint）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/upper_harnack043/certify.py \
  --upstream-root /Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds \
  --work-dir /tmp/upper_harnack043-replay \
  --output-dir /tmp/upper_harnack043-replay-output
```

仅新增本目录 `registration.json`、`certify.py`、`results.json`、本报告四个文件；
未改017/018/035或他人目录，未提交。下一步若做全域工作，需要提出有限覆盖及各中心的
值/梯度认证账本；本轮的全球 Hessian 定理提供可传递的区域接口，但一个中心不能替代该覆盖。
