# R003：统一原子删除 Lipschitz 命题的有理多尺度反例

2026-11-05。固定 Lebesgue 接收面积、全部正半径、轴平行中心正方形。只写本独占 scratch，不修改既有证明、轮次材料或公开 clone。

## 1. 结论和准确量词

记 \(H_\lambda(\mu)=\{M_\square\mu>\lambda\}\)，全部水平集均严格。

**定理。** 不存在有限绝对常数 c，使所有有限正原子输入及其中任一质量 ε 的原子都满足
\[
 |H_1(\mu)\setminus H_1(\mu-\varepsilon\delta_p)|\le c\varepsilon.
 \tag{1}
\]
反例甚至可同时满足：总质量1、全部位置和质量为有理数、所有原子在固定线段 \([0,1/2]\times\{0\}\)、被删原子在左端0。

更具体地，对每个整数 N≥1，有 N+1 个正原子的概率输入 \(\mu_N\)，其中 \(\varepsilon_N\downarrow0\) 是原点质量，且
\[
 |H_1(\mu_N)\setminus H_1(\mu_N-\varepsilon_N\delta_0)|
       \ge\frac N{16}\varepsilon_N.
 \tag{2}
\]
同一组有理数据还给明确非负、紧支撑、有限矩形分片的概率密度 \(f_N=g_N+b_N\)，\(\|b_N\|_1=\varepsilon_N\)，满足
\[
 |H_1(f_N)\setminus H_1(g_N)|\ge\frac N{32}\varepsilon_N.
 \tag{3}
\]
所以一般 L1 的正扰动也没有统一同阈值线性稳定性。该 L1 结论由直接盒见证证明，不靠不稳定的原子平滑极限。

本轮重点是有限原子：每个固定配置的等号水平集面积为零（有限开／闭 receiver 之差只在边界），仍然出现统一于 N 的删除放大。一般 L1 的 strict-level 连续性早已可被等号平台打破；本稿 L1 例的额外性质是新增证书区在旧极大函数严格小于1的区域。

这些删除比值无界，不意味着全局 weak-(1,1) 常数无穷。所给证书面积下界仍趋于0；本轮不改 \(C_2\) 的数值区间。

## 2. 有理概率族：每个前缀都临界

令
\[
 q_N=4^{-N}/16,\qquad \alpha_N=1+q_N,\qquad
 a_i=\frac{2^{i-N-1}}{\alpha_N}\quad(1\le i\le N),\quad a_0=0,
\]
\[
 \varepsilon=\varepsilon_N=1-\alpha_N^{-2},\qquad
 m_i=4(a_i^2-a_{i-1}^2),\qquad
 \nu_N=\sum_{i=1}^Nm_i\delta_{(a_i,0)},\qquad
 \mu_N=\nu_N+\varepsilon\delta_0.
 \tag{4}
\]
全部数据有理且严格为正，\(a_{i+1}=2a_i\)、\(a_N<1/2\)。望远镜质量式为
\[
 \sum_{j\le i}m_j=4a_i^2,\qquad
 \nu_N(\mathbb R^2)=4a_N^2=\alpha_N^{-2}=1-\varepsilon,
 \quad\mu_N(\mathbb R^2)=1.
 \tag{5}
\]
同时
\[
 \frac\varepsilon{a_1^2}
 =\frac{q_N(2+q_N)}{4^{-N}}=\frac{2+q_N}{16}<1.
 \tag{6}
\]
因此所有 i 都有 \(\varepsilon\le a_i^2\)。没有先用一个未归一化族再含糊重标质量；(4)本身已经是概率输入。

## 3. 旧 full halo 在负半平面为空

取任意接收中心 \((u,v)\) 且 u≤0，任意 r>0。

若 \(|v|>r\)，该盒不见任何旧原子。否则，其横坐标区间的左端 \(u-r<0\)，所捕获的正原子若非空必为一个前缀 \(1,\ldots,j\)，且 \(r\ge a_j-u\)。于是
\[
 \nu_N(Q((u,v),r))=4a_j^2\le4r^2.
 \tag{7}
\]
当 u<0 时不等式严格。空前缀质量为零。故
\[
 H_1(\nu_N)\cap\{u\le0\}=\varnothing.
 \tag{8}
\]
这个证明检查了全部 r，不只某个候选盒族；垂直坐标与横向捕获使用同一个半径。

也可从精确子集 receiver 式核对：任意旧非空子集的最大指标为 j，其质量≤\(4a_j^2\)，外接盒最右端为 \(a_j\)，故 receiver 左端≥\(a_j-\sqrt{4a_j^2}/2=0\)。非前缀子集不会生成隐藏的左侧 halo。

## 4. 互不重叠的新增条带与严格有理见证

定义
\[
 t_i=\frac\varepsilon{16a_i},\qquad \rho_i=a_i+t_i,
\]
\[
 S_i=(-t_i,0)\times\{v:a_i/2<|v|<a_i\}.
 \tag{9}
\]
因为相邻 a 加倍，所有纵向环带互不相交，故各 \(S_i\) 两两不交。

对任意 \((u,v)\in S_i\)，同一个正方形 \(Q((u,v),\rho_i)\) 捕获原点及全部旧前缀 \(1,\ldots,i\)：横向最远距离小于 \(a_i+t_i\)，垂直距离小于 \(a_i\)。其平均严格超过1，因为
\[
 \begin{aligned}
 4\rho_i^2
 &=4a_i^2+\frac\varepsilon2+\frac{\varepsilon^2}{64a_i^2}
 \le4a_i^2+\frac{33}{64}\varepsilon
 <4a_i^2+\varepsilon.
 \end{aligned}
 \tag{10}
\]
所有见证半径都为有理数；最少严格余量为 \(31\varepsilon/64\)。若盒还见到其他质量，只会加强见证。

结合(8)，每条带都位于真正删除损失内。其面积为
\[
 |S_i|=t_i a_i=\varepsilon/16.
\]
求和即得(2)。这是真二维面积的累积：每一尺度的横向宽度约 ε/a，纵向高度约 a；把问题当成一维区间会丢掉后一个因子。

旧固定 N 删除引理并未被反驳。对这类左端点／共线配置，含被删原点的坐标盒只有 N+1 个，旧单盒一-Lipschitz证明给 \((N+1)\varepsilon\) 上界；(2)表明其依原子数的线性增长在数量级上无法完全消除。

## 5. 明确的 L1 正扰动反例

仍使用(4)的 a、m、ε，取
\[
 h=t_N/4,\qquad
 B_i=[a_i,a_i+h]\times[-h,h],\qquad
 g_N=\sum_{i=1}^N\frac{m_i}{2h^2}\mathbf1_{B_i},
\]
\[
 b_N=\frac\varepsilon{4h^2}\mathbf1_{[-h,h]^2},\qquad f_N=g_N+b_N.
 \tag{11}
\]
全部系数和矩形端点有理，\(\|g_N\|_1=1-\varepsilon\)、\(\|b_N\|_1=\varepsilon\)、\(\|f_N\|_1=1\)。支持都在固定 \([-1,1]^2\) 中。事实上 h<1/128，统一平移 \((1/4,1/2)\) 后所有这些函数的支持都在 \([0,1]^2\)。

### 5.1 旧密度的全部尺度排除

令横向累计质量
\[
 F(s)=\int_{\{x\le s\}}g_N(x,y)\,dx\,dy.
\]
s≤0时为零。s>0时，令 j 为满足 \(a_j\le s\) 的最大指标（无则质量零），只有已经开始的块可能贡献，故
\[
 F(s)\le\sum_{i\le j}m_i=4a_j^2\le4s^2.
 \tag{12}
\]
不需要块之间不交，也没有把块内质量假装成端点原子。

若接收中心 u≤0，则任意 r 都有
\[
 \int_{Q((u,v),r)}g_N\le F(u+r)\le4r^2,
\]
其中 u+r≤0时左边为零。于是 \(H_1(g_N)\) 在左半平面仍为空。

更精确地，令旧总质量 \(M_0=1-\varepsilon\)。对固定 u=−d<0，平均还至多 \(M_0/(4r^2)\)，从而
\[
 M_\square g_N(u,v)\le\sup_{r>d}
 \min\{(1-d/r)^2,M_0/(4r^2)\}
 =\frac{M_0}{(2d+\sqrt{M_0})^2}<1.
 \tag{12a}
\]
两项分别单调增、减，其交点为 \(r=d+\sqrt{M_0}/2\)。同一上界也适用于旧原子族。故下面的新增条带没有使用一个旧的等号平台。

### 5.2 同一有理盒捕获真实块与 bump

取缩半条带
\[
 T_i=(-t_i/2,0)\times\{v:a_i/2<|v|<a_i\}.
 \tag{13}
\]
它们互不交。由于 \(t_N\le t_i\)，有 h≤\(t_i/4\)。对任意 \((u,v)\in T_i\)，半径仍取 \(\rho_i=a_i+t_i\)。

全部 j≤i 的源块横向最远距离≤\(a_i+h+t_i/2<\rho_i\)，垂直最远距离≤\(a_i+h<\rho_i\)。origin bump 的两个方向也满足这些界。因此该同一个正方形完整捕获前缀块质量 \(4a_i^2\) 和 bump 质量 ε。式(10)给严格平均>1。

每层 \(|T_i|=\varepsilon/32\)，全部处于新 halo 与旧 halo 的差中，得(3)。这里没有诉诸“平滑后 halo 差自动收敛”，更不需要从两个分别 Fatou 的面积不等式相减。

## 6. 至少 ε log(1/ε) 的同阈值障碍

对 q≤1/64，直接比较
\[
 q\le\frac{q(2+q)}{(1+q)^2}=\varepsilon_N\le2q,
 \qquad\varepsilon_N\sim4^{-N}/8.
 \tag{14}
\]
因此 \(\log(1/\varepsilon_N)\le(N+2)\log4\)。N≥2时，
\[
 \begin{aligned}
 |H_1(\mu_N)\setminus H_1(\nu_N)|
 &\ge\frac{\varepsilon_N\log(1/\varepsilon_N)}{32\log4},\\
 |H_1(f_N)\setminus H_1(g_N)|
 &\ge\frac{\varepsilon_N\log(1/\varepsilon_N)}{64\log4}.
 \end{aligned}
 \tag{15}
\]
这只是沿明确原子序列的下界，不判定有限原子类的最优统一模量。对一般 L1，等号平台本已排除全类的同严格阈值统一连续性；因此(15)不能读成一般 L1 的 sharp 阶，更不能据此建议一个统一 εlog(1/ε) 上界。

## 7. 正确的统一替代：保留阈值松弛层

**命题（单原子删除）。** 对任意有限正测度 ν、任意 p、ε>0、λ>0及 \(0<\theta<1\)，
\[
 H_\lambda(\nu+\varepsilon\delta_p)
 \subset H_{(1-\theta)\lambda}(\nu)\cup H_{\theta\lambda}(\varepsilon\delta_p),
\]
\[
 \boxed{\lambda\left|H_\lambda(\nu+\varepsilon\delta_p)
                  \setminus H_{(1-\theta)\lambda}(\nu)\right|
       \le\varepsilon/\theta.}
 \tag{16}
\]

证明：对一个不在第一个 halo 的接收点，任一原严格见证都使 ν 平均≤\((1-\theta)\lambda\)，所以它必含 p，并使原子平均>\(\theta\lambda\)。仍是同一个正方形与半径。单原子的 strict halo 是以 p 为中心、半径 \(\sqrt{\varepsilon/(4\theta\lambda)}\) 的开正方形，其面积恰为 \(\varepsilon/(\theta\lambda)\)，得式(16)。

因而真正同阈值删除差满足
\[
 \lambda|H_\lambda(\nu+\varepsilon\delta_p)\setminus H_\lambda(\nu)|
 \le\frac\varepsilon\theta+
 \lambda|H_{(1-\theta)\lambda}(\nu)\setminus H_\lambda(\nu)|.
 \tag{17}
\]
第二项是必须另控的 near-critical layer，不能由全局 weak 界擅自改成 \(C\theta\|\nu\|\)。这是替代估计的真实缺口。

若删除的是一般正测度 σ（或 L1 bump）且总质量 ε，同一见证论证给(16)的集合包含式，把第二项换成 \(H_{\theta\lambda}(\sigma)\)。采用已独立证明的非锐 weak 界9，就得面积误差 \(9\varepsilon/(\theta\lambda)\)。单原子版本的系数1来自其精确 halo，不是把已知全局常数重包装成进展。

## 8. 对最小反例能说和不能说什么

固定 K≥1。若 µ 是一个原子数最少的 K-反例、总质量1，且其实际超额为
\(\gamma=|H_1(\mu)|-K>0\)，删去任一质量 ε 原子所得 ν 在**全部阈值**均满足 K 界。否则由缩放便有更少原子的 K-反例。

由(16)，
\[
 K+\gamma\le\frac{K(1-\varepsilon)}{1-\theta}+\frac\varepsilon\theta.
\]
Cauchy 不等式及取 \(\theta=\sqrt\varepsilon/(\sqrt{K(1-\varepsilon)}+\sqrt\varepsilon)\) 给右侧关于 θ 的最小值
\((\sqrt{K(1-\varepsilon)}+\sqrt\varepsilon)^2\)，所以
\[
 \gamma\le2\sqrt{K\varepsilon(1-\varepsilon)}-(K-1)\varepsilon
       \le2\sqrt{K\varepsilon},\quad
 \varepsilon\ge\gamma^2/(4K).
 \tag{18}
\]
这是给定这个反例**自身超额 γ**的条件性质量下界。最小 K-反例的超额可能任意小，不能把另外一个更大比值输入的超额移给它。因此(18)不提供全局统一 N 近似，也不闭合任何候选锐上界。

## 9. 本轮分类与计算范围

| 类别 | 准确内容 |
|---|---|
| 已证 | 有理概率原子族(4)反驳绝对删除c；全尺度旧 halo 排除与有理严格新见证。 |
| 已证 | 有理矩形 L1 概率族(11)反驳一般正扰动的统一线性同阈值稳定性。 |
| 已证 | 沿该序列至少 εlog(1/ε) 的损失下界；阈值松弛(16)–(17)。 |
| 仅条件推论 | 最小反例自身超额 γ 对应的质量下界(18)。 |
| 未宣称 | 最优同阈值模量、全局原子数尾界、锐weak常数、极值达到性。 |

已实际顺序运行本机标准库 Fraction 脚本 `check_family.py`，原stdout完整保存在同目录。N=1、2、3、8、16、32、64的质量归一化、严格见证余量、source块捕获及两组strip精确并面积全部通过；N≤8另共266个非空子集检查旧receiver不越左半面。`sample-N8.json`保留一个完全有理的输入与证书。

没有调用或占用全核积分计算，没有把有限样本作为任意N证明。算术只核输入和证书，全部尺度旧halo排除及无界比值由解析论证承担。本轮没有把任何接收证书说成完整halo面积，也未声称一般L1的最优连续性模量。
