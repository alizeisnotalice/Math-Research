# upper_multibasis056：16 个合法 body 基的分象限区域证书

## 结论

新的多基曲率公式与 81 点闭盒判据成立。只读冻结的 052 batch0002，在其中 **162 个原失败 coarse-root 父盒中，162 个全部通过新证书，0 个仍未通过**。这 162 个 anchor 均严格位于各自闭盒内部。最小新 log 余量具有准确有理下界

\[
 \min_B m_{\mathrm{new}}(B)\ge\frac{293491}{1000000}>0.
\]

完整的 162 条 margin 证书位于 `results.json`，包括 anchor、盒端点、保存的 Delta 和梯度有向区间、旧余量、新余量及最差检查点。其最小新余量准确为

\[
 \frac{33983968010587124693717903513538128433237870479571267959365256633789859952279}
 {115792089237316195423570985008687907853269984665640564039457584007913129639936}.
\]

唯一实际批 PID **53688**，0.6074 秒，单进程单线程，峰值 RSS 288309248 bytes。检查 2592 个旧角点和 13122 个新分片顶点，全部采用 Fraction；16 个基的身份用 SymPy。**新 H、gradient、Green、kernel、log、exp、Gamma、Arb 调用均为 0**。

这是固定 162 个父盒的区域证书，不是全域 g 或 Cn 结论。本轮没有修改 052 的 state、cache、queue、controller、STOP，也没有重跑当前树。对后续批边界是否采用或合并该机制，仍由主任务负责。

## 1. 冻结输入、资格与读取范围

只选取
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover052/batch0002/state.json`
的 `splits`；没有取任何更新 batch。该文件运行前后 SHA-256 都是
`b79c134fdcec283c120d464cf78f131c669c905ac5cdacda309a34f07a9fb29b`。

同目录 `registration.json` 的 SHA-256 为
`585a5de78b8f159ca03f138013894041dba350d144c3217b7bbda09fd5c7c9e9`。
准确重建并核对状态 fingerprint：
`d2067c9ddfe6baba298ebe7bc9292bab91f26630774a49bfb69cb5088b2c9229`。

固定模型是 alpha=1 的五维 16 条等权独立 body Cauchy 方向，
\(H_\epsilon(x)=\int_\epsilon^\infty p_t(x)dt\)，epsilon=1/64，目标 g0=1/10000；数字中心数据只有 x2,x3,x4,x5 四个物理梯度，所有区域均在 x1=1。

采用的既有理论原件：

- `/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_harnack043/report.md`，SHA `a298ca56ba980fce4fe2542929e03befb97b528afd3d390ed9d15651da68322a`。
- 独审 `/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/audit_harnack047/report.md`，SHA `0b680c1bd7e4aea30c027fe32cc3e51436f7dda6dd8b74ddf521c2d74bb8ca16`。
- 原盒/证书算法 `Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover046/arb_cover.py`，SHA `df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028`；仅读源码，不导入或调用。

状态恰有 162 个 split，全部 node path 为空，即原 coarse root，不含其后代。状态内已经保存这些失败测试的完整 `Delta_lower` 与 `loggradient_directed`，因此不需要重新求任何超越函数。其来源仍是冻结 046 算法对准确 H 下端点作 Arb log 并取向下端点；本轮沿用该认证，不把展示小数当输入。

对每条记录，本轮另外核对：缓存五元区间有序、H 下端点严格大于目标、四个保存 log-gradient 区间包住所有准确导数/H 端点比、anchor key 一致、anchor 在盒内、原 16 点失败 margin 与保存值完全相同。生产核函数及其归一化资格沿已锁定的 043/046/047/052 输入；本轮不重新生产这些区间。

## 2. 同一个 H 的 16 个基不等式

043/047 已证明，对从 16 条方向中选出的任意可逆五列基 B、任意 x0,u，有

\[
 \log H_\epsilon(x_0+u)\ge
 \log H_\epsilon(x_0)+\nabla\log H_\epsilon(x_0)\cdot u
 -\frac{\|B^{-1}u\|_2^2}{\epsilon^2}. \tag{1}
\]

这是正 Cauchy 卷积与有界 score 方差给出的全球 log-Hessian 下界的积分形式，不依赖空间投影间距，也不需要剩余重尾变量的矩。

取 043 的固定基 N，其列为

\[
 (1,1,1,1,1),\ (1,-1,-1,-1,1),\ (1,-1,-1,1,-1),\
 (1,-1,1,-1,-1),\ (1,1,-1,-1,-1).
\]

对每个 s in {−1,1}^4，令 D_s=diag(1,s)、N_s=D_s N。尾坐标改号置换 16 条 body 方向，所以 N_s 的五列仍是同一模型中的五条不同方向；它可逆且 |det N_s|=48。以这些五个变量重新分组只改变正卷积证明中所选的基和剩余变量，**不改变 H_epsilon、中心 H 值或中心梯度**。因此 (1) 的 16 个版本可同时用于同一个 H。

令 u=(0,v)。因为 N_s^{-1}=N^{-1}D_s，限制二次型是

\[
 \|N_s^{-1}(0,v)\|_2^2=v^TQ_sv,\qquad
 Q_s=I_4/4-ss^T/36. \tag{2}
\]

每个 Q_s 特征值均为 1/4（三重）与 5/36（一重），严格正定。这里的基只取预先指定的 16 个改号版本，没有另做优化或搜索。

## 3. 逐点最大与准确曲率

对所有 s 的合法下界取最大，不需要对最大函数求导。共同中心项保持不变；二次罚项取最小：

\[
 \min_s v^TQ_sv
 =\frac{\|v\|_2^2}{4}-\frac{\max_s(s\cdot v)^2}{36}
 =\boxed{\psi(v):=\frac{\|v\|_2^2}{4}-\frac{\|v\|_1^2}{36}}. \tag{3}
\]

最后一个等号来自三角不等式 |s dot v|<=sum|v_j|，并由 s_j=sign(v_j) 达到；为零的坐标任选正负都不影响等号。故对所有 v，

\[
 \boxed{\log H_\epsilon(x_0+(0,v))\ge
 \log H_\epsilon(x_0)+g\cdot v-\psi(v)/\epsilon^2,}
 \quad g=\nabla_{x_{2:5}}\log H_\epsilon(x_0). \tag{4}
\]

它不要求 v 的符号预先固定，也不要求 anchor 的新核函数数据。原 Q=I/4−11^T/36 的罚项不小于 psi；准确改进量为

\[
 v^TQv-\psi(v)
 =\frac{\|v\|_1^2-(\sum_jv_j)^2}{36}
 =\frac19\sum_{i<j,\,v_iv_j<0}|v_iv_j|\ge0. \tag{5}
\]

这解释了混合正负位移方向上的收益。所有同号位移上可能完全没有收益，不能宣称每个点都有统一严格改善。

由 \(\|v\|_2^2\le\|v\|_1^2\le4\|v\|_2^2\)，还有

\[
 \frac5{36}\|v\|_2^2\le\psi(v)\le\frac29\|v\|_2^2. \tag{6}
\]

所以轴上与坐标变号处没有负曲率罚项或发散；v=0 时 psi=0，各基一致。分段表达式可能不可微，但 H_epsilon 本身仍是平滑正函数；(4) 是逐点不等式，不从拼接函数的 Hessian 推导 H 的新 Hessian。

## 4. 为什么需要分象限，而不能泛用原 16 角点

−psi 不是全球凹函数。例如在 (±1,5,0,0) 处取平均，两端值均为 −11/2，中点 (0,5,0,0) 的值为 −50/9，低于平均值，差为 −1/18。

更直接地，在包含原点的盒
\([-1,1]\times[-5,5]\times[-1,1]^2\) 上，16 个外角点的 −psi 都等于 −47/9，但轴向检查点 (0,5,0,0) 的值是 −50/9。若取零线性项与 Delta=16/3，外角点 margin 全是 1/9>0，轴向 margin 却是 −2/9<0。这是对“只查外角点即可认证该下包围函数”的代数反例，不是对实际 H 模型作反例声明。按任意正尺度缩小即可得到任意小的同类几何例。

正确做法是对每个闭正交象限分别处理。在符号为 s 的象限内，sum|v_j|=s dot v，所以 psi(v)=v^TQ_s v。对固定真实 g，(4) 的右边在该象限是严格凹二次式。

若梯度只有区间 g_j in [g_j,L,g_j,U]，安全线性下界是

\[
 L(v)=\sum_j\min(g_{j,L}v_j,g_{j,U}v_j).
\]

在每个闭象限内，各坐标所取端点固定，因此 L 是线性的；在 v_j=0 上两种端点都给零，连续拼接。令

\[
 \Phi(v)=\Delta_L+L(v)-\psi(v)/\epsilon^2,
 \qquad \Delta_L\le\log(H_L/g_0).
\]

Phi 在每个闭象限内是凹二次式，且所有 v 都满足 log(H_epsilon(x0+(0,v))/g0)>=Phi(v)。

## 5. 81 点闭盒定理及退化边界

设实际盒为 B=product[lo_j,hi_j]，anchor tail a 位于 B。写 l_j=lo_j−a_j<=0、u_j=hi_j−a_j>=0。每个坐标区间分成 [l_j,0] 与 [0,u_j]，四维得到至多 16 个闭象限矩形。每一矩形上的凹函数 Phi 最小值可以在其顶点取到，全部顶点的并集恰是

\[
 \mathcal V_B=\prod_{j=1}^4\{l_j,0,u_j\},\qquad |\mathcal V_B|\le3^4=81. \tag{7}
\]

所以准确充分条件为

\[
 \boxed{\min_{v\in\mathcal V_B}\Phi(v)>0\quad\Longrightarrow\quad
 H_\epsilon(1,x)>g_0\ \text{for every }x\in B.} \tag{8}
\]

等号 >=0 也可给非严格 >=g0；本轮只接受严格正余量。若 anchor 在边界、坐标宽度为零或多个检查点重合，去重后取更少点仍覆盖每个闭分片的顶点。零坐标的轴、边、面均已包括；不可省略。若 anchor 不在盒内，则不应直接套“每坐标三点”的这个版本，应另按实际相交象限处理。本轮首先逐盒准确确认 anchor 在盒内，实际上全部严格在内，因此每盒恰为 81 个不同点。

该证书在整个原盒上成立，即便盒内跨越原始投影重合面、S1/S2 或排序边界。正 H_epsilon 的全球定理覆盖这些位置；完整 Green 可能为无穷并不损害 G>=H_epsilon 的下界传递。

## 6. 唯一有限批的实际证据

登记后只读取指定 state 的 162 个 `splits`。原 root 的四坐标索引由已核对摘要的 roots 表取得，各闭盒边长为 1/32。逐盒使用保存的同一 Delta 和同一四个 gradient interval，对以下两者作准确比较：

1. 原 Q 的 16 个外角点；其最小值和最差 bits 与原失败证书逐条完全一致，均 <=0。
2. 新 psi 的每坐标“低端、零、高端”笛卡尔积；全部 13122 点都以 Fraction 计算，逐点核对 (6) 及 psi<=v^TQv。

所有 162 个新最小余量均为正，且逐盒均不小于旧最小余量。原点/边界检查点没有被抽样忽略；本批最差点最终都落在外角点，这是运行结果，不是省略其余 65 点的通用许可。

缓存 five-tuple 仅被读取。保存 gradient 区间对缓存 H/gradient 的准确区间除法包围逐条通过；没有用缓存中点。Delta 的已有向下认证沿用冻结来源；不存在本轮新增超越值或数值误差估计。

本轮最小新余量的精确分数和每条证书都在结果 JSON，正文的 293491/1000000 来自该分数准确向下取整，非展示小数转证书。失败测试本来只说明旧曲率包围未通过，不是 H<g0；新方法通过与它不矛盾。

已完成 16 基身份检查、162 个 anchor 包含性检查、2592 个旧点复核、13122 个新点计算。未选择新 anchor，未扫描更多父盒或空间参数，未读取新 batch 做增量扩展。

## 7. 可复核文件与锁定

`registration.json` 固定输入字节 hash、状态 fingerprint、源码 hash、162-box 范围及 120 秒可恢复时限。`requalify.py` 不导入 046/052 模块，也不包含 kernel/H/gradient 求值；其唯一数学依赖是 SymPy 的 16 个小矩阵身份及 Python Fraction。

每完成 16 盒写入可恢复前缀；实际一次完成，未触发 120 秒停止或续算。最终 checkpoint 仅保留运行元数据，不重复完整 margin 清单。峰值内存来自只读解析 52 MB 冻结 state；没有新大型缓存或复制整个 state。

锁后复核请使用新输出目录：

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/upper_multibasis056/requalify.py \
--output-dir /absolute/path/to/new-056-output
```

完整交付仅为本报告、登记、最小验证器、162 证书结果、实际收据和 checkpoint；各 hash 见 manifest.json。旧 047 及其他锁保持不变。没有提交、推送或 live queue 更改；162 个区域可供独立审计及之后经授权的批边界整合，本轮不作整合动作。
