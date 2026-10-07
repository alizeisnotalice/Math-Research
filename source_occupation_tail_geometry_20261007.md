# 真 centered-cube winner 的 deterministic 来源势：交集合同、截断尾与当前压力范围

日期 2026-10-07；沿用指定 6.1-sol high。专属 prefix source_occupation_tail_geometry_20261007；没有改主账或旧实验。使用 E01 的同一初始来源/占用质量范围守卫，已读其 SKILL、method、provenance、cube-interface；这里只用直接正 Fubini 和真实 cube 几何，不引用占用 LP 强对偶或免费停止定理。

**终态判断：没有证明平方根二阶/截断尾合同，也没有发现反驳它的真实来源反例。** 已取得的是精确联合核、已知 n 级粗界，以及有记录的三轮原空间探路。33 个数值阈值记录均非空，但全部 I1 置信下界为 0；这些数据不能支持任何维数阶数。shared-U 的随机 count 坏尾与本文 deterministic S 不同，未被当作其反例。

## 1. 检索范围与真实合同

已检索本目录 occupation/来源二阶/尾/列负载，阅读 [source_optimal_acceptance_coupling](source_optimal_acceptance_coupling_20261007.md)、[acceptance_count_projection](acceptance_count_projection_20261007.md) 的 coupling 边界、[gated_occupancy_bridge](gated_occupancy_bridge.md) 的 square 校准范围，以及下界总表 B38/B41/B62 与 A13/A14 范围。现有 source-only balanced coupling 使最小 overflow 等于 integral(S-K)+dmu，但它没有降低 S；shared-field/Palm count 反例也未构造同一个 deterministic S 的大尾。未找到现有笔记中的原 centered-cube 二阶反例；这不是全体文献不存在反例的声明。

固定完整非负有限来源 mu，W=mu(R^n)>0，有限物理节点 a<=L1<...<LJ<=2a。h_L=L^{-n}1_{[-L/2,L/2]^n}。完整响应 U_j=h_Lj*mu，M=max_j U_j，E={M>tau}，Aj 是 E 上的真正最大赢家分区，数学 ties 取最小 index。保持完整来源，定义

\[
 g(x)=\frac{\tau}{M(x)}1_E(x),\quad
 S(y)=\sum_j\int_{A_j}g(x)h_{L_j}(x-y)dx. \tag{1}
\]

本研究对象没有任意 receiver/source mask；不能另选有利 Aj 或先分组件再分别取 max。本文的有限 J 是原有限目录合同本身，**不是**连续 [a,2a] 最大值的数值认证。实际原 FIRST/history 的 joint kernel 若要归约到 (1)，还需其正共同参考及原门的接口；本次纯 hard 测试不声称已经做完该归约。

## 2. 精确一阶、二阶 cube 交集身份

对 x in E，令 Qx=x+[-R(x)/2,R(x)/2]^n，m_x=mu(Qx)，R=L_winner。因为 M=m_x/R^n，

\[
 S(y)=\tau\int_E\frac{1_{Q_x}(y)}{m_x}dx,
 \quad I_1:=\int S\,d\mu=\tau|E|. \tag{2}
\]

同一完整来源的二阶准确为

\[
 I_2:=\int S(y)^2\mu(dy)
 =\tau^2\int_{E\times E}
 \frac{\mu(Q_x\cap Q_{x'})}{m_xm_{x'}}dx\,dx'. \tag{3}
\]

这是完整真实两接收点 cube 交集质量；用 Lebesgue intersection volume 替代 mu(QcapQ') 需要另一个 density 假设。正方形/Gram 身份不自行支付其大小。对一般源，质量可以集中在一个很薄的交集中，单用 Lebesgue 重叠体积不能限制联合质量。

真实 finite winner 的全部行约束是

\[
 \frac{\mu(x+Q_{L_k})}{L_k^n}\le\frac{m_x}{R(x)^n}=M(x),
 \qquad k=1,\ldots,J. \tag{4}
\]

因此对同一 x 的每个合法 L_k 有质量上界 m_x(L_k/R)^n。它保留了整个 finite 参数未来比较，没有把捕获后的来源概率当成无条件均匀。两个不同中心的 cube 并非该行的嵌套族；把 Qx 或 Qx' 包进另一个中心的扩大 cube，扩大 sidelength 可能不在原目录或物理窗口中。不能免费调用 (4)。原 soft fullfuture 行帽也不是 (3) 固定来源 y 的列重叠费。

现有正几何粗 bound mu(QcapQ')<=min(m_x,m_x') 仍含未知的双 receiver 积分。真正欠项是针对这些由同一 mu 诱导、满足 (4) 的实际 Qx，控制 (3) 或直接控制 (1) 的高来源层。不是证明 positive joint kernel 或一阶守恒即可闭合。

## 3. 目前可证明的 n 级费，与更弱的截断合同

因 0<=g<=1、Aj 不交，

\[
 S(y)\le\int\sup_{a\le L\le2a}h_L(x-y)dx
 =C_n:=1+n\log2. \tag{5}
\]

内 cube a 的质量为 1；外 shell 用 full sidelength r=2||x-y||infty 的 Jacobian n r^{n-1}dr，积分 n integral_a^{2a}dr/r。因此 I2<=C_n I1，而 K>=C_n 时 deterministic tail 为解析零。这是已有 fixed-window 强 envelope 的 n 级费，**不是**新 sqrt(n) 合同。区间有界不表示节点复杂度低。

待证平方合同是 I2<=B_n I1，目标 B_n=O(sqrt(n) polylog(n))。它确实足够：对 K>0，(s-K)+<=s²/(4K)，从而

\[
 \int(S-K)_+d\mu\le\frac{B_n}{4K}I_1. \tag{6}
\]

取 K=B_n/(4delta)，delta<1，再由 I1<=KW+integral(S-K)+ 得 I1<=K W/(1-delta)。平方合同强于只在一个 cutoff 的尾合同；若只知道 tail<=delta I1，加上 (5)，至多由 S²<=KS+C_n(S-K)+ 推出 I2<=(K+C_n delta)I1，constant delta 仍允许 n 级平方费。因此失败的强平方路线若将来出现，不会自动否定 sqrt(n) cutoff 的来源尾。

## 4. 一个必要边界：点态 sqrt(n) 不成立，但不是 source-tail 反例

真实 cube 而非 sensor toy 的简单模型能逼近 (5)。a=1，mu=tau 1_{[-2,2]^n}dx+epsilon delta0，取几何节点 R_j=2^{(j-1)/(J-1)}，epsilon>0。对 x in [-1,1]^n，各候选 cube 完全处于背景盒中；包含 0 的最小节点唯一最大，因为响应 tau+epsilon/R^n 随 R 严格下降。外面没有合法节点包含 0，所有背景平均<=tau，故 E=[-1,1]^n（边界对体积无关）。于是

\[
 S(0)=\frac{\tau}{\tau+\epsilon}
 +\sum_{j=2}^J\frac{\tau(R_j^n-R_{j-1}^n)}{\tau R_j^n+\epsilon}. \tag{7}
\]

epsilon/tau->0 时趋 1+(J-1)(1-2^{-n/(J-1)})，J/n->infinity 时为 1+nlog2+o(n)。例如 epsilon<=tau/100 的固定真实源与 J-1=n² 已给点态 order n。这没有把 winner 指定为任意 gate。固定 n,J 后，以窄正盒替换 delta0，除 finite arrival shell 的零测边界外，响应/赢家及来源窄盒平均势按 dominated convergence 收敛到此模型；所以同一边界可 L1 化。

但 I1=tau2^n，而该 atom 的 I2/I1 贡献最多 (epsilon/tau) C_n² 2^{-n}；该来源点的 overflow/I1 贡献也至多 (epsilon/tau) C_n2^{-n}。这个高 S 来源极轻，不是 (6) 或所求二阶的反例。这里没有计算背景其它来源的完整 I2，不登记强平方失败。它只提醒：目标必须保持 mu 权，不能先要求每一个 y 有 sqrt(n) 列帽；随机 field 的坏 count 同样不能代换 deterministic S。

## 5. 已注册三轮真实原空间输入与估计

registration 在执行前落盘，所有 source/alpha/K 固定，无搜索后调参。三轮 n4/16/64，J9/17/33，layers2/4/8；来源样本数64/128/256，每个来源两个条件独立 inner receiver replicas，各128/256/256点。另固定 n16、J33，仅测试 B41-inspired layers2 与8。不同轮和 same-n cases 用新 seed；同 case 各阈值复用 common random receiver 响应。新 seed 是这些输入上的独立随机重复，并非冻结同一候选的三次验证，也不是 independent asymptotic confirmation。

完整来源均质量 1，先正混合再计算 full winner：

* capacity_weighted_grid：B38-inspired 的受限 positive product grid mixture；每个 grid 的一维 q=2(16n)k+1 点，k=1,...,layers，正交替 masses1/3、归一化 2q-1，部分层 phase=.25。未求解 B38 的一般 capacity LP，不代表任意网格最优权重。
* overlap_resolution_mixture：B41-inspired 的完整 n 维 finite tensor grids 正混合，k=1,...,layers，正权 proportional l+1；全部 coincident atoms 的物理质量在响应中自动相加。其 resolution 不是总表 B41 的 q_r=2^r，未实现完整 B41/A13/A14 下界资格。
* concentric_box_mixture：B62-inspired 的完整 n 维均匀盒有限 log-width mixture，width 从 n^{-1/2} 到 n^{1/2}，positive weights proportional1/(l+1)。它不是 B62 的连续 log-scale integral，也没有代替 A13/A14 的抖动、Sidon 标签、峰宽补偿或同步增长参数。

finite grid 即使 q^n 极大，所有实际 atom 仍存在，响应由每坐标的闭 interval 整数 count/正权 count 乘积精确表示，不枚举或抽样截源。盒源使用真实 intersection 各坐标长度乘积。浮点边界和 response 比较依然不是精确实数认证；代码有 near-boundary、near-max、near-threshold 屏幕，不把 heuristic screen 当 interval 证明。finite J 用所有节点 argmax，最小 index 处理计算中的精确 ties；附近潜在 winner 的 contribution 范围另存。数学上的真正 exact tie convention 要由稳定性或精确算术另认证。

预设 tau：两个 grid families 是 alpha(32n)^{-n}，盒 family 是 alpha2^{-n}，alpha=.1,.5,2。K=.5sqrt(n),sqrt(n),2sqrt(n)。same-n 改 layers 还改变最大 k、混合权和完整输入；它保持 n，但不是隔离只有一个 nuisance parameter 的实验。

### receiver Lebesgue 与 source 平方的正确抽样

对真实 source point y，proposal density

\[
 q_y(x)=\frac{\sup_{1\le L\le2}h_L(x-y)}{C_n}
\]

是真正 receiver 密度。core cube1 的概率1/C_n；其它概率 nlog2/C_n，抽 log r uniformly in [0,log2]，cone boundary omega 在2n面均匀、其它坐标均匀[-1/2,1/2]，x=y+r omega。shell Lebesgue Jacobian 为 n r^{n-1}dr dnu=r^n dt dnu（t=nlogr），所以 q_y=r^{-n}/C_n。估计变量 g h_R/q_y 在 [0,C_n]；source 没被冒充 receiver 体积。构造自身的 radial distance 已知，避免 x-y 相消用于这个 capture；但其它 grid 原子的响应 closed 边界仍必须单独检查。

先按完整 mu 采独立 y，然后分别估 S_A(y)、S_B(y)；I2 估计为 mean S_A S_B，conditional independence 给无偏平方估计，避免 mean S_hat² 的正偏差。I1 估 mean(S_A+S_B)/2。尾 plugin mean(S_hat-K)+ 为有偏凸估计，不能称无偏或解析零。保存 source y/component、两个 S replicas、near-response 候选贡献范围、所有 receiver M/index/events、输入 hash、seeds、numpy/Python version。

### 有限样本区间范围

registration 固定 33 records。inner variables 界 C_n，平方 cross-replicate cluster variable 界 C_n²；global union 用 .02 覆盖 sampled-source inner mean bands，.03 覆盖 outer I1、I2、三种真实 source-tail 的 finite-record Hoeffding。tail 用 sampled source mean 的上下 band，再加 true outer-function concentration。不同 record 可共享随机数，union 无需它们独立。

这些区间针对 ideal exact finite-source model；floating-point response/near screen 未作区间认证。宽 bound 常常 vacuous，不能以 zero plugin 认定 tail 为零或小。source inner/outer 误差不被 quadrature 节点增加消除。脚本只对自己的新 progress receipt 更新运行状态，没有覆盖旧数据。

## 6. 完成数据与判断

exec session **84840 正常 exit 0**，耗时46.3437秒；三轮9cases加两个same-n cases，共11 profiles、33记录，**33/33都有观测超阈事件**。12项 explicit scalar closed-atom response 比较与一个 exact box plateau tie 检查通过；这些是实现守卫，不是二阶 theorem。

|n/比较|capacity grid 的 I2/I1 范围|resolution mixture 的范围|box mixture 的范围|
|---|---:|---:|---:|
|4|.1014–.5079|.1144–.5383|.2035–.4307|
|16|.1295–.6937|.1238–.5423|.2681–.2966|
|64|.1176–.5747|.1306–.6419|.1975–.2056|
|same n16, layers2|—|.1193–1.0854|—|
|same n16, layers8|—|.07593–.5330|—|

最大观测 source S estimate 为2.0999（n64 resolution mixture, alpha2）。K=sqrt(n) 与2sqrt(n) 的全部 plugin tail 都为数值0；K=.5sqrt(n) 最大 plugin ratio约.06381（n4 capacity, alpha2）。**全部 I1 confidence 下界为0，全部 source-tail ratio 上界无界。** 因而没有通过任何非平凡尾合同，亦不能支持 I2/I1 的维数阶。source-y 未采中的稀有高 S 来源仍可能主导总体。

共有2867个 near-max-response observations；near atom/box boundary、constructed capture boundary 与 near threshold 计数在本批均为0。这些0只是 tolerance screen 的观测统计，不是闭边界不存在或 arithmetic error 为零的认证。11个profile hashes全部复核。最大变化出现在 same n16 layers2 vs8 的 alpha2：1.0854 vs.07593；这反映输入/稀有事件敏感度，也可能含 sampling error，不解释为层数规律。

这些 restricted families 与 fixed thresholds 没有实现总表中 sharp lower 输入，不能以三族非空替代“真正 critical rare events 已被测试”。全批保持封存；本轮不追调 alpha。

## 7. 下轮需要的 critical 参数映射与 tilting，尚未启动

一个可核的入口是单个**无权 interior tensor lattice**，step=1/k，半宽 A=16n，q=2Ak+1，L 满足 kL=m+theta、0<theta<1。在 interior receiver 的 uniform Lebesgue grid phase 下，每坐标 count C_i=m或m+1，高 count 的概率theta。令 B_n 为高 count 个数（此处符号非预算常数），则

\[
 \log\frac{U_L}{(2A)^{-n}}
 =n\log\frac{2Ak}{q}
 +B_n\log(1+1/m)-n\log(1+\theta/m). \tag{8}
\]

因此指定 deviation p>theta 的 critical 阈值映射为

\[
 \log\alpha_n
 =n\log(2Ak/q)+n[p\log(1+1/m)-\log(1+\theta/m)]. \tag{9}
\]

其 uniform-phase receiver 概率具有 binomial upper-tail 稀有度；当 n(p-theta)² 增大时，未经倾斜的固定样本可能全漏。登记的 log alpha=O(1)、k<=8 和有限 J 不保证覆盖任何总表的 moderate-deviation critical regime，更未覆盖 A13/A14 的完整输入。对 weighted grids 或正 mixture，式 (8) 仅是一个 component 的映射；**完整 M 必须重新用整份 mixture 计算**，不能 componentwise取max或仅保留超阈组件。

后续若要针对性压力，应先选择总表中确有资格的完整下界 input/scale/threshold，或明确只选式 (8) 可重建的较弱 lattice 子模型。冻结 n,m,theta,p 和 tau，再为每个来源 cone 的 free-coordinate count 制定有完全 Radon–Nikodym 权的 receiver tilting，保留 base proposal 正支撑；finite winner所有节点仍重算。若 source S 高层也稀有，outer source proposal必须相对完整 mu 合法加权，不能只偏抽获胜组件。训练选择 tilt 参数与独立新 seed 验证分开；不得看到本批后调 tau直到大比值，再称 independent验证。

连续 winner 压力需改为完整合法 arrival/节点包；本批 J9/17/33 可能漏窄的真实 winner 区间。对 tensor grid，每坐标实际 arrival sidelength 可解析列出，候选并集约 O(nk)，可在下一登记中重建连续 [1,2] winner；当前没有做这一步。对 A13/A14 必须另外保留其真实 winner 包、输入近似误差和共同 source。没有免费的 wide-window 或节点复杂度结论。

## 8. 收口

本轮主数学问题仍未解：是否能在 (3) 的实际同源 cube 交集上取得 sqrt(n) polylog(n) fee，或绕过平方以证明 (6) 的 deterministic source-tail。已知 n 级 (5)、表示身份、数值零或其它 coupling 的 count-tail 都不支付此目标。

新文件 registration/script/results/progress、11 profiles 与统一 receipt 保存于同目录同prefix；receipt 记录 hashes、session终态、scope。下轮只留下有根据的 critical映射/合法 tilting 要求，没有自动启动新搜索或任何旧实验。
