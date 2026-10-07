# 双平移完整来源 pooling：三轮精确 cell/log 守卫终态

2026-10-07；gated_radial_energy_audit，用户指定 6.1-sol high。只新增本 prefix，不改原解析稿或主账。使用已读 L03 的精确残差、向外区间、完整输入和冻结执行流程；本代码重新实现 log 区间与原atom cell响应，没有导入旧 oracle/probe。

**终态：一次成功运行，41413 项精确检查 PASS，三轮 n=1,2,3。** 所有完整来源、全部有限尺度、真实赢家和完整 halos均保留。严格区间认证 amplitude gain、cross-winner regret、合同残差均为正，而三轮 Φν−Φμ均为负。这验证原§4合同以及保留阈值/cross-regret的必要性；不是 nearmax、连续窗口、actual geom 或 weak-counter 认证。

## 1. 注册与完整原问题

执行前冻结 [registration](nearmax_translation_pool_guard_20261007_registration.json) 和 [script](nearmax_translation_pool_guard_20261007.py)。被核解析稿 nearmax_smoothing_obstacle_20261007.md SHA：05afeb387e63a4bc8d5102570e258d2171bde78b157283ac998473f3fa3aad2a。

每轮原完整来源
\[
 \mu=\tfrac13\delta_0+\tfrac23\delta_{(3/2,\ldots,3/2)},\quad W=1,\quad
 J=\{1,2\},\quad \tau=\frac3{8\,2^n},\quad z=(1/4,\ldots,1/4).
\]
两份整来源 μ±=T±zμ，各质量1；ν=(μ++μ−)/2也是质量1。原 Φ、两份平移 Φ、pooling Φ均来自同一完整有限尺度族，不重选输入或阈值。tie选择同一 smallest R，零响应也选择R=1。

对 μ、μ+、μ− 的每个 atom y、每个原 R∈J，列出所有坐标超平面 x_i=y_i±R/2，取它们的笛卡尔开 cells。ν的全部 atoms已包含于这组端点。每个 cell 内每个 atom/cube incidence恒定，故用有理中点**精确表示整 cell 响应**，不是抽样 Lebesgue。所有 faces Lebesgue零测；bounding box外所有原/平移/pooled响应为0。各cell volume为完整边长乘积，联合atom位置仍保留，不把各轴来源独立化。

每个完整来源的两个固定尺度响应积分都精确为1；φ±的合并 log coefficient maps与原φ逐项相等。这独立核验 full halo/source-once 与平移不变性。

## 2. 每 cell 的真实赢家及两份完整 cross-regret

定义 A=M+、B=M−、原真赢家 R±，并取
\[
 d_1=B-h_{R_+}*\mu_-,\quad d_2=A-h_{R_-}*\mu_+,\quad
 d=\min(d_1,d_2),\quad \beta=d/(A+B),
\]
A+B=0时β=0。逐cell精确检查 d_i≥0、d≤min(A,B)、0≤β≤1/2，以及真实 Mν 满足
\[
 \tfrac12(A+B)(1-\beta)\le M_\nu\le\tfrac12(A+B).
\]
ν的固定query响应严格等于两完整平移来源响应的半和。所有 M/winner均用 Fraction，未先把新winner冻结成旧winner。

记 T0=(A+B)/2、F0={T0>τ}、E±={M±>τ}，只在两个E的交集上累计
\[
 {\rm gain}=\tau\int_{E_+\cap E_-}
 \log\cosh[(\log A-\log B)/2]dx.
\]
其精确正log参数是 (A+B)^2/(4AB)，系数 τ·cellvolume/2。其余两费为
\[
 {\rm crossing}=\tfrac{\tau\log2}{2}|E_+\triangle E_-|,\qquad
 {\rm regret}=\tau\int_{F_0}-\log(1-\beta)dx.
\]
阈值始终strict，未删除plateau或把F0换成E交集。

## 3. 严格 log 区间与独立的正残差重组

新代码将任意正有理q按二进制精确分解 q=2^k r，1≤r<2。令 w=(r−1)/(r+1)∈[0,1/3)，用
\[
 \log r=2\sum_{j=0}^{99}\frac{w^{2j+1}}{2j+1}+\mathcal R,\quad
 0\le\mathcal R\le\frac{2w^{201}}{201(1-w^2)}.
\]
ln2用同一100项式、w=1/3。k<0时正确交换 ln2 lower/upper；有符号log系数按其符号选择外端点。完全Fraction，先合并相同log参数和取消零系数。输出256bit向外dyadic端点及原有理endpoint的integer-byte SHA；float字段仅诊断，不作证书。

另独立把每cell完整合同残差重组为一个正参数的log。设
\[
 B_x=\begin{cases}(A+B)^2/(4AB),&A>\tau,\ B>\tau,\\1,&\text{其它},\end{cases}
\]
\[
 r_x=
 \frac{\max(1,M_\nu/\tau)^2}
      {\max(1,A/\tau)\max(1,B/\tau)}
 \frac{2^{1_{E_+\triangle E_-}}}
      {B_x(1-\beta)^{2\,1_{F_0}}}.
\]
每cell以Fraction直接核 r_x≥1。故完整残差也等于
\[
 \Phi_\nu-\Phi_\mu-{\rm gain}+{\rm crossing}+{\rm regret}
 =\sum_{\rm cells}\frac{\tau\,{\rm volume}}2\log r_x\ge0.
\]
这个正重组保留所有原字段，不是删除regret以后换一个合同。直接四项log区间和正重组log区间也检查交叠。三轮的直接残差区间lower均严格正，没有 interval crosses zero 或未决符号；注册的 unresolved政策未被触发。

## 4. 三轮结果

[results](nearmax_translation_pool_guard_20261007_results.json)、[receipt](nearmax_translation_pool_guard_20261007_receipt.json)，以及逐轮完整cells：
[1D](nearmax_translation_pool_guard_20261007_n1_cells.json.gz)、
[2D](nearmax_translation_pool_guard_20261007_n2_cells.json.gz)、
[3D](nearmax_translation_pool_guard_20261007_n3_cells.json.gz)。

|n|cells|exact checks|交集体积|对称差体积|F0中β>0的体积|
|---:|---:|---:|---:|---:|---:|
|1|14|214|5/2|1|2|
|2|196|2762|11/4|9/2|3|
|3|2744|38434|29/8|43/4|7/2|

此外3个冻结script/proof/维数完整性检查，共41413。所有原fixed-kernel full-space质量积分精确等1。总耗时0.2353784999868367秒，执行一次，无随机种子，0旧oracle调用，无首轮失败或修复。

下表为严格区间的可读近似；符号依有理区间认证，不依这些小数。

|n|Φν−Φμ|amplitude gain|crossing cost|regret cost|合同残差|
|---:|---:|---:|---:|---:|---:|
|1|−0.136015583|0.008401140|0.064982548|0.128488688|0.049054513|
|2|−0.098695291|0.029103482|0.146210733|0.072879024|0.091290984|
|3|−0.064970313|0.037898722|0.174640598|0.020816384|0.092587947|

gain/regret/对称差都非零。完整来源Pooling的Φ下降不被正gain否定：它被真实cross-regret与阈值项抵销。没有从这些固定输入估计unknown ε/K，也没有据有限步输入声称nearmax。

## 5. 冻结 Hash 与范围

- script SHA：8be03dc2151573f985d986e42fbc7213c45b3e69f43791fc75c2911516645835。
- registration SHA：2d445aaeff18c080d8145fab8febe8382a855a5cc87af139872328efbffbff9f。
- results SHA：e1be0d9906cc3116433a206af58714ad74385a531213afb2cb0109f23c1d73de。
- receipt SHA：ac25bbcd24f27550a5831dd39cfff3f71e500493e2eefbe00ab7bdea82c22ce1。

这是新全空间有限family pooling合同守卫；没有 continuous帽误传、receiver MC、结构近优类推、actual FIRST/history样本或新weak反例。原一般√n目标保持未解；数值只给原§4保留全部项的解析合同可复核收据。

