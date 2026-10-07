# Local-clock entropy barrier：新算术守卫终态

2026-10-07。本 prefix 独立新增；使用并读取 D04/L03 的 SKILL、method、cube-interface、provenance，按同一条件律与零质量约定处理标签熵，按有理向外区间处理对数。所有生成元和对数余项在本报告自证，不借未核专题文献。未改他人 proof/主账，未运行旧 max 或 oracle。

## 1. 核验的精确公式与范围

有限正 posterior p，sum p=1，非负 rates a，abar=sum a_i p_i。Shannon H(p)=−sum p_i ln p_i。采用质量倾斜 Q 的无跳漂移
dot p_i=p_i(a_i−abar)，以及删除 i 的速率 a_i(1−p_i)；此速率不是原未倾斜死亡速率 a_i。若 0<p_i<1，删除后的剩余 posterior 为 p_j/(1−p_i)。

连续熵漂移为 −sum a_i p_i ln p_i−abar H。删除熵差为
[p_i H+p_i ln p_i]/(1−p_i)+ln(1−p_i)。把它乘真正的质量倾斜速率再相加，所有 ln p_j 的有理系数恰好抵消，余下

Q H=sum a_i(1−p_i)ln(1−p_i)，
D=−Q H=sum a_i[−(1−p_i)ln(1−p_i)]≤abar。

在 p_i=1 时删除速率零，既不评价空 posterior，也不评价 ln0；剩余项按 0ln0=0。a_i=0 时对应项精确零。未保留的零 posterior 标签贡献为零，不参与删除条件化。

标量界可自证：对 0≤p<1，由 −ln(1−p)≤p/(1−p)，得 0≤−(1−p)ln(1−p)≤p，p=1 端点延拓亦成立。这是解析全域证明；有限守卫不替代它。实际 accepted-history law、时间补偿及停时工作量由 local_clock_entropy_barrier_20261007.md 和 LC 解析证明承担，本批不认证它们。

## 2. 冻结的三轮与全部数值保存

代码与 registration 冻结后，新 main 只执行一次，工具 exit0；无失败尝试。三轮原权重依次为 (1,2,4)、(1,2,3,5,11)、(1,1,2,3,5,8,13,55)。遍历每个非空原标签子集，只重归一化质量、不重排标签；各用 rates：

- all1：a_i=1；
- alternating：a_i=i mod2，原标签从0开始；
- graded：a_i=(i+1)/(原N+1)。

共293子集×3=879行。每行保存全部 posterior、rates、abar、H/D/slack 的严格区间、每个合法删除的后验归一化和 entropy log 系数、无跳及跳跃 ln p_j 系数、剩余 ln(1−p_i) 系数。formal ln p_j 的抵消在区间评价前逐项 Fraction 验证，不凭近似小数取消。

|原N|子集数|行数|检查数|p=1行|D精确零行|
|---:|---:|---:|---:|---:|---:|
|3|7|21|240|9|10|
|5|31|93|1410|15|19|
|8|255|765|16089|24|35|

非退化行区间严格给 0<D≤abar，并保存正 slack 的下界；零退化 D=[0,0] 精确处理。**879行是纯 finite posterior 标量／生成元代数守卫，没有 source 坐标、实际时间 history 或累计 work，不称实际 clock 轨迹，也没有验证某个 policy 的可预见资格。**

全部879行存于 local_clock_entropy_guard_20261007_posterior_rows.json.gz，确定性 gzip。SHA256：
6cd936a95f4920518da07314db2e8a251f9572cb9057de7ac7fc7bf75b779cbd。

## 3. 对数区间的自含证书

每个正有理 q 精确写成 2^k m，1≤m<2。令 z=(m−1)/(m+1)，则 0≤z<1/3；ln2 使用 z=1/3。用100项正 atanh 级数

ln m=2 sum_{j=0}^{99} z^(2j+1)/(2j+1)+remainder，
0≤remainder≤2z^201/[201(1−z²)]。

因此先得到全 Fraction 闭区间；加 k ln2 时按 k 正负交换 ln2 上下界，最后向外取256-bit dyadic。q=1直接返回[0,0]。乘负 entropy 系数时交换端点，全计算不依赖 Decimal/float。区间跨零不会算 PASS；退化项先以精确零判定，避免用窄区间冒充等号。

此程序是本次独立新代码，没有调用旧 log/max oracle。所报小数只作阅读诊断，不参与证书。

## 4. 三组完整压缩 product-source 几何

另取 n=4,16,64，a=1,b=2,d=1/(8n)，完整来源

μ=(1/n^n) sum_{k∈{0,…,n−1}^n} δ_(dk)，W=1。

这三组不是879个标量状态的空间实现：它们是独立的完整压缩来源几何守卫。无需枚举 n^n 点，保存每维所有 k 的精确位置 dk 和半开 d-grid 的唯一 cell index k。逐坐标 floor(dk/d)=k，任意两个不同多指标有某个不同坐标，故来源点和来源格均不同。每个 atom 质量1/n^n，全部原来源总质量恰为1。

最大坐标 d(n−1)<1/8，所以同一 receiver0 的闭 Q_b 捕获全部来源，甚至 Q_a 也捕获。它足以核真实 co-capture clique 的几何资格；没有重新评价任何 maximal winner。

初始 posterior 等权，解析精确 H=ln(n^n)=n ln n。分别对 ln N 和 n ln n 求区间并检验相交，只是算术一致性；**区间相交不是两个真值相等的证明**，相等性来自精确 N=n^n 与 log 幂恒等式。

|n|几何／熵检查数|H诊断近似|
|---:|---:|---:|
|4|15|5.54517744448|
|16|39|44.3614195558|
|64|135|266.168517335|

N 的精确整数、全压缩 source 表述、atom weight、全部一维 cell 坐标与熵严格区间均保存于 results。没有拟合 n log n 的阶数，没有 giant N 枚举，也没有认证 nearmax、high-m 或 actual geom/history 资格。

## 5. 成功终态收据

**17938 PASS** = posterior 17739项 + product189项 + 冻结代码/skill9项 + 总数量1项。耗时0.36784137500217184秒，一次成功执行、失败0、旧 oracle0。只计成功终态，不叠加注册或文件阅读为数学检查。

- script SHA256：
  e2303b4a2ce34fd607faea047fd237e83ef4aca56e8dc6f499e7f30a4e1594fd
- registration：
  26e256cab81fb01fc8cb8af0ad4b155fa6eaa3a2d74cb1fab14c04a91735c3de
- results：
  a686495d5c1969dd4f6b5eef04d67b58a6ee303f161ed983cd442e3735d1c418
- receipt：
  e844b8ffcb4c987197a76f8f7dd4d6afa68add91f575772211f1332baf401260

Parent 已只读核验本批代码、注册与收据的算术范围，无需重跑。合法真实 process 与其熵停时下界仍依解析证明；本报告只是数学恒等式与压缩来源几何检查，不宣称一般 geom 余项闭合。

