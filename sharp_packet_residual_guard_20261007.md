# Sharp packet residual envelope：三轮标量与独立 sharpness 补充终态

2026-10-07。本 prefix 独立新增，未改他人文件、主账或旧 oracle。L03 SKILL/method/cube-interface/provenance 已读并使用严格有理区间；J01 同组必要文件新读，遵守从真实律推出生成元的限制，不由本批自由标量参数制造实际过程。以下全域标量及积分解释自含，不引用未核外部定理。两批各先冻结 code/registration，再只执行一次；均 exit0、失败0、旧 oracle0。

## 1. 全域解析解释与本批范围

固定 v≥0、u≥max(1,v)、0≤p≤v/u，定义
h=min(ln u,−ln(1−p))−p；p=1 时取 min=ln u，不评价 ln0。
准确的标量包络为

sup h=G(v)=ln(1+v)−v/(1+v)，
在 u=1+v、p=v/(1+v) 达到。

自含证明：记 p_c=1−1/u。对 p≤p_c，h=−ln(1−p)−p，其导数 p/(1−p)≥0；对 p≥p_c，h=ln u−p，导数为−1。因此固定 u 的最大值在 p=min(v/u,p_c)。当 u≤1+v，最大值为 ln u−1+1/u，导数 (u−1)/u²≥0；当 u≥1+v，最大值为 f(v/u)，f(q)=−ln(1−q)−q 单调递增，而 v/u 随 u 递减。两支在 u=1+v 给同一个 G(v)。v=0 时 p=0、h=G=0。

h 可以严格为负，不能将它无声地截成正值。此全域证明独立于有限参数 guard；本批不作 global numerical optimization 或全域数值认证，也不提供 source、真实 winner/history 或 actual-law 资格。

## 2. 主批冻结参数与结果

三轮 v：

- small：1/64、1/4、1；
- moderate：2、8、32；
- extreme：1/1024、1024、65536。

每个 v 取 u=max(1,v)、v+1、2(v+1)，每个 u 取 p=0、pmax/2、pmax、min(pmax,1−1/u)，pmax=v/u。全部108行保留重复标签，不删掉零点、p1 或负值。另保存 v0、u1/2 的两个工具端点。

分支使用 Fraction 精确比较 u(1−p)≥1，不通过近似 log 排序。h 与 G 的 log 项按其精确有理 argument 合并；饱和点同一 log argument 的系数及有理项抵消为零，才认证等号。没有将区间相交当相等证书。

|轮|参数行|检查数|h零行|h负行|h正行|符号饱和行|p1行|
|---|---:|---:|---:|---:|---:|---:|---:|
|small|36|213|12|6|18|6|1|
|moderate|36|222|9|1|26|6|3|
|extreme|36|239|10|2|24|6|2|

**主批688 PASS**：三轮674项、冻结 code/skill9项、v0端点4项、总数量1项。elapsed0.01882341699092649秒。全部 log 的符号系数、h/G/gap 区间、分支及符号均保存 results；18个饱和行 gap=[0,0] 是精确 symbolic 证书，其余行 gap 下界严格正。9个负 h 保持原符号，6个 p1 从未求 ln0。

对数证书为100项正 atanh 级数、精确 power2 range reduction、余项2z^201/[201(1−z²)]，最终向外256-bit dyadic。负系数交换端点；log argument1直接为零。普通小数只作诊断。

主批还核验
G'(v)=1/(1+v)−1/(1+v)²=v/(1+v)²
的有理恒等式，并对已指定正积分
G(v)=∫_0^v t/(1+t)²dt=v²∫_0^1 s/(1+vs)²ds
保存精确有限 Darboux 外包。分割含 v*j/4、由1/4起的 powers2、0/v 和落在区间内的临界点1；每段 extrema 由端点及临界点取得，非浮点 quadrature。G 的严格区间位于独立 Darboux 外包内仅为一致性；积分相等性由解析 primitive 与精确替换证明。

## 3. 独立 supplement：真正混合 log kernel

主批已冻结后 root 指定另一个有用混合恒等式，故单独注册、未更改或重跑主批。取
k(t)=2t/(1+t)³。
有精确 primitive/tail：

∫_0^v k(t)dt=v²/(1+v)²，
∫_v^∞ k(t)dt=(1+2v)/(1+v)²；
∫_0^v k(t)/t dt=1−1/(1+v)²，
∫_v^∞ k(t)/t dt=1/(1+v)²。

故 k 与 k/t 的全积分各为1，前者是 mixture 质量，后者是 weak-threshold 变换时可能用到的权，二者不混淆。对 v>0，用非负 Tonelli 展开 log(v/t)=∫_t^v ds/s，

∫_0^∞ k(t)log_+(v/t)dt
=∫_0^v [s²/(1+s)²]ds/s
=G(v)。

v0 时 log_+(0/t) 定义为0（t>0），积分为零。此推导不借未定义端点 log0。补充在原9个 v 加 v0 的10记录上，只用 Fraction 核 primitive 导数、正性、CDF/tail 和各全权1；没有巨大 quadrature，也不把有限 primitive 检查当全域积分证明。

## 4. 独立 supplement：完整真实来源 sharpness

偶维 n=4,16,64，a1、b2、w1、τ=2^(−n/2)，完整来源

ν=τ1_{[-b,b]^n}dx+wδ0，
W=τ(2b)^n+w。

真实背景全部付入 W，background rate0 只是该 sharpness 解释，不免费丢质量。令 V=w/τ、m=V/b^n=1/V。对 receiver x∈Q(0,b)、原任意 R∈[a,b]，逐坐标 |x_i|+R/2≤b，所以 query 完全落在背景支撑。其背景响应为τ。原 packet 捕获响应为 H_B=w/R^n，真实唯一 winner 为 R=max(a,2||x||∞)。因此 M=τ+H_B，在此 receiver 区域 u=1+v、p=v/(1+v)；区外 packet 不能被原尺度捕获，响应≤τ，所以完整 E={M>τ}=Q(0,b)（closed-face 零测边界照原约定）。

完整 residual 的 receiver 积分为 τ∫_E G(H_B/τ)dx。方体径向变量 t=R^n 的 volume Jacobian 给

residual/w
=(τ/w)[a^n G(V/a^n)+∫_{a^n}^{b^n}G(V/t)dt]
=ln(1+m)/m−1/(1+V)，  a=w=1。

primitive 为 t ln(1+V/t)，其导数正是 G(V/t)。故公式保留 inner core 和全部背景，没有把 receiver 只取子域冒充完整 E。随着 n 增大，V→∞、m→0，此 ratio 解析趋1。这说明相对于 packet 质量的 sharp residual 不能免费变成一个趋零常数，不是近极值或一般 weak 型下界。

补充保存全 source 压缩描述、背景质量、完整 W、receiver box包含、radial core/end volume、严格 logratio 区间及解析有理上下包：
1−m/2−1/(1+V) < ratio < 1−1/(1+V)。
三轮 ratio 诊断近似为：

|n|ratio诊断|
|---:|---:|
|4|0.692574205256839|
|16|0.994160895824617|
|64|0.999999999650754|

严格区间见 supplement results；这些小数不参与证书。**补充72 PASS**：冻结 code1项、三组 source每组7项、10个 kernel记录每组5项。一成功执行、失败0，elapsed0.004056292003951967秒，旧 oracle0。主批与补充分别计688与72，不合称一次执行。

这个补充有完整正来源与解析真实 cube winner，但没有认证 nearmax、actual FIRST/global order 或一般 geom 支付。很大的背景质量保留在 W，因此也不能把 ratio≈1 解读为 residual/W≈1。

## 5. 文件收据

主批：
- code 5e7d32e6535946e600363828fbd01d35b8007dff586b895021bfc501feb1d22c
- registration d5bf63104289b1c3bab97e94dc76c4d25c0562c53a0a72ce67548443661a2545
- results 892c8bd45e8f45cf309d03f3f74a45f5234dbf3b3b94ccb3fe0f43fb6c418c90
- receipt e0ee676dbba5a16cacd3fe9910326688cb8e4f015717109e177ba59a878703bb

独立 supplement：
- code b1d28bd62c0eb1ef095ce7c8a26550212b08d8bb76024cfb3239f71ac2613cef
- registration 207283e15c0aa2e438dca174e28ef74d21f889737a17ef60630852d57442ae54
- results f73d4c9bcb6707d00bc696f199e8801663255726f2d6e9b0d463bde60490e604
- receipt d168cc1d29591695ad7a3177ba353d2aa09f1271da5cc2e5610df6c752143df2

以上仅 scalar 参数、精确积分 primitive 与特定完整来源 sharpness 的守卫。没有新增来源一次空间上界，不宣称主 geom 目标闭合。

