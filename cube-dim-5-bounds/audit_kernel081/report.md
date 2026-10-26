# audit_kernel081：070闭合覆盖、新核独立值审与全域g资格

**本轮全部前提通过，不再是pending_values。** 从冻结077源出发，五批070新增1318个anchor、1001个通过叶和317个split，最终97960个passed、45600个split、52360个root、queue/inflight均为0。旧178722缓存记录逐对象保留，1318个新增H及四个物理导数共**6590个分量**全部经不同算法、不同算术路径独立包围，新的512-bit准确dyadic区间**严格位于**原256-bit生产球中。

在已锁历史局部资格、075及同alpha1的029分析定理前提下，闭合覆盖已具备并在本报告中接成

\[
G(x)>\frac1{10000}\|x\|_\infty^{-4}\quad(x\ne0),\qquad G(0)=+\infty.
\]

因此同模型canonical Green在闭单位cube上的essential下界可取 **g0=1/10000**。主任务保留最终接纳责任。**没有更新C5/Cn，没有重算、联合或更新I，没有把这一g结论直接写成常数突破。** 后续不同target须另行资格。

## 1. 冻结范围及真实运行

唯一旧源为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_framecover077/final/state.json`，SHA
`c55918ec39ef62629d53d20e1aaed7a5faf454d52212c5897bd16fd15056ea9b`。

070最终full state为
`/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover070/final/state.json`，SHA
`8727af764ea088b1c30522f0204d6630ad6e0a4d804c7698c5ed2d646cedf343`。
070冻结manifest锁为
`5f8caa2fe34775f1dfec518f484870828054966d10fc13575b28a211a8c5fdd2`，最后overlay锁为
`19fffc35d12b56e96db3f392e9bf08266543e22076a4dce1878a95ea7f862d9a`。

用户已取消程序wall硬上限。本轮所有程序均无alarm/deadline，但输入集合和每次点数有限。没有计算旧178722个点；没有派发生产新核、MATLAB、提交或推送；只新增audit_kernel081。旧062/069/078共20文件锁在所有阶段保持不变。

|阶段|实际PID/范围|时间与状态|
|---|---|---|
|最终结构/全部角点|71965|8.787308667秒，exit0，RSS270860288 bytes|
|独立数值|12个有限进程：首1点、10批128点、末37点|累计408.378293080秒；最长40.870454333秒/进程；全部exit0；峰值RSS197115904 bytes|
|最终联合资格|72817|1.229006375秒，exit0，RSS85606400 bytes|

数值PID依次为70949、70966、71134、71696、71991、72065、72186、72275、72321、72384、72572、72610。全部成功点立即原子保存，下一有限调用从准确成功前缀后开始，**没有重跑成功前缀**。每个数值、最终结构及联合进程运行后实际OS核验均已不存在。五个生产master和60个worker的65个PID也独立核验OS不存在，所有池close/join/exit0收据与真实返回点的PID一致。

## 2. 数值独立性：有理时间积分，直接极点权，整数log和π

固定模型为alpha1、16条**未单位化且权1**的body Cauchy方向、Fourier逆因子(2π)^−5，epsilon=1/64。共享的解析输入只有已经认证的017真fan及其同模型归一化、正epsilon时间resolvent身份；没有调用046、070、039/044代码或其DD函数。

017表绝对路径在登记中，SHA为
`8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9`，280个整数vertex、3504个真cone、分母144。每个新点单独准确构造共同分母D及整数投影k，重新核全部280投影非零、互异，最小整数gap与raw记录一致。每点实际遍历全部3504锥，全批 **4618272锥**。

在一个锥的五个不同投影k_j上，令

\[
w_j=\prod_{m\ne j}(k_j-k_m)^{-1},\qquad a=144D\epsilon=A/B,\quad M_j=A^2+B^2k_j^2.
\]

对正实q，直接部分分式为
`∏(q−ik_m)^−1=Σw_j/(q−ik_j)`，其中i^4=1使权为上述实有理数。先在[a,T]积分，再令T趋于无穷；由Σw_j=0消去共同logT，实部恰为`−Σw_j log(a²+k_j²)/2`。把a写成A/B时共同logB同样准确抵消。故H为

\[
\frac{D^4}{4608\pi^5}
\sum_{\mathrm{cone}}\sum_j-\frac{\det_{\rm num}\,w_j}{2}\log M_j.
\]

这是一条直接有理时间积分/极点权路径，**没有DD商差递推**。对任一物理尾坐标，固定D令s_j=D(vertex_j)_tail，

\[
\sigma_j=\sum_{m\ne j}\frac{s_j-s_m}{k_j-k_m},\qquad
\partial H=\frac{D^4}{4608\pi^5}\sum_{\mathrm{cone},j}\det_{\rm num}w_j
\left[\frac{\sigma_j}{2}\log M_j-\frac{B^2s_jk_j}{M_j}\right].
\]

代码不做数值微分，也不复用作者jet传播。若d_j=∏(k_j−k_m)，则σ_j的分子准确写成`T_j=Σ(s_j−s_m)(d_j/(k_j−k_m))`，log系数为`det*T_j/(2d_j²)`，全部通过整数积和有向除法求值。四个**物理**导数均走此路径，不作basis坐标变换。

数值算术只用Python整数与Fraction，**不导入flint、Arb、mpmath或任何超越函数库**：

- 用2的整数幂把任意正有理数约化到u∈[1,2)，令y=(u−1)/(u+1)∈[0,1/3]，以282项`2Σ y^(2j+1)/(2j+1)`包log u。
- 每一步在2^−768固定网格向外floor/ceil。准确尾项上界为`9/[4(2N+1)3^(2N+1)]`；log2亦从y=1/3用相同带显式尾项级数产生。
- π用Machin式`16atan(1/5)−4atan(1/239)`，各216项交错有理级数加下一项余项界。恒等式由`(5+i)^4(239−i)=114244(1+i)`及正确角度范围给出，不依赖外部π常数。
- 所有有理系数、π^5归一化与正负导数区间乘法均准确有向；最后向外取512-bit dyadic端点，再逐分量验证严格内含于生产球。

全部 **1318个H+全部5272个物理导数**通过，未仅选固定一分量或若干采样点。候选最小新缓存H点为第1311号（0-based，batch5 ordinal1），坐标
`[1,522241/524288,4177921/4194304,33292289/33554432,2097153/268435456]`，其H及四导数的独立端点均保存在qualification_result，其他1317点亦有全部五分量端点。

新生产H下端最小准确值为

\[
\frac{50383906330716966665052244156296703307497946433492289435430235345070422896065}
{474284397516047136454946754595585670566993857190463750305618264096412179005177856}
>\frac{53115711}{500000000000}>\frac1{10000}.
\]

独立更窄包围严格落在这些生产区间内部，因而本轮并非只重复“all_cones/quality”标签。

## 3. 每个新点的资格、1318个证书及3959事件

每个raw point checkpoint先核point_manifest中的字节SHA，再核task哈希、point/record/cache一致性、真实worker PID、complete/partial=null、precision_index=0、bits256、一次3504锥完整attempt及一次segment。准确重建全部candidate中心/微扰顺序、非零/互异280投影guard；没有对零/重复投影调用普通泛化DD。

原球五元区间有序，H下端严格高于1/10000，实际宽度逐项满足相对2^−40质量条件。对四个保存log-gradient区间，核H和导数端点的全部交叉不等式，正确处理导数正负与区间除法。Delta另外用本审计的768-bit整数log下端核`Delta_saved≤log(H_lower/target)`；没有沿用作者512-bit Arb log复核作为证据。

物理box由root及path每隔四位的二进制子串独立还原；source和当前split的axis、两个闭child由该几何一致确定。对每个新certificate准确清分母后检查16角点

\[
\Delta_L+\sum_j\min(g_{j,L}v_j,g_{j,U}v_j)-256\sum_jv_j^2.
\]

全部 **21088角点**的准确minimum、worst_corner_bits、通过/失败符号一致。1001通过、317几何未证后split；没有把一次几何未证解释成H反例，或靠升精度消去几何余量。新通过叶最小准确余量为

\[
\frac{102579975995791144455165715583734087316776101019968069452076155995857145398164825}
{1942668892225729070919461906823518906642406839052139521251812409738904285205208498176}
>\frac{52803633}{10^{12}}>0.
\]

五批从各自initial_overlay独立检查全部dispatch/received/merge/finish事件。规范event哈希链、序号、ordinal派发/接收/合并唯一性、有限注册前缀、raw→返回record→新cache→certificate完整关联均核对；再由事件与自己的闭父子几何构造整个post-overlay对象，与冻结state**逐字段完全相等**，非只比摘要。累计 **3959=3×1318+5事件**。

每批counterexample=null、faults=[]、inflight/resume为空；首counter为零，累计newstarts/newcomplete/segment_calls及cache_hits=0准确。父overlay路径/SHA、run/code/math/source/cache_ref锁及12worker身份都核对。初始684点与冻结077困难清单顺序相同，后续child按准确父H下端/失败margin/root/path顺序处理。

|批|新点|新通过|新split|剩余queue|
|---|---:|---:|---:|---:|
|1|684|474|210|420|
|2|420|333|87|174|
|3|174|158|16|32|
|4|32|28|4|8|
|5|8|8|0|0|
|合计|1318|1001|317|0|

## 4. 源保留、精确树与最终schema

对旧077及070最终full state同步流式解析，直接比较对象：旧 **178722 cache、96959 passed、45283 split**完整前缀与对象逐项不变；新 **1318 cache、1001 passed、317 split**与小overlay逐项相等，cache新旧key互不相交；其余 **50个旧顶层值**全相等，只新增task070/status070/fingerprint070/conversion070四键并把queue置空。

fingerprint070准确绑定run_registration SHA及最后overlay SHA；conversion070的源/overlay完整路径及SHA、新计数、数学标签、closed-cover候选与false全域/C5生产claim逐键相等。保留旧status/fingerprint等源元数据不冒充当前状态。

源及每批frontier在每个root上转成共同深度D的dyadic整数区间，排序后准确从0连续铺到2^D；全部52360 root prefix-free且Kraft和严格为1。活split集合同时准确等于全部叶path真前缀；新child为物理最长边（最低坐标优先）闭二分，公共边界没有遗漏。

最终
`97960=52360+45600`，queue=0；亦有`684+317−1001=0`。所以已认证叶构成整个root覆盖，而不是passed与待办共同覆盖。旧叶局部资格由锁定078/069/062继承，没有重新计算旧178722个值。

## 5. 从闭ordered face到全域G：只在全部值审后消费

075报告/manifest锁为`3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225` / `1d5e270eb1664d714ec80f0b371420d40028b24c3943196b9ab41a0b3de1f9d1`。其正Cauchy纤维—时间混合证明给全空间H_epsilon正、有限、C²及full-frame曲率；区间min-affine项本身全局凹，减正二次仍凹，因此16角点控制整个闭盒。S1/S2、零/重投影节点、符号边界和相邻闭盒边界均由该全空间定理支付，不需要在这些点重新求普通DD。

根列表准确为全部`31≥r2≥r3≥r4≥r5≥0`整数四元组。对任意`1≥y2≥y3≥y4≥y5≥0`，取`rj=min(31,floor(32*yj))`；单调性保持顺序且闭cell包含y，包括yj=1。因此根覆盖完整**闭ordered face**。queue=0和所有旧/新叶严格资格给
`H_(1/64)(1,y)>1/10000`，从而完整`G(1,y)≥H_(1/64)(1,y)>1/10000`。

16个未单位化body方向的无向line集合在全部signed coordinate permutations下不变：联合程序对5个坐标改号与4个相邻坐标交换生成元逐一作准确整数集合核验。因Cauchy变量对称、同权独立，完整p_t与G具有这些对称。任何cube面点都能经该群变到上述ordered face，故全部闭cube面均有同一G下界。

同alpha1的029报告锁
`298f0de7277ea6a7a516b85344c5189562bf7026f3f8de43d4d3f304ff6fed81`，manifest锁
`9246272cbe732b603365f3f2ad3210a0dba19b011c8b3a7f18557273bc82a2c6`。已完整读过该原件，模型归一化与本轮一致，不使用alternative_alpha的1.5版本。由其完整Green扩展值齐次性，对任意非零x取r=||x||∞、y=x/r：

\[
G(x)=r^{-4}G(y)>\frac{r^{-4}}{10000}.
\]

这里缩放的是**完整G**，没有把H_epsilon当作齐次函数。对cube内r≤1即G(x)>1/10000；原点由p_t(0)=t^−5p_1(0)>0的正时间缩放得G(0)=∞。奇异集合允许G=∞，不破坏非负积分或下界。一般L1、所有正尺度的传递资格由同一029给出，但本轮没有消费I去更新任何Cn/C5数值。

联合程序在6590个严格内含、每个raw/代码/输入锁、闭合树及理论锁全部通过之后，才把状态从pending_values改为qualified。它没有再求H。以上是有限数值审计加锁定分析证明的组合，未声称测度论或Python整数实现已形式化验证；主保留最终证明及接纳责任。

## 6. 文件、锁与复现范围

必要白名单由manifest逐项列出：登记；common/arith/numeric/structural/qualify五个标准库模块；数值完整前缀、结构/角点完整结果；联合登记/结果；本报告及manifest。没有复制全量source/cache/raw点文件。numeric_result约5.2MB保存全部6590个独立dyadic区间，结构结果约1.2MB保存全部1318准确证书核验。

关键结果锁：

|文件|SHA-256|
|---|---|
|registration.json|`8ef301297b5b2469a94a6dc3da8499b66557940dd07dd21c31397f672107cb3c`|
|numeric_result.json|`7260d2fa8b2859341130a1ecbf69d532c8bb29aa1d06880346f7dc9b90db6b68`|
|structural_result.json|`88ea6601332abdc41c5c68124a0897f92276fce8591f66551ec3cd1a044f6db7`|
|qualification_registration.json|`a26e084b2dc5e6284b5112dd1a6e15839cb4044eed5d23228a759ee6eea82316`|
|qualification_result.json|`e13c430251175ca47a89912d1a74c58733ba49c3783620af84d46681eb8dc59d`|

标准库Python即可，无python-flint依赖。复现请在新的独立输出副本保护本轮锁：`structural.py --run`重做结构/Delta/角点；数值从空前缀用`numeric.py --count 64`或`--count 128`逐次推进，每次成功前缀不重跑；`qualify.py`按完整联合登记锁作最后只读值核验与理论接入，并只写该副本的联合结果。不要在当前已锁输出上重复调用numeric脚本，避免新增空调用收据改变结果hash。
