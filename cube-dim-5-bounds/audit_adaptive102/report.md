# audit_adaptive102：095 首1024父盒的独立资格审计

**登记范围完整通过，未发现阻断095完整I上证的数学或实现gap。** 13406个新子盒、978638个保存端点所定义的正band、全部新range/mean/chord、1024个完整父替换前缀及最终空间账本均通过本机独立核验。0新H/log/exp/Arb计算，没有重跑978606个新H。

\[
\boxed{I\le54702517/10^{12}<6010267/10^{11}.}
\]

内部I上包约5.470251657598456×10^-5，合法g=1/10000配对cost上端约40.8996707458；准确32门为false。正式C5上界仍32，本报告没有新Cn结论，也不把上包未过门当成真实C5的下界。

## 1. 冻结范围与来源

095最终manifest `ad4e78f0e484d217c1d14bd13f585cd8f85240039b13915ffbd16b96b8d2f433`，report `e1222db0269ec75dd4257768deaa9a4c7f2657acc7bbcee121e0541e0b5d8d69`。本次核全部34件095冻结文件，并核092/082/080/088/094/091/085/076来源manifest全部成员hash。014的source092_primary_acceptance原字节及其frozen_inputs逐项核对；eps512、已接受g100、旧Icap和完整正面slot均一致。

102计算前登记SHA `1177eedacbe4564478941dcc08020d3273d8b1b5f7e4a71fa029184bab5ceda4`，check.py SHA `703d310b77a787af56200c5eb5a9a36a095e1d378624746a137047415e6695b5`。固定13406 child、1024 parent prefix、8个新ratio、至多12个预定义selector的PF样本；8worker、无wallcap，有限子阶段checkpoint。实际运行没有失败、没有恢复；结束后临时checkpoint汇总至run_receipt并清理。

旧094已接受的091全参数比值证明、092初图和旧89个ratio证书直接消费。075/085/091所给正fiber、面内25/24 Hessian、25/48 Taylor和chain ratio数学前提保持不变。原017 fan与原始Arb H生产仍是明确共享输入，不能把本轮下游有理审计或有限PF样本称为全部H数值独立重生产。

## 2. 数学与算法审阅

实际审读095 preflight/prepare_inputs/kernel/run、机制、登记和作者两条verifier源码，不执行其数学核验程序。同步四轴二分后取canonical sort的机制正确：粗盒ordered index ii的16个细盒索引为2ii+bits；坐标置换与该分割相容，canonical child的唯一父为floor(child/2)。因此完整子orbit并集等于完整父orbit，子权和=16父权，体积相同。单轴二分不得套用此式，095未如此使用。

range以同一真实anchor/box的h∞和正band左端a消费091的全线段比值RU；mean使用真实μ与M2而非居中假设。L/U/M向正确方向量化，z^(5/4)的递增凸弦支付盒平均，不能反向用Jensen。095只在所有孩子均成功后安装整组；父self约束仍保留，effU=min(selfU,children sum)、effL=max(selfL,children sum)。未完成父保持旧付款，历史数值bestcover与physical frontier分开。此规则独立于实际完成顺序。

kernel保持body16未单位化、alpha1、Fourier逆因子(2π)^-5、eps=1/512、73端点和同一017整数ray/144归一化。D仅为点坐标共同分母，3504真锥全部聚合再形成H球；不是eps齐次换算，也没有新jet。尾段以H∞=0支付完整[最后端点,∞)。本轮只审这些接口和固定源身份，不重新证明旧017 fan或一般L1传递。

## 3. 全部新geometry、ratio与mean

从092全部量化记录构造56139个旧self，再用080 roots与082两批children拓扑独立重建初图，与094已接受的冻结092 state逐核心字段一致；未重算任何旧chord/H。初始55819个physical leaves的准确weighted self gap排序得到完全相同的1024父。

对每父展开全部有序坐标置换及16个二分模式，独立比较细格索引集合。累计260080个ordered子盒恰等于13406个canonical child各orbit的互不重复并集；唯一父、权重、level及原先不存在的新key全部通过。新child levels为L6/L7/L8，与登记相同。

全部13406 raw逐SHA读取，核point、indices、parent、weight、geometry、bounds与冻结jobs/certificates一致。从实际有理点重建共同分母D和280个投影，每点全部投影非零且互异；anchor严格位于真实满盒内部，u1=0。h∞、μ、四维体积均由真实端点准确重建。

97项catalog中旧89项逐字段等于092表，且消费094已有代数资格，不重新求值。仅新增8个x用不同缩放公式

\[
R(x)=\frac{(\sqrt{16+x^2}+x)^2(\sqrt{256+x^2}+x)^8}{16^9}
\]

和224-bit整数平方根包围，全部落在生产192-bit区间内；真实全部新盒所需x准确被覆盖。新8项及完整外端保存在ratio_results。

L6、L7的offset与094已接受几何逐端点完全相同，分别从旧key `6:62,61,61,61`、`7:127,126,125,124`消费准确M2，未重跑旧moment基线。只有L8是真正新增offset几何；独立编写分段CDF多项式积分

  M2=R²−∫0^R 2r·Π_j[min(r,−lo_j)+min(r,hi_j)]/(hi_j−lo_j) dr

给准确值

  35446092762170780176585335155032095147 / 13937965749081639463459823920405225941237760，

与生产一致。三个几何值逐child精确匹配，不能仅凭level推断anchor相同。

每child的73个保存H端点满足0<lower≤upper；相邻差和最后无限尾共978638个band全部严格正。以2^-256逐项有理向外累加新range和mean，确定2^-160准确bin；13406个L/U/M均与生产精确相等，0模糊bin、0Fraction回退、0全局cap使用。全部新凸弦使用li^5/ui^5的纯整数四次根核方向，不调用Arb或幂函数；13406个积分下/上端均通过。

## 4. 真实point/time注册表

扫描冻结089 base180040缓存及batch8 overlay4744，核其源SHA及epsilon=1/64。与全部13406新点准确比较，恰命中32个H64；每个物理点、time=1/64、values[0]上下端与095的H64_reuse_values、jobs及raw原始cache origin完全相同。复用H64不被误当H128，也不从g105目标失败推断H球失效或g105合格。

全部新点与旧56139个H73点交空，与pending9的12个原始point key交空。新增要求的H128对账也实际通过：095全部13406个point/time=1/128与098 batch0001冻结194条记录交空，194条epsilon逐项恰为1/128。绑定098原point_manifest SHA `293622876e07bb998227bad5b1a2f14af56a9727d17e3242023a43ca77d7d06c`；103追加overlap说明SHA `11b66fe5ca5fcdde12b7f3ccc443c7a49501d7ed21455a4ec35d3f34f7e4bc8d`仅是来源索引，本轮从两份原点集重新准确求交，没有仅接受作者结论，也没有改095锁。

## 5. 11个新stencil的不同路径

预登记12个selector去重后得到11个样本：首/末key、最小normalized投影gap前2、最大D^4Σ|stored weights|前2、到120个body两方向span并集的准确平方距离最小前2、H512保存下端最小/上端最大、H64复用首/末。最大系数变化量只是选择诊断，不冒充真实log项抵消比。选择和全部指标见sample_selection。

每样本全部3504锥，共38544锥，独立展开整数多项式P(X)=Π(X−node)，Horner计算P'(node)，用−det/(2P'(node))构造partial-fraction系数后折叠node²。此计算路径不同于生产直接乘分母和作者递归DD。保存全部folded weights逐个准确一致，0/2次矩消失、4次det矩恒等式均通过。没有重求任何H或log；其余13395个新点的完整stencil仍共享作者生产/全量递归读回资格。

## 6. 全1024前缀及最终I

冻结事件hash链1024项全部核验。每事件只安装完整父children；所有未替换父仍保持初始self且children为空。本轮对受影响root的整棵子树重新递归，再从52360个roots重求全局总和，区别于生产沿祖先路径累加delta；每项global sum、root delta、volume=1及未完成父付款完全一致。最终69545 nodes的核心字段与冻结state逐项相同，68201 physical frontier与冻结state及portable frontier清单一致。

共同分母N=2^160，最终

  LH=38413887677479752724438204644362255442485903/N，
  UH=75063710980782909270112604401565474168288679/N。

LH只属于截时∫H^(5/4)下证。076完整eps512短时范数E逐字同源；全F直接积分，不另加064tube。准确核B>E、(B−E)^5≥UH^4、Iupper^4≥B^5、costupper^4≥(5B/g)^5及完整有理32门margin。简单cap54702517/10^12严格优于014接受的6010267/10^11；32门false。本轮不据cost显示下降更改正式Cn。

## 7. 实际执行与交付

本机Python标准库整数/Fraction主PID94282、session34528，exit0，69.721909916秒含加载；ledger阶段26.928083542秒、registry1.462575208秒。全部13406新raw与新数学对象一次实际运行，无失败、无恢复；0新H/log/exp/Arb，0MATLAB。作者978606新H生产值明确共享，32旧H64只读精确复用。

登记、协议、源码、8项ratio结果、矩结果、真实registry结果、全部1024前缀trace、有限PF选择/结果、总结果与汇总收据是本次必要白名单。约31MB逐raw诊断checkpoint已汇总必要统计/hash后清理；不把可重生成的大诊断复制进仓库。恢复运行没有在本次验证，收口结果存在时程序拒绝重复成功基线。

复核依据为check.py和registration绑定的同SHA输入；原raw大文件及缓存保留在producer外scratch，不声称仅精简本目录就能重建原H数值。旧094/087及其它旧锁最后复核不变。仅新增本目录，未提交、暂存、推送或修改生产代码；主保留终审、一般L1传递及正式常数整合责任。
