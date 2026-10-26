# A062：高端正事件的完整Green预算

对原body32、rho=1/64、q=6/5及全部31个unsigned tubes的mean_B>3/4高端union，作者证明与两个不同算术路径给出严格预算

    ||G||_(L^(6/5)(H))
      < 27779670841319/1000000000000000000
      = 0.000027779670841319.

H可进一步扣除主空间分区中的cap、low与mid；删除这些集合只会降低范数。本结果支付全部正时间和两个真实参数事件，允许31条管任意交叠。该预算约为条件性全face阈值4/48750的0.338564738379，但未支付middle长时、far或global g，**不是全域I或新C6界**。

主流程已审读并放行理论；本包的纯有理保存检查也已PASS。该检查由本构造者编写，属于不同算术路径复核，不能代替下一步未参与构造的代理独审。最终数学证明与总账集成仍由主流程负责。

## 新机制与准确范围

[theorem.md](theorem.md) 保留执行前的完整证明，固定新阈值b=1/4。它与旧042帽预算的数值批不同，也不使用060的ones-first Phi表达式。

在真实高端管内，每个坐标>=47/64。对补事件|a_ones|<1/4，剩余参数所需位置y=x-a_ones*ones每个坐标>=d=31/64。这给出全部剩余方向的有效下界

    |y|>=sqrt6*d,
    dist(y,span v)>=sqrt(10/3)*d,
    dist(y,span(v,w))>=sqrt2*d,

其中v,w均不是ones。证明无需枚举：正坐标的线距离平方从d*ones向坐标正向增加不降；两body-line span只含原来的两条body lines，故删去ones后任一pair plane至少有一个mixed-sign支持块，两个相反符号坐标提供平方距离>=2d²。

真实事件p_t=P_t+B_t分别为|a_ones|>=b和其补集，其他31个参数始终不受限制。P使用原Q5/coarea，得到

    P_t<=D*t^(-4)/(b²+t²), D=R1/(32pi^6).

B在31列上正删参三次。先在受限事件内提走三个距离分母，再保留f_t(a)把正a积分扩大至全实轴，恢复原Q3峰。因此全部26970个ordered triples仍由已接受的

    Crem=44743079/(2970pi^6)

支付，分母为(6d²+992t²)((10/3)d²+960t²)(2d²+928t²)。不能颠倒这一步骤或在扩张后的a上继续声称原gap成立。

真实face Jacobian为sqrt(km)，四维正规球的maxcone比例为1/k，高端mean区间长度1/4。由正volume union得到唯一整体上体积

    V=(pi²/8)*rho^4*sum_(k=1)^5 binom(5,k-1)*sqrt((6-k)/k).

这是允许clipping与交叠的体积上式，非实际体积等式。P短时使用P<=p取得完整密度的正face切片L1上界t/pi，不假设P自身反射对称，也不按31条管重复支付同一个P分支。

只取一次松式解析最优epsilon=((D/b²)*pi*V)^(1/5)。P完整长尾采用准确的

    1/(3b²epsilon³)-1/(b⁴epsilon)+atan(b/epsilon)/b⁵,

B则使用三个正二次因子的闭式积分覆盖0<t<infinity。重复正根通过连续性延拓，没有小根差除法。两个P时间区间与整个B事件经Minkowski合成，没有遗漏长尾或补事件。

## 保存结果

以下为展示值；认证使用results.json的完整向外dyadic端点和checked_results.json中的准确有理比较。

|项|展示值|
|---|---:|
|唯一高端union体积上式V|2.1041299455419554e-6|
|唯一解析epsilon|0.002346804812033683|
|P短时范数上预算|1.3274306818353634e-5|
|P从epsilon到infinity的完整尾预算|5.1608857464946701e-6|
|B整个正时间轴范数上预算|9.3444782764687584e-6|
|完整高端预算式|2.7779670841317063e-5|

纯Fraction路径独立包围了相同固定公式，其总upper经向上整理为

    27779670841317063046273263/10^30,

严格小于本包声明的27779670841319/10^18。最终声明也被准确比较为严格小于4/48750；后者只作为条件性整体预算参照，不能据此宣布全face总账完成。

## 门控、运行与180秒限制

原registration、theorem与main_release未改变。主release SHA为56350004c9fae96e85a1bfe62f65f8fc2a974b6b3a5a1412a3a275c8afeb9ff1；execution_binding在2026-10-26T09:44:23.961312Z、任何科学执行之前补登记producer、checker和run_guard的SHA。

严格执行已登记9个生产配置与1个纯检查配置，没有额外配置、参数修复、失败或重试。4个结构/符号配置核几何常数、一般b的tail及optimizer身份、三正二次因子积分及重复根身份、真实体积Jacobian与31权重。5个Arb768配置只计算本轮固定V/常数、epsilon/Pshort、Plong、Bfull和合成上预算。

|程序|PID|tool chunk|退出|受守卫墙钟|程序内计时|
|---|---:|---|---:|---:|---:|
|唯一producer|64682|93502a|0|0.378558417秒|0.067766041秒|
|唯一pure checker|64698|fede0d|0|0.130960208秒|0.104525334秒|

两个程序都由run_guard.py以新的专属进程组启动，并实际设置communicate(timeout=180)。若超时，只对该组SIGKILL并wait/reap，保存收据和已写检查点，禁止自动重跑；本次均正常完成并已回收，无活句柄。理论准备和主等待不被误称为一个180秒计算程序。

纯检查不调用Arb。它采用固定Machin120/40项有理交错夹逼pi、固定96项小参数atan(epsilon/b)夹逼，再用atan(b/epsilon)=pi/2-atan(epsilon/b)核完整尾。分数幂采用60位十进制网格的整数根，每个上下端独立核对准确幂不等式，累计26次端点检查。它从准确R1、Crem与固定几何参数独立包围三个预算，未使用060失效的“根与其幂同相对外扩后互相包含”逻辑。

旧M积分、旧042帽预算、旧060预算、kernel/native、网格、全方向/triple枚举、MATLAB、递归agent、stage/commit/push均为0。没有环境安装。所有固定来源在运行前后保持原SHA。

## 文件与后续独审

results.json保存生产区间，checked_results.json保存纯有理夹逼结论，run_receipt及两个process_receipt保存真实执行与180秒守卫信息。preparation_record/preparation_manifest和theorem开头的“尚未执行”是原准备阶段历史快照；当前状态以本报告、最终manifest及执行收据为准，未改写历史记录。

原执行入口为

```sh
python -B run_guard.py producer
python -B run_guard.py checker
```

冻结目录禁止重跑；run_guard拒绝覆盖既有process receipt。独审应先只读证明、准确输入和保存端点。若另有授权重放，须使用新输出目录及独立登记，保持相同物理参数和来源SHA，并维持每程序180秒限制。源绝对路径记录本次实际只读位置；异机需恢复相同哈希文件后显式创建新的路径绑定，不能修改本次生产release或收据。

下一步由未参与构造的代理审核真实事件、全部remaining-subspace gap、正确恢复Q3的顺序、真实volume union和完整时间支付，再由主流程纳入不遗漏空间账本。本包不提前声明全I、global g或Cn突破。
