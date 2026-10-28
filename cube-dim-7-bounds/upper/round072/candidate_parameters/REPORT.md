# R072 · L=4096的条件参数检查与有限支付改进

## 结果身份与缺口

给定基础候选的42项共同算术gates及完整cap/tau/energy/scalar支付均通过。单进程运行PID99908，墙时约0.002秒；这只是新的精确有理单候选/12行有限表，不是旧网格或生产扫描。

**目前不能据此声称任何新的C_n界。** 缺少已核准且适用于独立非同分布数组的Berry–Esseen定理及常数C≤1；R070第二profile引理在原运行时只有作者传来的接口、尚无最终可绑定源hash；本tuple还需要完整operator适配审。下文所有正态保留率及Gamma结论保持这些条件。R066已接受的原机制不能自动提供本tuple的新保留率。

原始verification.json诚实保存“H2 source pending”“BE source missing”，不追溯改写。既有R066/relative-coordinate源仅为公式和域来源，不是新tuple的算子证书。

## 1. 基础候选：准确预算

固定L0=4096，s=3/2，qmin=2/3，qmax=91/100，sigma=19/200，chi=9/200，x=(181/200)L，delta=2^(-19)，rlo=1+delta，rhi=273/200−delta，h=1/(512L)，theta=exp(−chi*x)，k=nh theta/rhi。

固定eBaseline=2^(-16)、eJ=2^(-20)、eN=2^(-19)、eRad=epsD=2^(-20)、epsA=1/65536、epsRet=1/3、mL=311/1000。

P0=1−eBaseline−6/512，Ulo=P0(1−eN)/(1+eJ)，Uhi=(1+eN)/[P0(1−eJ)]，M=2(1+eRad)Uhi。

基础epsS=3/40，epsT=1/4096给

b=2(1−epsRet)(1−eRad)Ulo−epsS−epsT，
Gamma=mL*b/M
 =15757510577752772683186277517/82529650777166176415068979200
 >19/100。

更强的可支付条件也准确通过：(1−epsD)Gamma>gamma=19/100。不是从显示小数判断该不等式。

新的七项odd alternating log1p上包给price=log(200/181)<999/10000。完整根式的剩余正项按Taylor8上包为

R8(L)=2/L²+8!/(sigma L)^8+2sqrtL*8!/(sigma L/2)^8。

L0处R8约1.19212·10^(-7)<1/10000。因gamma−999/10000−9/100=1/10000，基础display c=9/100的完整scalar支付通过，不能只比较gamma扣减。

## 2. 两个profile及正常近似的精确条件接口

所有新tuple几何/浓缩检查使用N(L)=L^8/8!≤n；N与N/L递增。相对几何为R=sq−h、g=(s−1+epsilon)q−h，epsilon=1/(xn²)。准确基点检查保证h<delta/2、epsilon qmax<2h/3、active R>1、rho≥4r、g/R≤1/3及所有normalizer/count/radial-width正域。

Relative-coordinate lemma实际还给|∂q logJ|<3/rho；matching宽4h/3及n/rho≤(3/2)L给product bridge≤6Lh=6/512。因此本P0使用较锐的已证局部常数，而非R066为展示保留的256/rho粗常数。第一profile的8/rho导数给shift≤16xh≤x/100000；individual J精度与负ell_pin导数方向亦按新tuple单独核验。

目标集合需要六项坏事件：baseline logJ、两端count、原cap、第一profile和第二profile，各按epsA/8预算，总计6epsA/8<epsA。新20项even alternating target log lower实际给

(1−epsA)log(rhi/rlo)>311/1000。

最后一项第二profile的解析输入在原运行时尚待R070最终源：作者传来h≤2^(-20)、n≥2000、x/n≤1/1000时，baseline H2和的mean在[.00528,.005295]x，variance≤n/(125²rho²)，故双边[.00527,.0053]坏率≤5760000/n；匹配shift≤xh/5。新基点满足这些finite guards，shift/x<10^(-6)。最终source-face方差必须用真实中心化式

Var(V|face j)=Σ_(i!=j)[B_i/J_i−(A_i/J_i)^2]，

且删除任意source pin必须支付B_j≤1/(125rho)，不能用target pinned坐标的半界。传来的R070接口在这些guard下给

(13/2500)x≤Var(V|face j)≤(531/100000)x。

在最终源/独审到达前，这一行是明确的条件输入，不由数值gates提升为定理。方差范围只用于每个captured source face，不能套给face mixture。

第一profile给每face统一mean上包

E V≤183/1000+[2211/(50000(1−eJ))]x=:mbar(x)。

局部Qpin≤Q5=34187/187500<183/1000，所以每个free坐标phi∈[0,183/1000]；pinned贡献确定，不计variance。在L0，准确检查

(chi*x−mbar)^2≥(3/5)^2(531/100000)x，chi*x−mbar>0，
(13/2500)x≥19。

于是条件于该方差范围，standardized threshold z≥3/5。第一比值随x增大而增：其非恒定部分为a sqrtx−b/sqrtx，a,b>0；variance lower也随x增大，故这两个条件覆盖全尾。

**待补的BE假设**是：对这组独立free坐标，适用误差≤C·ΣE|X_i−EX_i|³/variance^(3/2)，并且C≤1。由有界幅度，Σ三阶中心绝对矩≤(183/1000)variance，所以若此假设成立，误差≤(183/1000)/sqrt19<21/500。这里不是把尚未有合格文献的BE结论当作已证事实。

正态CDF本身的下包不需要BE假设。正积分恒等式

∫_0^1 t^4(1−t)^4/(1+t²)dt=22/7−pi>0

给pi<22/7，而(199/500)^2<7/44。结合exp(−u²/2)≥1−u²/2，得到

Phi(3/5)>1/2+(199/500)[3/5−(3/5)^3/6]
 =90559/125000。

减21/500后为85309/125000>2/3；因此在BE及H2条件真正获证后可实现epsRet=1/3。整个论证先逐face、再混合其保留概率，不把独立variance错误地赋给混合律。

## 3. 正项S64支付及无限尾

记S64(z)=Σ_(j=0)^64 z^j/j!≤exp(z)，所有计算为完整正项有理前缀，没有浮点指数或有限采样外推。

Csrc=(qmax−qmin)/qmin，eta=9/182。N≥10000时准确eta−chi−eta/N>1/225，故在实际k归一化下

cap/k≤(Csrc*rhi/h)/S64((181/45000)L)。

准确基点upper约0.07312124504561815<3/40。宽band的Csrc*rhi因子被完整保留，不能无依据用旧nd，也不能把此因子丢为1后仍声称3/40通过。

完整选择阈值tau=sqrt[A(B−gamma n)]与theta共同给

tau/k≤(rhi sqrtL/h)/S64((271/40000)L)。

基点upper约0.000162590938306846<1/4096。完整root正域来自旧B≥T−4n/L，L≥64、x∈[L/2,L]；A≤ne^(−sigma L)、B≤nL。剩余fullroot使用第1节R8，不是遗漏scalar price后的原扣减。

对于任意正项prefix，mu(z)=zS64'(z)/S64(z)是权重z^j/j!下j的mean，dmu/dlogz=Var(j)≥0。新准确基点分别核mu((181/45000)L0)>1、mu((271/40000)L0)>3/2，所以L/S64(aL)与L^(3/2)/S64(bL)在全部L≥L0严格递减。这证明cap/tau全尾，不把前缀当成只在基点有效的exp近似。其余finite guards为明确负powers；scalar R8按L^(-2)、L^(-8)、L^(-15/2)递减。

## 4. 登记的12行支付表及改进候选

只改变epsS和epsT，几何、h、两个profile及正态条件接口全部保持不变，故没有自动外推未证H2域。表为epsS∈{.075,.074,.073,.072}、epsT∈{1/4096,1/5000,1/6000}的12行，每行均保存原pass/fail，没有删除失败结果。

结果：6行条件算术PASS；6行FAIL，失败均是epsS=.073或.072小于准确cap/S64 upper。所有三种tau预算通过。

将每行(1−epsD)Gamma向下取到1/100000有理格点，登记表中的最高安全展示系数由以下一行达到：

epsS=37/500，epsT=1/5000，gamma=19109/100000，
c=9109/100000。

准确(1−epsD)Gamma>gamma，且gamma−999/10000−c=1/10000>R8(L0)。epsT=1/6000同样达到这个rounded c；此处保留1/5000一行，以给tau更宽余量。未声称此表外最优、连续全局最优或真正Cn改善。

该条件c=.09109严格大于基础.09；改进只来自本固定几何下更细的支付预算与gamma余量使用，不代表新算子机制已经获证。

## 5. 明确的最小余量、来源与停止条件

MARGINS.json给准确Fraction最小余量，并明确以1−required/available归一化。它是在指定critical catalogue和四项支付检查中的最小值，不是对未登记条件或所有参数取极小。

基础候选的catalogue最小值来自target mass lower，约0.0004704414490342264，准确>47/100000。改进候选最小值来自energy budget减所选gamma，准确为

19283255242054881602967286603559/
2067102182483004369477860757363515175
 >9/1000000。

两者的完整exact ratios、负的失败行余量和code/receipt身份均已保存。显示小数仅供阅读，所有pass/fail与选点来自Fraction。

本轮没有资格跳过BE文献/定理核准、H2最终源锁定和新tuple完整operator/history/generalL1/allscale适配。完成这些之前，报告及每行receipt均维持conditional身份，不进入Cn账本。若后续H2/BE任一接口不符，必须撤销对应条件结论而非沿用数值保留率。没有MATLAB、Git、旧baseline、生产重跑或n7工作。
