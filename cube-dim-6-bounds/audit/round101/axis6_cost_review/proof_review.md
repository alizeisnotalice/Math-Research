# A101：等权六轴 Green 正部上界构造族的独立证明审查

结论：接受实际核质量成本 `66<U<12000` 与完整真实face积分
`1/100000<I<1/300`。这是上界构造族的必要障碍，**绝不是C6的下界**。
用户已停止C6下界推进，本审只审查上界候选方案是否可能让其质量成本小于64。
不由该结果推出C6>=66，也不排除其他核、修正、不同轴权或对同模型的其他上界机制。

## 1. 同一真实模型与global corner floor

六条单位坐标轴、同一Cauchy时间t的实际密度是
`p_t(x)=prod_i t/[pi(t²+x_i²)]`，不是body32、augmented38或alpha3/2。
每因子对|x_i|单调减，所以对完整[-1,1]^6以及每个t>0，
`p_t(x)>=p_t(1,...,1)`。正积分给G(x)>=g，包括允许G=+infinity的点。
没有通过只查有限corner样本推出全域floor。

corner处的Beta公式是已有精确handoff，本审不重新评价旧参数点。
Gamma半整数递推的有理系数为
`(1/2)*(15/8)*(3/4)/120=3/512`，带一个pi，故
`g=3/(512pi^5)`。这里1/2来自t²换元，Fourier/密度的pi^-6未混用。

corner邻域所有坐标都离0有正距离，small-t统一O(t^6)、large-t统一O(t^-6)；
可积支配给G连续。cube内正体积邻域逼近corner，因此inf与essential inf都恰等于g。
这补的是旧数值pilot没有证明的全域最小性质，不是重算它的65536 Sobol点、96个时间节点或corner数值。

有k个非零坐标时small-t恰为正倍数t^(2k−6)，故真实无限Green集合是support size<=2。
这些坐标轴/二平面及交叠保留在后面的全空间分布与积分中。

## 2. 完整face、全部时间的有限上包

q=6/5、I=int_[0,1]^5 G(1,y)^qdy、B=I^(5/6)。
对正时间truncations做Minkowski后单调收敛，得到
`B<=int c_t(1)[int_0^1 c_t(s)^q ds]^(5/q)dt`。
这不是点wise kernel upper的错误Lq交换。

0<t<=1时，以真实半线概率质量1/2和sup c_t=1/(pi t)插值：
`int_0^1 c_t^q<=1/[2(pi t)^(1/5)]`。
取25/6次幂并乘c_t(1)<=t/pi，得到
`2^(-25/6)pi^(-11/6)t^(1/6)`；0端可积，积分为6/7倍常数。
没有额外乘/除32的orthant因子，五个一维积分中各自的1/2已经在此指数内。
t>=1时全密度<=pi^-6t^-6，face体积1，尾积分1/(5pi^6)。
因此所有时间和所有奇异邻域都已支付。

只读E098已认证的准确有理pi区间给157/50<pi<22/7。
`(157/50)^11>8^6`、25/6>4及pi>3严格推出
`B<3/448+1/3645=11383/1632960<1/128`。
这是真实整个face范数的上包，不是旧QMC或单个分区预算。
结合准确g，`6B/g<8(22/7)^5`；独立有理幂残差
`12000^5-[8(22/7)^5]^6>0`给真实族成本U<12000。

## 3. 真实mean/Jensen下包，方向没有倒置

face体积恰1，Holder/Jensen给B>=A=int_FG。真实product与正Tonelli给
`A=int_0^infinity t/[pi(1+t²)] * [atan(1/t)/pi]^5 dt`。
固定t∈[1/4,1/2]，1/t>=2>sqrt3，所以主支atan(1/t)>pi/3；窗口外贡献非负。
由此 `A>log(20/17)/(486pi)`。
对x=3/17，`d[log(1+x)-x+x²/2]/dx=x²/(1+x)>0`，
所以log(20/17)>93/578。它是解析严格lower，不是拟合log数值。

成本下包必须先用**同一准确g**消去pi：
`6A/g > (512/243) pi^4 log(20/17)`。
然后才能用pi的lower157/50和log的lower93/578，得到
`R_L=301356307696/9144140625`。
不能拿I的一个lower再除以g的lower来伪造cost lower。
本审独立核保存完整有理残差`R_L^6-66^5>0`，所以真实U>66>64。
这个下端是**实际核质量成本**的下端，来源是真实A和B，不是某个upper budget太大。

另用pi upper22/7，得到A>217/2059992；准确比较
`(217/2059992)^6>(1/100000)^5`，及
`(1/128)^6<(1/300)^5`，给完整真实I的粗严格区间。
这些I与U不等式均已由独立Fraction整数幂实际核验；旧750.7412显示没有进入任何决定。

## 4. 全空间概率源与六维frame接口逐项迁移

axes的symbol满足|xi|<=sum|xi_i|<=sqrt6|xi|，L2生成元域是H1。
热核上界常数是C_axis=pi^-6，不能继承body32的1/(160pi^6)。
实际正、保质量、自伴、order-one scaling和n=6>1满足001的一般finite positive spanning-Cauchy接口。

对测试phi，坐标二差分给Lphi∈L1∩L∞∩C0；
`||(P_tphi-phi)/t-Lphi||∞→0`由integrable Fourier multiplier给出。
`int G|Lphi|<=int min(||Lphi||∞,C_axis t^-6||Lphi||1)dt<infinity`，
故全空间半群导数配对合法，LG=-delta0。没有从singular planes的补集延拓点值方程。

有限I及cube-radial正积分直接给R_c=(G-c)+∈L1；G∈L1_loc也由transient短长时拆分成立。
Tonelli给`G=P_tG+U_t` a.e.，正部凸性与1-Lipschitz给`P_tR_c>=R_c-U_t`。
用非负测试、除t并极限，得`mu_c=LR_c+delta0>=0`为正Radon测度。
大cutoff由||Lchi_R||∞=O(1/R)和R_c∈L1回收总质量1。
写h=min(G,c)有mu=-Lh；小cutoff的L1缩放eps^5及h<=c给mu({0})=0。
整个论证保留全部真实无限Green层，不省略任何坐标二平面及交叠。

001的障碍逼近在axes只需把热核常数换成C_axis：离散势上包是
`F+(6/5)(C_axis M)^(1/6)F^(5/6)`；频率尾由psi~|xi|控制，
局部紧性、u的L1回收和大cutoff质量守恒的维数指数仍为6和5。
所以对bounded compact f>=0，有`f=sigma-Lu`、u>=0在L1∩L∞∩H1，
0<=sigma<=kappa、sigma=kappa在同一个D={u>0}、int sigma=int f。
这里没有移用body32的具体数值常数或不同维的源结论。

## 5. 同一坏集控制所有common radii与general L1

对0<c<g，`K_c=R_c/[64(g-c)]`非负even且>=1_cube/64。
order-one scaling给
`K_(c,r)*Lu=(mu_(c,r)*u-u)/[64(g-c)r]`。
在同一个D补集u=0，右边非负，因此`K_(c,r)*f<=kappa int K_c`。
先合并正有理r的可数null sets；bounded f与L1 dilation continuity给
`||K_r*f-K_s*f||∞<=||f||∞||K_r-K_s||1→0`，故同时扩到所有正实r。
最后单调截断|f|给一般L1，且kappa|D|<=||f||1。

cube kernel是(2r)^-6 1_[-r,r]^6，六轴在这里共用同一r。
product-Cauchy semigroup没有变成六个自由尺度的strong maximal operator，
也没有使用六个一维maximal operators的composition。

## 6. 所有可用spatial cutoffs的精确质量成本

G(rx)=r^-5G(x)。l-infinity径向坐标的Jacobian是r^5；
64个orthants与6个最大坐标face共给384I，ties为零测集。
对每方向的正radial积分，直接得到
`int R_c=384*(5/6)c^-1/5 I=320c^-1/5 I`，因此
`int K_c=5I c^-1/5/(g-c)`。
本审核384→320→5，没有遗漏face、sign或体积64。

准确essential floor使64(g-c)正好是此正部核可用的最大常数normalizer。
c>=g时essential inf R_c=0，不能以有限正常数使它dominate cube；
c<=0时核不可积。因此这族可用的空间截断正是0<c<g。
端点质量均发散。设s=c/g，第五幂的可变因子是5^5/[s(1-s)^5]，
分母导数为(1-s)^4(1-6s)，唯一最大在s=1/6。
独立多项式系数/有理数检查验证该factorization与最小factor第五幂6^6。
所以真实理想族成本恰为
`U=6^(6/5)I/g^(6/5)=(6B/g)^(6/5)`，下端>66排除全部此类截断的质量成本改进64。

这是固定等权axes6 Green空间正部的核质量上界机制的失败证书。
`C6<=U`与U>66绝不推出C6>66；也不排除其他核、修正、weights、
time-truncation或更细的其他上界机制。当前C6上界64保持。
