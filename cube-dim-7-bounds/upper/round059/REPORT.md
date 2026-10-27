# 第059轮：宽来源带构造将完整原算子成本推进到T−n/10（eventual）

## 已证结论和准确完整成本

最终证明源fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f、完整独审c8aff39bf68749d25f9978cef11be82ccb7907ec7e72c5678ac66b06a65a33d4已FULL ACCEPT。主核完整机制及最终4hunk差异后接受：存在本固定tuple自己的N，对所有integer n≥N，在Lebesgue背景、中心轴平行立方体、各坐标一个共同半径、一般L1输入与全部正尺度上，

    Cn≤Phi59=F(A_n(t),B_n(t)−(1637/12500)n)<T−n/10，
    T=nlogn−nloglogn，t=n/[(97/100)logn]。

未认证显式N，不能把R056的log n≥2^20或任何有限采样点移来激活n/10界。R056的显式T−n/40继续在log n≥2^20使用；所有旧界依其自己的真实域取min。公开首项系数没有新断言，新结果仍leading nlogn系数1、内部未发表。

固定d=391/1905，vmin=(283/2032)d、vmax=(63/64)d，qmin=32429/40640、qmax=3760307/3870960，Rbase=122311/121920；target endpoints为qmin+d+2^-20、qmax+d−2^-20。sigma=3/100、chi=147/10000、h=2^-48/(Ld)、theta=exp(−chi*x)、x=(97/100)L。epsJ/N/S/T/D=2^-20、epsRet=1/4096、epsA=1/65536，baseline tolerance2^-20，P0=1099510579195/1099511627776。选择干净deficit gamma=1637/12500，严格低于完整准确gamma_budget；c_raw=121289/1212500，对1/10保留严格准确slack39/1212500。参数来自完整分数源件，未从展示小数重造。

## 宽来源带机制与必须重证的部分

完整sourceband上先定义G及0<bsoft=min(1,theta/G)≤1，再取P=D∩sourceband，目标A是trimmed annulus和good events交D，最终K=1A1P1matching1capture bsoft。A/P径向不交，inactive source列0。对active source，trim确保R=q+d−h≥Rbase>1，而非将wholeband的Rbase−h错误说成Rbase；全matching区间宽2h、baseline q_b=r−d+h恢复R=r。局部一维lemma使用已重新证明的d≤1/3域与64 a.e/AC，individualJ和productJ分开控制。

由于本d<1/4，Kw≤n≤4nd不成立，新增Kw=1+Bin(n−1,2d/r)好事件，Cheb坏率≤1/(4nd²)。baseline、Kw、两端count和Dcap共五坏事件uniform o(1)，全半径/每face先fulltarget law再交D，固定epsA后共同大N保留目标mass≥159477n/1000000。两端sourceface count单调夹逼统一所有q，取消pin半因子，Ulo*n≤U≤Uhi*n。

每captured signed sourceface先独立按capture/(ell J)归一：phi=-logell的准确free矩R Qpin²/rho、2R Qpin³/(3rho)。pinlog<1/9，meanrate<293/20000、varrate<1/900，chi−293/20000=1/20000，逐face Cheb坏率≤40000²/(900x)→0后混合概率，未误用mixture variance、未在Dcap后假独立。准确column因Ni/(2ell)≤1付2k，M=2Uhi/qmin>2。

宽sourceband全mass不一定≤nd，明确改为Csrc*n、Csrc=(qmax−qmin)/qmin。sourcecap/k≤(Csrc/h)exp(−[(1−1/n)eta−chi]x)，eta=vmin/(2qmax)，准确eta−chi=499871/37603070000>1/100000。tau/k的另一指数sigma/2−chi(1−sigma)=741/1000000>0。两者均为polynomial(L)*exp(−positive*L)→0；每个guard对固定tuple统一所有半径/面，有限多guard共一个未知N，不给自适应family共同N。

同一真实history部分保留：每a_alg≥8n^4，source真实ri≥rho>4，实际target rj<3而均先于target；center displacement由zeta−(beta+2)/(a_alg+1)>0支付。每个b≤1，实际row只付原同一target一次eligibility tau，无inverse factor。lambda_a=ca lambda在两support精确成立，energy ca²、degrees/targetmass ca、ca≥1只改善缺额。原一般L1 witness有限cover和内正则性保留全部正共同半径，仅对已统一成立的标量界取a_alg极限，不要求selector/history收敛。n≥5奇异子空间、低维线外连续和小时间可积性都未作为迁移假设。

最后付完整F：B≤T+3n/97+1，A+2sqrt(AB)≤nexp(−sigma L)+2nsqrtL exp(−sigma L/2)=o(n)，在准确正slack39/1212500内得到严格T−n/10。不是以R058理想c_model约0.10或单个budget提升直接宣布Cn改进。

## 完整新旧比较及来源身份

旧Phi56的B正buffer必须保留：由I_n(t)≤2t/(n−1)及−log(1−7/100)≥7/100，Phi56≥B56−gamma56*n≥T−n[gamma56−7/100+4/L]。括号在L=2^20准确小于1/10且随后递减。因此Phi59<T−n/10<Phi56在两个真实operator域交集成立；不能使用过粗的T−gamma56 n与展示1/10比较。其余旧R053/R055在R056已做完整根式排序，新界只在自己的未知N后进一步改善；较小n仍取各自有效界。comparison/receipt实际新增3个准确gate核余量与下包，不重跑旧基线。

作者原64项atanh检查的certificate291fc41...和certify332cd3...逐字保留，其当时dependency记录6d174...。当前最终lemma97430...数学内容相同，只来源身份header/末节变动，roots已核可逆恢复；作者全文重读当前lemma并给DEPENDENCY_ADAPTATION逐项绑定。不得把原receipt追改成对最终lemma的新运行。最终4hunk REVISION.patch只补P0/P/K定义、fullband b延拓与provenance，独审反向恢复06a3原SHA相同，无估计/参数变化。生成时candidate/pending字样保留为历史身份，当前结论由最终FULL ACCEPT激活。

## 实际运行与最小复核材料

作者新单tuple64项有理atanh证书PASS；Terra另48项有理interval单点证书，Native Terra用不同Fraction实现独立核其conditional scalar；当时均只证明条件标量，现在完整机制独审后才接入Cn。Sol完整独审另以40-even/41-odd交错log1p(不用作者helper)完成18个Fraction和3个SymPy新checks，完整一般L1/allscale/finitehistory和所有uniform eventual guards全部通过。主读完整源、最终定义补丁，新增3个比较gate，不重跑这些已通过程序或R058 grid。

主另实际调用共享MATLAB新脚本验证两个全参数count符号恒等式，并算理想模型局部gradient；code/receipt/tool输出分列primary_local_math，workspace/path未清空重置。该运行早于最终理论锁审，既不是新的生产grid也不证明Cn或最优性。MATLAB原脚本保留当时绝对outputpath，仅作运行来源；可移植最小算术核验另见comparison及各最终Fraction代码，不将建议重跑命令冒称实际新运行。

只归档最终证明、准确原证书、必要差异/来源适配、最小代码、独审与实际receipt/hash。旧06a整份源/重复论文/private草稿不复制，patch可逆恢复。未知来源未跟踪同名2文件不暂存/改删。当前paper仍R054冻结20页、范围至R055，R060正完整整合R056显式n40与本轮eventualn10；未经新全文审不替换论文。R061分层与R063显式门槛仍新任务，不提前激活其猜测。goal active，原资料只读，lower/n7专项及旧暂停项目保持原状态。
