# 第064轮：乘法匹配与有界目标径向校准，完整原算子上界T−n/5

最终源17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6、FULL ACCEPT独审443699173985a45738e535fcf7e0b2407c67f06b5d55832e59f361e1793a55af锁定。主全文重推并核对后接受：存在本固定tuple自己的N，所有integer n≥N有

    Cn≤Phi64=F(A_n(t),B_n(t)−(61/200)n)<T−n/5，
    t=n/[(909/1000)logn]，T=nlogn−nloglogn。

无显式N，不能迁用R063L33或R056L20。当前已量化域仍为R063 log n≥2^33的n/10、R056 log n≥2^20的n/40。Lebesgue、中心轴平行立方体、一个共同半径、一般L1输入和全部正尺度保持。内部未发表，首项nlogn系数仍1，不声称公开纪录或最优。

## 新乘法几何及两种有界权

固定s3/2、q∈[2/3,909/1000]、target r∈[1+2^-30,2727/2000−2^-30]，sigma91/1000、chi1/20、h2^-64/L、theta=e^(−chi*x)、x909L/1000。改matching为|r−s q|<h；参考R=sq−h≤r、a0=(1−zeta*t)q、g=R−a0。active matches由trim给R>1，首层inactiveedge可R<1但其column0。相对g/R≤1/3由zeta*t*q≤h/s精确付出，g绝对可>.333，不能套旧absolute d≤1/3域。

SOURCE权bsoft=min(1,theta/G_R)在完整sourceband先定义，Dcap删除随后支付；TARGET新权a_rad(w)=r/rhi在targetannulus且≤1，外延0。最终K明确1A1P1matching1capture*a_rad*bsoft，两个权都≤1保真实eligibility。baseline q_b=(r+h)/s使R=r，整个qmatching宽2h/s。

新local R'=3/2、g'=1/2+epsilon不同于旧dconst：稳定Psi/J公式新推a.e导数<3/rho(安全用256/rho)并AC，free/pin baseline矩重新按R≥1、g≤R/2、rho≥4r核，不将oldsmall-d范围越用。g/R≤1/3和r<1.5保证这些常数，captureonset有真实kink只作a.e/单侧/AC，不称处处平滑。无条件allw导数使新A无需Kw事件。

sourceface Nsum随q升而ellpin' =−(1−epsilon)h/(2R²)为负，故Nsum/(2nellpin)递增，endpoint准确均值及两浓缩事件控制整个matching区间。baseline+两count+Dcap四个fulltargetface坏事件uniform o1，保留λ(A)≥(1−epsA)nlog(rhi/rlo)>mL*n。逐点J→1和productJ由a.e/AC分别控制，不能用n个individual误差粗乘替代产品。

## 同一核的完整径向行列预算

源逆权全face系数仍U=(prodJ)/(2ellpin)ΣNi/Ji，Ulo*n≤U≤Uhi*n。精确source logwidth是log((r+h)/(r−h))；乘目标r/rhi后，width位于2h/rhi*(1±e)。另一方面，完整target-face系数的每坐标Ti/(2ell_i)≤1，给nG_r≤nG_R；目标a_rad使(dr/r)*(r/rhi)=dr/rhi精确取消Jacobian，ordinary targetmatching r宽≤2h。因此column≤2nh theta/rhi=2k，k=nh theta/rhi；row按同一个k下界2(1−epsRet)(1−e)Ulo*k，上界M*k，M=2(1+e)Uhi>2。没有漏s或rhi、没有独立平均column预算。

fullsource phi=−logell的准确矩不变；Qpin≤log6/5，用五项odd Taylor上包Q5=34187/187500<1/5，R/q≤s给sourceface mean≤1/5+(999/20000)x、var≤x/150。x≥8000给pin/x≤1/40000，chi−999/20000=1/20000，逐face Cheb坏率≤40000²/(150x)→0，之后混合概率。不能把mixture variance当独立和，也不在Dcap后假独立。

宽sourcecap用Csrc=(qmax−qmin)/qmin=727/2000，删除未加权mass除k有必要因子Csrc*rhi/h。准确eta−chi=91/1818−1/20=1/18180>0，tau/k指数sigma/2−chi(1−sigma)=1/20000>0，二者固定幂乘正指数衰减。所有geometry/targetprobability/retention/cap/tau条件uniform于全半径/面，有限固定条件共一个未知N，不从有限probe推无限尾。

## 准确完整成本、真实history及新旧排序

所有误差与baseline/radialwidth tolerances为e1/1024；P0=1−e−2^-50、Ulo=P0(1−e)/(1+e)、Uhi=(1+e)/(P0(1−e))，mL=1239/4000，b=2(1−e)^2Ulo−2e。原准确64项atanh证书给(1−epsD)mL*b/M>61/200。Phi64扣减gamma=61/200；完整scalar价格先付91/909，raw37249/181800，对1/5保准确slack889/181800。A、cross-root及1/n的全部正余项除n为o1，付清后严格得T−n/5，不只模型价较好。

source ri≥rho>4、actualtarget rj<3，原samehistory中真实source先于同一target且center guard对a_alg≥8n^4成立；sourceb和targeta均≤1，actualweighted row≤a_rad(w)*tau≤tau。lambda_a=ca lambda在A/P精确成立，energy/degrees/mass按ca²/ca/ca缩放、ca≥1使缺额不减，保gamma n。原finitecover/真实f witnesses及内正则性给一般L1/allpositivecommonradii，最后只取统一标量a_alg极限，不假设selectors/history收敛，不使用n≥5奇异子空间、低维线外连续或小时间核可积性迁移。

旧R061的完整cost下包保留B正buffer：Phi61≥T−n[429/2000−3/50+4/L]，L≥2^20括号<1/5且递减；旧Phi59≥T−n[1637/12500+4/L]也严格高于T−n/5。因此Phi64<T−n/5<Phi61/Phi59在各自真实operator域交集成立，未知N不因scalar排序变显式。这不是Cn下界、也不是两个展示扣减排序替代完整F。主新增3个对应准确gate实际PASS，不重跑旧完整验证程序。

## 审核、实际运行及归档

作者64项log新单tuple准确证书PASS；roots新relative lemma及SOURCE_ADAPTATION冻结c60，star局部完整row/column/count/actualμ与ca接口另7个新Fraction例PASS，full Sol独审另用40/41交错log1p、41个不同Fraction/SymPy检查实际PASS并通读全部operator链。主全文核数学与冻结身份，无待修。所有pending/candidate字样保留生成时状态，只有最终FULL ACCEPT与本报告激活。

R065另有真实15worker230400有限状态及1024 direct/256 baseline核、初scalar误抄tau符号的真实错误receipt及修正版；数值原件和scalar/facts审计已单独提交87904a9de0961cba253397fe9021fccbc57769df，不能以它们代替以上全scope证明。它的保守mL3097/10000不同于此作者1239/4000，来源身份分开。没有重跑passedbaseline、生产grid、n7或下界，未用共享MATLAB。

精简归档只最终证明/原准确code与receipt、relative基础、独立radialdegrees及full审/hash，不复制重复R065大raw或阶段draft。当前24页R060论文冻结至R056/R059，本轮新机制及R061/R063后续结果先见各自报告，后续统一论文修订。goal active；R066新增target first-profile仍独审，R067数值结束不代表理论接受；旧资料只读、旧暂停项目及lower/n7专项保持原状态。
