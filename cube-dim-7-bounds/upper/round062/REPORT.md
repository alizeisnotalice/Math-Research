# 第062轮：不交分层的本机前沿探索与128层准确条件证书

本轮围绕新的层间费用计算，不重复旧R058全band最坏值网格。将source qband分成固定J层，理想模型逐层计算qlo/qhi加权target logmass；每层用自己的mu/eta兼容条件，费用只以max mu进入同一global scalar价格。模型目标为sum_j log((qhi_j+d)/(qlo_j+d))*qlo_j/qhi_j−2max_j mu_j，mu_j=((qlo_j+d)/qlo_j)[−log(1−d/(2(qlo_j+d)))]²，eta_j=(1−qhi_j)/(2qhi_j)。模型中没有付全概率/finiteh/finitehistory等真实费用，因此数值提升不等于Cn提升；R061另给完整分层算子证明和独审，二者身份分开。

## 固定128层准确构造

先对d=7/25、q从18/25至24/25等分128层、target每层trim2^-30，给完全有理区间证书。epsJ/N/S/T/D=2^-20、epsRet1/4096、epsA1/65536、P0=1−2^-20−1280*2^-48，Ulo/Uhi按原逆列公式，M_j=2Uhi/qlo_j。每层48项正atanh log及明确几何余项，向外取到2^-200，记录所有phi_upper、mu_upper、chi_j=ceil(10^6 mu_upper)/10^6+1/10000、m_j/b_j/degree/gamma_j。全部eta_j>chi_j、sigma/2−(1−sigma)chi_j>0(sigma3/50)，完整total gamma_lower约.2145684777。准确原件SHA b889e35ef909bd979937752c67410bf3913eb76cdf7dc6c720556e093d4c0711。

后补独立gate用同一原件精确确认total gamma>429/2000，保留原证书/run不追改；clean scalar 429/2000−3/47=14163/94000，超过3/20为63/94000。Native Terra不导producer helper，独立Fraction逐层重建128层、log interval与2^-200 round方向（包括微小负数）、mu/chi ceiling/两rate/全部budget，并实际PASS。此处仅接受条件性标量证书；真正Cn结论须经R061的single-history直接和完整定理激活。

## 两个有限15进程生产批次

pilot固定注册4096个组合、245760个分层项、128动态块；实际全部完成0异常，3098组合满足模型局部兼容。网格内最高d1297/4650、qhi747/775、J128，c_model约.154296979575。墙时.784s、worker CPU7.85s；仅0.5s一次采样15worker/1051.2%CPU/454320KiB RSS，不推称持续满核。

frontier_v2独立新注册：256个d中点在[53/200,59/200]、256个qhi中点在[479/500,121/125]，J256或512，共131072新组合、50331648分层项、256动态块；Fraction比对与pilot重复0。实际全部完成0异常，84808满足模型局部兼容。数值最高d14669/51200、qhi246503/256000、J512，c_model约.1555301880665；最小eta−mu约2.53382e−6。对这个已选点另做48项有理log及2^-160向外取整的512层准确证书，Native Terra逐层独立核全部512行和方向通过。它只证明该点模型为正/兼容，不是连续或全网格最优证明，更不是operator。

v2真实墙时136.02s、worker CPU1846.22s，平均约13.57核；15worker采样峰1458.8%=14.588核，RSS合计峰458416KiB≈448MiB。该136秒有意义数学负载已结束，没有以worker数量冒称持续满峰或满内存。注册/summary/monitor分别保留。原chunk没有源码hash字段，正确链条是注册锁源码、summary锁各原chunkhash；后来source-to-chunk manifest另标post-run，不当原raw运行字段。两batch数据和128/512层证书、postrungate互不混作一次运行。

## 审计、复核与作用边界

Luna对两组128+256原chunk逐个hash/计数/覆盖/retained top合并/代码guard和全部资源采样核对，FACTS ACCEPT。Native Terra单点两套准确layer证书另独立审算，代码和receipt分开。每块保留top32，summary合并top64；这些decimal排序仍条件于producer真实执行和其数值精度/局部保留正确，未独立重扫所有省略参数。边界/有限grid和理想model限制保留，不声称最优。

两组原UTF8 chunk无损压在一份gzip，完整固定candidate和原sourcehash保留，50m中间项不存盘、不上传大表、缓存或重复裸chunk。包约1.65MB，producer先在干净临时目录实际运行新portable summarychecker exit0；主复制后再核retained hashes/counts/top与原证书字段，不重跑生产。这种归档验证不是新的grid数学证据。原运行和postrun archive/checker身份分别标明。

公开与内部比较：本数值轮无独立Cn新界，R063完整显式n/10仍按自身域使用；新R061分层理论独审在单独报告激活，不能用此模型值替代它。一般L1、全正共同半径和真实history要求由理论完整证明承担。新型更大源带的模型成绩不得称公开发表、全球最优、leading系数小于1或下界进展。任务active，lower/n7专项及旧暂停项目不恢复。
