# 第058轮：更宽来源带的2097152项有限模型搜索与准确单点区间证书

新机制问题是放宽sourceband而非重复R052旧tuple搜索。令vmin=a*d、vmax=b*d、qmin=1−vmax、qmax=1−vmin、Rmin=qmin+d，phi=−log(1−d/(2Rmin))，mu=(Rmin/qmin)phi²，eta=vmin/(2qmax)。理想inverse-row模型的费用节省定义为

    c_model=log((qmax+d)/(qmin+d))*qmin/qmax−2mu。

这里忽略了实际goodset丢失、individual/productJ误差、有限clip保留、正tau、sourcecap余量、trim、所有有限history及operator守卫。因此即便c_model增大，也不能称Cn改善。mu<eta只是这个理想clip/sourcecap指数可兼容的条件；真实固定tuple需重新证明全部量。

有限注册i,j,k=0..127，d=1/5+(2/15)i/127，a=1/16+(3/8)j/127，b=9/16+(27/64)k/127。15动态workers实际完成128chunk、每块16384项，共2097152不同有理参数，80位mpmath两log计算，无异常。报有1469673项mu<eta，其中1414566项c_model>0；这些符号及列表排序为数值判断，未逐项区间化。每块保留top32、全部计数/PID/耗时/源hash及完成身份；128原始UTF8收据无损压缩，可独立合并4096条retained记录与top64。

有限网格最佳点d=391/1905、a=283/2032、b=63/64，展示c_model约0.1017216430123904761，mu约0.0146412731，eta−mu约0.00007202025。b处于所扫上边界，不能称连续或全局最优。保留top64在每chunk至多5项，少于top32；这支持保留策略，但完整排序仍依赖producer确实执行每块且其80位排序正确，省略的所有tuple并未独立重算。

该单点另外以48项有理atanh展开和几何余项给准确上下界：对y=(x−1)/(x+1)∈(0,1)，log x位于2sum_{j=0}^{47}y^(2j+1)/(2j+1)与其加2y^97/[97(1−y²)]之间。源参数从准确分数原件读入，严格得mu_upper<eta与c_lower>0，不从展示小数重造。该certificate证明完整理想标量式，仍没有真实覆盖族的参数域/一般L1/allpositive-scale结论。

实际墙时3.93s、workerCPU合计50.43s，平均约12.83核。批次很短，只有0.5s一次ps样本：15workers、1052.6%=10.53核、453456KiB≈443MiB RSS；没有5s峰样本，不根据worker数声称持续满核。后补portable摘要checker只验证保留hash、计数/合并和单点有理证书，不重跑2097152项。Native Terra独立不同Fraction实现审单点、128收据和top合并，另登记五个新角点(含被排除c<0点)；没有用producer helper冒充独立审计。

代码、注册、准确certificate、每块收据、summary、monitor、独立审阅/原始运行receipt及来源hash精简归档；没有全tuple大表、缓存或阶段稿。公开与内部两种比较均无新Cn。R055/R053等完整上界仍用各自域和全部根式取min；R059拟从本候选建真实宽来源带证明，未闭合前不能把c_model约0.10写进Cn账本。七维专项和下界保持暂停/停止，首项系数小于1没有证明。
