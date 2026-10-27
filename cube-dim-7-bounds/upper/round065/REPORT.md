# 第065轮：乘法匹配相对核的新有限诊断与错误符号审计

本轮为R064候选的独立计算证据，不能以有限passed点激活Cn。新的几何路径是s=3/2、R=sq−h、a0=(1−epsilon)q、g=R−a0，故R'=3/2、g'=1/2+epsilon；与R057固定absolute d导数不同，不重复旧通过网格。它允许g绝对大于1/3，却控制相对g/R≤1/3，验证目标为allw的a.e/AC导数安全界256/rho。真正全域相对lemma、径向row/column和原算子结论各由R064完整独立审承担。

## 条件性准确单点

准确s3/2、qmin2/3、qmax909/1000，trim2^-30，targetrlo1+trim、rhi2727/2000−trim，sigma91/1000、chi1/20，所有e=1/1024，P0=1−e−2^-50。Ulo=P0(1−e)/(1+e)、Uhi=(1+e)/(P0(1−e))，M=2(1+e)Uhi，rowb=2(1−e)^2Ulo−2e，保守mL=3097/10000。48项有理log加几何余项证gamma_budget=(1−e)mL rowb/M>61/200，Qbar18233/100000>log6/5、mu=sQbar²/(1−e)<999/20000、varrate=(2s/3)Qbar³/(1−e)<1/150；eta−chi=1/18180、sigma/2−chi(1−sigma)=1/20000。raw61/200−91/909=37249/181800，对1/5余量889/181800。

初次scalar源错误将tau里的(1−sigma)抄成(1+sigma)，产生假失败；原错误source及真实JSON(receipt中的tau=false、与错误source hash配对)完整保留，修正code/certificate是另一身份。不是数学候选排除，也不把修正结果追写成初次运行。Native Terra独立40项有理atanh与余项重建全部scalar通过，不导producer helper。R064作者另用mL1239/4000与Q5而非此mL/Qbar，两种保守证书可同时真，不能强行称同一run或相同gamma_budget。

## 新15进程有限导数批

注册128动态chunk×1800状态，共230400有限状态，覆盖新q/rho/h、pin与captureonset单侧邻域；closed160十进制位，独立piecewise source积分/单侧差分256十进制位(大于512bit)，不是512十进制位。另有1024直接积分/差分与256个baseline均值核，实际128chunk全部完成0异常。最大观察|rho*d_q logJ|约.7794256218<256，最大缩放独立值/derivative误差约9.20e−181/6.36e−73、baseline约1.37e−258；均为有限观察值，不是区间全域误差或最大值认证。captureonset须按a.e/单侧/AC解释，不因有限差分通过就称处处可微。

wall19.7459s、workerCPU247.9498s，平均约12.56核；15worker真实采样峰1448.8%=14.488核，RSS合计峰1910224KiB≈1.82GiB。峰值是采样，非全时持续，也未宣称满内存。q和epsilon是局部核数值参数约束，未将这批有限kernel例子说成某个巨大n的完整operator实验。

## 审计与复核范围

Luna逐个核128原chunk/hash与gzip原UTF8一致、全count、source/registration/summary及实际资源，FACTS ACCEPT；错误符号的真实receipt身份也实际核实。Native Terra只另审scalar精确式，不代审几何。原source+register+summary+monitor、128chunk无损gzip和两版scalar、独立audit/registration/receipt均精简保留。新portablechecker在干净temp exit0，主复制后核retainedbytes/counts/准确字段，不重跑积分/旧baseline。归档check不是新生产证据。

公开与内部两种比较均无本数值轮独立Cn界；R064完整同history finitea、generalL1/allpositivecommonradii与完整F证明通过后才给原算子意义。模型价好、有限导数无反例、单点准确证书都不能代替这些全scope要求；不称发表、全球最优、leading系数小于1、七维或下界进展。所有old项目只读和原暂停状态保持，goal active。
