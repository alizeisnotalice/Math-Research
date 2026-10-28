# R080：加权帽补突破原覆盖损失

FULL ACCEPT：所有整数n且自然L=log n≥4096，Lebesgue中心轴平行共同半径立方体极大算子的弱(1,1)常数满足

    Cn≤Phi80=F(A_n(t),B_n(t)−(39029/100000)n)<T−(5473/20000)n，
    t=n/[(89/100)L]，T=nlogn−nloglogn，F(A,B)=(sqrt A+sqrt B)^2。

机制是直接估计源帽补的加权捕获积分：原count事件保证nu*n个可用坐标，任意sourcepin至多删1，预帽补inverse product law给指数小概率；bsoft≤theta/G在匹配宽积分中消掉theta，得到capLoss/k≤M exp[−(nu−1/n)eta*x/(1+eJ)]。旧unweighted估计的全band宽因子和eta−chi限制不再需要。一般覆盖构造改进，不涉及七维专项。

固定s3/2、qmin2/3、qmax247/250、sigma11/100、chi51/1000、D1536，Ret999999/1000000、mL196691/500000，gamma39029/100000、price5827/50000、c5473/20000，完整标量余量1/10000。全机制、六target条件、真sourceface法则、无限维数尾、同history finitea、一般L1和全尺度完整证明见weighted_threshold/PROOF.md；主审核及独审均接受。

准确比较：Phi79−Phi80>2073n/20000在L4096全域成立（相同A,B,t）。通过保留原B价格的正下包，Phi80<Phi78在L≥8192成立，Phi80<Phi69在L≥2^30成立；没有把旧界门槛迁入新域或只比模型成本。其余已接受界仍按各自真实域取min。公开比较沿用已有文献核查，不宣称公开纪录或全局最优；本轮是内部严格改进。

实际新运行：terra先登记10点条件批，首次成本变量残留导致错误FAIL并单独保留；登记修正后同批9点PASS，qmax.989因真实weighted cap超额被排除，.988固定候选通过53项有理门。sol另用N32/log12/S32/精确sqrt4096核18组新守卫；独审另用N16/log20/P40敏感算术并全文审核。主新增R079同参数完整成本精确比较。所有程序输出、失败与修正身份分离，原件hash不改，旧通过基线未重跑。

归档只选固定候选、一个失败边界、必要复核代码/输入、证明及独审，原全批与初始错误批保留本地。最小复现：python3 weighted_threshold/verify.py --compute-dir compute（在本轮目录）；或python3 compute/case.py 988/1000。复现会产生新输出，请在副本执行以保留已归档运行身份。论文不改，研究继续，下界暂停，七维专项停止。
