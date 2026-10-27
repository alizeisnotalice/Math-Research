# 第052轮：7372800个准确有理tuple的有限搜索与保留证书

本轮本身是当时conditional的有限参数模型，不把目标改善自动变成Cn上界。R051完整generic机制及R053/R055原算子门槛已另轮接受后，对满足域的指定tuple才能激活Cn结论；原注册、实际run、源code和results的conditional身份均未倒写。后续激活身份单列ACTIVATION.md，主可从完整数学源核，不从本轮排名推断算子定理。

有限列表a=A/8,A=32..256，rho=t/16,t=1/2/4/8/16/32/64/128，s=2a+rho，d=k/2^20,k=1..4096，共7372800 tuple。P0=3579/4096，epsJ/S/T/D1/1024、epsA1/64。准确目标是gamma−sigma/(1−sigma)，所有预算按R051公式生成。M=4先用整个registered区间的单调端点分数不等式证明；其余公共分母由符号恒等式重导，生产ranking纯整数交叉乘，非浮点排位。

1800动态块每个固定A与512个k、遍历8个rho，producer实际计数3975632合法预算tuple、3397168严格d-guard拒绝；sigma/bracket无拒绝。这里合法并非说每tuple都给正c。最高保留值a73/8,s293/16,d1023/2^20，c_raw=6468305130094850733787980147/370780099533215320770594052505600。原全结果仅336KiB，作为sole准确bundle直接保存，未再压缩或从小数复造。

完整保留64个global候选及15worker各32个证书；producer的另一个Fraction实现实际重构这些证书并验证。独立Terra审阅者另从原公式重构全部top64、域、严格c_raw>2^-16、总数及guard的225项公式，并用SymPy证明公共分母恒等式、在6个新边界点做独立Fraction核。producer自检不是独立作者，两类receipt分别归档。

每块仅保32本来不能一概保证global64；本数据有max同块入榜21<32，使给定每块完整且正确保留top32时，可逻辑排除被截断项进入global64。但原每块完整结果和独立执行ledger未保存，所以不能独立证明所有块确实无遗漏执行一次或正确截断。完整性与finite-grid最优仍条件于producer的冻结运行；该局限不影响一个准确可行tuple依R051得到Cn，也不证明连续参数最优。

实际15worker，wall0.4428145s、累计CPU4.198576s，首0.5s采样前结束，monitor为空，不能声称瞬时满核；相较R047已改动态细队列。登记的全有限列表完成停止，不为占CPU额外重跑。7.37m全枚举只执行一次，后续portable_check及干净临时结构smoke仅复核保留证书/计数，不再枚举，不更改原receipt。完整bundle、原run、注册/依赖hash、代码和最小portable checker均保留。

source原report为独审引用hash故原字节保存为R052_REPORT_producer_original.md；新的R052_REPORT.md只澄清复核命令措辞，未改数学或实际run。公开与内部Cn比较由R051/53/55完整rootcost和真实域负责，本轮不产生独立的新Cn界，不借失败或被拒tuple声称数学不可行。下界及七维专项未启动，旧资料和项目状态不改。
