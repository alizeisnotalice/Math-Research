# 122：保留062构造的CEH预算回收

四个固定配置与唯一saved-only核验均已实际通过，现封存停止。回收只来自H的真实补事件；原kernel、空间分区、rho、事件阈值、gap、Pshort、Plong、解析tau/epsilon、体积及C/E预算全部不变。

## 结果及准确账目

旧062/111 paidCEH准确I上端为

    996880729065386184584374670023140717152990867507023323/200000000000000000000000000000000000000000000000000000000000

新CEH准确I上端为

    25773099890807522001737635999358247103027591282385388967249739077652164733383373634675967546583/5444517870735015415413993718908291383296000000000000000000000000000000000000000000000000000000000000

严格回收的完整有理量为

    delta = 1364574831131708916327094813880161730534735223049043357480170240058603124182097300539055516457/5444517870735015415413993718908291383296000000000000000000000000000000000000000000000000000000000000

展示值为旧CEH约4.98440364532693e−6，新CEH约4.73377083200356e−6，delta约2.50632813323375e−7。证书与所有比较使用上面的完整有理数及results.json完整端点，不从这些小数构造。

H的变化如下；C、E准确端点仍直接读取原证书。

|量|旧保存预算|新预算展示值|
|---|---:|---:|
|H补事件完整时间范数|9.34447827646876e−6|7.63113762289608e−6|
|H完整Green范数|2.77796708413171e−5|2.60663301877444e−5|
|H账户I上包|3.40778757108703e−6|3.15715475776365e−6|

完整111 middle没有重算，其原far充分额度两端同时加同一个exact delta，更新为展示区间

    [4.7877381304874166e−6, 5.0702288908180828e−6]。

完整原始有理两端在certificate.json和results.json cfg4中。该额度仍条件于全域g>=1/65000，是充分分配，不是actual far积分、下界或全局S。

## 回收机制与适用域

062的真实参数事件先给出逐点恒等式 p_t=P_actual,t+B_actual,t，B_actual是原a_ones落在|a|<1/4的非负subdensity，其余31参数无限制。因此在任何删除、triple union、空间union与gap放大之前，已知B_actual<=p_t。

在原high union，每个坐标>=47/64，旧受限事件gap d=31/64给出B_actual<=Bbar。Bbar是放大的26970正triple上包，绝不声称Bbar<=p_t。独立的108原law full6D峰给p_t<=216/(2^30π^6)t^−6，因此合法新包是

    B_actual <= min{Bbar, 216/(2^30π^6)t^−6}。

这是对真正补事件的两个支配式取小，不是给deleted law生搬BL。它覆盖所有正时间、全部其他大参数与cancellation事件，以及原H允许的所有tube交叠。

两条上包的比值在u=t²上严格递增，存在唯一交点T，实际保存区间约T=0.01586005200804845。T由原界强制确定，不是扫描或新选的tau。192次精确有理cubic bisection及平方端点检查支付交点。新积分为0到T的三resolvent有限primitive，加T到infinity的完整rank6尾。Sympy准确核了partial fractions符号；Arb512在两个新端点计算atan primitive，双端有理包围保存完整时间积分约0.41054631172846717。旧完整B积分、P积分、volume与epsilon未重算。

把这个新B积分只乘原保存V^(5/6)，再加原保存Pshort/Plong得到新H范数。严格比较通过FULL旧H公式范数下端与旧H账户I上端两道检查。C/E/H最终相加用059真实优先分区的互斥物理集合；旧证书可覆盖更大的raw union，不要求这些raw union互斥。

## 直接投影替换诊断

另一个配置准确比较了033的1个rank5、3个rank4、5个rank3真full-projection峰与108通用峰。全部9类的BL/033比值平方都>1；直接替换均更弱。相同正系数与M(A,C)单调性意味着不能从这些直接rank替换回收C/E/H。

这些比较在原law恢复后才合法：P的大参数line coarea所删ones本身投影为0；C/H补事件在拉出受限gap并保留同time Cauchy因子延拓a后，才恢复真原law Q3。E的037 Q4/Q3中已删除方向都在所投影span内。未恢复的31列deleted law不满足108完整tight-frame条件，未用BL替换。新rank6峰直接用于B_actual<=p，未插入rank5或rank3密度槽。

## 实际执行与复核

主全文数学/source gate后的parent_release SHA为0320ec76272cdf92bf2c30616de62c7072089a3ec838649195bd646aa1e77e59，随后才运行一次四配置。tool chunk5dc759、exit0，科学PID73840、observer73832。完整科学子进程从spawn前到最终receipt、解释器退出与reap后，wall0.2143243330065161s，CPU0.211136s；observer自身的写收据/退出不计入该child范围。最多一科学CPU，无时间上限、无worker/native，未占用其他任务的生产队列。

唯一保存核验为tool chunk053882、exit0、PASS_SAVED_ONLY。它不导入producer、不调用旧程序，精确检查九类比值、cubic signs、平方端点、保存primitive与双端积分公式、第五根不等式、旧C/E/H账及delta同加。复核入口：

    /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_saved.py

31项sourcepins结束时未变，包括原062十项原件pins。0旧passed网格/峰/M/full-B/P/volume/epsilon重放，0native、MATLAB、其他frame、下界算子研究、新LLM代理、stage或push。一般L1与全部正共同尺度仍依赖同一001/027正部源接口；actual far、globalg和完整S仍未闭合，不能称Cn突破。未授权新批次，现停止。
