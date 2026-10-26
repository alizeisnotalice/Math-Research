# A101：axes6上界构造族成本独审

**PASS_INDEPENDENT_ACTUAL_AXIS6_UPPER_FAMILY_COST_OBSTRUCTION**。
接受真实equal-axis Green空间正部族的严格成本`66<U<12000`和完整真实积分
`1/100000<I<1/300`，并接受已限定修复的v2 guard在本次固定执行中的默认硬fuse设计与真实收据绑定。

这是**上界构造族的必要成本障碍**，不是upper预算太粗而失败，也**不是C6下界**。
真实global floor是`g=3/(512pi^5)`，族内唯一最优cutoff为`c=g/6`。
全空间正源、general L1及所有common radii迁移的详细审查见`proof_review.md`。
用户已停止C6下界推进，本审没有推进任何lower-bound项目。

接受范围仅为六条单位坐标轴、共同时间、alpha1、同一个G的空间正部与既定核质量机制。
不排任意其他axis kernel/修正/weights/time-truncation，也不从U>66推C6>=66。
原C6上界64保持；无需为缩小已明确非竞争的宽上端再扩大计算。

## 冻结来源与数学核验

作者目录：`/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round101/axis6_cost_preflight`。
final manifest SHA `72672d8892ce141d3ddb248c9e319bd1f5d33c059a1cd3199f6b1112ea60a02a`；
theorem SHA `3b963c3449e23078687554828dc9874eac3e188595aea5c340fdf682bfe41e5a`；
results SHA `4ab4229489fd73a62fec9c80829b9ae9e539d3c947f1fbf6e0e8e93b101d982e`。
全部20个artifact、22个readonly source pins及9个execution pins逐项核对。
原科学child SHA `2aaf588f0fe28cff165f61a609a27aab7a5f6bb37056d0747dae7ecf04533340`在修复前后不变。

独立最小程序没有import作者或旧数学程序，只读保存原件，用标准库Fraction、整数幂和static AST：

- 直接绑定已有corner公式与E098准确pi证书，检查Beta/Gamma有理系数；0旧corner数值评价。
- 核cube-radial的384、320、5归一化，准确多项式导数符号及唯一c/g=1/6。
- 核固定窗口primitive参数20/17、严格log下端93/578、真实mean/Jensen的相关pi消去。
- 独立重算保存幂残差，确认真实族成本>66、<12000及完整真实I两端，未使用旧750.74显示。
- 核限定guard diff、主74bf96放行、原失败/修复收据与fuse/CPU/PG限定、旧新源hash。

本审在proof review中重新确认corner附近连续性，故global inf与essential inf一致；
保留support<=2的全部真实无限Green集合，不只在其补集声称源方程。
axes heat-kernel常数是pi^-6，谱比较是|xi|到sqrt6|xi|，未混body32的常数。
全time upper在0和infinity都可积，完整face I有限，general-L1/allpositive-common-radii量词没有弱化。

## 原CPU setup错误及明确v2修复

原guard39093在Popen前设CPU hard2，child39094继承后不能升hard28，
在check_exact.py第14行、数学try之前抛ValueError。原配置0、断言0、stdout空、exit1，
child已reap；这不是timeout或数学反例。原源、原失败receipt和traceback均保留。

主在执行前静读并核22+9pins，chunk74bf96明确放行一次v2。
独立逐字确认diff只包含：spawn前soft2/hard30，spawn后parent hard2，及exclusive新receipt名。
child数学/参数、177/150秒默认fuse、160秒wait和5秒own-group cleanup不变。
这不是自动重试。原source和first failure没被重写成成功。

v2 guard39912与child39913真实exit0，3个科学配置、77断言；
数学CPU0.002046秒、child OS CPU0.026661秒、guard CPU0.029646秒。
guard内部0.032273792秒止于最后receipt写入前；
封存的原工具观察另以Date.now跨完整await测得0.091秒tool roundtrip，
包括receipt/stdout/os._exit/reaping及transport，不把它当成纯数学CPU或内部单调计时。

原首tool_receipt把raw wait字段/描述叫“complete”，与其实际约1e−5秒含义不符；
本审按最终报告/manifest的明确限定，把它只作raw tool wait，首尝试独立完整roundtrip未记录。
首内部0.043823875秒亦止于最终receipt前，不冒充完整过程值。
v2新roundtrip0.091秒独立保留，不用它回填首失败的缺失计时。

本审检查的是**完整外guard及child源码**：两者均SIG_DFL，177/150秒timer没有取消、重设或Python异常handler；
父guard的最后receipt/stdout/flush/os._exit仍在硬fuse内，cleanup等待有界。
这与099/100的取消总alarm情况不同。本次固定源码满足所审180秒完整程序限制设计；
没有运行一个新timeout故障试验，也不声称该设计能超越正常OS调度/内核终止能力。
仅对四个原任务PID作一次ps，均无行；不按名字扫描、不发信号，不据此声称全部孙进程已reap。

## 独审实际执行与交付

仅1个新saved/exact组一次PASS：原chunk `503818`，PID42189，exit0；
内部wall含结果保存0.027183166秒、CPU0.015578秒。
本审程序默认SIG_DFL硬179秒、CPU硬15秒，timer一直到最后receipt及exit，不cancel。
上限3组，剩余2组未用，之后停止数学。

0作者3配置/原guard/旧65536Sobol/96time/旧corner参数数值/核/MATLAB重跑；
0递归分派、stage/commit/push、原件写入。全部写入仅在本clone的`audit/round101/axis6_cost_review`。
090、099、100维持冻结，没有额外lower-bound工作。
最小代码、准确结果、真实receipt、登记、完整proof review和report由manifest白名单交主集成。
