# M04 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-55ddb3291a49ef12 lines254–258/399–402/547–552是semi-analytic/Borel dominatinghyperplane选择；P-a594f1076246b512 ThmA/B lines217–252引用externalPettis选择器，Thm4.3 lines330–393是标量可测 cwk(X) HKP 可积性，非一般最大点选择。

## 证明骨架与已检查推导

给可测图像/Carathéodory与非空紧值→按准确maximumtheorem取值函数可测性/最大点映射的图像→选择器；若还需integral须L¹/Pettis各自conditions。

本次实际检查：K(t){0,1},φ=(t−.5)u，u0左/u1右Borelmaximizer；非LebmeasurableA对应唯一选择器1Ac不可测。

## 正反例与失败处理

正例：令 T=[0,1]、K(t)={0,1}，并最大化 (t-1/2)u。取 u(t)=0（t<=1/2）、u(t)=1（t>1/2），这是达到最大值的可测选择。

条件缺失例：取 [0,1] 的非 Lebesgue 可测子集 A，并令 K(t)={0}（t∈A）、K(t)={1}（t∉A）。唯一选择器不可测，说明不能省略多值映射的可测性假设。

## 中心立方体迁移推断

只有在积分中需要把逐点最优解组成随机输入或参数依赖见证时，才调用可测选择器。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

逐点存在、tightcompact、pointwiseLipschitz selection不自动general measurablemaximizer；本源external theorem未独立证明。 另实际回读P-a594原文lines260–393：Definition4.1区分范数HK与弱HKP；Theorem4.3(ii)打印HK而证明只给HKP，外引[39]待核前不把该证明称为强Henstock等价。

逐点存在、光滑有限维选择定理或值集紧本身，不足以在缺少图像可测性等条件时得到可测选择器。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
