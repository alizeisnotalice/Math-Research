# N02 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-8cdc5c07bc1a147e lines96–108/159–161是discretePL/BBL，lattice与非退化假设不同；continuous边缘PL需明确连续版本，不能将discrete转换无条件使用。

## 证明骨架与已检查推导

写h((1−λ)x+λy)≥geo度量/pmean of f,g→非负、可测、可积→准确dimensionadjustedmean→integralbound；边际对数凹性用该sourceversion逐层积分。

本次实际检查：λ.5 f=g=h1[0,1]前提及integral1；f1[0,1],g1[3,4],h0 在x.5,y3.5中点2 violates0≥1。

## 正反例与失败处理

正例：取 lambda=1/2 且 f=g=h=1_[0,1]。只要 f(x)g(y)>0，就有 (x+y)/2∈[0,1]，从而 h((x+y)/2)>=sqrt(f(x)g(y))；三个积分均为 1。

条件缺失例：取 f=1_[0,1]、g=1_[3,4]、h恒为0、lambda=1/2。令 x=1/2、y=7/2，则中点为2，h(2)=0<1，故 Prékopa--Leindler 前提不成立。

## 中心立方体迁移推断

用边际分布保持对数凹性，将多维几何问题降到低维密度或截面计算。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

等号不自动稳定性，pdomain及zero-conventions须跟版本，cube截面/weights不由labellogconcave认证。

插值不等式方向或 p-凹约定写反会使结论反向；仅有“对数凹”标签不足以在可积性不满足时推出边际公式。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
