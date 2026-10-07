# N03 蒸馏与局部审计

审计者：GPT-6.1 Sol。本记录读取真实证据卡并执行所述有限演练；full_read只记录源阅读与覆盖，不是独立证明认证。来源行号指转换原文，准确卡/页码由provenance回查。

## 原文结论与适用量词

P-b1a1dc66a1923470 Thm3.1 lines423–441/proof553–691：固定局部 对 Newton 多面体实非退化的实解析相位 f0=gradf0=0，截断支撑足够小；λ−1/df log^{k−1}非整1/df、整时log^k，上界非渐近等价。Lemma2.3 lines273–395单项式models另含振幅权重。

## 证明骨架与已检查推导

局部toricmonomialization→Jacobianpowers与faceindices→非零横向导数→局部vdC/sublevel预算→有限坐标图求和；toric/Greenblatt依赖未独立证明。

本次实际检查：∫0¹e^{iλx}dx=(e^{iλ}−1)/(iλ)，bound2/λ；x²驻点0一阶下界失效，localchangeu=√λx给λ−1/2尺度。

## 正反例与失败处理

正例：对 lambda>0，令 I(lambda)=∫_0^1 exp(i lambda x)dx。相位导数恒为1，且 |I(lambda)|=|e^(i lambda)-1|/lambda<=2/lambda。

条件缺失例：对相位 phi(x)=x^2 在 [-1,1] 上声称一阶导数法给出 O(lambda^-1)；但 phi'(0)=0，均匀导数下界失效，驻相尺度是 lambda^-1/2。

## 中心立方体迁移推断

用振荡衰减和显式几何误差控制有限尺度分解后的残差，并逐项展示可求和余项。

迁移前逐项给输入映射、工具前提为何成立、界的方向、n/scale常数和连接引理。公开原文没有自动提供这些接口。

## 待验证命题

constantCφ 固定相位与截断，无uniformdimensioncoefficient；half-density/Kaufman私有definition未给，指数对A/B旧源另待核。 Newton版本还须确认相位非恒零，Newton距离与主面确有定义；不能对空多面体机械套指数。

相位导数为零处不能直接使用一阶 van der Corput 界。遗漏端点、Jacobian 或驻点区域可能使余项主导最终估计。

行为结果只认证记录中的有限输入；部分演练通过不表示高级定理全部证明或一般cube主张验收。旧文逐段审计保持局部scope，不认证同段其余结论。
