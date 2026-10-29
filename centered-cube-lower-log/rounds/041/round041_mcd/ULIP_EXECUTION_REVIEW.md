# factor6最终执行的独立算术复核

结论：通过。使用自己实现的log-odds Cantelli/Markov和factor6指数公式，未导入root的upper_math或uniform_upper，未取新样本、未重建source、未反演任何CI。输入/统计均为同一40轮数据；FP与原mean/BE/r0/Lipschitz资格仍保留。

复核脚本 `uniform_offline_check.py`，实际输出 `uniform_offline_check.json`。FINAL_RESULTS的SHA256为 `6ac83c0fe19269df55c6520d2ffed86ee7dc9769f8e643dfb15e0b202310d00b`，FINAL_PLAN为 `80e8650d5469ad6c502000d458d98947c4f2c4d7ffdce4707f49ff6a5c5170b7`。

## 全量核对

- 对4R、3模型、N256、J233、21阈值，共15031296行重算per-scale variance/3 Cantelli/Markov与uniform6 minimum，再作J union。
- 指数严格先向上截到−700，再乘R（fixed0乘1）；S0单独处理，positive gap时给0，strict等号使用冻结实现保守的1。
- floor调用总10584119次，逐R为2723299、2626946、2664834、2569040，与FINAL全部一致。
- 从保存的旧lower和logD恢复B并计算φ(B)=3B/(1+2B)。逐点sharp lower≤uniform6 upper≤variance-only upper≤old upper，所有贡献帽均无正违反。
- 与保存factor2 upper比较最大正差4.2633e−14，属于独立浮点算式的roundoff容差。向上floor本身允许极微小正差；其解析额外conditional上界≤J·terms·exp(−700)，乘λW后仍可单独计入。故没有要求不可能的bitwise factor6≤factor2。

独立重算FINAL的1008个profile统计row（old/sharp lower、old/sharp upper，每项mean/SE/ESS/max），最大绝对差4.8317e−13。factor2与variance-only诊断均值也吻合。84个paired mean/描述SE及simultaneous endpoint减法最大差8.8818e−16；全部按相同h、相同counts行配对。这里仅核对保存CI端点的减法，没有重新验证CI公式或覆盖率，CI交另代理。

## 数据身份与边界

FINAL_PLAN的四个input NPZ hash、receipt、aggregate/uniform源码、uniform证明及method-uniform CI证明hash均匹配。逐R直接比较旧40与base_replay的三条保存stream hash，以及全部旧payload键/数组，均相同；upper inputs的h和logD与旧40逐位相同。未保存的原始count/source stream字节没有在此重新计算；其同协议身份由已接受observer重建及hash链承载。

C固定4.05，source r0为10^−6，receipt记录的所有已检查r约≤2.512e−7。uniform6依赖连续uniform隐藏坐标与整bin Lipschitz帽；不能只凭振幅套用。当前source h=δ/6使同一C·log product surrogate有效，数学独审已另完成。

向上floor只防这一指数分支因exp underflow变零，不代表normal CDF、Markov、线性代数或完整FP获得interval认证。报告与代码均维持该限制。

## 同h读数

在各joint sharp-lower经验最佳h处，下面的prior upper仍明显更大：

|R|h|joint sharp lower均值|prior uniform6 upper均值|
|---:|---:|---:|---:|
|4|4|7.680632|10.140977|
|8|4|5.137640|7.116614|
|16|3.5|4.370863|6.393516|
|32|3|4.589723|6.561977|

全部四R的同h paired MC正差数和保存endpoint正差数均为0，与FINAL一致。这是当前有限grid上下包未分离，不是证明gated来源无潜在收益，也不是完整weak排序。

实际复核耗时.371秒，Python3.15/NumPy2.5.4。没有更改alpha、冻结文件或生成新CI/MC。本次交付仅为算术/身份/摘要复核。
