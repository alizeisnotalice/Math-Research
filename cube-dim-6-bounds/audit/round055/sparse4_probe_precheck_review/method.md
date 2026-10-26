# B055 独立预检方法与边界

这份文件先于实际结构检查准备。C055 稳定源码和输入的 SHA 锁定前，不执行新数学检查。主代理独占生产 full 数值放行；本审计不调用作者 writer、编译 binary 或任何完整 G/H。

## 依赖与完整性

F050 已独立验证 21 父锥细分及 13229 个完整有序锥、478 条原始射线。这里不重做 active-system ledger 或 fan 覆盖。绑定其已审 primitive determinants，再核生产输入的每条射线、每个 simplex 和每个权重，足以检出删项、换 fan 或错误 scale。

F052 已独立验证一般正满秩六维 frame 的有限部分公式。生产 scalar G 必须首先通过该新 frame 的真 S2 检查；原 body32 资格不能继承。公式只针对标量，不延伸到零节点梯度。

## 六个实际检查配置

1. **frame。** 独立组合生成 32 个六支撑方向与 120 个四支撑方向，后者幅度 1/80。以分母 80 的整数列与冻结 frame 全比对。检查 Gram 与六个 coordinate Cauchy rates。取五个相邻置换和第一个坐标的符号翻转，逐列验证 mod 全局正负的加权 frame 不变；这六个生成元生成完整 signed-permutation 群，故所有群作用都保留真实符号 ψ，不必执行 46080×152 的重复验证。
2. **真 S2。** 对全 11476 个 unordered pairs 用独立整数 Gram 行列式比
   \[
   \operatorname{dist}(x,\operatorname{span}(a,b))^2
   =\frac{\det\operatorname{Gram}(a,b,x)}{\det\operatorname{Gram}(a,b)}
   \]
   精确计算。方向的非零缩放不改变 span，可用未加权整数代表、将点分母准确放在最终比中。使用一般 Bareiss 行列式，区别于作者 normal-equations 展开。两个行列式计两个 matrix operations，共 22952；pair 阳性即包含所有 line 的资格。额外 line 残差只作标量参考。任何 pair 为零硬拒绝 Green finite part。
3. **群和轨道。** 独立对固定点枚举全部 720×64=46080 个作用，保留零坐标的全部符号稳定子。检查每像计数，global-negative 无自配对，代表通过“首个非零坐标正”选取，不能只检查首坐标。核 zero-coordinate 在六个位置的分布以及 orbit×multiplicity 恒等式。
4. **完整输入与 scale。** 每原始射线的真 rate 是 152 个带幅度绝对线性形式之和。若 q=p/d，R=d r、Q=p，则 b=(R·x)/Q=(r·x)/q。逐列缩放的 determinant 乘以各 d 的乘积，故
   \[
   \frac{|\det(R_1,\ldots,R_6)|}{\prod Q_i}
   =\frac{|\det(r_1,\ldots,r_6)|}{\prod q_i}.
   \]
   逐条核 478 rays 与 13229 simplices，继承 F050 primitive determinant，不做新 13k matrix determinants。六个 rates 不能以 int64 求乘积；weight 要任意精度或逐次 Arb 除法。
5. **标量语义与归一。** G 的 sign-value functional 无 π 时归一为 1/(128π5)；若 core 已把每个非零 sign value 乘 π/2，则应与 H 一同除以 (2π)^6=64π6，不能再次乘 π/2 或重复除 π。H 使用真实 atan(b/T) normalized confluent jets。所有节点排序、相等性和正负必须准确；零节点和重复节点不跳过，分母区间含零要失败而非容差合并。固定同一正 T、Arb512、仅 scalar。
6. **driver。** 审完整 source/input/binary/runtime 绑定、预算守卫、逐字段完整逻辑项与实际 tuple 数、记录出口和失败语义。对于给定点，正确计数应逐字段为 13229×46080 个逻辑项，global-neg 压缩后为 13229×1920 个实际 tuple，每代表 weight24。一对 G/H 共享一次 tuple traversal 可以，但仍有两个 scalar 字段；不得将实际 traversal 乘24后冒称新算了全逻辑项，也不得把两字段合成一个未区分的物理标量值。

## 必要条件结论边界

唯一正时间和唯一点须从稳定注册读取，不引入新参数。K=G−H 的单点下端超过 T/4096 仅通过必要点约束；它不证明全域 m_T。若严格上端不超过 T/4096，才可排除此固定 frame/T 的小于64目标；T/4000 对应62.5目标，二者不可混淆。任何结果均不直接改变 C6。

本审计只提供 launch 前的代码与结构意见，不赋予生产实际运行的 PASS；主代理仍须独立放行，生产输出之后还需另行审计。
