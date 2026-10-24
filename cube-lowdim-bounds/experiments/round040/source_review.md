# 第40轮独立静态源码审查

最终GO：已全文审signed_once.py、PLAN及wrapper，无阻断，未运行任何模型或积分。

- signed_once.py SHA：`0c6531e82741ec785b4a0f66ba84a4630871c73cf84d9ba797858d589af2babf`。
- PLAN SHA：`62d5a208e1cc7a25713dd05388de5011e362d1813ac2285d29b1b07d6c493fe9`。
- wrapper SHA：`a1158395ca46585e80c18ea7b8782590c3736422876f2e94ba35fd9d562afa77`。

四输入hash、node661身份、原深度/local/physical、基线参数及整盒正外区且两坐标<R均准确断言。
已追踪原保存链，旧Psi Bernstein是未扣余项/未乘c的T1；AS及原行仅Fraction读取核对，不重算。
模块导入后设192bit，再构造Kern。唯一新数学调用是近角点及交换两次_Hint(...,Hab)，无F0/梯度或完整Psi函数调用。
局部acb/实例Hint代理原样转交参数、callback和analytic契约，finally还原；真实积分与Hint独立记录且各最多2次。
首非finite/异常/预算即停，无重试、提精度、换盒或其他积分器。返回事件先completed并保存，再guard，不重复计已返回层。
两个Hab绝对球上端经man_exp转准确分数，Muv求和再乘准确正半宽du*dv；其单侧合法性见analytic_audit.md。
新行同基底逐项组合，保持旧AS余项，索引/区间顺序/旧行身份全部核查；strict>0才accepted，不使用浮点容差。
完整保留原绝对扣项、旧行和差；记录deduction是否减小，但该布尔量不替代接受判据。
60秒外进程组硬限含准备/导入/两积分/准确行/保存，已有result/signed_model/meta/log拒绝覆盖，子程序另拒旧calls。
0新AS/完整Psi/const/F0梯度/点/细分/T/LP/全域checker，预算没有扩大。
数学证明仅给单侧下界，不宣称双侧绝对误差；成功只限这个闭盒，不能称新上界或完成3.4目标。
