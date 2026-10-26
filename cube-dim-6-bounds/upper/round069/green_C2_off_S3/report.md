# A069：真 S3 外的 Green C2 与 zero/repeat 导数有限部接口

本轮给出解析证明：原 unit-weight body32，以及更一般有限非零正权 full-rank
六维 Cauchy frame 的 Green 函数，在**真 S3 补集上属于 C2**。据此可在严格
off-S3 资格下，把正 cutoff 的真实 copied DD6/DD7 导数取完整 Laurent 常数项，
合法处理零/重复投影。主线程保留最终证明审定与生产接口放行责任。

十项预登记配置一次通过，程序墙钟0.387030秒，PID73307，低于180秒硬限。
真实 tiny axis6 完整工作为六字段各64个 orthant 项，共 **384 tuple-field项**。
Body32 fullkernel、body32 全 triples 资格运行、旧 baseline/reference 重算、
MATLAB、delegation、retry、stage、commit、push 全为0。

## 解析结论与关键接口

取固定六列基 N，c=(1+||N^-T||F^2)^-1，S=sum_j||a_j||1，
R_delta=delta/(2S)。dist(x,S3)>=delta 时，任一表示x=sum a_jz_j中的大参数列
集合必须具有至少4的秩；若秩至多3，其和在真S3，余项范数却至多delta/2，矛盾。
这一步按秩贪心选独立四元组，适用于 dependent columns；没有假设随便选四个
大参数就独立。

四个尾密度给 t^4，二维正交商的峰为1/(2pi c^2 t^2)，因此

    p_t(x)<=M_delta t^2,
    M_delta=sum_(independent |J|=4) det(A_J^T A_J)^(-1/2)
                         /(2pi^5 c^2 R_delta^8).

同一个固定六basis的**完整**正混合独立给出

    |D_u p_t|<=C1(u)p_t/t,
    |D_uv p_t|<=C2(u,v)p_t/t^2,
    C1(u)=||N^-1u||1,
    C2(u,v)=||N^-1u||1||N^-1v||1+2 sum_i|(N^-1u)_i(N^-1v)_i|.

没有对任何 event indicator、rank selector 或 moving fiber 求导。两个估计
结合，并支付大时间t^-6-k尾部，给每个正距离邻域上的一致导数收敛及

    0<=G-H_epsilon<=M_delta epsilon^3/3,
    |D_uG-D_uH_epsilon|<=C1(u)M_delta epsilon^2/2,
    |D_uvG-D_uvH_epsilon|<=C2(u,v)M_delta epsilon.

这些不是样本推断，而是 `theorem.md` 的局部支配/一致收敛证明。Rational frame
可用binomial(m,4)d_A^4/[64 c^2 R_delta^8]作为完全有理的粗M上界。

## 有限部求导的准确顺序

先在epsilon>0对实际atan(b/epsilon)求 labelled node 导数：一阶插入一个copy，
二阶不同label插入两个copy，同label另有factor2。数值相同的不同label不误乘2。
再通过线性physical Jacobian组合，最后取完整fan的Laurent常数。

DD6最高负渠道为epsilon^-5，DD7最高为epsilon^-7。off-S3上的有限真实导数极限
迫使每个完整physical field的所有负渠道精确抵消。因此导数也使用常数functional
及归一1/(128pi^5)。每锥zero jets只是固定epsilon展开的系数，不能先把052人工
zero jets当可微函数再微分，也不能让有限常数或近零principal sums代替资格证明。

本轮不判断body32真S3上的实际gradient/Hessian是否存在或发散。严格off-S3 gate
失败时输出unsupported，不能解释为数学上的不可微定论。S3相关正moment障碍
仍保留，原positive-band g-cover路线继续有效。

## 实际新计算

资格检查：axis6点(1,2,3,4,0,0)的全部20个triples给最小距离平方1。含重复列和
dependent triples的8列frame，对两固定点各检查56个triples（其中8个rank2、
48个rank3）；合格点距离平方1/2，真S3点距离平方0并被导数gate拒绝。共132个
准确triple残差，未在body32上执行这一枚举。

新的formal DD6 tuple (0,0,0,0,0,0,1/3)给

    (pi/2)*729 -243/epsilon +9/epsilon^3 -3/(5epsilon^5)+O(epsilon).

新的formal DD7八零tuple准确为−1/(7epsilon^7)，常数0。它只是formal tuple，
不被报成任何不合格点的Green导数。两项均与独立normalized Hermite矩阵及实际
atan Laurent jets吻合。

tiny完整参考固定axis6、x=(1,2,3,4,0,0)。独立路线先对product-Cauchy密度真实
求导，再用SymPy准确时间partial fractions及

    integral_0^infinity (t^2+a^2)^(-m)dt
      =pi*binomial(2m-2,m-1)/(2^(2m-1)*a^(2m-1))

积分。另一条路线完整累加64orthants的DD5/6/7常数与principal channels。两条
路线在以下六字段精确一致：

| 字段 | 完整Gamma常数 | 真实Green字段 |
| --- | --- | --- |
| G | 16/315 | 1/(2520pi^5) |
| D1G | −164/4725 | −41/(151200pi^5) |
| D5G | 0 | 0 |
| D11G | 2414/70875 | 1207/(4536000pi^5) |
| D55G | −8/63 | −1/(1008pi^5) |
| D56G | 0 | 0 |

每字段的64个base tuple均含零和重复投影，全部保留；每个完整
epsilon^-1,-3,-5,-7系数均准确为0。D55直接使用
−2pi^-6 integral[1/((t^2+1)(t^2+4)(t^2+9)(t^2+16))]dt，验证了真正的零坐标
二阶导数及同label factor2。384项是实际tiny完整计算，不是0所有kernel工作，
也不是body32真实数值结果。

## 可执行未来资格与数值意义

`derivative_contract.py` 提供：

1. `qualify_off_S3`：绑定准确frame/point，按实际秩核全部triples，返回严格距离
   和一个固定六basis；dependent triples不跳过。
2. `local_tail_contract`：从上述合格结果和小于delta_point的radius，返回有理
   M上界、六个gradient尾界及6×6 Hessian尾界。可据此把matching H_epsilon
   jet扩展成full-G jet区间。输入若为finite band，仍须另付上端时间尾部。
3. `bind_derivative_request`：核对frame哈希、字段及外部complete-fan/
   multiplicity声明。非空proof reference不等于本函数已经认证外部fan证明。
4. `formal_dd_channels`/`formal_cone_field`：明确标为formal的局部系数工具，
   没有production body32 fan dispatcher；调用者必须先满足完整资格和fan证明。

这是数学合法域与误差预算的新增接口。当前粗常数可能非常大；尚无body32实际
full-G Hessian新数据、全face连续I证书、跨S3的全Green余项或全局g下界。匹配
模型、时间窗口、完整fan及精确符号/相等性仍是不可省略的契约。

## 冻结与复现

固定源的哈希在 `input_manifest.json`，前后均未变化。必要文件由最终
`manifest.json`白名单记录；`checkpoint.json`仅为本地重建性检查点，不进入
精简白名单。原运行命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_fixed.py
```

脚本拒绝覆盖已完成收据/结果；本批已封存停止。复现须用保留冻结来源的干净副本，
不修改本轮原始产物。没有任何追加fixture、新点或production kernel授权。
