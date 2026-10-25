# 四维 pure-body minimum：增量反例与薄弱接口审查

2026-10-25。**本次有限审查未发现反例或阻断性缺口；下列准确检查全部通过。** 已完整阅读作者锁定证明 u11_minimum_proof.md，SHA256 `5546ee2fd0c8dddab3b88d0db79143a846a81711f13c6967b4fa922fa9c4757e`。审查对象固定四维、八条未单位化 body 方向、真实 Green、目标 (21/25000)。不扩展 τ，不计算 I，不改变账本或最终弱型常数，不替代原 independent_audit。

原指定作者 checker 的 dd71… 版本在读取期间按作者预告删除纯 timing metadata；本脚本首次锁断言立即停止。已明确换锁到 `35489f9fa4c52a09154105ad949a8d3410f61bf226e8a24dd8d57359f0a3e060`。最终数学结果 hashes 在 n4_review_result.json 中，作者已确认字节冻结。没有执行作者会重写成果的 producer/checker，没有修改作者文件。

## 1. 准确代数恒等式：已核

本审查用逐项有理函数约分相加验证固定 star，未调用作者 common-denominator 构建脚本或重跑表示搜索。六类 orbit representatives A/B/C/D/E/F 分别涉及 48/48/12/48/12/12 项：

- A/C/E/F 的有理函数分子准确为零。
- B/D 与锁 star 字典交叉相乘后准确相等。
- 完整 proof 中新给出的紧凑 (S_B,S_D) 公式也与重构 star 交叉相乘准确相等。
- 三个既定 global (N_{17}/D_{20}) 与二十个变换 star 的全和准确恒等，不由点值拟合或有限采样支持。
- 作者三个 (P=[a^3N-kD](Ru)) 的全部稀疏系数与本审查重新展开的系数逐项相同，非零项数为 1739、1702、1702，全部非负。

这里 (k=(21/25000)(22/7)^3=3993/153125)。原始 1152 锥到真实 Fourier Green 的正常数与符号是共享输入，来源 u11_representation.md；本审查读了推导，但有限代数检查不是另一份从独立谱模型开始的重建。

## 2. 样本决定的二十项符号：整锥证明补足

作者 global producer 确实使用每 case 一个 sample 决定符号。单靠代码中的 sample 不足以证明全域，因此本审查另做整锥检查：对每条 ray (r) 和该 case 的四个生成 ray (R_j)，准确计算四个 (r\cdot R_j)。二十条 ray 的每组四系数均同号（允许零），且至少一个非零。因此对所有 (u_j>0)，

\[
\operatorname{sign}(r\cdot Ru)
\]

恒定，并与 producer 的 sample sign 一致。60 组系数完整保存在结果 JSON。

还可直接看出：B 的十二个 (e_i\pm e_j\) 在严格排序正 chamber 全为正。D 的八个 ((1,\pm1,\pm1,\pm1)) 中，仅

\[
a-b-c-d,qquad a-b-c+d
\]

可能变号，其余六个严格为正。对应 signs 为：

|case|条件|第一项|第二项|其余六项|
|---|---|---:|---:|---:|
|I|a−b−c>d|+|+|+|
|II|a−b−c<−d|−|−|+|
|III|−d<a−b−c<d|−|+|+|

所以 sample 的使用在有上述解析/准确检查后合法，不存在遗漏的内部 sign subchamber。边界零投影不在此赋 sign=0，后文用真实 G 延拓。

## 3. 标准 star 到全部 orbit 的等变性：已核

准确检查 1152 个带绝对 determinant 的锥 multiset 在三个相邻坐标交换及一个坐标变号下不变，同时 rates 不变。这四个操作生成整个 signed-permutation 群。对任意该群矩阵 T，

\[
S_{Tr}(x)=S_r(T^Tx)
\]

由锥一一对应、Jacobian/rates 不变及内积恒等式逐项成立。因而 B 的 24 个和 D 的 16 个 star 均由标准 B/D 公式覆盖，A/C/E/F 的零恒等式也覆盖全部对应 orbit。

对 T=−I，三个分母差同时变号，故 (S_{-r}(x)=-S_r(x))。与 (operatorname{sign}(x\cdot(-r))=-operatorname{sign}(x\cdot r)) 配合，每对 ±r 的乘积贡献相同；因此全和系数由 −1/32 变成二十条 line reps 的 −1/16，未漏乘二或改变符号。

作者的 B 置换替换恰为 (T^Tx=(x_i,\sigma x_j,x_k,x_l))；D 为 (T^Tx=\operatorname{diag}(r)x)。本审查按这个映射独立构造二十个式子并验证 global 恒等式。

## 4. 非负 cone 覆盖与 D 的内部正性：已核

取 (a\ge b\ge c\ge d\ge0)。三个 cone 的准确逆参数为

\[
\begin{array}{ll}
I:&(a-b-c-d, b-c, c-d, d),\\
II:&(a-b, b-c, b+c-a-d, d),\\
III:&(c-d, b-c, (-a+b+c+d)/2, (a-b-c+d)/2).
\end{array}
\]

对应列依次为 ABEF、EBCD、EBDF，determinant 绝对值为 1、1、2。本批用准确矩阵逆核对。将 (a-b-c) 与 ±d 比较，恰得三个闭 cone 覆盖整个排序 chamber，内部不交。不是从体积相同推断覆盖。

对各 case 的 D 每个线性因子，本批代入 (Ru)、按因子符号定向，确认系数全部非负且至少一个正；整体常数严格正。因此 (D(Ru)>0) 对 (u_j>0) 成立。结合准确 P 非负可合法除以 D：

\[
G(x)=N/(\pi^3D)\ge k/(\pi^3a^3).
\]

此式先在 generic 严格内部建立。cube 内 (a\le1)，且 Arb128 已核 (pi<22/7)，给出所需目标。没有把 (D=0) 处直接代成零或除以零。

## 5. D=0、原点及方向线：解析延拓足够，但不能只用 lsc

仅有 G 的下半连续性并不足以从稠密集合的下界推出边界下界；下半连续函数允许向下的孤立值。作者 proof §4 使用的是**真连续性在线外、真无穷在线上**，因此逻辑方向正确。本审查给出可独立审核的具体概率证明，避免这一接口只停留在文字引用。

令 (S) 为八条 body 直线，(x_0\notin S)，(d_0=\operatorname{dist}(x_0,S)>0)。在 (B(x_0,d_0/2)) 上距离至少 (d_0/2)。取 (delta=d_0/64)，写过程为八个独立尺度 t 的 Cauchy 标量沿 body 向量相加。若至多一个标量绝对值大于 δ，其余小向量和范数至多 (16\delta<d_0/2)，故这种事件的密度不能贡献该邻域。

于是只需取至少两个 large 标量的 28 对并集。任意一对 body 向量均可补成四维 basis；准确检查全部 28 对有非零补基 determinant，且最小非零值为 8。对 large 两个标量，Cauchy density ≤(t/(\pi\delta^2))；另两个 basis 标量 density ≤(1/(\pi t))。条件化其余四个标量、乘 Jacobian 后得到局部统一界

\[
p_t(x)\le \frac{28}{8\pi^4\delta^4}\qquad(t>0, x\in B(x_0,d_0/2)).
\]

正时间 p_t 连续来自绝对可积 Fourier 逆；在 (t\ge1)，再用 (p_t(x)\le p_t(0)=11/(3072\pi^4t^4))。前后两个可积控制给 (G=\int_0^\infty p_tdt) 在线外有限且真连续。这个论证在每个线外点给邻域统一控制，因而确实足够处理所有 (D=0) 的非线边界及任意 zero/repeated 投影点。

在线上 (x_0=s v)、(s\ne0)，选含 v 的 basis，把非 basis 四个 Cauchy 标量限制在 (O(t))（概率有固定正下界）。剩余 basis 密度中，沿 v 的标量在 (s+O(t)) 给 (gtrsim_s t)，另外三个在 (O(t)) 给 (gtrsim t^{-3})，所以 (p_t(x_0)\gtrsim_s t^{-2})，积分发散为 +∞。原点直接用准确 (p_t(0)\asymp t^{-4})。故方向线及原点无需依赖从外侧极限猜值。

generic 严格 cone 点避开有限投影/重合超平面，组成稠密集。任何闭 cone 边界点可从其内部 generic 点逼近；在线外用上面真连续性，在线上用直接 +∞，即可延拓整 cube 下界。这也解释了为何不要求每个 rational N/D 都在边界拥有显式约分后的可代入式。

## 6. 有限边界测试：通过，但只作公式核验

按预登记，每个 cone 的 15 个非空 ray-support 选一次准确点，共 45 次、31 个不同点。每个 D=0 的模式在对应 support stratum 上由非负线性因子的零系数模式确定。沿 (x(t)=x_0+t\sum_jR_j\) 计算 N/D 的准确首项阶数和系数，与独立全节点 residue evaluator 的准确有理数比较。

- 45 次全部一致。
- 14 次 D=0，其中 2 次是同一 body line 在两个 cone 的表示，准确内向极限为 +∞。
- 其余 D=0 点均有有限值，未被误判为奇线；例如 (C(e_1)=3/32)、(C(1,1,0,0)=1/24)、(C(1,1,1,0)=17/640)，其中 (C=\pi^3G)。
- 输入边界点使用原 cone ray 和，未归一化，某些点在 cube 外；这些测试只核等式/极限，未把它们的 raw 数值拿来比较 cube 下界。

所有层的下界来自第 5 节的解析论证，绝不来自这 45 个点。测试完整记录见 n4_review_result.json。

## 7. 复核与剩余审核责任

实际最终有限批耗时约 **3.588 秒**。最小复核命令，在本目录执行：

```text
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B n4_review_check.py
```

需要现有 SymPy、python-flint、NumPy、mpmath。脚本只读作者固定 geometry/结果，不运行其会写文件的生产程序；运行输出仅写本目录 n4_review_result.json。一次 JSON fmpz 序列化实现错误在数学断言全部通过后发生，已改成 int 并重跑同一集合，没有扩大测试范围。

审查边界：共享 1152 锥真实 Green 接口来自作者锁源；本增量审查专门补足 fixed sign、star 群等变性、global 恒等式、D 内正和边界合法性。原 independent_audit 仍负责最终验收。本结论不提供 I 上界、不单独宣称 (U_4<16)，不重跑已经独审通过的三维证书。
