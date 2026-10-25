# 第007轮单候选下界：F独立审计通过

**可登记严格新界：C6 > 103346206909/50000000000 = 2.06692413818。**

唯一候选为a2=(1000001/1000000)(u³+u⁴v³)/4，u、v分别是2、3的正六次根；其余五个间距保持原构造。此结论来自完整静态union、独立整数内缩及具体有限L1见证，不使用正导数或预测小数。它严格超过第007轮候选前冻结的完整准确旧上端和本阶段phase frozen完整上端。

上界仍64，双侧达标条件仍未满足，后三轮计数不开始。只属于同一构造的可靠新实例，不宣称新下界机制或全局最优。

## 1. 实际做过的独立验证

1. 阅读E的certify_candidate.py及其复用的36系数有理算术，核对固定参数、4096长度向量、radius=(Πell)^(1/6)/2、半胞元端点、包含剔除和静态slab union。E没有读取D的导数或dualunion系数，只对候选来源文件作hash；其普通数域重建与D的方向计算独立。
2. 阅读主inward_smallstep.py。复制到本目录inward_replay.py，只把ROOT替成主输入绝对路径，算法逐字不变；实际执行退出0，用时约0.324秒。127个非空内缩盒、96个保留盒及完整有理union下端与主一致。此回放不导入E的数域程序。
3. F新增纯Fraction程序独立重算有限源比例、全部完整准确比较，并用自行二分得到的10^90分母六次根区间重新评价E的36项最终多项式。所得新区间完全包含于E已保存区间；主/F内缩R严格低于E区间下端，符合几何内缩方向。

没有重复D扫描、没有重新执行E全部数域union、没有MATLAB或新kernel。E数域构造本身的执行证据由E保存；F的实际回放对象是独立整数内缩和转移比较。

## 2. 内缩几何、正裕量与并集的方向

主/F以10^50为根分母，逐个用整数六次幂包住2、3及m=Πell的正根。每一维间距[a_low,a_high]均由正系数单项式的单调性取得。奇数ell取中心块[-floor(ell/2),…,floor(ell/2)]，偶数取[1-ell/2,…,ell/2]。

每个候选接收盒用半径r_m=root_low(m)/2和最不利的间距端点收紧，再以10^35网格严格向内多移一格。代码对每个坐标的最左/最右原子位置逐项验证

|x_j-a_j k_j| ≤ r_m−margin，margin>0，且(2r_m)^6≤m。

凸区间性保证中间原子同时包含；六维取Cartesian积得到m个单位质量。内缩盒严格位于(0,a_j/2)各轴中，不存在半胞元边界的双重体积。整数union按所有端点切slab并递归合并余维，计算准确Lebesgue体积；包含剔除只删除已由另一盒覆盖的盒。

R_inward=64|U_inward|的完整分数见inward_smallstep_results.json。乘64来自六轴反射到一个完整胞元，而不是密度或体积中漏掉/重复间距乘积。每完整胞元只有一个单位原子，归一化比例就是胞元中的接收体积R。

## 3. 一个具体有限L1见证

令K=2^45=35184372088832。接收胞元的六个整数索引各在[-K,K]；源原子的索引各在[-K−2,K+2]。长度ell≤4的见证块及其反射只用[-2,2]索引，故collar恰取2足够。F独立枚举这四种一维长度与反射，确认最大绝对索引2。

在每个源原子位置放总质量1的均匀小方体密度

b_eta(x)=(2eta)^(-6)·1_[-eta,eta]^6(x)，eta=5×10^-36=1/(2×10^35)。

eta等于实际几何margin的一半。因此每个选定接收点的中心轴平行六维cube、以同一个r_m为六轴共同半径，仍完整包含其m个质量1 bump；平均值≥m/(2r_m)^6≥1。

取theta=(R_inward+target)/(2R_inward)，target=103346206909/50000000000。F纯有理核0<theta<1。故接收集合全部在{Mf>theta}中；源函数非负、有界、紧支撑，属于原问题允许的一般L1类。其质量为(2K+5)^6，接收体积为(2K+1)^6R_inward，于是

theta |{Mf>theta}|/||f||1 ≥ theta R_inward[(2K+1)/(2K+5)]^6。

实际准确分数的有理小数包含为

**[2.066924138180051561, 2.066924138180051562]。**

与target的严格正余量至少5.1561505693×10^-14。此为有限、可定义的L1函数见证，不只是无限原子测度的形式比例。

“全部正尺度”在下界侧无需证明每个尺度都表现坏：原极大算子对全部r>0取上确界，当然包含本证明所用的有限见证半径。六坐标始终共用单一r_m，没有改成强极大算子的独立坐标半径。

## 4. 完整理想数域值的转移与准确比较

E的理想半胞元union由有限个非空盒组成。取其开内部，每一点对所用见证原子有严格正包含裕量；逐步内缩盒并令bump半宽趋零，体积趋于理想union体积。再令theta↑1、K→∞，上述有限L1构造给C6≥R_ideal。因此E完整代数值可以作为该实例的完整准确记录；本轮的具体有限见证已足以推出所展示的严格目标，不需要依赖这个极限来弥补有限源失败。

F比较采用完整Fraction端点，实际逐一通过：

- 第007轮frozen_lower_before_candidate.json的旧当前完整上端；
- phase_20261025/frozen_thresholds.json的lower与internal_current_best完整上端；
- 候选前旧ledger的完整当前上端、R6-001/k1/k16历史上端；
- 继承R4完整上端、动态五维完整上端、已核公开下界区间上端。

新E区间下端、R_inward、具体finite_ratio及展示target均严格超过适用旧当前/phase上端。没有拿展示小数替代完整旧区间。

旧ledger比较输入已从SHA256=096aaff33ea186359046528cfccdde86367610f7fa081d0bcf7d9b273f69ad6b提取并冻结在frozen_comparands.json；主后续更新ledger不会令回放误拿新候选与自身比较。完整旧当前上端与phase冻结上端准确相等也已核。

## 5. 证据hash与重放

关键SHA256：

- 主整数内缩源码：570336112117aef31aa4a76192a22da10c227ad9e144b779549714c79f5e580d。
- F仅ROOT调整的回放源码：273d27fd77c332af6ecec3f1be413c324f5a30cebf97de4984e51328ee965f6f。
- F实际内缩结果：8031301028affb9e5111672d196f624c5bd77bb218c279f6639ce470f4368ac9。
- E完整候选结果：b5661c0eb23c264286085053b3c59ab9f992918dcb1dc4a79fdb01e1ddcb6b81。
- 第007轮候选前冻结文件：9af38aee2138b1c64c28dba7b8b61e827f33ad5191a10ddb0ebeb57af4c1091c。
- phase冻结文件：e0528774b908a5135e24337dacfd6c15c5cc3e156da4578cd694fdece5eecd6d。

E源码及其36项算术源码hash也由F实际核对，全部输入和本轮checker hash见comparison_transfer_results.json。回放命令：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-6-bounds/audit/round007/lower_candidate/inward_replay.py
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python cube-dim-6-bounds/audit/round007/lower_candidate/check_comparison_transfer.py
```

状态为PASS_INDEPENDENT_REPLAY_FULL_COMPARISON_AND_FINITE_L1。主可据此更新本阶段lower突破标记；最终账本、归档和推送由主负责。上界尚未突破，不能设置dual checkpoint或开始后三轮计数。
