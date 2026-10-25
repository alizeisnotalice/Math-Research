# CI4B：四维真实 Green 积分的连续上界证书

2026-10-25。固定 purebody、Lebesgue、centered axis-parallel cube、共同半径、一般 L1/全尺度。本文只认证积分；不搜索 g/tau、不独自更新 U4。完整 Green 全和表示及 g 的独立审计由原任务负责。

目标

\[
 I_4=\int_{[0,1]^3}G(1,y,z,w)^{4/3}\,dy\,dz\,dw\le19/100000
\]

由一个固定的 2440-cell 有理多项式证书实现。实际所得上界表达式约 `.0001841447010739849311214820756`，含真实 corner 帽；准确有理上端及最终 checker 状态见 `ci4b_check_results.json`。独立于 LP 的 checker 实际89.48秒通过，重核4,932,474个非零 majorant 系数全部非负；帽上界表达式为3.5270801566565500745092543062e-9。该表达式小于 `184145/10^9 < 19/100000`。完整表示独审 pending，因此这是依赖锁源表示正确性的积分证书，不能将其单独当作最终弱型常数已审结果。

新批预登记 SHA `a5696336c394920573905c491ad38c789a8d499831cf6865b36cb5ccdca68d53`。先固定10个代表 cell 实际全部通过，再等作者锁源后运行登记的2440 cells，没有增加分辨率、机制修订或回搜参数。旧 CI4 的351区域失败报告仍如实保留，没有重新包装为证书。

## 1. 锁源全和表示与人工极点

只读共享作者接口：

- `u11_minimum_global_interface.py` SHA `6c8d1a5b4912103df8105f7c7662ac2265b0978f91de00651be03acfb1790bae`；
- `u11_minimum_global_results.json` SHA `43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f`。

变量 a>=b>=c>=d>=0，三case为 I: a>=b+c+d，II: b+c>=a+d，III: 其余。接口给 H=pi^3 G=N17/D20，以及同阶齐次对 N20=a^3 N17、D20；在 face a=1 上 H=N20/D20。这里没有自行重做作者的20个真谱顶点/star消去，也没有以点拟合代替它。

旧逐 cone 评价在248个区域的人工极点**不再阻止积分认证**：新方法不对各 cone 或边界 D=0 逐项做除法，而直接核全和多项式不等式。所有2440 cell 均给出有限证书，无未结算 cell。全和 D 仍可在排序域边界有可消0/0；其内部正性和邻域不等式由下一节证明，不声称每个边界分母都非零。

## 2. 单 tetra 的准确仿射上包围与积分

令 tetra 顶点 x_i=(1,b_i,c_i,d_i)，i=0..3，且 barycentric v_i>=0、sum v_i=1。将 N20,D20 作准确仿射齐次替换 x=sum_i v_i x_i，得到次数20的 N(v),D(v)。checker准确验证 D 非零且所有系数非负，因此在 tetra 内部 v_i>0 有 D(v)>0。

每个 cell 保存四个非负有理 P_i，令 P(v)=sum_i P_i v_i。checker不运行任何优化，只准确重算

\[
 R(v)=P(v)D(v)-N(v)\sum_i v_i
\]

并验证其全部次数21系数非负。故在内部

\[
 0\le H=N/D\le P.
\]

其中 H>=0 来自真实保守 Cauchy 半群。幂函数 u^(4/3) 凸，故

\[
 H^{4/3}\le P^{4/3}\le\sum_i v_i P_i^{4/3}.
\]

由于 tetra 中每个 barycentric coordinate 的积分为 volume/4，得到准确连续积分预算

\[
 \int_{\rm tet}G^{4/3}
 \le\frac{|\rm tet|}{4\pi^4}\sum_{i=0}^3P_i^{4/3}.
\]

所有 Pi 和 volume 均是有理数，pi、正 cube root 及求和由 Arb 最终256bit向外包围。D=0的 tetra 边界是有限多面体的零3D测度边界，可以在 Lebesgue 积分中忽略；并没有删除其邻域，也没有将人工 node-zero 平面冒充真实奇线帽。

候选生成器中的浮点 HiGHS LP 只寻找 P_i，以 sum_i P_i 为目标。设某 degree21 系数行为 A_j·P>=B_j，A_j各项来自 D 的非负系数，s_j=sum A_j>0。把候选转成有理数后统一加

\[
 \Delta=\max\{0,\max_j(B_j-A_j\cdot P)/s_j\}
\]

到每个 P_i，就准确支付全部 LP/有理化误差；s_j=0且B_j>0会明确失败。本批没有这种失败。最终checker完全绕过 LP，只核保存的 Pi 和完整 R 系数。因此浮点 solver 成功、采样或点值从未承担证明职责。

## 3. 三物理 tetra 及无叠覆盖

正 face octant 因 b,c,d 置换对称，积分是排序域 T={0<=d<=c<=b<=1} 的6倍。三个物理 tetra 由作者 cone rays 除以各自第一坐标得到：

- I: (0,0,0),(1,0,0),(1/2,1/2,0),(1/3,1/3,1/3)，体积1/36；
- II: (1/2,1/2,0),(1,0,0),(1,1,0),(1,1,1)，体积1/12；
- III: (1/2,1/2,0),(1,0,0),(1,1,1),(1/3,1/3,1/3)，体积1/18。

这些不是由体积相等猜测覆盖：它们分别是三个物理条件 I/II/III 在排序域中的非负 cone 截面，准确非负参数表示就是锁源 U11 表示的三cone map；内部不交、并满 T。体积之和1/6仅作额外核对。

case I 固定做8轮准确最长边二分，得到256 cells。每次二分由一条边中点和两个保留的反顶点构成两 tetra；二者在包含其余两个顶点及该中点的公共平面两侧，内部不交、并满父 tetra。边长用有理平方长比较，ties按登记索引处理。

II/III唯一真实方向线 corner 是 q=(1,1,1)。在各根 tetra 中，以 q 为原点，另外三个顶点为 A,B,C，写 x-q=lambda1(A-q)+lambda2(B-q)+lambda3(C-q)，lambda_i>=0，t=sum lambda_i<=1。第 k 壳为 s/2<=t<=s，s=2^-k，k=0..27。

设 outer 顶点为 sA,sB,sC（此处相对q），inner 为 sA/2,sB/2,sC/2。壳的三个 tetra 准确是

1. (sA,sB,sC,sA/2)；
2. (sB,sC,sA/2,sB/2)；
3. (sC,sA/2,sB/2,sC/2)。

无叠覆盖由三个阈值区域给出，而非只比体积：

\[
 \begin{aligned}
 T_1 &: t+\lambda_1\ge s,\\
 T_2 &: t+\lambda_1\le s\le t+\lambda_1+\lambda_2,\\
 T_3 &: t+\lambda_1+\lambda_2\le s.
 \end{aligned}
\]

在 lambda>=0、s/2<=t<=s 中它们并满整个壳，且内部不交；逐一求四组边界交点即为上述四顶点。每个壳 tetra 再做 max(0,8−2k) 轮最长边二分。II/III合计2184 cells，加I共2440，未超过4096登记上限。

28个壳覆盖 t∈[2^-28,1]。准确覆盖体积为

\[
 |T_{\rm covered}|=\frac16-\frac5{36\,2^{84}}.
\]

每次二分及全域体积都由 Fraction 行列式再核，最后未覆盖体积正好是II/III两个角 tetra 的体积之和。

## 4. 真 corner 帽只付一次

II/III的三个非corner顶点距 q 均不超过 sqrt(2)。故剩余 t<2^-28 满足

\[
 |x-q|\le t\max(|A-q|,|B-q|,|C-q|)
 \le\sqrt2\,2^{-28}<2^{-27}=\epsilon.
\]

两case剩余角域的并包含于同一个真实方向线球帽；不分别多付两个完整球。主 030 已审 body 常数为

\[
 D=\frac1{48\pi^4},\quad B=\frac{180224}{\pi^3},\quad G(1,y,z,w)\le D+B/r,
 \quad r=|(y,z,w)-(1,1,1)|.
\]

来源 `proofs/round011_face_caps.md` SHA `b1564bf4c0cb4a6e2f72786ea5ffae30c922ad45c0cb9874aae6bfa55d4fa8de`，审计 `audit/030_audit.md` SHA `dcbc9f31561668e0b31922990b659786038f205cdb106230921ace53a72db640`。完整三维球径向预算

\[
 E_4(\epsilon)=2^{1/3}4\pi\left\{D^{4/3}\epsilon^3/3+
   (3/5)B^{4/3}\epsilon^{5/3}\right\}.
\]

在 [0,1]^3 内 q是corner，只占球的一个 octant，故 I4 的真实帽预算 E4(epsilon)/8。先把2440个排序cell预算乘6，再加此正octant帽**一次**：

\[
 I_4\le
 \frac6{4\pi^4}\sum_{\rm cells}|T_j|\sum_{i=0}^3P_{ji}^{4/3}
 +E_4(2^{-27})/8.
\]

该式正是checker实际向外包围的表达式。没有把全faces的64倍、排序域6倍、球octant的1/8混用。

## 5. 交付、运行、依赖与认证边界

- `ci4b_certificate.json`：2440个准确 Pi 及确定性几何标识；不是采样数据。
- `ci4b_certificate.py --check`：不调用/不导入 SciPy优化，重新生成确定几何，逐 cell 准确核 D 和 R 全部系数，并用256bit Arb重算完整积分上端。
- `ci4b_affine.py`：共享锁源接口、准确几何、候选构造；仅 `majorant` 被调用时才加载NumPy/SciPy。
- `ci4b_pilot_results.json`：固定10-cell实际先验诊断。

最小复现（在本clone，输入在作者upper clone和主clone只读）：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python \
 cube-finite-dim-bounds/upper/continuous_integral_n4/ci4b_certificate.py \
 --check cube-finite-dim-bounds/upper/continuous_integral_n4/ci4b_certificate.json \
 --save cube-finite-dim-bounds/upper/continuous_integral_n4/ci4b_check_results.json
```

生成器和最终checker共享Python Fraction、python-flint fmpq_mpoly/Arb算术内核及确定性几何代码；“不依赖LP”不冒称不同区间库或完全独立实现。全和表示作者源码、star/null-star身份和更早U11谱锥表示仍由原independent_audit核；本积分证书须与该来源审计一起通过。g>=21/25000、正部迁移、一般L1/全尺度、最终成本整合不由本批代审。

没有占用共享MATLAB、没有自行运行成本优化、没有提交推送、没有修改其他任务源码或总账。成功后交原independent_audit复核；只有主任务有权更新最终上界状态。
