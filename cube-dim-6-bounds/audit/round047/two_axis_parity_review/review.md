# 第047轮 E：两轴 parity 质量的一阶排除独立复核

**结论：通过，严格限于冻结 E008 spacing 上的 127-label 新子族。**

\[
R(\varepsilon)=R(0)+c|\varepsilon|+O(\varepsilon^2),\qquad c\le0.
\]

不能推成 c=0、严格负 c、二阶项符号、有限 ε 排除、整个质量族最优或 Cn 新界。两侧门槛不变。证明使用新的两轴质量与合并块，而不是 evenness、concavity 或 round039 的 25-label 静止删除结论。

## 冻结与实际预算

先登记 `registration.json`，实际 source manifest SHA-256 与指定值一致：

`e2495c43312ebc2c550d695f593c63cf5e3db8fb15bf4c7727bed0b621a2257f`。

其全部声明文件 hash 均通过。独立最小核没有导入或执行作者 writer，仅使用标准库整数、Fraction、集合、JSON 与 hash。

实际 PID **50382**，执行 chunk **eaaa28**，exit **0**，数学 CPU **0.007819 秒**。恰好 12 个显式关系/库存谓词；其中类范围谓词在继承准确 ranks 上遍历全 127-label metadata，不是 31 次新空间曝光计算。全部 1524 个面签名、460 个异速 pair 的 rank bookkeeping 数量均透明记录。没有新 union、ε 候选、baseline replay、root 近似、新数域符号判定或 MATLAB。登记至报告及 manifest 收尾约 5.5 分钟，上限 12 分钟；CPU 目标 2 秒以内。完成后停止。

## 1. 新模型的准确质量、四 phase 与共同半径

源位置固定为 a·k，质量为 `1+ε(-1)^(k₁+k₂)`，|ε|<1。period 为 (2a₁,2a₂,a₃,…,a₆)。四个源 phase 的常数质量系数均为 1，ε 系数为 (1,−1,−1,1)，因此周期质量恒为 4；本核只算系数，没有代入非零 ε。

中央块 J_ℓ 的起点仍为奇长 −⌊ℓ/2⌋、偶长 1−ℓ/2，定义 S(ℓ)=∑_{k∈J_ℓ}(−1)^k。直接有限系数求和给

\[
M_\chi(\ell,\varepsilon)
=\prod_h\ell_h\left[1+\varepsilon\chi\rho_\ell\right],\quad
\rho_\ell=\frac{S(\ell_1)S(\ell_2)}{\ell_1\ell_2}.
\]

任何一个前两轴长度为偶数时 ρ=0；其余 (1,1)、(1,3)/(3,1)、(3,3) 分别为 1、−1/3、1/9。本核对全部 127 标签与四 phase 的 ε 系数作了直接整数求和核对。

选定块仍只有一个六轴共用半径：

\[
r_\chi(\varepsilon)=\tfrac12 M_\chi(\varepsilon)^{1/6},\qquad
r'_\chi(0)=\chi\rho_\ell r_0/6.
\]

非 cell 下/上端速度分别 −r′、+r′。ρ≤1 且基线 m>0，使半径在零附近 analytic。局部反射 h_j↦−h_j 保持 parity；每个接收 phase 有 64 个半胞元反射像。四 phase 的 character 两正两负，故

\[
R(\varepsilon)=64[2V(\varepsilon)+2V(-\varepsilon)]/4
=32[V(\varepsilon)+V(-\varepsilon)].
\]

这里 evenness 只确定两侧对称，不说明导数为零。

## 2. 完整 31 组库存确实是新质量方向

用冻结 E008 的全部准确 endpoint ranks 与唯一右边界 tie 合并重建基线位置。127 个盒全部坐标仍严格正宽；因此下面的异向“基线已填满邻域”判断没有忽略零宽新生盒。

本核从新的 ρ 重新生成全部 1524 个 face occurrence 的速度签名，与作者 31 组输出逐标签、axis、rank、orientation、mass、ρ、速度完全一致。结果为：

- axis1 upper：质量 2、4、8、16、2304，共 5 组；
- axes3..6 lower：质量 4、8、16、2304，各四轴，共 16 组；
- axes2..6 upper：质量 108、162，各五轴，共 10 组。

共 31 组，确含 round039 原十个位置组之外的新 21 组。旧库存只用于核对“新增 21”这一集合差，不输入旧质量方向的速度或支配结论。所有 31 组都具有同一组内共同基线质量与共同法向 source length。

使用继承 ranks 对不同速度 pair 的五个横向区间作准确 max-lower/min-upper 比较，共得到 **460 对正五维交**；此批次没有零面积 pair。这里没有求任何五维或六维 union。所有 460 个 pair 的 C=max(A,B) **实际在冻结 127-label 集内**，共有 **32 个不同 C**，并逐 pair 核对 mC>m。

因此归档结论不是“C 的长度≤4 所以算在族里”的跳步：真实成员关系已单独核对。也没有假定 C 静止；其中 **48 个 pair 的 C 具有非零新 ρ**，确实会动。

## 3. 整个共同五维面为何有严格 collar

精确中央源链为

\[
J_1=\{0\}\subset J_2=\{0,1\}\subset
J_3=\{-1,0,1\}\subset J_4=\{-1,0,1,2\}.
\]

对同一异速组中不同标签 A、B，基线半径均为 r，因为质量相同。令 F 为其闭共同 facet 的横向交（法向为共同非 cell 面坐标）。任何 x∈F 的中心 cube 在半径 r 时都包含 A 与 B 的全部源原子。轴平行 cube 是坐标区间乘积，所以包含二者的每轴源端点，进而包含中央块的坐标 Cartesian hull，即 C=max(A,B) 的全部源原子。这里并没有把“两块原子的集合并”误当作 Cartesian 积；使用的是包含它们的 cube 的坐标凸性。

不同且等乘积的两个正整数向量不可能分量可比，故 mC>m，rC>r。更强地，对每个 C 源原子 y 与全部 x∈F，

\[
\|x-y\|_\infty\le r=r_C-\delta,\qquad \delta:=r_C-r>0.
\]

所以存在覆盖**整个闭 F**的统一相对 H collar：

\[
\bigl(F+(-\delta/2,\delta/2)^6\bigr)\cap H\subset B_C(0).
\]

若 F 碰到横向 cell 边界，不能宣称跨出 H 的绝对 R⁶ 开 collar；上式的 relative-H collar 正是需要的完整表述。法向面是非 cell，故在该方向也有真实两侧余量。

本核还独立用准确 ranks 直接核对了全部 460 对的几何结果：C 在法向的 lower/upper 严格夹住共同面；每一个**非 cell 横向端点**都被 C 严格超越；遇到横向 cell 端点时 C 取相同 cell 边界，符合 relative-H collar。因此没有只核法向或漏掉其他轴非 cell 端点。

## 4. C 会动时为何仍遮住一阶 coexposure

C 可能具有 ρC≠0，但 rC(ε) 在零附近连续且 rC(0)−r>0。A、B 的面与投影端点位移为 O(|ε|)，C 半径位移也是 O(|ε|)。因此选足够小的 |ε|，这些位移总量小于例如 δ/4，前节 collar 仍覆盖两个移动面的共同横向邻域。

每对的 δ 可以不同，但只有固定有限 460 对，可取统一正小邻域；无需给出数值 ε，也不能据此声称较大/全部 |ε|<1 的覆盖。多重交会若含两个不同速度参与者，就含上述某个被遮住的 pair；因此不会出现遗漏的正五维多面 coexposure。

这处理了 48 个 moving-C pair；不能用 round039 的静止 C 论证代替。

## 5. codimension、异向、固定 cell 与余项

正一阶 cusp 的局部机制是同向外包络的速度范围 `max(v)−min(v)>0`，乘正五维曝光面积。前节已遮住全部这类正面积异速 pair。剩余同向外片只有一个一阶速度，ε 与 −ε 的线性项相消；同速但高阶不同的 germs 仅影响 O(ε²)。

若同一横向内部片上同时有正宽盒从法向两侧结束，基线已覆盖双侧邻域，小扰动只能保留或开 gap，因此 even 一阶贡献非正。若为固定 H 的 cell 面，内侧基线覆盖、外侧始终被 clip，仍无外侧增益而可能有内缩损失，贡献非正。

高重交的 O(ε²) 处理不需要一般位置：取**所有轴全部基线端点层级**，包括非 cell 与 cell 面。先合并同轴同坐标超平面；这些高重合面仍属 codimension1，不能删除，已由全部 31 速度组处理。只有同时位于两个**不同坐标方向**超平面的 O(|ε|) 带中，才划入忽略集合。它是固定有限个两条坐标薄带的交，体积 O(ε²)。任意多盒在同一 edge/corner 相交只增加常数，不改变阶数；正交坐标方向避免了退化夹角。没有漏掉横向非 cell 端点。

外于双带的局部占据是一维法向问题。所有源位置固定、半径 analytic、clip 为有限 min/max，端点有统一 O(ε²) Taylor 余项；有限盒并集对端点局部 Lipschitz，使该误差仍为 O(ε²)。每侧 analytic 端点次序也在充分小穿孔区间稳定，故完整两侧一阶展开存在。

综上，真正的局部面分析给

\[
c=32[V'_+(0)-V'_-(0)]\le0,
\]

并非从周期总质量、evenness 或六次根 concavity 直接得出。

## 验收范围与停止

结论只排除该冻结 127-label 两轴 parity 新子族的**正一阶 cusp**。准确 c、严格负性、二阶项与任何非零 ε 改进均未计算；未排除原族之外的新标签、其他质量族或全局机制。无 Cn 新界，两侧门槛不变。

可审材料为本目录 `registration.json`、`minimal_review.py`、`review_results.json`、`review.md`、`manifest.json`。未暂存、提交、推送、扩招或改 ledger；已停止，上界工作继续优先。
