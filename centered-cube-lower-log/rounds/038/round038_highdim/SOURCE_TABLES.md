# Common-D source table module

模块 `source_tables.py` 独占提供一维来源/查询核；driver/proposal 由父线负责。本模块不建立 R×J 缓存，不抽高维样本。当前冻结 SHA256 为 `85c592ac80d54b376a43038af7395cfc757482b4851b996ffe4d72be2a2fa12d`。

## 来源与 API

`Source(R=4,w=4,v=8)` 固定 K255、n6704800、p=.5、M3、δ=1/(4Kv)、Q=6wKv=48960，query bins 为 `[q/Q,(q+1)/Q]`。属性包括 `n,K,Q,delta,pchi,bins_center,binwidth`，另保留 `atoms,labels,c,field` 供独审。只允许该组参数；平均接口允许 L=1 或 `1+δ+h<L<2`，其中 h=1/Q。

field_labels 函数逐字复制第34轮 `block_probe.py`，验证脚本通过 AST 取原函数再逐字核对。表 b 在 j<127 的 plus cells 使用 `(label_j+b) mod4`，其余 cells 使用 `(label_j−b) mod4`，U=(label+.5)/4。每格峰质量 p/K，全部 table 使用同一个原 D_M、c_m，box 都留格。

`evaluate(L,b)` 每次只返回当前 L,b 的：

- `mean_plus`, `mean_minus`：各 Q 数组，真实 class 核的 bin 均值近似。
- `relative_variance_class_max`：各 Q 数组，两类相对方差帽的最大值，均值分母使用计算均值减去下面的 error。
- `analytic_mean_error_abs`：对每一个 class 均值的解析中点积分误差帽。
- `floating_point_guard_abs=1e−13`：明确标记的经验浮点 guard，不是向外舍入证书。
- `mean_error_abs`：上述两项之和。driver 应同时保留它们的不同证据资格。
- `normalizer`：scalar `mean(mean_plus+mean_minus)`，供实际 tilt 概率归一化/C_Z 审核。
- `range_plus/minus`, `transition_count_plus/minus`：各 Q 诊断数组。

L=1 使用 class 质量的恒等式，range/variance严格为0；pchi 的浮点表示仍有 rounding guard。`point_class_kernel` 是独立于 bin-CDF 的真实 box 点核，仅供验证少量查询。

## 均值计算与解析误差

对1<L<2，设 ℓ=L−1，query x 的 extra window 长ℓ、中心 y=x+.5。class 核为一个完整周期的 class 质量，加 extra window 的 class 积分，再除L。共同 plus 质量

`pchi=a+(p/(2K))Σ_m c_m sin(2πma)`，`a=127/255`，

当前 pchi=.498047800888639，minus质量1−pchi。完整周期没有另乘π，也未将class质量设为.5。

peak bin 均值采用 source uniform δ 与 query uniform h 差分布的 trapezoid CDF。对于每个 extra-window 端点，仅把距离≤(δ+h)/2的 atom 作为 partial candidate，远离端点的完整 atoms 用 `searchsorted` 整数计数。每端最多一个 candidate，因为相邻atom间距≥1/(4K)>δ+h。partial区间不重复计入full内部端点；外部零贡献端点可保留。此算法避免四个巨大平方 primitive 的差除以h。

背景的常数部分在 query bin 内对 class/window 交界点分片；overlap为分片线性，中点积分在实数算术下精确。对 D 部分也按同一交界点分片，点处 D 的区间积分用稳定恒等式

`(c_m/2)[sin(2πm u)−sin(2πm l)] = c_m cos(πm(u+l)) sin(πm(u−l))`。

在一片内 extra-D 积分作为 y 的函数 F，有 `|F''|≤2 Dprime_bound`。各片中点积分总误差除bin宽后≤`Dprime_bound h²/12`；因此 class均值的解析绝对误差是

`Eanalytic = (p/(KL)) Dprime_bound h²/12`。

本系数公式给 `c_m²≤(23328/87)/(π⁴m⁴)`，因 sqrt(23328/87)<20、π>3，可以取 `πm c_m≤20/(3m)`。由此采用确定系数帽

`Dsup_bound=(20/3) H_M`，`Dprime_bound=(880/21)M`，

后一式用了π<22/7。这些不是查询网格上寻得的最大值。它们同时给背景下界 `1−p−(p/K)Dsup_bound=.4760348583…>0`。

L7/6时 Eanalytic=7.3451996032e−12；不包含普通浮点四则/三角函数舍入。单独列出的1e−13 guard未获interval认证，本模块的数学解析误差与实现roundoff资格不能合并宣称为严格计算证书。

## Range 与相对方差帽

记 N为该class在两个window端点的partial候选数，0≤N≤2。query在宽h bin中变化时，class常数背景overlap的振幅≤h；每个候选峰的质量振幅≤min(1,h/δ)；D-mask 积分的导数绝对值≤2Dsup_bound，故

`range ≤ [(1−p)h +(p/K)N min(1,h/δ)+(p/K)2Dsup_bound h]/L`。

mask跳变不会破坏该一阶Lipschitz帽：窗口积分的导数几乎处处是两端density-mask值之差；class边界处积分仍连续。没有把不同表、不同端点或不同核视为独立。

任何取值落在区间振幅R内的变量方差≤R²/4。用真实均值下界 `mean−error>0` 除，得到每class relative variance帽，再对两个class取max。解析数学上这一步成立；FP实现对error/coefficients的证据边界同上。R4各表内部得到自己的帽，driver可按分层Minkowski组合，勿先把不同表transition计数合并。

## 已运行的一维验证与纠正记录

`validate_source_tables.py` 在L7/6对每表冻结52个bins，覆盖预定6个cell的峰进/退出、class边界和uniform控制，分别用真实点核按breakpoints分片做24节点Gauss积分。它另对64个class/bin以Fraction积分完整piecewise-linear peak；33个有理count均值不是整数，确实含partial transition。

主验证为 `SOURCE_VALIDATION_normalized_reference.json`：4表 normalizer均为1.0，最高class relative-variance帽≤2.51131e−7；mean对reference最大差1.416e−12，小于解析帽；观测振幅/帽最大.997643，参考variance≤range²/4。Fraction count的最大误差折算class kernel为4.19e−15。类均值在全Q平均后，plus与pchi差≤5.56e−17，minus与1−pchi差≤1.12e−16。Gauss/振幅观测只是数值控制，Fraction认证仅覆盖有理peak部分。

首次确定验证耗时.9554秒，原 `SOURCE_VALIDATION.json` 和 `validate_source_tables_initial.py` 保留。发现reference权重直接以浮点两端差除精确h时，常数积分的权重和会有约1e−12的消减偏差；这与source算法无关。验证脚本改为归一Gauss权重，使常数积分为1后，按同一批确定bins重算，耗时.9380秒，输出新文件并记录原结果hash。source_tables源码和hash未改变，两次均无随机数或高维样本。

环境为 `/opt/homebrew/opt/python@3.15/bin/python3.15`、NumPy2.5.4。四表的单L evaluate实测约.05秒。重放验证须在新目录保留source模块；原输出存在性保护会拒绝覆盖。driver的高维PB、actual prior、proposal/RN、C_Z与统计区间均不属于本模块验证。
