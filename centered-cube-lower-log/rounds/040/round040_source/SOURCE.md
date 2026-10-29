# 第40轮 common-D count-gated source模块

`source_tables.py` 为第38轮冻结模块的独立副本扩展，原目录未改。冻结SHA256：

`3d49b03dab553e1d828638248a724ee5da9d5201b9eba57a7a8643fa76f3c616`。

支持 `Source(R=4/8/16/32,w=4,v=8)`，K255、B4、M3、p=.5、n6704800、δ=1/(4Kv)、Q=6wKv=48960不变。`evaluate(L,b)` 的字段、均值误差口径、分片积分和range算法与38一致，无R/J大缓存。table b和所有坐标共用同一个L。

## 相位、几何和R4兼容

GF标签仍逐字沿用round034的field_labels(8,0x11d)，原标签y_j∈{0,1,2,3}。plus为j<127。新表为

`U_plus=(y_j/4+(b+.5)/R) mod1`，

`U_minus=(y_j/4+(.5−b)/R) mod1`，

`a_bj=(j+U)/K`。

实现不用接近整数的浮点mod，而保存奇整数

`N_bj=(y_j R/2+1±2b) mod(2R)`，`U=N_bj/(2R)`。

R4时它逐点等于旧 `(y_j±b mod4+.5)/4`；本机已对原子、标签以及全部evaluate输出逐位核对，下文给范围。

每格内最小U与1−U均为1/(2R)。box的cell-unit半宽为1/(2wv)=1/64，所以R≤32都留格，R32允许恰好触格边。构造时用整数条件 `N_bj wv≥R`、`(2R−N_bj)wv≥R` 验证，而不是浮点严格大于。接触边界只有零Lebesgue质量，不改变cell/class质量。实际实现的普通浮点roundoff仍未interval认证，整数证书说的是定义中的精确位置。

在一个class内部，同表相邻cell的U之差是四分格差，中心间距≥1/(4K)。跨plus/minus边界不承诺这个更强间距，但算法各class独立搜索，足够保证宽δ+h的endpoint候选区间最多有一个atom。周期class内部间隙也已精确检查。

## 真实质量、共同η和误差资格

每表每格有相同峰质量p/K，box留格，背景始终是共同

`ρ_bg=1−p+(p/K)Σ_m πm c_m cos(2πmx)`。

给定系数公式可粗证πm c_m≤20/(3m)，故背景下界.4760348583…>0。表中每个完整一维密度非负、积分1。各表cell mass共同：

`q_j=1/K+(p/(2K))Σ c_m[sin(2πm(j+1)/K)−sin(2πmj/K)]`。

plus质量η在所有R、所有b下相同，约.498047800888639，不能替换为.5。count-gated来源只要使用源格决定的固定count分割，其联合cell质量是∏q_j，总质量1；实际prior仍是Bin(n,η)段质量。相位扩展没有改峰质量、n、v、D或c_m，也没有取得原匹配下界参数资格。

API仍返回 `mean_plus/minus`、`relative_variance_class_max`、`mean_error_abs`、`normalizer`，并返回分列的 `analytic_mean_error_abs`、`floating_point_guard_abs` 及range/transition诊断。完整period加extra弧的分解、trapezoid峰均值和class边界分片沿用38算法。中点D的解析均值帽是

`(p/(KL))Dprime_bound h²/12`，其中h=1/Q、Dprime_bound=(880/21)M。

振幅帽沿用

`[(1−p)h+(p/K)N_transition min(1,h/δ)+(p/K)2Dsup_bound h]/L`，

每class N_transition≤2、Dsup_bound=(20/3)H_M。相对方差用振幅²/4再除均值下端²。经验FP guard仍为1e−13，明确不等于严格舍入证书；解析误差也不包含普通BLAS/三角函数/四则的完整roundoff。

## 已执行验证

`validate_source.py` 是确定计算，无随机seed/MC。实际总耗时1.474秒，`SOURCE_VALIDATION.json` 保存源码hash、完整记录及Fraction peak证书。

- R4原子和标签与38逐位相同；4表×L=1,10/9,7/6,11/9的完整Q个evaluate值和全部scalar/diagnostic字段共16次逐位相同。
- R4/8/16/32全部60张表在L7/6分别计算完整Q均值/variance帽；质量与class均值残差最大2.23e−16。
- 共1560个class-bin真实点核reference检查，包含预定cell峰进/退出、class边界及uniform控制。采用按全部breakpoints分片的24点Gauss，权重归一以准确积分常数；最大均值差7.783e−13，小于解析帽7.3452e−12，振幅观测/帽最大.997643，参考variance≤range²/4。
- 对全部15300个box用Fraction验证每格保留质量恰为1，含R32边界touch；class内部最小周期gap≥1/(4K)。峰均匀测度μ的plus/minus质量分别127/255、128/255，加权密度中另乘p。
- 24个包含partial进入/退出的bin额外用Fraction积分piecewise-linear peak，折算kernel误差均低于1e−12。Fraction认证只针对有理peak几何/积分，Gauss与完整trig背景仍是数值检查。

运行解释器 `/opt/homebrew/opt/python@3.15/bin/python3.15`，NumPy2.5.4。验证只读导入原38模块并禁止生成其bytecode，未改冻结字节。

## Builder与streaming边界

`build_cache(R,grid=None)` 返回src、Fraction grid、features、norms、errors，字段与38一致，features预分配shape `(4,R*J,Q)`，逐表逐尺度填，没有list后stack双份内存。`grid=None`保持原生产请求的17点均匀共同尺度 `10/9+j/144,j=0,…,16`。可显式传grid；所有长度须Fraction、已排序去重、与Q对齐并在旧窗口内。

`exact_grid(src)` 提供完整233点自然网格union旧9点，供父线streaming pilot使用；`short_grid(src)`提供17点。`norms=mean(exp(feature0))` 与真正nominal tilt定义一致。errors含每表/尺度解析误差、FP guard、min total、loginflate、p_sum_error及C_Z/C_E；完整实现沿用38公式。

只执行 `build_cache.py` 的小型preflight，未调用evaluate、未分配features、未实际full build。预检的feature字节数如下，小返回arrays另加O(RJ)：

|R|J17|J233|
|---:|---:|---:|
|4|106536960|1460183040|
|8|213073920|2920366080|
|16|426147840|5840732160|
|32|852295680|11681464320|

R8/16/32的完整自然cache仅主数组已超过2GB。builder会在分配前拒绝超过2GB的返回arrays；这不是运行时RSS保证。父线采用streaming：先用选定anchor抽counts，逐表/尺度evaluate并累计充分统计和全部anchor RN，本模块不实现或替父driver审核这一步。

`CACHE_PREFLIGHT.json`记录全部八种组合及Darwin ru_maxrss bytes。preflight耗时.00370秒、峰值约34.4MB；它不是full-cache成本测量。没有任何大cache落盘。
