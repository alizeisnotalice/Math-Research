# audit_lower028_106：固定 period-3 位置路径完整 PWA 的独立代数读回

本次有限独审通过。对 028 已冻结的完整静态数据和 descriptor 一维响应，使用新写的稀疏 Fraction 多项式算术、独立五次根隔离及不同的 mask 聚合路径，核对全部 965 个准确闭事件、964 个 cell、91 个 descriptor 权重、87633 个保存响应的加权结果、963 处闭值连续性和全部残差符号。

仅 epsilon=0 的残差准确为零；其余 964 个闭事件均严格负。结合完整 PWA 资格和首 cell 的已锁定负斜率，得到同一个受限结论：

    R(0)=C_cert,
    R(epsilon)<C_cert，0<epsilon<a.

这只排除固定方向 z=(1,0,−1)、登记 1024 长度族的有限 functional。没有 C5 上界、一般 period-3 位置族全局最优性或一般输入排除。正 epsilon 上没有统一的正 gap；epsilon=a 只作为碰撞处的连续极限，不作为三个互异原子的候选。主代理保留最终证明责任。

本次没有重建 456976 个静态 cells，没有重跑 963 个 union、旧 C 或 025 导数。原始静态 mask 权重与保存的一维 descriptor 响应是明确复用的生产依赖；本报告不把独立算术读回夸称为独立几何重建。

## 1. 冻结来源与本轮登记

生产目录：

    /Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-continue/
    cube-dim-5-bounds/lower_round028/

主提供的四个关键 SHA 全部匹配：

| 文件 | SHA256 |
|---|---|
| archive_whitelist.json | 4de82a528cbfd32d19fef8535c154e721775d89c7ec48d129c0a6bdf2c1df832 |
| pwa/archive_whitelist.json | 281d231c51737e626cf29d0eefb8d3c6e68da0f4e092f445a40a5a8c20273192 |
| pwa/window_result.json | 81469a0008732d591b2028734c49e3ac9584922a896275397246f61649f29ef3 |
| pwa/REPORT.md | 8c1eaab72aca0f2cba1e53e7649b00ac723dedc9f6916d31eec52a01f9a876b1 |

准备白名单列出 33 个文件，追加白名单列出 57 个文件；加上白名单自身、明确引用的 C/025 几何和导数、算术来源及本 clone 旧锁，共登记 132 个只读文件。实际运行前后逐项 SHA 相等，白名单文件字节数也一致。没有修改 028、027、101、068、075、083 或其他研究者目录。

本轮登记 SHA：

    97e4a1de8a413075461322a02e7a54cd0078e4607215b65c715a028b587298f7

审计源码 `check.py` SHA：

    194ea00767236b40a696fe639ac59611c0c3d37a71153e33429fdf6674e97379

登记先于唯一审计执行。所有数学检查在一个 nice10、单线程进程内完成；没有程序 wallcap，但对象数在登记中固定。运行没有失败、重试或扩张。

## 2. 全窗几何资格：不能沿用局部 51 形状

令 u=2^(1/5)、v=3^(1/5)，

    a=1+u^3/2, s=(1+u)/2,
    p_j=j*a+epsilon*z_(j mod 3), z=(1,0,−1).

其余四轴为 s 倍整数，原子单位质量。一个周期三条相邻 gap 为 `a−epsilon,a−epsilon,a+2epsilon`。所以在 `0<=epsilon<a`，首轴严格有序、互异，五维 Cartesian 索引给不同原子；平移周期为 `(3a,s,s,s,s)`，质量为 3。

对长度 ell∈{1,2,3,4}^5，N=∏ell，共同半径为 `r=N^(1/5)/2`，与 epsilon 无关，全部五轴使用同一个 r。尾轴 canonical 连续索引在奇长时以 0 为中心、偶长时以 1/2 为中心；在 `[0,s/2]` 中，它逐点支配该长度的其他平移。尾坐标因此固定，首轴捕获区间确为

    [p_(j+ell1−1)−r, p_j+r] ∩ [0,3a].                     (1)

这里依赖真实有序性，不能在碰撞之后不加检查地延伸 first/last 原子公式。

安全 starts 在全窗成立。若 j<=−2，则 j=−2 时 p_j=−2a，而 j<=−3 时 `p_j<=j*a+epsilon<=(j+1)*a<=−2a`；由 r<=2<2a，receiver 上端在 0 左侧。若 j>=5，最后索引 k=j+ell1−1>=ell1+4，且

    p_k>=(k*a−epsilon)>(k−1)*a>=(ell1+3)*a.

ell1=1 时 `r<=u^3<a`，因为 u^3<2；ell1>=2 时 `r<=2<ell1*a`。故首轴 receiver 下端严格大于 3a。于是 j∈{−1,...,4} 完整覆盖全窗，不是小位移局部估计。

本程序从全部 1024 个长度重新构造固定尾轴区间及精确半径，核对旧几何表每个 family ID、长度、N、radius 和 tail box。得到 213 个正尾形状、811 个固定空或零宽形状，与冻结表逐项一致。由此重新生成全部 6144 条 raw 会计，逐项等于 `capture_forms.json`：

| 类型 | raw 数 |
|---|---:|
| 固定 tail 空/零宽 | 4866 |
| ell1=4 的全窗跨度排除 | 540 |
| 全窗保留 | 738 |
| 合计 | 6144 |

ell1=4 的首末原子相隔一个完整周期，距离恒为 3a，而 `2r<=4<3a`，故 90 个正尾形状的这 540 项全窗首轴为空。其余 123 个形状进入事件；并未套用 025 局部 51 形状名单。长度 3 在大 epsilon 下的出生得到保留，例如 (3,1,1,1,1)、j=0 的首轴宽为 `v−2a+2epsilon`，在 epsilon=a−v/2 后打开。把该盒因原点为空而全窗删除会是实质错误，028 此处没有这样做。

本轮没有对这些 1024 个形状求五维 union；上述是有限 endpoint/radius 资格检查。

## 3. 静态 section 压缩与真实归一化

对固定尾坐标 y∈[0,s/2]^4，静态 mask 保存所有尾盒覆盖 y 的 shape 身份。同一 ell1 内，若 N 增大，r 增大，则对任意 start j 和任意 epsilon，(1) 的 lower 下降、upper 上升。最大 N 对该 ell1 的其他区间形成包含。

因此，在每个静态 section 内只保留 ell1=1、2、3 各自最大 N，是准确等式。该最大 N 来自一个实际属于该 section 的 shape；没有把不同 tail 盒的半径或质量拼成新见证。ell1=4 已全窗排除。

本次独立流读现有 `full_tail_histogram.json.gz` 的全部 139843 条 mask 权重，不重建静态 cells。新路径按每个 ell1 的 N 递减组构造 bit masks，通过组交得到 descriptor；没有调用作者逐 set-bit 遍历函数。检查了 mask 唯一性、每个 bit 都属于全部 213 个资格形状，再逐系数聚合。

准确得到原 91 个 descriptor；其全部权重与 `descriptors.json.gz` 和 `pwa/interface.json.gz` 逐项一致，每项经新根隔离严格正，权重总和严格等于

    (s/2)^4.                                                (2)

又从每个 descriptor 的 N 与 ell1 独立生成最多 18 对首轴 raw form 索引，逐项等于 interface 的映射。没有只核对 descriptor 标签或总权重而跳过逐项数据。

四个尾轴的反射给出 16 个互不交内部副本；首轴使用完整 3a 周期，周期质量为 3。所以

    R(epsilon)=(16/3)*Σ_d W_d*L_d(epsilon),                  (3)

其中 W_d 是物理四维体积，L_d 是该 descriptor 的首轴区间 union 长度。没有额外胞体积除数。对不相交 section 求和是准确体积，不是把重叠 family 区间长度直接相加。

## 4. 279 forms 与完整 965 事件

本程序由全窗保留的 738 个 raw 捕获区间及固定域面 0、3a 重新构造 affine forms 集合，精确等于保存的 279 个不同 `(constant,slope)`，斜率为 −1、0、1。raw-to-form 映射也已在第 3 节逐项检查。

对全部 `279*278/2=38781` 对 independently 计算不同斜率的交点

    epsilon=(d−c)/(m−n).

只保留准确落在 [0,a] 的交点，加上 0、a，得到的集合与全部 965 个保存事件完全相等。另核对每个保存事件的来源 pair 等式、事件编号、严格递增顺序和有理包围。不是只验证已保存的事件各自合法；还检查没有漏掉任何 pair 交点。

相邻事件之间没有 affine form 交换次序。裁切 max/min、区间空/正宽、重合端点、区间交叠和 union 连通分量因此固定。每个 descriptor 的长度以及 (3) 都在开 cell 内准确 affine。965 事件是安全超集，不宣称它们全部是最终 R 的有效拐点。

## 5. 生产一维响应的源码审阅与独立性边界

已手审 `pwa/exact_cells.py` 的两条一维路径：

1. 按准确端点次序裁切、删除非正宽区间、排序并 merge 连通分量，沿最终左右端面传播常数项与整数斜率。
2. 取所有保留区间端点 cuts，逐相邻 cut 判断是否被任一区间完整覆盖，累加一次端点差。

在端点排序准确、原始 form 映射正确的前提下，两种算法都给出准确区间 union。touching 的相等端点不产生长度差；重复区间不会重复计数；零宽出生在完整 arrangement 事件处发生。裁切域面也是 forms 的成员。未发现这些算法的数学接口错误。

但两种生产算法共享同一排序、同一静态权重和同一 Fraction/Q(u,v) 内核；生产的 `both_engines` 标签本身不能算独立代理验收。106 没有重跑它们，也没有逐个从坐标重算 87633 个 descriptor union。106 的运行是：

- 对 963 个新增 cell，在其保存的内部 dyadic 上用新根包围独立严格排序全部 279 个 form；所有相邻包围不交，排序 SHA 与保存值相同。
- 逐项读回全部 87633 个 descriptor 的保存 affine 响应，用新卷积独立加权，得到每个 cell 的 A、B；全部与保存的总体 A、B 完全相等。
- 不把该操作描述为新的 union 求值或有限参数采样。

因此这里保留两个明确的生产依赖：原始 139843 个静态 mask 权重的几何正确性由已冻结 fine-cell/pair-factorization 证书提供；每个保存 descriptor response 的几何正确性由冻结的生产双算法及其执行记录提供。本次核对了它们的完整来源锁、数学接口、映射及所有后续代数，没有独立重建这两层。它们是本次授权复用的证据边界，不以“hash 相等”替代本报告的其余数学审阅。

## 6. 新写算术与独立根隔离

审计没有 import 生产者 `field`、`full_static`、`verify_pairs` 或任何其他数学模块。只用 Python 标准库 Fraction、整数和新写的稀疏多项式：

    Σ_(0<=i,j<5) c_ij*u^i*v^j,
    u^5=2, v^5=3.

乘法直接作稀疏 monomial 卷积，按五次关系分别约简两个指数；加法及标量运算直接作用于 Fraction。所有检查输入的指数范围和有理系数格式均受验证。逐系数零给出准确恒等式；对非零差值必须通过下述根区间分离符号，不依赖浮点容差或显示小数。

重新以 S=2^256 作整数二分，分别用恰 256 步隔离 2^(1/5)、3^(1/5)。验证相邻整数 l、l+1 满足

    l^5 < p*S^5 < (l+1)^5, p=2 或 3.

由此得到宽 2^-256 的正根区间；对正 monomial 分别取两个根的 lower/upper 幂，根据系数符号加到相应界。没有复用旧 10^60 根端点或调用生产者区间求值函数。保存事件和 cell 残差区间都包含本次较细的独立区间。

各有限数组的 equality、连续性和权重使用准确系数，不由区间近似等于零推出。符号检查中没有出现未分离项。

## 7. 全部 cells、闭值和首 cell 的正确继承

Cell 0 的常数项逐系数等于锁定完整 C，斜率逐系数等于已锁定 `lower_round025/plus_derivative.json` 的 D025。106 用新根隔离核验 D025<0，并准确重算 `D025*epsilon_1` 等于首 contact 残差。这里仅消费旧 C/D 的证书，没有重跑旧 C 或其局部 union。

第一 cell 的继承合法依赖于第 4 节完整全局 arrangement：从 0 到首正事件之间没有任何 raw topology 改变，因此已经知道的局部 C+D025*epsilon 身份延伸至整个首 cell。这个论证不把局部 51 形状删减沿用到更远的窗。

其余 963 个 cell 的新 A、B 来自上一节准确加权读回。对每一 cell，106 重新计算两端残差

    Delta=A+B*epsilon_event−C,

逐系数与保存值比较。相邻 963 对闭值全部准确相等；每 cell 的最大端点选择也由新算术核对。16 个压缩包和 receipt 的范围连续、不交、总 prefix=964；cell ID、registration SHA、artifact SHA、source/interface SHA、previous-receipt 链及实际记录的对象数一致。963 个新 cell 对应 87633 个新 descriptor 响应，cell 0 明确单独继承，没有遗漏或重复计入一个 cell。

端点公式随 epsilon affine，max/min 裁切 Lipschitz；一维区间由空集出生时宽度趋于零，有限 union 长度仍连续。静态权重固定，故开 cell 的 affine 身份唯一连续延伸到闭事件。每个 affine cell 的最大值在两端达到。这里不是用一个内部点的值外推，也不是只检查若干参数点。

独立符号结果为：

    closed events = 965,
    exact zero residual indices = [0],
    strict negative residual count = 964,
    strict positive or unresolved count = 0.

首 cell 的负斜率处理紧邻零点的正 epsilon；其他 cell 两端严格负，内部也是严格负。由此支持固定有限 functional 的整个 `0<epsilon<a` 排除。非初始有限事件集合的最大残差仍在 event 1，与冻结结果一致；其接触点 gap 不能被称为所有正 epsilon 上的统一 gap，因为 R(epsilon)→C 当 epsilon↓0。

## 8. 碰撞端点与算子结论范围

epsilon=a 时三个相位位置在同一周期内碰撞，原“不同索引对应不同单位原子”的证明不再适用。106 将最后事件仅作为 (1)–(3) 的连续极限，重算其残差并验证严格低于 C；未把标签重叠的原子按不同物理点重新登记，也没有由此制造一个有限候选。

对所有合法 epsilon<a，几何输入的互异性、周期质量和共同半径正确。若未来某个不同方向或不同有限族给出改进，仍需严格 receiver 收缩、共同半径余量、有限 collar 和同质量平滑 L1 传递。当前只是同一有限见证 functional 不超过既有 C，不能从中推出真实周期输入的全部尺度都不超过 C，更不能推出一般 L1 算子的上界。

本次未发现所审全窗公式的反例或未解释的符号缺口；上述生产静态/union-response 复用边界是刻意保留的独立性限制，而非声称已经取得全新几何证明。

## 9. 实际运行、结果与复现

唯一数学进程 PID **96551**，nice=10，六个线程环境变量均为 1，exit 0。记录落盘前耗时 **17.67457820798154 秒**；不把它称为包含最终写文件的完整 wall。`ru_maxrss` 原始记录为 **91340800**，本机 macOS 单位为字节。没有后台进程、失败重启或额外计算批。

| 实际核验对象 | 完成量 |
|---|---:|
| 全长度 tail/radius 资格 | 1024 |
| 全窗 raw 会计 | 6144 |
| form pairs 完整事件比较 | 38781 |
| 存量静态 mask 读回 | 139843 |
| 独立聚合并核对的 descriptor 权重 | 91 |
| 新 cell 的 form 排序位置 | 268677 |
| 保存 descriptor 响应加权 | 87633 |
| 总 cells / 相邻准确连续等式 | 964 / 963 |
| 闭值符号 | 965 |
| 重建静态 cells / 新 union / 旧 baseline | 0 / 0 / 0 |

`results.json` 保存新根端点、准确残差序列的 SHA 摘要、首个非初始 contact 的独立残差包围、碰撞 limit 包围、运行记录及依赖边界。`checkpoint.json` 是最后完成阶段及 prefix 的记录；最终通过状态以 `results.json` 为准。来源完整性在执行前后都核对，而不只在开头核对一次。

复现入口如下，但程序故意拒绝覆盖已有 results；任何复现应复制到另行登记的验证目录，不日常重跑已通过 batch：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 \
  nice -n 10 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B \
  cube-dim-5-bounds/audit_lower028_106/check.py
```

交付状态：`passed_independent_algebra_readback_reusing_static_and_union_responses`。没有提交、推送、MATLAB 调用或正式 C5 账本更改。
