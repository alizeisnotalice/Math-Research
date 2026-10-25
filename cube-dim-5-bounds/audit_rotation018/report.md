# 018 · 固定旋转格 level 1 全尺度完整性的有限补审

2026-10-26。**有限补审通过，采用明确的生产证据依赖边界。** 固定精确间距 `a=1+u³/2`、`s=(1+u)/2`、`u⁵=2`，固定旋转 `cos=4/5,sin=3/5`，阈值仅为 1。尺度补全的数学论证、已保存证书的有限不等式、一个独立 side 4 profile 和最终有理汇总一致。

本轮没有独立重建四份全局 histogram，也没有重新运行 22 个区间的全部平面 arrangements。因此下面的完整数值 enclosure 仍依赖锁定的生产 histogram、案例、几何及 interval-profile 证据；本报告不是对所有全局几何输出的独立重算证书。范围不扩张为一般 C5 算子上界、其他角度/间距的排除或本输入其他 threshold 的结论。017 冻结四件没有修改。

## 实际工作与资源

只读使用 `lower_round016/delivery007/` 的 15 件精简交付和原有完整 C 引用；不需要旧 50/52/68 文件包。独立源码 `check.py` 不导入生产模块，标准库单进程运行：PID 24951，nice 10，Python thread 1，六个数值库线程变量均为 1；总时间 0.062217 秒、峰 RSS 28,803,072 字节。运行前生成 `registration.json`。实际结果与输入哈希在 `results.json`。

只重新构建一个必要 profile：T=4 的平面最大计数。没有重跑旧约 3M 个全局 union cells，没有启动 MATLAB/NumPy，没有提交。

## T≥8 的 Minkowski 计数界正确

令 T 为 cube 边长。三未旋轴的格点数逐轴不超过 `T/a+1,T/s+1,T/s+1`。旋转平面按 s 缩放并旋回整数格，捕获区域是边长 L=T/s 的旋转方形 Q。每个被捕获整数点加上 `[-1/2,1/2]²`，所得单位方块内部两两不交，全部落在 `Q+[-1/2,1/2]²` 中。

旋转方形的 x、y 投影长度均为 `(4/5+3/5)L=7L/5`。先加水平单位段再加垂直单位段，面积增加依次为 y 投影长度与扩张后的 x 投影长度，故总面积恰为

```
L²+(14/5)L+1.
```

于是平均值不超过

```
F(T)=(1+a/T)(1+s/T)²(1+(14/5)s/T+s²/T²)/(a*s⁴).
```

正因子各随 T 严格下降。独立整数第五根二分重建 2^-80 的 u 包围，回算 a/s 包围，取分子全部 upper、分母全部 lower，所得有理 F(8) upper 与保存证书逐位相等，展示为 **0.9341795728870501<1**。因此所有 T≥8 不可能有 strict level 1 witness。

## 22 个 [4,8] 区间的尺度逻辑与有限不等式

平面格点投影均为 `(s/5)Z`。窗端点是原子坐标±T/2，且被接收胞元边界 0、5s 裁剪；任意端点顺序变化只可能在 `T=j*s/5`。域边界的事件 `T=2j*s/5` 亦包含在该集合中。三未旋轴最大计数的变化只可能在 `T=k*a` 或 `T=k*s`，后者已包含在投影事件集合。

独立枚举上述事件，加 4、8，精确生成 23 个有理包围互不相交的临界值；表达式系数、顺序、包围端点全部与 setup 一致，故 22 个开区间没有遗漏。每个开区间中，整个 arrangement 的端点顺序及 1D 最大计数固定，平面最大 overlap 也固定。

对全部 22 份**生产保存** profile，独立重算其内部 side 包围、1D 最大计数以及总质量乘积；检查保存 inner/outer 平面最大值相等、两生产算法一致标志及 setup 哈希。逐区间实际验证

```
max_total_mass < (left_endpoint_lower)^5.
```

这里全部为严格小于，比排除所需的 ≤ 更强。22 项精确有理余量已存结果。**平面 maximum 数字在这些 22 个内部 side 上没有全部独立重建**；它们仍是 `radius_certificate.json` 绑定的生产结果，不能把本轮的有限乘法检查写成 22 次独立 arrangement 核验。

临界 T 自身不会漏掉 strict witness：若 N>T⁵，把同中心 cube 稍放大，仍有 N>T'⁵，原原子保持被包含，且 T' 落在右侧开区间；T=8 由 F(8) 处理。因此开区间排除足以排除所有 T≥4，包括临界值。

## side 4 的必要独立 profile 与 T<4 的实际计数

独立直接枚举 `(m,n)∈[-20,20]²`，按物理系数 `((4m−3n)/5,(3m+4n)/5)∈[-5,10)²` 筛出 225 个邻胞原子，并与唯一生产原子清单逐项比对。T≤8 时 r≤4<5s，所以这 225 原子完全覆盖所需接收区：更远胞元至少一坐标距接收域 ≥5s>4。

T=4 重新构建 inner/outer rational atom windows；每侧 2601 个 cell，用直接整数中点逐点计数，独立得到最大值 **15、15**，与 setup 一致。由于 4 不是任何 `j*s/5`，不存在只在 opposite endpoint 接触处新增、却不出现于邻近 cell 内部的 overlap；闭方形最大值也为 15。

计数随 T 单调，故 T<4 的平面实际计数 n≤15。由 a>7/4、s>1，未旋轴实际计数上界为 `(3,4,4)`。因此任何 strict witness 的质量 N≤`3*4*4*15=720<4⁵`。

在一个 1D 间距 b 的格上，长度 l 的连续块中心属于 bZ（l 奇）或 `b(Z+1/2)`（l 偶）。接收 x∈[0,b/2] 时，最近中心分别是 0、b/2；块端点的最远距离等于固定半跨度加中心距离。因此规范块不会比任何实际捕获块更远，可逐轴替换而不破坏同一半径包含。

一个实际 strict witness 满足 T⁵<N；把边长放大到 N^(1/5)，三个规范块和平面至少 n 点仍被同一个 cube 捕获。所需长度≤4、n≤15 已在 l≤4,n≤16 库内；新边长仍<4。故所有 T<4 的 strict witnesses 已纳入有限 family。代码核对全部 256 `(K,n)` 案覆盖，以及 26 个正面积案恰与 geometry bundle 的 26 个 family 对应。**256 案的每个零/非零平面判定没有在本轮重新计算**，依赖 `case_evidence.json` 及其绑定来源。

反向包含在体积意义上成立：有限固定原子见证集所形成的全原子交盒，其内部可稍缩小 equality radius，得到平均严格>1；有限盒的边界是零测集。因此有限 union 与此固定输入完整 strict level 1 集合的体积相同，允许用 inner/outer 数值 enclosure 包围。

## 四份聚合 histogram 和最终 rollup 的依赖范围

生产输入的关系为：

1. `rotation_atoms.json` 给固定物理原子；`case_evidence.json` 保存 256 案、零案和真实源哈希。
2. `geometry_bundle.json` 给按固定 K/n 编号的 26 个 family，以及 inner/outer 未旋盒和携带固定原子 ID 的平面矩形。
3. `aggregate_histograms.json` 给生产从这些几何计算的 plane inner/outer、unrot inner/outer 四份 family-mask measure。其 68 个原分块内容没有包含在精简归档，仅保留来源哈希、x 区间及运行信息；可按需用 replay 分块重建。
4. `result.json` 绑定上述输入以及半径 setup/certificate；`manifest.json` 绑定精简交付文件和完整 C 的已有引用。

本轮核对 manifest 每件哈希、result 的全部输入哈希，及 C 引用哈希。通过从 geometry bundle 提取 x 端点，核对四组生产分块记录索引和边界连续、无缝、无重叠覆盖完整 x 范围（plane 每侧 32 块、unrot 每侧 2 块）。**这不认证保存的每个块内容哈希所指向的计数，也不从几何重建聚合 measures**；不能用“分块索引覆盖通过”代替该重建。

在明确把四份已保存 histogram 当作生产证据输入后，独立以纯整数/有理运算计算

```
V = sum_{plane_mask & unrot_mask != 0} plane_measure*unrot_measure,
ratio = 8*V/(25*(2^40)^5).
```

该相交条件保证同一 family 的平面与三轴几何配对，没有跨 K/n 拼接。两个精确端点与 `result.json` 完全一致；同一次有限算术也精确恢复 K1/K8 和 K1/K2/K8 的历史签名比值。展示结果为

```
[1.73995838416010979469, 1.73995838418300349718].
```

上端点与完整 C 的下端点作严格有理比较，确实较小。本轮据此给出的是：**生产 enclosure 的最终算术与尺度补全论证通过这次有限补审，数值结论依赖锁定的生产全局直方图和未独立重算的 profile/case 输出。** 不宣称已经独立重建完整 all-scale enclosure，更不宣称一般 C5 上界或 threshold 以外的结论。

## 复现和白名单

```
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 \
 nice -n 10 python3 cube-dim-5-bounds/audit_rotation018/check.py \
 --project-root /path/to/cube-dim-5-bounds
```

必要本轮文件仅 `report.md`、`check.py`、`registration.json`、`results.json`；生产输入沿用 delivery007 的独立 15 件白名单与完整 C 引用。源码没有绝对生产路径，不使用旧大包、不改 017、不提交。
