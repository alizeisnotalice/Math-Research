# 028追加：固定period3位置路径的完整exact-PWA证书

对准备阶段固定的 `z=(1,0,−1)`、`0≤epsilon<a1` 路径，已完成全部965闭事件、964开cell的准确有限functional认证。结论为

```
R(0)=C_cert,
R(epsilon)<C_cert  对所有0<epsilon<a1.
```

闭窗口连续扩展的最大值准确等于原完整C，仅在epsilon=0达到。没有新赢家、有限epsilon认证或Cn改进。本结论只排除该**固定方向、登记1024长度有限functional**的非局部窗口，不是一般period3位置模型、一般输入或全尺度上界。

本文件是新授权后的追加证书；原028准备报告、源码、静态数据、965事件及34文件白名单保持冻结，不修改其哈希。

## 完整有限输入和未重算的部分

全部213正tail形状的139843个static masks、全部6144 raw会计、全窗ell1=4排除及exact91-descriptor压缩，均沿用准备材料。第一轴原子在epsilon<a1有序、互异，period mass为3，共同半径 `r=N^(1/5)/2` 恒定，归一化为16V/3。没有套用旧51形状局部过滤。

cell0 `[0,epsilon_1]` 已由准备阶段的完整PWA资格和原C/025导数证书锁定。本追加直接继承该cell的A=C、B=D025，没有重跑025/026、旧C基线、static histogram或cell0的union。

因此逻辑inventory为964×91=87724个descriptor-cell对象；其中91个cell0对象由原证书身份替代，实际新计算963个cell、87633个descriptor dual-union对象。每个新cell仅计算一次。

## 为什么一个内部数足以生成整个cell的准确公式

准备材料收集全部合法raw的capture端点及固定0、3a1裁切面，准确去重为279个 `c+m*epsilon` 仿射式，m∈{−1,0,1}。965事件是这些式两两相等在闭窗口中的完整安全超集。91-descriptor的端点均是它们的子集；相同ell1选最大N的身份已经证明对整个窗口准确，未随epsilon重选或近似。

在相邻事件之间，没有不同仿射端点相交，所以全端点次序固定。max/min裁切、区间正宽/空集、重复端点身份、区间交叠和并集连通分量身份随之固定。每个descriptor的首轴并长度准确为

```
L_d(epsilon)=alpha_d+beta_d*epsilon,
R(epsilon)=A_cell+B_cell*epsilon
         =(16/3)*Σ_d W_d*(alpha_d+beta_d*epsilon).
```

程序在每个开cell选择一个简单严格内部dyadic，只用于确定端点排序；不运行参数点体积估计或插值。对279个端点，先取Q(u,v)的向外有理包围，再加上准确的m*q。所有相邻值区间必须严格分离才接受排序；实际963个新cell全部通过该精确检查，不使用浮点tol。

随后沿所选端面直接传播常数系数及整数速度。对每个descriptor分别运行：

1. 一维区间排序merge，按并分量端点差生成alpha、beta。
2. 独立端点cut-coverage mask，按被覆盖小区间端点差生成alpha、beta。

两个算法在所有87633个新对象的两系数逐项相同。算法共享已认证端点排序与精确Fraction/Q(u,v)算术，不声称独立底层算术或独立代理审计。

零宽出生会在raw lower=upper的完整事件处发生；在开cell正宽判定按固定排序完成。重复端点/区间按相同global form身份合并，覆盖体积只计一次；裁切域面也在global forms中。较大epsilon新出生的长度3盒继续保留，没有因为基线为空而删除。

## 闭事件连续性与准确C比较

所有raw端点是epsilon的仿射函数，max/min裁切是Lipschitz，区间由空集打开时宽度连续趋0；有限区间union长度仍Lipschitz。tail sections和权重固定，所以R在整个闭窗口有唯一连续扩展。故开cell生成的仿射公式准确延伸到两闭端点，而非仅在内部数处有效。

全部963处相邻cell在共同事件的值，均在Q(u,v)基系数中严格相等。每个cell的仿射最大值在端点达到；因此965个闭事件的完整比较覆盖整个窗口。

比较直接减去 `lower_round003/a2_equals_a3_results.json` 的完整C系数，再做有理根隔离：

- event0/epsilon0残差的全部系数准确为0；
- 其余964个闭事件残差的隔离上端全部严格负；
- 没有符号未分离项，也没有正残差。

所以全闭窗口max准确等于C；对正epsilon，首cell有已锁定负斜率，其余cell两端均严格负，因此每一点都严格小于C。

这里不存在正epsilon上的统一正gap：epsilon→0+时R仍趋向C。除event0之外，**有限闭事件集合**中的最高值在event1，展示为 `1.9938368324333606`，其与C之差约−0.00047962913512514751。这个contact gap不能误称为所有epsilon>0的统一gap。

epsilon=a1的碰撞连续limit展示约1.2818155721024915，严格低于C。它仅用于upper/连续扩展，不按3个互异单位原子登记witness，也没有以碰撞参数作新输入认证。

## 有限批、恢复与实际资源

先单独登记16个新cell1–16（cell0继承）。PID87548用0.3430秒、峰值33357824字节完成1456个dual对象；准确最大contact为event1，剩947个cell。随后15批逐批先登记实际prefix、源码、interface及state hash，每批≤64个新cell，最后一批51个。

每个cell完成后原子保存完整A/B与91个descriptor响应；cell内每16个descriptor有准确inflight prefix。已经完成的cell只读复用，代码拒绝重算；发生新正证据会写候选LOCK并停止批次优先认证。实际没有fault、恢复重算、正LOCK或wallcap。

| cell范围 | PID | 秒 |
|---|---:|---:|
| 1–16 | 87548 | 0.3430 |
| 17–80 | 88178 | 2.1846 |
| 81–144 | 88207 | 2.7356 |
| 145–208 | 88231 | 2.2843 |
| 209–272 | 88259 | 2.2613 |
| 273–336 | 88277 | 2.2625 |
| 337–400 | 88298 | 2.2873 |
| 401–464 | 88316 | 2.2840 |
| 465–528 | 88345 | 2.1990 |
| 529–592 | 88362 | 2.1522 |
| 593–656 | 88390 | 2.0769 |
| 657–720 | 88420 | 2.0418 |
| 721–784 | 88446 | 2.0872 |
| 785–848 | 88465 | 2.0553 |
| 849–912 | 88482 | 2.2469 |
| 913–963 | 88504 | 1.6708 |

963个新cell累计33.1726秒，单进程峰值41746432字节。所有数学进程nice10、单线程，六个数值线程环境变量均为1；不设wall上限。固定有限inventory完成即止，没有参数网格或无限搜索。

## 可直接审核的最终archive

每批压缩cell artifact包含准确A/B、两闭残差、内部排序证书及全部91响应，16批共1429442字节。它们是scratch individual cells的必要完整副本，未来审核不需要外scratch或重算union。cell0另存为继承证书。

最终 `audit_archive.py` 以这些archive为输入：重新核对所有raw form映射、每cell279个端点的严格向外排序、87633个保存响应的加权A/B、965闭值、963相邻连续性等式、逐批registration/receipt链和完整prefix。PID89488用12.7060秒、峰值50905088字节完成，没有新union或参数点计算。

审核命令：

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 nice -n 10 python3 cube-dim-5-bounds/lower_round028/pwa/audit_archive.py
```

`window_result.json` 保存完整C等式、965闭事件符号分类、最大值、最接近的非初始contact、碰撞limit及全部实际PID资源。`archive_whitelist.json` 只列本追加必要文件，并引用已冻结的028准备白名单；不把964份重复scratch文件归档。

未修改028准备包、025/026或其他研究者文件，无MATLAB、提交或推送。结论只对固定z路径和登记finitefunctional成立，未推广到其他大幅period3位置方向、其他长度/输入或Cn上界。主代理保留最终proof、L1和collar审阅。
