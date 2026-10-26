# 083：body-center双coset源接口及一个明确三mixed子族的准确排除

本批完成并停止。未发现同一六维075参数body-center模型的既有准确研究，
因此不是去重actual0。本批证明了新两coset模型/共同半径接口，并在一个
从准确体积ceiling导出的权重 `w=1/16` 上认证三个具体mixed见证。

对“旧127族U + 这三个mixed盒及其反射/平移”这个**明确有限子族**，完整
归一化比值的保守上界约 **1.9695369181296791**，严格低于R075完整
理想下端。因此这个子族在w=1/16不能改善。没有计算其实际union，没有
排除其他mixed盒、其他w或整个两coset模型，也没有新C6。

## 去重的确切范围

只读检查当前lower-clean、主cube-dim-6及可用cube6-lower clones中的
satellite/coset/body-center相关源码/报告，以及039、047、019。
039/047在原整数晶格上改变一轴或两轴parity质量，源位置不变，且在008
间距分析特定一阶项；它们不包含新的半整数coset，不能迁移其排除结论。
旧n2/n3/n4 satellite脚本属于其他维数、自由相位/其他间距，不能当成同一
六维075参数模型已付。六维旧003只有spacing屏及明确“未排satellite”说明。

这不是声称全仓库外部历史中绝不存在相关想法；本批检索的具体上下文源
和SHA在manifest中。没有重跑任何旧探索、科学程序或union。

## 新源、distinct性、胞质量与反射

正根u^6=2,v^6=3，固定075的间距

```
a1=1+u^4/2, a2=...=a6=a=(1+u)/2.
```

新周期正源为

```
unit: diag(a)*n，n∈Z^6，每点质量1；
satellite: diag(a)*(n+(1/2,...,1/2))，每点质量w>0有理。
```

两个orbit逐坐标一为整数、一为半整数，严格distinct。半开原胞
`prod[0,aj)`恰有一个unit点和一个body-center点，真实胞质量为 **1+w**。
坐标反射将satellite整数指标`n`变为`−n−1`，仍在同orbit且质量不变；
unit指标变为`−n`。周期平移和各坐标反射均保持整个新源。

保留旧U仍合法，因为新增正质量只帮助原unit witnesses。H的64个反射
内部不相交并铺满一个receiver基本胞，故有限见证族的比值为

```
R_new=64*|U_new∩H|/(1+w).
```

不能把(1+w)错当体积密度再额外除胞体积，也不能把各body-center边界副本
当作64个不同胞质量。反射索引公式在实际所选sat块上逐轴准确核对。

## 先验ceiling与唯一几何权重

**仅在level λ=1、当前临界半径 `(m+kw)^(1/6)/2` 目标下**，任何receiver
集合至多填满一个物理周期胞，因此

```
R_level1 <= V_cell/(1+w), V_cell=prod aj.
```

自然等权w=1时，准确ceiling上端展示为1.2071733198746126，严格低于
R075完整理想下端。本批仅据此排除该w=1、λ=1目标；这不排该双coset输入
在其他λ的weak ratio，更不排整个模型或C6。

令`h=V_cell_lower/R075_upper−1>0`，实际h约0.1646543377820849。
只选不超过h/2的最大二进分数，得唯一几何权重

```
w=1/16，w<=h/2<2w.
```

这是从准确ceiling余量导出的单个值，不是grid或优化。共考虑2个w：
w=1仅做上述ceiling；w=1/16才实例化三个mixed盒。没有增权重或盒网格。

## 最小两orbit见证接口与三个有依据实例

每个key给unit和satellite各一个非空整数Cartesian块的lo/hi指标。
分别计数m、k，质量严格为`M=m+k*w`，共同六维半径为

```
rcrit=M^(1/6)/2.
```

在坐标j，unit物理极值为`aj*[unit_lo,unit_hi]`，satellite为
`aj*[sat_lo+1/2,sat_hi+1/2]`。必须先取两者真正min/max，再构中心区间
`[source_max-r, source_min+r]∩[0,aj/2]`。不能只把旧r改大。
Cartesian块的每个极值都实际取到，所有源原子每个坐标都在这些极值之间，
所以逐轴核极值充分且必要；这不是遗漏内部源点的抽样检查。

三个实例分别为：

|key|unit块/旧参考|卫星选择|m|k|M(w=1/16)|
|---|---|---|---:|---:|---:|
|minimal_interior|all2，原始54|每个unit坐标相邻点之间的bodycenter|64|1|1025/16|
|max_relative_satellite_count_inherited|`(3,4,4,4,4,4)`，原始126|内部半整数块|3072|486|24819/8|
|unit0_plus_single_bodycenter|unit0，原始0|单个`(1/2,...,1/2)`|1|1|17/16|

前两个使用旧中心unit starts q，satellite整数starts同q、结束q+ell−2，
所以k=product(ell−1)。单位块的极值包住全部选中satellite，扩大共同r
合法，但已实际核真实两orbit极值。第一个是最小interior块；第二个由
原127中逐标签准确最大化k/m选出，并非任意增加块长度。
第三个真正source hull为`[0,aj/2]`，不等于旧unit0的零宽源hull，专门
验证“只增半径”会错误的情形。不是satellite-only盒的体积加法。

通用根接口使用整数有理六次幂包围，D=10^60，逐个准确验证

```
rlo^6 <= M/64 < rhi^6.
```

用向外spacing区间及rlo构严格向内盒，统一源余量eta=10^-35；同时用rhi
构真实mixed盒的有理外包。三个内盒均正宽。每个坐标检查两个轨道的
四个真实极值、两个接收端点及双向距离：144端点/极值节点，共288个
准确有理源不等式。共同r在六坐标相同，极值计数、反射索引及distinct性
全部闭合。所选源及反射的整数指标collar安全上界为2，但本批没有由此
自动生成有限L1证书。

## 完整可靠的新增体积上界，不是actual union

三个mixed真盒记B_i，其外包为O_i。各自旧参考内盒I_i分别来自冻结075
真实原始54、126、0；只读准确端点后向内舍入，故I_i⊂旧U。没有重跑旧127
union或旧科学程序。

逐个计算准确有理数

```
delta_i_upper=|O_i|-|O_i∩I_i|.
```

由于`B_i\U ⊂ O_i\I_i`，有合法保守界

```
|U∪B1∪B2∪B3|-|U| <= sum delta_i_upper.
```

各项的旧U及mixed彼此重叠只会让这个上界更宽；没有把上界说成实际新增
体积。相交部分均用完整6维闭盒交的逐轴max/min和准确有理乘积计算。

三个delta_i_upper展示分别约
`0.000007789296822`、`0.0000277133332622`、`0.000271018474322`，
总和约`0.0003065211044060701`。要超过R075完整理想上端，新增half-cell
体积的安全充分门槛为

```
((1+w)*R075_upper-R075_lower)/64 ≈ 0.002024429321123824.
```

更直接且用于排除的完整子族上界是

```
R_three_upper=(R075_upper+64*sum delta_i_upper)/(1+w)
             ≈1.9695369181296791
             < R075_lower.
```

比较全部用完整Fraction，负差约−0.10347870670111646。因此**仅对旧U加
这三个mixed盒、w=1/16、λ=1的指定有限子族**，不能改善已经准确排除。
没有新union下端或可晋升候选；不把此结论外推至其他mixed blocks、
其他卫星采集位置/质量、其他λ或整个body-center模型。

## 剩余缺口与可认证接口

最小接口已给实际源key、m/k/M、共同r包围、真实两orbit极值、内外捕获盒、
严格源余量、反射规则、胞质量和正确归一化。若继续此新模型，需先提出
能越过本批保守损失的**不同mixed family或参数机制**，再对oldU及全部新增
盒的真实union增量做完整认证；当前三个盒已准确排除，不能反复跑它们。

旧a2仿射g没有继承或参与这次模型计算。新模型源位置、质量和半径均不同，
039/047一阶结论、079旧源局部极大亦不自动迁移。

## 实际执行与封存

唯一主执行PID86492，exit0；4配置=一个ceiling分类+三个新mixed实例，
144新源极值节点。CPU0.006708秒，保存前墙钟0.006847秒；180秒alarm
覆盖exclusive save/打印cleanup。0union、0旧程序重放、0失败或重试。

主脚本已有results即拒执行。saved-only checker已PASS，只核保存Fraction
根包围、源数/极值、盒体积及保守上界，0新几何配置、0union、0写文件。

```sh
python3 -B check_saved.py
```

必要六文件为registration、preflight、results、check_saved、report、manifest。
新源/输入及去重上下文hash和收据齐备。没有stage/commit/push、MATLAB或
递归分派；所有旧轮次未改动。当前任务在有限子族排除处停止。
