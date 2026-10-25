# 第八轮 E：给定准确端点次序区间的并集多项式接口

本件只验证 D 交付的有限 affine 盒族及次序区间，不寻找事件，不优化参数。先运行 8 个有理解析玩具；真实盒族最多接收一个，且仅在 D 交付后运行。新候选的比较对象是 **E007 完整理想 R 上端**，不是显示值或较早 frozen 门槛。下界推广仍待主代理的有限 L¹ 与独立复核。

## 1. 输入定义

系数域为正实嵌入的 Q(u,v)，u⁶=2、v⁶=3。每个字段用 36 个 monomial 系数表示（索引 `i,j` 表示 u^i v^j，0≤i,j<6），可用非零系数 dictionary 或长度 36 的 dense list。每个系数必须是精确有理数字符串或整数；JSON 浮点和布尔值均拒绝。

有限盒的每个端点是准确仿射式 `c0+c1*t`。输入 schema：

```json
{
  "dimension": 6,
  "interval": ["0", "1"],
  "endpoints": {
    "x0": {"c0": {}, "c1": {}},
    "x1": {"c0": {"0,0": "1"}, "c1": {"0,0": "1/2"}}
  },
  "boxes": [ [["low0", "high0"], ["low1", "high1"], "...六轴..."] ],
  "order_certificate": [ [["相同仿射端点ID组"], ["下一组"]], "...每轴一项..." ],
  "normalization_factor": "64",
  "candidate_parameter": "1",
  "static_check_parameter": "1/100"
}
```

上例的 `...` 仅说明结构，不是可运行 JSON。`boxes` 每坐标恰有 lower/upper 两个已定义 ID；也支持 flat 的 2d 个 ID。参数范围是有理闭区间 [L,U]，0≤L<U；所有几何排序在开区间 (L,U) 内严格稳定，闭区间边界可有事件相接。静态核对参数必须严格位于内部。

每个坐标的 `order_certificate` 是该轴实际使用的端点 ID 的准确 partition。组内全部端点必须是同一个仿射函数，即 c0、c1 均准确相同。不同组按内部从小到大排序。不是同一仿射函数的端点不能因为小数接近而放进同组。

这个接口假设输入端点确实表示所研究盒的真实几何端点。若原几何包含 min/max，D 必须先在当前区间证明所选分支固定，或再拆区间；不能把尚会换分支的函数冒充一个仿射端点。E 只验证给定 affine 数据与 total-order 证书，不负责替 D 寻找事件或决定 clip 分支。

## 2. total-order 证书的充分性

对于两个相邻组的代表端点，令差 g(t)=high(t)−low(t)。程序准确验证

\[
g(L)\ge0,\quad g(U)\ge0,\quad \max(g(L),g(U))>0.
\]

因为 g 是仿射函数，它在内部是两个端值的严格凸组合，所以 g(t)>0 对全部 t∈(L,U) 成立。相邻验证经传递性给出全部端点的准确内部总次序。若两端都等于零，则两个 affine 函数恒同，必须用准确同组表示，不能把它们当作两个严格排序组。若任一端差为负，证书被拒绝；符号区间不能分离也停止，不用容差合并。

这一步允许事件恰在 L 或 U；不允许次序事件位于区间内部。事件寻找完全由输入方完成。

## 3. 并集体积为何是至多 d 次多项式

在开区间内，每轴各组的整数 rank 固定。一个盒若 lower rank 大于 upper rank，则该轴始终倒置，盒始终为空；若两 rank 相同，则该轴恒零宽，盒为零体积。其余盒在全部内部参数时均为正宽。出生在边界的盒只要内部 rank 正宽便保留，不能因为 L 处宽度为零而删除。

沿首轴全部端点 ranks 切片。两个相邻空间切点之间的条带内，活跃盒的投影族由整数 rank 包含关系准确确定，与 t 无关。条带宽度是一个 F[t] 仿射多项式。递归切片各余维，末端返回零维体积 1。每个完整 cell 的体积是 d 个仿射宽度的乘积，因此次数≤d；有限不交 cell 相加仍次数≤d。d≤6，即得至多六次多项式。

程序用普通多项式卷积保留全部次数，**没有**把 t² 截断成零，也没有使用上一轮的一阶 jets 近似。系数的加/乘在已审 36 系数算术中准确完成。不同盒重叠由条带投影并集处理，不能直接把 raw 盒体积相加。

以上等式先在开区间成立。有限盒并集体积对端点局部 Lipschitz（第七轮已独立证明），所有端点在闭区间连续，因此多项式等式延续到 L 和 U。即使边界发生排序相接、薄层出生或盒子消失，也不影响这个延续。

## 4. 独立静态交叉核对

同一输入的一个小内部有理参数 t₀（默认 L+(U−L)/100）被单独代入所有端点。静态路径重新按数域端点准确排序，重新过滤空盒，并采用直接 F-coordinate 的切片/末维区间合并。它不复用输入 rank 总次序来决定静态顺序，也不调用 D 的静态程序。

要求这个独立静态体积与多项式在 t₀ 的值的全部 36 系数完全相同。只比较显示小数或包围重叠不算通过。

对于真实候选，程序读取 `--comparison-results` 指定的 E007 结果文件完整 `exact_interval[1]`，在候选参数代入并乘准确有理归一化系数后，以完整有理包围的下端作严格比较。不会只超越 `2.06692413818` 就登记新界。

## 5. 有限资源、异常处理与输入追溯

一次运行总计上限：10⁶ 个递归状态、1200 秒、10000 个盒、100000 个端点。参数多项式与独立静态核对共用节点预算。真实模式最多一个输入盒族，无事件循环或优化器。

输入、E007 comparison 文件、共用算术源与本源码均保存 SHA-256。数域根隔离分母为 10⁸⁰，正根包含由整数次幂不等式核验。字段不合法、次序证书不符、内部交叉、符号不分离、静态不匹配或超预算均输出 `stopped`，不输出新 C6 界。

## 6. 已实际运行的 8 个玩具

| 例 | 盒族/机制 | 预期体积多项式或结果 |
|---|---|---|
| 1 | 两区间交叠且重叠随 t 变 | 2 |
| 2 | 两个恒相同的二维盒 `[0,1+t]×[0,1]` | 1+t |
| 3 | 正常盒，加倒置空盒、恒零宽盒 | 1 |
| 4 | 与旧区间不交的 `[2,2+t]` 在左边界出生 | 1+t |
| 5 | 二维两薄宽度 t、2t 同时出生 | 2t² |
| 6 | 六维 `[0,1+t]^6` | 1+6t+15t²+20t³+15t⁴+6t⁵+t⁶ |
| 7 | `[0,1+t]`、`[2,3]`，t∈[0,1] 在右边界相接 | 2+t |
| 8 | 同上，错误扩大到 [0,2]，t=1 为内部 order 事件 | 准确拒绝 |

全部 8 例通过（例 8 正确拒绝），共 37 个节点，约 0.008 秒。每个接受例都有独立有理静态并集与准确多项式交叉核对；完整输出在 `toy_results.json`。该 toy 批次没有运行真实候选或 D 的事件寻找；随后唯一真实盒族的运行另见第 7 节。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cube-dim-6-bounds/audit/lower008/parametric_union/parametric_union.py --selfcheck \
  --arithmetic /path/to/cube-dim-6-bounds/lower/round001_rectangular_certificate.py
```

真实模式使用 `--input ... --comparison-results .../audit/lower007/candidate_number_field/candidate_results.json`；全部输出限定于本件目录。MATLAB、代理扩招、Git 暂存、提交、推送及 ledger 更新均未执行。

## 7. 唯一真实首事件盒族：已实际通过

D 随后交付 `lower/round008/spacing_chamber/chamber_interface.json`，SHA-256：

`99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933`。

其参数为

\[
s(t)=s_{007}+t(s_*-s_{007}),\quad 0\le t\le1,\qquad
s_*=u^2v^2/3+u^2v^4/6,
\]

其他五个 spacing 固定，单位质量普通晶格，归一化系数 64。E 只使用交付的 affine 盒端点和次序证书，没有读取 D 的静态答案系数，没有寻找事件，也没有决定新的参数步长。

实际输入为 127 个 raw 盒、59 个已命名端点。E 准确核验 120 个相邻 order gap；唯一边界相等发生在 t=1、第二坐标的相邻组 11/12，内部不存在次序相等。127 个盒在内部均正宽，几何去重为 117 个。程序准确算出并集多项式的次数实际上为 **1**，完整体积/比值多项式的两个 36 系数向量都在 `candidate_results.json`。

在 t=1/100 的独立静态 F-coordinate 并集与多项式逐 36 系数完全相同。另实际准确比较了比值多项式常数项与 E007 的全部 36 系数，完全相同；因此这里确实从上一条完整理想 witness R007 起步。

候选 t=1 的完整 36 系数与完整有理上下区间也已输出，显示值为

\[
R_*=2.06739733326738152771911793722176824278660056639758\ldots
\]

准确有理下端严格大于 E007 完整理想区间上端，认证差值约为

`0.00047319508586851229325029519124730576618766531605815547`。

新上下区间宽度约 `3.81216185627e-74`。比较输入 E007 结果文件的 SHA-256 是

`b5661c0eb23c264286085053b3c59ab9f992918dcb1dc4a79fdb01e1ddcb6b81`。

真实运行 PID **21185**、exec session **65216**，exit **0**，耗时 **10.675937458 秒**，参数并集与静态交叉核对共 **7159** 个递归状态。真实候选恰好一个，未达到节点或时间上限；进程已结束。真实结果文件 SHA-256：

`6849d3774ab3db1524969715970bbb3aecd07de4467faee5d77b69775cf5ed3d`。

真实运行与最终 toy 批次共用源码 SHA-256：

`743bd799f4783a75075f841dc6cb6451252502c099f214d288a28aaa76b9915c`。

复现命令：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  cube-dim-6-bounds/audit/lower008/parametric_union/parametric_union.py \
  --arithmetic /Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round001_rectangular_certificate.py \
  --input /Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-astra/cube-dim-6-bounds/lower/round008/spacing_chamber/chamber_interface.json \
  --comparison-results /Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-round005-endpointfix/cube-dim-6-bounds/audit/lower007/candidate_number_field/candidate_results.json
```

状态仍为 **pending 真实源几何、主代理有限 L¹ 与 F 独立复核**，程序保持 `new_C6_bound=false`。本件证明交付 affine 盒族的区间多项式和动态门槛比较，不提前替代完整 C6 下界验收。
