# 第37轮：外区父格 (52,1) 的 const-only 基线准入

状态：主侧已批准有限范围；export_cover.py 静态实现，待源码审读后由主侧唯一运行。当前没有执行新模型、点、细分或 LP。当前已证界仍为 B1；本阶段不产生新上界。

## 范围与目的

只处理物理 h=1/32 父 knot 单元 (52,1)，固定 B1 准确核 SHA
`21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`，α=β=47/40、R=33/20、N=80。物理盒为 u∈[13/8,53/32]、v∈[1/32,1/16]，其 r_min=53/32>R，距径向截断界至少1/160，严格离两轴。

此前有限 LP 没有为这一外区父格列入连续安全行；第34轮全域核验实际发现失败，第35轮点诊断确认其中有真实 AK 负值。现先检查已证基线能否通过常数型充分行，以便后续为该区域补行。本轮不修改候选、不重跑全域 checker。round13 三旧父格的1017叶及 round14 T模型均不读取、不重算。

采用公共 GRID 的18×18=324根闭盒，GRID固定为0,1/256,1/128,1/64,1/32,1/16,1/8,1/4,3/8,1/2,5/8,3/4,7/8,15/16,31/32,63/64,127/128,255/256,1。根深度0；不接纳则准确中点四分，最大深度4。不改变 GRID，不作非中点种子切分。

## 公式与每盒账本

复用 pin SHA 的 independent hl_core/hl_check，所有模块导入后设 Arb 192 bit，再重建 MC3/MC3T，才构造 Spec/AKChecker。与旧实现共享这些区间原语及公式，不声称独立重推。

每个唯一 (p,q,x0,x1,y0,y1) 盒只允许：

1. 一次 `AS_model(52,1,x0,x1,y0,y1)`，保存 S1 的4×4多项式显示值、16个 Bernstein 准确向外端点和余项准确上端 E。调用已有 chord 节点分支，不假设节点处有界 Hessian。非有限量或负余项立即停止。
2. 一次 `APsi_lower_const(u0,v0,u1,v1)`，参数顺序严格如此；保存准确常数下端 Lψ。每个生成盒都断言 u0>0、v0>0、u0+v0>R。API在此使用外区单调下界和右上角闭式点值。
3. 对16项准确核查 `ρ=c1*Lψ+BernS_lower−E>0`，均严格正才接受。保存线性列 c_column=Lψ、s_column=BernS_lower−E、baseline_slack=ρ。Arb有限 exact lower/upper 通过 `man_exp()` 转准确 dyadic Fraction；显示字符串仅诊断，不参与接受判定。

本轮 **0 Taylor、0 T方向模型、0 LP、0额外点值**。不调用隐含重复模型的 subcell_lower/check_knot_cell。父格虽全在径向外区，但 u 跨过坐标值 R；不能直接套用 APsi_out_taylor 的整坐标分支条件，本轮完全不调用该接口。

## 硬上限、首失败即停

唯一300秒外部 hard，由主侧 run_once.py 执行。内部时钟包括导入、核 SHA、构造对象、计算、转换有理数与写产物，并在盒循环、组件调用前后和转换循环内检查。最多1500个实际 started 盒，AS≤1500、常数≤1500，总组件≤3000；started/completed/failed/interrupted 分列，组件已返回但后续时钟终止时仍记 completed。

**首个深度4仍未接纳盒即停止**：先保存该完整模型与16行、把该盒列 unresolved、保留所有其余 pending，再退出循环。不继续处理其他根，不自动 Taylor、不提高精度、不加深、不删盒、不放宽严格正门槛。非有限、异常、计数上限或时间上限同样保存终态并停止。剩余队列可能包含大块未访问域，不称为已覆盖。

同盒组件键拒绝重复调用。弹出盒后立即 checkpoint in_flight；完成组件并转换后保存 partial 模型/hash；最终 node 模型另存不可变文件。中断保留当前阶段及已有部分模型。硬 kill 时以最近原子 checkpoint 为证，不能把未完成调用冒充完成。

## 输出与成功条件

- `result.json`：状态、源码/输入 SHA、限制、实际调用和盒账本、complete、接受叶/未决数量；状态遵从实际结果。
- `leaf_cover.json`：沿用 round13 schema，含全部324根、nodes及中点子节点、准确 local/physical盒、depth、accepted leaves、unresolved、pending、in_flight；每接受叶含16条准确 c/s 行及模型 SHA。元数据注明 strict_exterior_const_only。
- `models.local/node_NNNNN.json`：各完成盒（接受、分裂或深度失败）的模型、准确 Bernstein 端点、余项与行；`.partial.json` 是阶段 checkpoint，允许同一盒推进时更新，节点保存其当前 hash，不把它当不可变最终模型。

只有无 pending/unresolved/in_flight、全部根由接受叶覆盖、接受叶 local 准确面积和为1且所有准确行严格正，才标 geometry_complete。独立树/矩形覆盖核验由主侧另做；它不新增模型。失败只表示该有限 const-only 准入未闭合，不能据负充分下界否定已证基线的真实 AK。

准入若成功，也仅得到这个父格的基线两列连续行。以后12个 T方向模型、与旧三父格合并以及新增充分结构约束 c≥s*c1，均须另行冻结预算；本轮不实施。
