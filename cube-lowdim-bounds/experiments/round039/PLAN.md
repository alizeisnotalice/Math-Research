# 第39轮：继承父格 (52,1) 的未访问树，补齐连续安全行

状态：范围已批准；extend_cover.py静态实现，待源码审读后由主侧唯一300秒硬包装运行。当前未运行新模型。不续跑或修改第37/38轮原件，不重跑全域核验。

## 固定来源与初始状态

继承 round18/leaf_cover.json（SHA `1b545178105b1b344557e1b0540f0ce5af126555df62cd3cd021ca95556ab6b4`）完整324根、385个已处理盒、311个const接受叶、231个pending及其原始深度。继承 round19/new_model.json（SHA `1cb0a6baed805ce17162179988ba6004acc49df49ef811c8c19a53bc1f7f3298`）的node615多项式准确16行，与round18原AS模型（SHA `ca863a29c6cd9fc479c1e198ae8c6af9f470d12b655f7e2a9120e55323dc4800`）交叉核对。新副本把唯一旧unresolved node615标为已补证接受，保留previous_status及双来源，不回写旧失败原件。

初始接受叶312个；只处理原231 pending及其必要中点后代。原栈次序、节点ID、父子和root ID、local/physical端点不变；深度始终相对最初GRID根，不能在本轮重置。385个旧AS、全部385个旧const、原311接受叶和node615的新Psi均不重算。每个旧已完成节点必须有最终model_file并逐个核SHA，最终文件缺失即停止，不用重算补齐。旧partial仅是本地过程checkpoint，归档不含它们；保留带前缀路径/原SHA并标historical_partial_optional，不读取、不要求存在，也不以它们替代最终模型。

固定B1 SHA `21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`，α=β=47/40、R=33/20、h=1/32、N=80、c1>0。仅父格(52,1)；全盒u0,v0>0、u0+v0≥53/32>R。父格会跨坐标u=R，采用第38轮已审C1零延拓及分片有界Hessian引理，不调用有整坐标分支限制的independent APsi_out_taylor。

## 新盒仅两个模型

每个新开始的唯一盒：一次independent AS_model，再一次修正first Kern._outside_taylor。固定independent hl_core SHA `c3d864b1368ad8150fd3366bf77f687a9a3cdd49b924a20d576aa678aed42ebf`、hl_check SHA `0a5ee85e9b6c250494aaf88abce55bb28a4ba5fc46075c531eeefaa0cc2888f8`；修正hlcert SHA `4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05`。所有模块导入后设置Arb192bit、重建MC3/MC3T，才创建两套核对象。旧模型只读取准确端点；新AS与新Psi共享python-flint，不能称完全独立区间实现。

Psi返回未乘c、已经升至双次数(3,3)的4×4 Bernstein数组，禁止重复升阶或乘c。两模型须对应相同物理端点及u/v索引，准确行：

`c_column = PsiBern_lower − Psi_remainder_upper`

`s_column = ASBern_lower − AS_remainder_upper`

`baseline_slack = c1*c_column+s_column`。

所有有限Arb端点通过exact lower/upper.man_exp转Fraction；不靠显示字符串、float或分开取两个最小值决定通过。每盒16项全部严格>0才接受，否则depth<4时准确中点四分。首个depth4仍不接纳，先保存模型/行和pending，立即停止整个阶段，不继续其他根。

0新APsi_lower_const、0独立点值、0T、0LP、0全域checker；不fallback、不换积分器、精度或模型，不改门槛。非有限、负余项或任何异常同样立即停止。

## 实际调用预算与效率

最多1500个新started盒；新AS≤1500、新Psi≤1500，每个Psi最多12次_Hint请求及12次真实acb.integral，总各≤18000。L≤0时Hint返回0而没有积分，必须分列真实计数。每盒组件键拒绝重复，已有处理盒不能回入新计算。

300秒外部hard含全部导入、核SHA、读取与核对继承文件、构造、模型、Bernstein转换、有理行及持久化；内部盒循环和逐项转换检查时间，每层调用前后也检查。主侧run_once.py拥有独立进程组并负责外部硬终止；本脚本不另启子进程。

AS、Psi、Hint、实际积分各有started/completed/failed/interrupted；completed表示有限数值返回，模型是否已转成并保存准确数据另列AS_models_saved/Psi_models_saved。记录每个Hint的(a,b,kind)、node_id与积分归属。代理原样转交原函数/回调/参数，完成后还原模块局部acb及实例Hint，不改公共源码或flint全局类。

调用事件只追加calls.jsonl并flush，result.json仅保存计数、在途栈、最后事件与摘要，不反复重写增长的calls数组。cover按AS组件、Psi组件和完成盒原子checkpoint；阶段.partial模型推进可更新，final模型不可变。超时后只能用最后checkpoint作为证据，不将未完成积分称完成；未访问域完整保留。

## 输出schema及成功界限

- `leaf_cover.json`：继承全部原树，追加新节点。旧const叶的model_file改为`../round18/models.local/node_N.json`、SHA不变，origin=round18。node615保留原AS模型路径和SHA，再加supplemental_Psi_file=`../round19/new_model.json`及SHA，mode=outside_poly_reused，rows读取第38轮准确行。新叶mode=outside_poly，model_file指向本轮models.local。旧parent/node/leaf元数据只为新文件引用与provenance调整，不改变旧geometry/数值行。
- 新完整模型含AS={bernstein:[index,lower,upper]×16,remainder_upper,power_balls}，Psi={bernstein同式,remainder_lower/upper,unscaled_by_c:true,Bernstein_degree:[3,3]}及rows16、准确几何和源码SHA。未接纳分裂父盒也保存，不隐去其计算。
- `result.json`：新阶段账本与历史账本分开、终态、complete、接受/未决数、来源和产物SHA；`calls.jsonl`是独立append-only层级调用事件。新cover把旧Taylor_calls=0改为historical_round18_Taylor_calls，避免它被误读为本轮外区Psi模型数；新调用以result.counts.Psi为准。

只有无pending/unresolved/in-flight、全部接受行严格正、面积准确和1，才由exporter标geometry_complete；完整分区/树/行来源另由主侧Fraction核验，不以面积相等替代无缝覆盖证明。首深度失败或预算耗尽不得标成功，负充分下界不等于真实基线AK负值。

成功仅表示父格(52,1)的连续基线c/s行完整取得。与旧三父格合并、12个T方向新模型、c≥s*c1充分结构约束及有限LP另立预算。本轮不是新的C2上界，目标3.4仍未完成。
