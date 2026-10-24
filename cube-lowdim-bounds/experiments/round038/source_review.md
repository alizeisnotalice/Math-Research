# 第38轮固定单盒Psi模型比较：独立静态审查

**最终 gate：GO。** 已全文审`compare_once.py`、PLAN及`run_once.py`，并实际读取核对最终SHA：

| 文件 | SHA256 |
|---|---|
| compare_once.py | `67d0cfd66ec099681c02b58418b1974c62800b3f0345ac8bf642b7d74f71a44e` |
| PLAN.md | `b138a76aa2ba12d260a19fbc03267a676f1aff33aa1870f0cc81a29b90a4e6d5` |
| run_once.py | `a2a152b3be0b8623f5f13d032da99cb2d4ac55c93a9c6cd91d19d7240f43424e` |

该审查未导入数学模块或运行模型、积分、点、LP；运行许可仅为已批准的node615单次模型实验，不预言accepted。

## 输入及模型对应

修正hlcert、B1、旧node615完整模型和旧cover四项准确hash固定。先核node/depth/local/physical/模型hash/基线身份，再准确核h映射及整盒正坐标、u1,v1<R、u0+v0>R。没有换盒或访问新点的选择逻辑。

旧AS仅使用保存的16准确Bernstein区间与非负余项；逐项核索引、区间顺序、原const列及准确行身份，不调用AS_model或旧const函数。固定模块导入后重设192bit，才构造Kern；不复用96bit时已构造的核常数。

唯一数值模型调用`_outside_taylor`使用物理端点。其返回4×4已升阶、未乘c的Psi数组，源码不再升阶或额外乘h/c。每个(i,j)与旧AS同物理盒/同方向/同次数(3,3)，准确新行是`c1*(Bpsi_lower-Epsi_upper)+oldBS_lower-oldES_upper`。Epsi必须有限且非负；所有新端点经man_exp准确转换，显示字符串不参与接受。16行strict>0才accepted，保存旧行及准确差，失败下界不作为真实AK负点。

数学余项、Hbb边界项和近角Hessian审查见`../round18/next_model_audit.md`。本次同坐标分支，不借跨R解析草案扩大实验域。

## 三层调用账本及代理

`recorded`分别限制Psi_model≤1、Hint≤12、acb_integral≤12，调用前先保存started/in_flight及trace，有限返回才completed。返回后时钟超限时已返回层仍为completed，外层未返回调用记interrupted，不双计；非finite或异常立即向外传播并停止。

模块局部AcbProxy的构造器原样委托原acb类，integral原样传递callback、args、kwargs和选项；不改flint全局类。当前Kern实例的Hint代理原样转交a,b,kind，保存Hint_id供真实积分归属。只有代理进行日志/guard，不改原解析性标志或被积式。finally还原模块acb和实例Hint；失败路径没有换积分器、提精度、重试或第二模型。

正常源路径的12次请求分别为中心6次H/Ha/Hb及近角6次Haa/Hbb_int/Hab；固定盒中L均正，各实际积分由AcbProxy独立记账。Hbb额外边界幂不算积分。不能由模型数推报积分数，实际终态以三层counts和trace为准。

## 时限、输出及结论范围

外层60秒从启动前起算，含导入、输入hash、构造、所有积分、准确转换及保存，独占新session进程组超时杀组并等待，终态清理该组。内部guard不冒称能抢占单个积分。已有result/new_model/meta/log拒绝运行，子程序也拒绝旧calls；不覆盖证据。

new_model保存16准确Psi包围、余项、来源及组合行，result保存实际counts/accepted/最小准确行，calls记录三层started/terminal。硬kill时只依最近checkpoint说明进行中的调用，不能冒称完成。

本轮0新AS、0新const、0独立点值、0细分、0T、0LP、0全域检查。成功仅表示这个固定小盒的新充分条件通过，不是父格完整覆盖、LP可行性或新上界；失败亦不否定已证B1。
