# 第29轮运行前独立静态审计

结论：最终 PLAN、Python driver、MATLAB batch 及外包装器未发现运行阻断，可由主侧执行既定唯一一次发现。此次只读源码、依赖定义及文件 hash；没有导入研究模块、构造数值模型、评价点、求 LP 或执行成本证书。

实际核对冻结 SHA256：

- `search.py`：`341b32476884efc96dcb85a418db43d6147e6a72c1b22d9821d6492f596ecfc2`。
- `batch_solver.m`：`d529b0e56a9dc8b6f44a41cbffec2f8251ec7e18083d84fa52ca12050a3c2279`。
- `run_once.py`：`9d52c59f5d8fcdc01ba301e90971bf367358a486183bed53a31ad86c16fd7d15`。

## 数学模型与列索引

固定算子 alpha=47/40、R=33/20、h=1/32，beta 仅依次 alpha、6/5、5/4。K首列用真实截断 `Psi_beta`，以 expm1 稳定实现内区并在外区置0；A首列用 `chi(r)*A_alpha(Psi_alpha)`，其中 `chi=(beta/alpha)*min(r,R)^(alpha-beta)`，没有将 chi 移入 A 或乘到样条列。该保守列的解析证明及 `c1/mu` 可行基线已见 `convex_bound_review.md`。缓存是 float，不能把它的“保守”模型身份当向外区间。

128列为新绝对 c、绝对 s、30粗、96细。实际抽列索引 `[1] + [2+COARSE.index(ij)] + range(192,288)` 与旧192/288布局及原 `global_S` 定义吻合：旧第0列是 c1 乘奇异列，不能当新绝对 c 列；代码正确另建此列。粗类按固定 `i+j<=9` 抽取，细类保留既有96类顺序，canonical/不重复/`i+j<=48` 均断言。旧 `direction_order` 和前192数值列还作完整一致性检查。

固定点池 K26604、A24970，共51574行，点序与缓存来源先核 hash；通过 `Geometry.__new__` 装载既有 B/L，不调用构造器或旧闭式。细化仅为准确系数算术，沿用已核 mask，逐粗类核对完整 D4 半质量。粗项展开和细项均支撑 L1 半径至多13/8<R。外区 K 非负结构只适用于 r>=R，非局部生成元仍须联合认证。

每组目标用原 K1 全池值的一半 slack 加强，先拒绝任何负 baseline slack，不 clamp、删行或修 rhs。标准模型使用 `-values*x <= -strengthened_targets`，方向正确。beta=alpha 控制记录新旧奇异列及基线差，但差记录本身不声称严格相等。

## 成本、复查与出口

目标首系数为 `beta R^(2-beta)/(2-beta)`；S1取准确有符号半质量的 float，粗/细方向分别用 H²/h² 及实际 D4重数，不取绝对值，也不重复计轨道。完整展开成本另作相同质量公式回算。保证可行替换本身会提高奇异成本，已在独立凸性 review 说明，不能由基线支配宣称改进。

每 beta 最多1浮点全池批；只有原/加强/raw门槛及成本 gain 均过才进入该 beta 唯一1次参数有理化。先将128参数按2^-50有理化，再用 Fraction 细化方向准确构造 T 和 `a=s*a1+T`，没有逐完整系数独立舍入。c>0、s>=0、T支撑和结构身份再次断言。完整回代使用原 B/L 与展开系数，而非仅相信128列模型；原/加强/矩阵三门槛固定 K1e-9、A1e-8，非有限值拒绝。A回代仍是保守模型，不是 beta 真实生成元的直接值。

两套成本与值的差及 argmax 坐标/尺度均记录；差没有单独设为数学严格误差界。最终资格依据完整回代门槛和完整成本。三组结束后只输出合格者中成本最低的一个候选，或无候选；不会自动认证、修补、重算、导出第二个。candidate 含 alpha、beta、c/a 及准确 `structure={s,T,baseline_sha256}`，与成本封闭接口一致。有限点容差不保证全域非负，准确成本也不保证 C1/C2/C4。

## 配额、IPC 与退出

本轮最多154722个变换后奇异列条目、6批309444个完整保守模型标量值，0新 B/L/奇异原语、0新T列。基线缓存算术、仿射矩阵乘法及 MATLAB 已有模型残差另列计数，未混作新闭式评价。每项先记 started，完成才增 completed；超时可能留下未完成 checkpoint，须结合外层终态。

唯一专属 MATLAB batch，至多3次 linprog，各 MaxTime<=20秒为 soft，不声明 basis warm start或独立HiGHS内核。输入矩阵/目标/盒及点列序在请求前保存，模型文件不可变，各 beta 独立编号。request、started、result、ready 绑定 iteration/input SHA；ready 对结果文件另核 SHA，原子文件发布避免半写入。非optimal拒绝该beta并按固定顺序继续，工具/输入错误全停，无 fallback/重启。

coordinator 断言自己为 wrapper 专属进程组首进程；batch继承该组并核对 PGID。finally发stop，等待后仅终止自有child；wrapper在正常或异常退出后清理其创建的整个专属组，不触及共享 GUI。180秒外层 hard 包括 imports、MATLAB启动与清理；内层计时另列 import耗时，不将内部表观时间冒充全轮耗时。运行结果与清理是否实际完成仍需终态日志确认。
