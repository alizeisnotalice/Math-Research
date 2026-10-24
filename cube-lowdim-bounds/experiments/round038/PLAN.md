# 第38轮：固定 node615 的一次外区 Psi 模型比较

状态：主侧批准静态准备；尚未执行。仅源码审读通过后，由主侧唯一60秒硬包装运行，不续跑第37轮，不改变 node615 或精度，不尝试第二模型。

固定输入：第37轮 node615，local=(27/256,7/64,35/128,9/32)，physical=(13339/8192,3335/2048,163/4096,41/1024)。整盒严格 r>R 且0<u0<u1<R、0<v0<v1<R；避免本轮跨坐标 R 接入。几何由 local/physical和 h=1/32 精确交叉核查。

输入 SHA 固定：

- 修正 first hlcert：`4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05`。
- B1：`21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`。
- 原 node615模型：`ca863a29c6cd9fc479c1e198ae8c6af9f470d12b655f7e2a9120e55323dc4800`。
- 原 cover：`1b545178105b1b344557e1b0540f0ce5af126555df62cd3cd021ca95556ab6b4`。

全部先核SHA，核cover node与模型准确坐标/基线一致，核原16条 const 行的 Fraction 算术；不重算旧 AS、const 或旧函数值。原 AS 只用准确 Bernstein lower/upper及E，不重parse显示幂系数。

导入固定 hlcert 后将 ctx.prec重设192 bit，才构造 Kern。唯一新数值调用为 `_outside_taylor(u0,u1,v0,v1)`：其返回未乘c、已升为双次数(3,3)的4×4 Bernstein数组和二阶余项，不再次升阶。解析适用性承接第37轮 next_model.md，零延拓/跨分支证明待独立审读；本次同分支盒不依赖跨 R 用法。

代理只统计调用，完整转交原 `_Hint` 参数、原 acb构造参数、原 integral函数/回调/积分参数，不改数值公式、精度或积分选项。用模块局部 callable acb代理，不修改flint全局类；用当前Kern实例的_Hint代理，不改源文件。完成或失败后还原。

最多1次Psi模型、12次_Hint和12次真实acb.integral，三层各有started/completed/failed/interrupted账本；Hint记录a,b,kind及integral归属。同盒两坐标均<R，预计12_Hint均实际积分；以实际日志为准，不由预计值冒充完成。首非有限或异常立即停止、无重试；模型返回None也为失败。completed定义为有限返回，非有限记failed。

总60秒外部hard含导入、SHA/旧JSON读取、构造、积分、准确行与写文件，内部每层调用前后及逐行转换检查时钟。一次积分内部可能由外包装硬终止；此时仅最近原子checkpoint为证，in_flight不称完成。不运行旧checker、独立点值、AS、const、细分、T或LP，所有此类计数固定0。

新16条准确行是 `c_column=lower(Bpsi)-upper(Epsi)`、`s_column=oldASlower-oldES`，基线slack=c1*c_column+s_column；c1>0。保存原const slack、差与最小索引，只有全16项严格>0才accepted。零或负充分下界只表示模型未接纳，不当作真实AK反例。

输出：`new_model.json`保存16项Psi准确端点、余项、来源SHA、与原AS/const逐项组合行；`result.json`保存实际终态、accepted、计数/时间/最小准确slack与产物SHA；`calls.jsonl`保存每层开始/终止事件。不得覆盖既有产物或重复执行。

若通过，仅说明该固定闭盒的充分条件闭合，不是完整父格覆盖或新界；若失败即保存停止，不新增点/盒/模型。后续使用c≥s*c1结构约束及连续LP仍须另立预算。
