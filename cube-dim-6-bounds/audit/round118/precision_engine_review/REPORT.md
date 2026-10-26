# 118 精度接口独审

结论：**128-first 的软件正确性门通过，可供主协调者在下一轮采用有界逐点回退接口。**本轮不改117、不部署新 wrapper。四个样本的数值宽度与计时仅是有限校准，不证明所有点都在128位成功，也不保证整体 speedup。

唯一保存 checker（a5a3a3，PID69606，exit0）通过全部40个完整 dyadic 包含：band、六个物理梯度、p_a、p_T、Euler。原512 raw/receipt、准确点和时间窗、来源 hash、新128 stdout及拓扑计数逐项相符；全群967680逻辑项/483840配对求值保持。源码严格等于旧C016仅将 PREC=512 参数化的替换，patch全文一致。编译/链接收据均exit0，128 binary hash正确；当前13个实际include头文件逐字等于封存archive。本独审0 kernel、0编译、0 binary launch。

## 正确性与失效边界

沿用已封存C016公式及pinned Arb/FLINT合同：准确有理输入先以 outward Arb 表示，所有算术、atan、Hermite差商、端点密度和物理链式汇总均继续 outward。因此在合法精度128/256/512和合法数学域下，PREC只改变宽度及能否完成，**不改变外包含义**。该论证依赖原C016数学公式与可信库/工具链，四个相同算法的512参考不是独立新公式证明。

Rat约分、排序、等重节点、零节点和cache相等性均为准确整数判断。比较的int128乘积可容纳int64交叉乘；fan有界ray与六项点积也在int128范围，转换回int64有拒绝门。精度下降不会把不同节点错误认作重复。不同准确节点若Arb差含零就抛错，不跳项；重复节点使用解析Hermite jets。a>0、a<T保证atan合并式分母aT+b²正且无分支跳跃，零/重节点同样合法。

失败处理必须连同wrapper看：输入/范围/缺失fan/ABI/节点分离/非有限输出有异常门，native异常exit2；输出过程中异常可能留下部分stdout，必须凭exit0、完整严格JSON和全部字段原子接受。fan的21-cell/范围检查不证明它是真实fan，wrapper必须核准确fan SHA及source/binary/library pins。Euler含0只是必要一致性检查，不能单独认证六个梯度；旧C016明确记录的全群梯度独立交叉缺口没有被本轮关闭。

`finite/Euler failure` 将非有限与Euler排0合成同一错误，无法准确归因为“精度不足”。不得对此消息静默自动回退；应fatal/quarantine，或在下一轮另行登记的接口中给出不改变数学公式的类型化错误。Euler排0、pins/schema/group/input不一致、整数溢出均是完整性/输入失败，不是升精度即可修复的证据。

## 下一轮准确回退合同

见 `fallback_contract.json`：同一冻结point/window/fan依次128→256→512，每级至多一次、每点最多3次真实native，先前attempt与预算不清零。只对已分类的精度分离失败、明确数值不可用或有效但不满足预登记质量指标的输出升精度；结构/来源/整数/Euler错误停止。接受一个完整同次jet，不混拼字段、不从短小数重造、不为通过而截断区间。512后仍不能满足指标则保持unresolved，不升级到无限精度或静默新增点。

所有级均须验证exit0、完整有限dyadic、status/precision/count、准确输入绑定、Euler含0、consumer预登记的证书质量条件。若有已保存512参考，可另要求全部10字段包含该完整参考；这是一项校准合同，非一般Arb外包正确性的必要条件。没有该参考的新点必须用预登记的width/目标决策条件，不能凭“128返回了JSON”就通过。升级不保证区间宽度逐级嵌套或最终可用。

已识别的精度守卫可在无完整JSON时exit非0：保存错误、丢弃全部partial字段、wait/reap后才按合同升级。其他未分类非0停止。与可信同输入参考仅“不嵌套但重叠”可按校准合同升级；出现不相交则是完整性冲突，不能靠升精度掩盖。三个精度binary须事先逐一批准并hash绑定，不在逐点路径静默编译。

## 执行与范围

作者v2原证据：一次128编译加post-link，四次新native，256与旧512重跑均0。128记录kernel时间约2.91–4.24s；旧512来自别轮独立运行且无可比outer时间，本报告不计算公平speedup。v1保留未运行，不混入第四个原v1样本。

本独审唯一exclusive started，完整await工具往返0.138s；保存结果后/receipt前wall0.036560542s、CPU0.030622s。无子进程，正确ps仅PID69606无行，不声称未知后代全部reap。用户程序无时限，source/状态/旧收据不改。作者manifest_v2=42715472a8a4b2ba9785a2e6098b6d4981b6c406b10a50f92691d55681f2bbcf，final receipt及全部来源在本审输出中绑定。主保留最终数学证明、release与集成责任；无新C6界、下界工作、MATLAB、扩招、push或新方向。
