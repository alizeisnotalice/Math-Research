# 第五轮唯一运行的结果审计

2026-10-24。**本轮在初始addRows接口处停止，0次LP、0候选、无新上界。没有证明核族不可行，也没有在本轮修复并重跑。** 第一轮已证上界B1继续保留。当前结果只表明已构建缓存与基线检查，以及一次API调用未返回代码要求的kOk。

已只读result.json、discovery.run.json、discovery.log、point_manifest.json及冻结源码。另用本机数学Python对保存计数、整数乘积、hash和文件存在性作有限核对，全部断言通过；没有重新评价点、构列或运行LP。

## 终态与范围

外层status=returned、exit0，墙钟2.8253314169996884秒；driver内部2.593880042004457秒，status=failed_or_rejected，reason=`Rejected: add_rows_not_kOk`。exit0是异常被捕获并保存结果后的正常进程退出，不代表发现成功。300秒预算未耗尽，最多6次LP的配额没有补用或重试。

源码启动前、运行后及当前hash均为169b1e2d33684e7250da45dbf482d8c2faff9d60d59ad9fefbd481209f600407，source_unchanged=true；实际日志hash与run记录491e4267b6d224eea821aa054e9bb5befa0075b7ce0dc3e6f7c8fa78c2d6568d一致。reviewed helper、K1基线、诊断坐标、refinement结果及B/L/奇异公式模块的实际文件hash均与result依赖记录相同。

LP_started=LP_completed=0，solves为空，candidate_exported=false。last_model.mat、candidate.json、matlab_result.json均不存在。未运行准确候选成本或全域checker；MATLAB本轮只进行了选项能力及CodeAnalyzer检查（主侧记录），没有执行有限模型linprog求解。不得将此前CodeAnalyzer通过写成MATLAB解通过。

## 已完成的算术账本

固定验证池实际K=26604、A=24970，合计51574；初始入行索引数为K16270、A14698，合计30968。point_manifest实际长度与结果一致；本审计没有重新生成点池。

| 类别 | 保存started=completed | 独立整数关系 |
|---|---:|---|
| B缓存标量条目 | 20423304 | 2×(161+37)×51574 |
| LB缓存标量条目 | 9888120 | 2×(161+37)×24970 |
| 奇异基列条目 | 51574 | K+A |
| 完整基线构建条目 | 51574 | K+A |
| 基线残差检查条目 | 51574 | 1批×(K+A) |
| 方向缓存条目 | 9902208 | (190粗类+2整体方向)×51574 |

所有上述started与completed相同，direction_column_in_progress=null。仅有1个基线残差批，0次浮点LP解/有理候选残差批；affine_matrix_rows_checked=0。192方向列已构建，未激活细类，距离288方向上限尚有配额。方向条目9902208低于17280000，残差51574低于480000。B/L底层原语条目、缓存方向条目及残差算术属于不同工作单位，不能把这些数字相加当作新的不同核点数量。

保存的K1原尺度基线检查passed，K与A的最大违约均0、bad_indices为空，最低值均0（验证池包含可为0的区域；这不是D内K的最小值）。因此本轮没有观测到基线模型不符；此有限池基线结果也不替代已存在的第一轮全域证书。

## API失败含义与下一步接口要求

源码在初次append_rows调用solver.addRows后，check_status只接受HighsStatus.kOk，其他状态统一抛`add_rows_not_kOk`。原始返回枚举、求解器详细日志与该调用后的模型行数没有保存，因此**无法从本轮证据区分kWarning与kError，亦无法断言矩阵已全部载入或被全部拒绝**。没有HiGHS.run调用，不存在本轮LP的Optimal、Infeasible或TimeLimit结论。phase仍显示cache_A是最后更新的粗阶段标签，真正失败位置由reason和源码确定在初始入行；不能据旧phase否定已完成缓存计数。

下一轮应把原始API状态、操作名称、载入前后行列数及必要诊断写入结果，先保存可复现的有限矩阵/元数据，再进行求解器装载；若容许warning，需先确认警告含义、模型完整性及未缩放残差，而非简单删除检查。接口修正属于下一次新预算运行，不得回写或称本轮已修复成功。当前失败不否定四倍恒等式、Tomar机制或更大核族达到3.4的可能；也没有提供向3.4推进的优化成本证据。
