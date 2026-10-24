# 第32轮局部叶覆盖导出静态审计

结论：冻结导出器及300秒包装器无数学或预算阻断，可执行唯一导出；独立覆盖核验器另在下节记录schema联调。此次只读源码、计划及hash，**没有构建模型、评价点、细分实验、Arb运算或LP**。

实际hash：`export_cover.py` 为 `114298a0fc80a9b995ec9cbfb170abcade2e995122ec69f6b7a7b377bb0481ea`；导出wrapper为 `ed577bef26424b9c525c98f05aedefaa956ed29baefdf6e96eccdfcace51843f`。本轮新导出几何是将旧成功基线局部覆盖转换为可复用输入，不宣称逐字节恢复历史未保存坐标的旧叶，也不冒充全域重跑。

## 调用与接受条件

固定三个父格、公开19个GRID端点产生972根闭盒，各根depth0，中点四分至多depth4，最多1500个实际started盒。基线、历史报告、两个公共源码均pin；192位在导入后设置并重建MC3/MC3T，随后仅构造一个基线Spec/AKChecker，缓存不跨不同核/方向。所有盒在内区并离轴/原点/截断线，节点接端点由公共弦分支处理。

每盒 `call_keys` 保证AS_model恰一次、APsi_lower_const恰一次，只有const未通过时才调用至多一次APsi_taylor。常数API物理端点按其特殊顺序u0,v0,u1,v1正确传入。未间接调用会重复AS的subcell_lower/DFS入口。各组件started/completed/failed/interrupted独立计数，盒上限与组件总上限区分明确；4500是三组件调用上限，不是1500函数调用。

AS与Taylor的Bernstein球端点通过finite/exact+man_exp准确变成Fraction，统一余项取准确上端且要求非负。const模式16行是 `c1*Lpsi+(BernS1_lower-ES)>0`；Taylor模式是 `c1*(BernPsi_lower-EPsi)+(BernS1_lower-ES)>0`。**只有一套模式全部16行严格正才接受整叶**；代码没有逐系数取两种模型之max后混拼，也没有把等于0当严格正。signed T尚未计算，当前只有基线c/s两列。

组件返回但随后到时的调用记completed，盒仍可记interrupted，不掩盖实际调用；组件输出非有限或转换失败时不会接受。每盒弹出即保存in-flight，在pending之外保留其node及stage；预算/异常保留剩余队列、unresolved和在途对象。完成盒的模型不可变文件包含实际AS/const/必要Taylor及selected mode/16行，叶引用其hash。未接受父盒模型也保留并计数。

## 覆盖与终态

精确中点产生四个无缝子盒，父子ID、rootID、depth保存，几何key拒绝重复。这里须区分PLAN与实际实现：**导出器终态只检查每父格接受面积恰1，没有另行逐split执行完整Fraction树分区核验**。面积单独不是覆盖证明；随后主侧独立slab检查提供不依赖递归树的完整覆盖证据。pending尚存、depth4未决或在途时不能geometry_complete。准备阶段失败可能没有完整根清单，仍只能解释为incomplete，不能虚称未访问域很少。

包装器创建专属进程组，300秒外硬包含import/历史输入/保存与阻塞模型；超时SIGKILL、无自动重试或续跑，保留源前后/loghash。单次工具exit0不等于geometry_complete，终态必须另审。没有T模型、候选优化、LP或新C2界。

## 独立Fraction核验器

`verify_cover.py` 的x-slab方法不同于递归树构造：以所有准确x边界切片，每片y区间必须准确从0相接到1，结合闭矩形端点、无正面积重叠及面积1证明三父格完整覆盖。它检查固定baseline/hash和alpha/R/h、local/physical一致、模型hash及parent/local/mode/rows一致；逐行从已存端点/余项重组c/s，核 `c1*c+s==saved>0`，各叶16个Bernstein索引唯一完整。这是纯Fraction/文件核对，不是重新证明超越函数模型。

运行前审读发现过mode大小写联调问题：导出器用 `Taylor`，核验器初版只接受 `taylor`。最终补丁已实际复读，核验器改收 `const/Taylor`，并显式要求pending/unresolved/in_flight空，余项非负、Bernstein端点有序；该问题已解决。90秒外层包装器限制最多1500叶/24000行的独立核验，拒绝旧输出且不重试。

最终核验器SHA256实际为 `4eceaf7a2338205bc57ef2a8021a4028414345ea9b18473a6a78952e7be8191a`；其90秒wrapper为 `13b104c4ab6a64141ec0d880f4cb3a4b54b073ad8e40906c99c82b8e77e84e57`。静态审无剩余阻断，实际终态留给后续 `independent_result_audit.md`，本文件不把未运行核验器当已通过证据。
