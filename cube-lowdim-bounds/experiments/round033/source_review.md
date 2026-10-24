# 第33轮方向 builder / 单次LP静态审计

最终结论：builder、solve、MATLAB及阶段wrapper均已完整静态审读，未见剩余阻断，已通知主侧可启动唯一build及在其complete后启动既定唯一solve。本审计只读源码、计划和hash，没有运行模型、点、矩阵乘法、LP或准确成本。

## 已审求解链

实际SHA256：`solve.py` 为 `a0773677506bb595ea58a92678c6d941c504e02ac9df804e35f3826096d5d190`，`solver.m` 为 `15a713982a8028918b09165e76bf62f9f9e9b5b0087eff6b56d97911b0856327`，`run_stage.py` 为 `4d2f9ff5fb5d4f28da5fc49cebc0a604928c33411be7429610d3d8d3d5eaec51`。

26变量顺序为绝对c、绝对s、12个t、12个z，其他114个旧T自由度不入模。c列使用未乘c1的旧singular缓存，s列为旧288矩阵第1列整体S1；T固定抽列 `[2,3,4,5,198,199,201,202,222,244,267,268]`，并与builder输出方向/行文件hash绑定。最终准确方向映射、质量及抽列语义仍须结合builder审读，不能只信列号。

51574点行保留原目标及half baseline slack；局部16272行的RHS为已存正baseline_slack/2，绝对变量不需另减一次baseline或翻成负halfslack。24条abs辅助行是t-z<=0及-t-z<=0，方向正确。总矩阵67870×26；成本使用完整奇异半质量、有符号S1半质量及各T准确半质量，z成本0。

MATLAB仅一次linprog、dual-simplex-highs、20秒MaxTime为soft，源注明共享HiGHS。输入MAT文件先保存hash，MATLAB结果回传实际文件hash，Python核同；LP_requested与实际calls_started/completed已分开，started紧邻linprog前并在options成功后写。非optimal、无finite26向量拒绝，无fallback、换求解器或第二次求解。

候选检查只作一次2^-50参数有理化，z准确取abs(t)，c/s/t准确盒界核对；全部16272局部行以Fraction无容差检查>=预定halfslack，记录坏行数、最小分数及leaf/index定位。随后51574缓存矩阵行分别按K1e-9/A1e-8核原/加强目标，矩阵乘积非finite拒绝。此为**旧缓存模型行筛查**，未用最终完整a重新执行B/L回代，报告不得称全核独立回代；这些点数值也不构成全域证书。

只有局部准确行、点池门槛与成本gain均过才准确展开 `a=s*a1+sum(tD)` 导出一个候选，不对完整a独立舍入，不修补取整，不自动与基线混合。后续准确成本及全域认证仍必需。

## 已解决的问题与包装边界

初版曾将K/A合并使用1e-8且未绑定build输出hash；现已分别恢复K1e-9/A1e-8并核direction/rows两个hash。另已补MATLAB inputhash、requested/started口径、非finite乘积拒绝及最坏准确行定位，均已局部复读。没有为此执行测试或矩阵算术。

wrapper分别给build300秒、solve180秒外进程组硬限，涵盖imports/准备/MATLAB启动与读写；拒绝已有证据，单次运行不重试。solve断言自己为wrapper-owned组首进程，MATLAB默认继承该组；wrapper正常/异常终态均仅SIGKILL自己新建的组，不触碰共享GUI或他人任务。实际退出及清理仍需终态证据，不能从静态代码声称已完成。

## 最终builder审读

`build_directions.py` 实际SHA256为 `9d48f92efb82220033fbceec6137ed71af02efe402333471805b74406a6fba63`。准确pin第29轮原参数与顺序，从保存坐标振幅排top12，固定4粗/8细及旧列号；选择不是按质量贡献排序，也不构成最优性或独立正性证明。

4粗方向调用已核四倍refinement展开至共同h，8细为单位D4类。准确轨道半质量与原spacing半质量核同，非零索引和<=48；点列oldcolumn由同一名称/原顺序反查，没有只把粗索引乘4或额外乘16。方向文件保存准确映射、half_mass、bound、来源索引及旧参数。

旧1017叶cover及独立Fraction核验均pin，16272条c/s行准确复用并检查正基线slack，不调用baseline AS/Psi或重构geometry。每方向c=0的独立Spec/AKChecker正确使用其自身mlist/缓存；不跨方向按行号错用phi表。192位在导入后设置并重建MC3/MC3T，再构造Spec。固定叶仍在单knot内，端点触节点由公共chord模型处理。

每方向每叶恰一次AS_model，入口key及12204 started上限约束。Arb端点通过已pin helper的finite/exact/man_exp准确转Q；统一余项要求非负，Bernstein中点与半径准确计算，penalty=rad+E_upper>=0。追加至旧c/s行时核leaf ID及Bernstein索引逐项一致，12方向按唯一固定顺序追加，符号与solve中的负penalty一致。

模型按方向流式写JSONL并每项flush；输出保存power/Bern球、准确端点/中点/radius/penalty及Spec hash。完成方向文件保留hash；部分方向与in-flight在异常/超时终态诚实记录，不成为完整LP输入。只有12204模型完成且每行12个mid/penalty齐全才置complete；最终directions/rows各自atomic写入并将两hash写至build_result，和solve读取schema一致。

本结论为运行前静态审计，不将预算、假定性能或待生成产物写成实际成功。终态仍须核真实counts、源前后hash和完整输出；局部16272行不代替其余域认证。
