# 第六轮可观测装载接口独立审查

2026-10-24。最终静态结论：完整审读修订driver及与第24轮差异，局部复读最终regex/快照修正后无启动阻断，已通知主侧可执行唯一既定预算run。本代理不运行新核点、LP、API probe或MATLAB求解，不修改作者源码。上一轮真实结果为0LP的接口停止，不能将本轮修订称作原运行已经修复成功。

## 方案边界

本轮保持α=47/40、R=33/20、h=1/32/N80、190个H=1/8 D4粗类与两整体权重、最多96细类、同一固定K26604/A24970池和激活规则；保持300秒发现外硬限、6次LP各20秒软限、8残差批/480000、方向条目17280000、唯一有理化出口。缓存重建因第24轮未持久化而必要，仍属于本轮真实新增工作，必须重新计数，不假称读取旧缓存或0评价。MATLAB对最后已有intended模型一次30秒共享HiGHS交叉求解，另列≤60000已有行，不产生核点或新核候选。

数学模型与第24轮相同：粗层支撑可到5/2、细层声明窗口41/16。沿用已审准确四倍细化及完整D4质量、真实K/AK无削减、完整细系数残差。实际点池需与既有准确manifest一致，不能借接口修复增加点相位、核族或预算。

## 必须核实的装载证据

首次addRows及后续addRows/addCols前先保存完整原尺度矩阵/rhs、scaled矩阵/rhs、row scales、列顺序/准确行坐标、变量界/目标及基线常数，并保存MATLAB标准模型。这样求解器装载失败后也有可交叉检查的有限模型。文件保存耗时属于300秒，后续增量前同样先存证。

getLp比较对象必须是当前全部scaled intended模型，不是raw模型。先读真实format，再正确转CSR；未知/partitioned格式如果没有明确解释应停。所有维度、rhs上下界、变量上下界、目标及每个系数都要比对，不能只比nnz或范数。getLp的数值“精确相同”指存储的IEEE浮点数据逐项相同，并不使原闭式点行成为数学精确证书。

kError即停。kWarning只在本次API日志明确指向小项删除、实际模型唯一改变恰为模≤读取到的small_matrix_value的输入项置0，其余全等时可继续。仅看到历史日志有过small warning不足以解释当前warning；未知警告/额外变化/读回失败均停。kOk也应保存实际模型核对结果。API实际枚举、前后维度、日志片段及模型diff文件需保存。

对于自由正负变量，删除小系数不保持LP解集，也不是简单的保守强化。可接受它作为有记录的数值发现近似；须保存逐行Σ|ΔA|·max(|lb|,|ub|)盒扰动界，若声称此界严格应由IEEE数转Fraction准确求和，不能以普通浮点和冒称向外界。此误差界不用于放宽原raw残差或K1e−9/A1e−8容差。最终原始完整行和独立完整细展开全池gate不变，才保护数值候选出口。

主侧另有固定3次addRows/getLp的接口probe；它是独立预算记录，非核点或LP。本代理不重复；即使其观察到warning原因，也不能倒推上轮未保存的API enum。

## 最终源码与实际实现

最终search.py SHA256为b5d9f5de7259b5c430ef6f0cb6e7e2c75310a41e2731497739e877d2e54cefd5，实际读文件hash匹配；wrapper沿用已审300秒进程组硬截止版本，hash48caf9325e8a51d3e96266f1a3fdd22640e29a393f2c42ff2a777e7a72efc84e。本审查未导入driver、未新增API调用或数学评价。

数学族、coarse/fine准确映射、D4半质量、真实生成元列及原完整展开残差/成本门槛承接第24轮；diff未改变这些部分。重新构造的准确point_manifest与上一轮文件逐字节hash比较，且实际K26604/A24970、初始K16270/A14698逐个断言，不借本轮修订增加点。B/L、奇异项、初始192方向及baseline缓存持久化到cache.local并记hash/字节数，后续新增方向也保存；声明这是新重建工作，非旧缓存复用。中断前未写完的文件不冒称完成。

LoadAudit读取实际稀疏format并分别按CSR/CSC构造，未知或partitioned格式拒绝，规范去重/消零排序后比当前完整scaled模型。shape及cost、column lower/upper、row lower/upper逐项IEEE相等；系数delta中只接受actual=0且source非零、|source|≤实际读取threshold、delta=−source的删除。其余变化或error日志即停。每次API记录真实enum与当前callback片段，prev_deleted用于区分本次新增删除；warning条数和所报阈值必须与本次新增删除及实际选项一致。kOk若有新删除或warning同样拒绝。最多11次装载读回，不增加核点/LP调用。

删除项按Q.from_float准确转有理，以当前有限变量盒max(|lb|,|ub|)加权得到每行误差和最大行误差，并保存全部删除条目的hex/准确模。它严格包住“所存IEEE模型”在该变量盒的scaled行变化，不是原点态函数的严格区间，也不用于放松原始候选门槛。由自由变量符号可知删项不保原LP解集，driver明确记为numerical approximation；最终raw活跃行和完整细展开全池复查仍使用未删除系数。

save_intended为每次装载保存唯一interface_N_raw_input/scaled_input/input_vectors/row_keys及独立_pre_api.mat，其hash保存在API记录。solved=True只刷新last_model并另写唯一_solved.mat，不重写历史input NPZ或_pre_api.mat，末尾last_model另记finalhash。发现并修正了最初同前缀覆盖快照的风险。按主侧明确范围，初始化addVars/changeCost记录enum；完整先存证承诺限定于首addRows和后续addRows/addCols，不声称首次任意API前均已有矩阵。MATLAB共享HiGHS、一次30秒soft-limit的含义与schema承接前轮，若Python装载停止仍可对已保存intended模型交叉求解。

## 日志回调与warning格式的静态源核对

已读本机安装的highspy/highs.py：构造函数enableCallbacks，cbLogging.subscribe启动kCallbackLogging，Python分发把消息交给event.message。本机未附完整C++源码，因此另静态查阅相同版本官方tag；没有为此重跑3例probe。

官方[HiGHS v1.15.1 HighsIO.cpp](https://raw.githubusercontent.com/ERGO-Code/HiGHS/v1.15.1/highs/io/HighsIO.cpp)的highsLogUser在output_flag启用后，将callback与console/file独立判断。因此当前output_flag=True、active callback下，log_to_console=False不关闭回调；原文附带warning/error级别前缀。

官方[HighsMatrixUtils.cpp](https://raw.githubusercontent.com/ERGO-Code/HiGHS/v1.15.1/highs/util/HighsMatrixUtils.cpp)对模≤small_matrix_value的项删除，非零删除发warning；实际格式使用单数`|value|`。本审计发现初稿regex只匹配复数`|values|`会误停，作者已改为`values?`，并保存log_type与原文，当前API warning不再借历史消息解释。匹配数目/阈值及实际delta仍是继续条件，未知消息依然停止。此源核对不是对原第24轮未知enum作事后推定。

## 启动结论及证据限制

最终剩余阻断为0，已给主侧明确启动结论。允许一次批准的300秒/最多6×20秒探索，不能据静态审查预报其将有optimal解、可认证候选或成本<3.4。若本轮又停，应保留实际装载状态/差异及矩阵，以便区分工程障碍与数学约束；不得补跑同轮或称未完成API已修复成功。新上界仍需实际准确成本及两全域证书，不以warning处理或MATLAB有限模型一致替代。
