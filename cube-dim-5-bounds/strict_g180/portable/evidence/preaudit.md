# 180：第二轮 G 证书审计准备

状态：已完成源码和实际输入锁；尚未运行任何本轮数学 reader。仅拥有本目录，独立 clone，无提交或推送。审计与 178 的 I 全量 reader 分工，178 的 P8/down192/166 四次幂字段核验只读继承，未重放。

固定登记为 `plan.json`（SHA256 `a80327c7c0dd81f0bcab3a0095e1530ccfd5a005d823f789ee0cd7cccf04eac9`），代码为 `reader180.py`（SHA256 `d80c831a253fa3f97bd961da3b83bceb60b7bf4effe4caf5dd1d65be0f8d4900`）。80 个实际文件锁覆盖 151 原 typed jets、全部 155 源包与 q journals、保存新 q/余量、9 个原 state、理论源及正式主接纳证书。准备期间发现 clone 的已跟踪快照不含尚未提交的 166 主接纳文件，已明确锁定原 main 的实际接纳文件；保留原输入包并等值复用，未运行数学。

固定 15 个新分片按 canonical ID 的 hash 分配，与旧 canonical_index%15 分配不同。有限合同是 107530 个实际参与源、3812 个 H64 质量支撑、109446 个原 jet 回连、111201 个闭叶、52360 个有序根、58841 个分裂；不得删去任何合同外或数值未分离对象。Delta 合同为绝对值不超过 16，1024bit 有向整数算术、512 阶 Taylor。超范围或失败保留输出并报告未覆盖，不能解释成对象无关。无程序 wallcap、自动重跑、新 H/DD/梯度/网格、Astra 或下界工作。

## 解析接口的独立阅读

已读 075 的固定纤维正产品 Cauchy 表示、coarea 常数1024、固定可逆五列的处处有限性、时间 t^-5 尾及有界 score 支配；微分不要求 Cauchy 绝对矩。BB^T=16I 和固定 lift B^T/16 给固定时间二次罚项 1/(16t²)。129 以固定 anchor 的概率 p_t(c)dt/Hε(c) 使用指数 Jensen，因而得到对所有实位移成立的同一 anchor 支撑，不把 anchor 后验沿空间线段冻结为 Hessian。

153/154 的 Gaussian 混合给 t^5 p_t(c) 非减。用 t^-5dt 权下反向 Chebyshev，完整尾、前 band 和后 band 的 t^-2 均值分别为 2/(3ε²)、7/(10ε²)、1/(6ε²)。因此 Qε≤min(2A/3,A(7/10−8r/15))，其中 A=1/(16ε²)、r=H2ε/Hε。合法同点质量端点 L≤H2ε、Hε≤U 给 r≥L/U，递减分支保持上界；可再取已资格 oldq 的最小值。EPS64 无 H32 质量时 L=0；EPS128 的非零 L 只来自真正覆盖该 anchor 的闭 H64/g105 支撑。不能用 G 或 H128 下界倒推出 H64，也不能要求新 q≥256。

每个证书的 Delta 和四个物理方向梯度都来自原同一 Hε anchor。独立 reader 从正 H 区间和四个物理导数区间作准确区间除法。如果该独立区间比保存 gamma 更宽而包含门失败，应诊断原相关性/算术路径，不能 clip，也不能直接宣称 C5 无效。

新 Delta 检查使用 exp(d)≤Hlower/target，完全不调用历史 log 数学模块。d>0 的每个正 Taylor 项向上取 dyadic 整数，余项不超过 term512·d/513/(1−d/514)；d<0 则先对 exp(−d) 的正 Taylor 有向下端，取倒数获得 exp(d) 的上端；d=0 精确比较。公共 log(23/22) 使用一次反向 Horner 的有理 atanh 多项式及正几何尾，另核同一源 Delta115=Delta110−log_U。

角点函数为 Delta+Σ min(gammaL_j v_j,gammaU_j v_j)−qΣv_j²，全局凹，包含盒外 anchor、负 Delta、闭边、零/重复投影和子空间交叠。先构造每个坐标的两个端点因子，再枚举同一 anchor 的全部16个和；不在角点之间拼接赢家。q≤oldq 和新余量≥保存旧正余量保证原旧目标域保留。完整森林使用兄弟成对自底向上归并至每一个有序根，无旧 Kraft 求和代码。每个已排序 y 的夹断 floor(32y_j) 给其中一个有序根，尾坐标置换对称性覆盖 F 的全部点，包括 y_j=1 的上面。

所有闭叶只使用 EPS64/EPS128，故 H128≥各叶 Hε，在 F 上得到 log(H128/g0)≥m。只有完成全量 reader 后才能把保存的最小 m 与真实 witness 列为本轮资格。正 Taylor P8(m)≤exp(m)，178 已实际核过保存的九项 recurrence 及 g192≤g0P8(m)<g192+2^-192；本轮只绑定其代码、登记、真实结果。H128 floor 仅在 F；先用 G≥H128，再用完整 G 的−4齐次性及 cube 对称扩到 Q。共同半径、Lebesgue、一般 L1 与所有正尺度的解析传递归主全文审计负责。

## 独立性和恢复边界

本轮没有新的 rigorous H enclosure engine。所有实际使用原 H/导数数值区间，包括97112个 H64_inherited052，仍是共同依赖；原 raw hash、模型/epsilon/point、保存区间和新证书运算将被重新读取。不能把新增 reader 的完成标签代替对原 enclosure 来源和完整解析证明的审阅。

typed 标签的缺字段不是由字符串自认证。80锁中的 `upper_global151/data/documents-0001.jsonl.gz` SHA256 为 `63ee68cb41c36f33460305fd68b8be4cf4b9d5b95a7ca463fc1250671ac57d80`。其 `document_id=6cb840b9c1f9690792578b15d553216c0ea2061e044c2702e582f810f97b63b4` 行包含原件 `/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_continue123/binding.py` 的全部 UTF-8 字节，已实际解码并核 SHA 相等。该文件明确声明本模型字符串、绑定 `upper_mixed098/epsilon_kernel.py` 及其 SHA `128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0`。后者原件已只读，固定 epsilon=1/128、3504锥、四个D缩放物理斜率以及相同 Fourier 因子；没有执行。

同一锁定 documents 包还含 046 登记（SHA `a5d9e1d670ff578f73ab690440e135f46e1c68e8e50c61c472a52cc702ee8bcf`），其 model/epsilon/API 字段明确未单位化等权16body、epsilon=1/64、value加四个物理尾导数和全部3504项；046 原 kernel（SHA `df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028`）明确默认 epsilon=1/64 和物理导数尺度。缺 model/epsilon 的三个来源为 H64_inherited052（97112）、H64_new070（105）、H64_new089（2426）；其固定原 state、cache/new_cache 字段和 epsilon 分别由 metadata.numeric_stages 锁定。081 报告 SHA `76276ff907887dd68d984aaabc7a91677f6c92f607d9e51c3e4d2d749fbd15ce` 与 090 报告 SHA `fb83d5017e162bf593d199036939f88a4b78bccb0bdd8f5f9667a133636731ff` 均在该包内，已实际阅读同一 epsilon/模型、物理导数及 raw→task→cache 的资格链；081明确旧178722值不曾独立重算。

old098的194个raw缺 model但有 epsilon；对应100报告 SHA `3f9a3fed53baec4cacb1872616aca03979b7e3e52b87a1043454292a73534f23` 也在该包并已读，明确1/64与1/128两个typed外层桶、194个H128及其970物理分量的独立部分分式审计，以及 H64(2x)×16、物理导数×32的正确变换。其余H128_new107/123/126/137/1679的raw本身均有 model与epsilon字段。它们的文献资格链均由numeric_stages.numeric_qualification_documents引用同一包的真实原件。以上补充属于实际源文本阅读与哈希闭包，不是新增H真值计算。

每256个对象落盘带 hash 链的 chunks 并原子更新 checkpoint。故障保留已经持久化的前缀以及当前待写记录，不自动启动第二次 reader；任何继续需读取已冻结前缀并另登记授权，当前命令拒绝覆盖既有输出。没有声称已经自动执行恢复。

正式主界仍为 U166≈28.936078340172045，便携上端28936079/1000000。下界停止、L_ref未超、没有双侧突破时点、后三轮0/3；本轮审计完成也不表示整个研究任务已完成。
