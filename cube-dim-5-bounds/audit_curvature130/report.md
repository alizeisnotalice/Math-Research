# audit_curvature130：11330个新曲率单锚点证书的独立整数审计

**局部组件独立资格通过：11330个q、181280个实际闭角点全部与128作者的准确结果一致，5903个child通过、5427个child继续待证。** 当前前沿为104482 passed、57549 splits、5427 pending，覆盖52360根；待证叶分布于1369根。完整g110尚未成立，全局S记0，正式 \(C_{\rm cert}\le C_5\le32\) 不变。主代理保留最终证明接纳责任。

本审计独立实现整数有理数、实际盒角点和前沿结构检查，没有导入作者candidate checker。129理论全部冻结只读；没有H、物理梯度、log、kernel、Arb或旧数值基线求值，没有新split，没有调用MATLAB、Astra或推进已停下界，没有提交或推送。

## 1. 先审源码，再登记执行

本目录先生成准备登记，锁28个实际来源和三份独立源码。主全文读完源码后，授权先运行唯一structure gate，通过后再运行唯一11330q/181280角点批。主授权收据路径及SHA：

`Math-Research-cube-dim-5/cube-dim-5-bounds/round014/primary_audit130_authorization.json`

`b1d1c9aa1bb08087abc4cb0db6ab47cf269534b897fb24f953c8c942666b868c`。

准备登记的 `new_q_and_corners_authorized=false` 保持原字节，未事后修改。结构门通过后，另生成全新 `numeric_run_plan.json`，把结构receipt、独立candidate包、123终态、作者实际结果及123全部最终成功材料锁定，再执行数值门。

| 固定文件 | SHA256 |
|---|---|
| `preparation_registration.json` | `62ded086193321555429617cefdbc17c74bb1fa34e0d35c713a03983ec2b1f9d` |
| `numeric_run_plan.json` | `8df339dbfe7b902f3adebc53e6f9cc27018c283b68948f58348f3c924e8b0088` |
| `integer_core.py` | `47bad8dc97899eda598a801e1a383f3b55dc445ea6b36375f5c73fb93e6c405c` |
| `prepare_structure.py` | `46a4975f2e8ea02bf75c89249cb64671afa55f695623e41367d1389f23692f6e` |
| `audit_integer_corners.py` | `0c85c2807d02822492c71569dc2cdf2dfa3901d9405d89c014b8c42ff7ea7123` |

实际两次数学入口均只执行一次，正常退出，无恢复或成功前缀重算：

| 门 | PID | 脚本记录秒数 | 实际退出 |
|---|---:|---:|---:|
| structure/ancestry/closed-root | 40146 | 36.776655624998966 | 0 |
| 11330q / 181280整数闭角点 | 40426 | 4.905621832993347 | 0 |

运行后实际ps检查确认两个PID均已不存在；没有pool或后台计算。三份数学代码及61条准备/数值输入锁全部再次核hash未变，详见 `final_qualification.json`。

## 2. H64祖叶、H128 jet和完整来源依赖

固定同一五维Lebesgue、alpha1、16条未单位化body方向、权1、共同时间、Fourier逆因子 \((2\pi)^{-5}\) 的模型。目标 \(g=11/100000\)，新组件用 \(H_{1/128}\)；尾质量下界为同锚点的 \(H_{1/64}\ge21/200000\)。

独立重读096完整g105 compact packet的52个chunk，逐chunk核来源hash，核其与主已接受g105收据绑定。105的全部4973待证源，经source mapper的node/epsilon/准确box/源leaf、certificate和record hash，逐一对应原packet中同node的已认证 **H64整叶**。没有从 \(G\ge g105\) 推断H64下界。107的4096个H128锚点、123的6237个H128锚点均落在这些准确原H64祖叶中。

123的5665个失败parent各为两个准确child。5659个child的parent锚点在child外，全部仍在合法原H64祖叶中；129已证明的固定锚点全空间式允许这些外锚点。每个child只使用直接parent的同一H128 record、同一保存Delta和同一保存log-gradient四区间；不按角点切换锚点。

关键源锁：

| 来源 | SHA256 |
|---|---|
| 096完整g105 packet manifest | `c53cef7f2626e68042b73d4bb5e80387841eb810a8fbcc1c511f9ce1d8b4780b` |
| 096完整packet checker receipt | `5ad237758ac62c860e43e1977056b24e0ed9b366854dbd275dfb5aa3139c701d` |
| 123 batch2终态 | `6c1d4c6e075053c427d797798b5ccfc20e62e46a46151511e01a330b38213301` |
| 128作者实际results | `765200e4186f5f519691d1516fec85bf1aaafa10c54f4ff87b33d7aa0ae36d67` |
| 128作者实际registration | `a701dfa9615f7992ec204aad33ceb71bfd5d3fbf29d0a2e40b82f2aea04127e2` |
| 128作者冻结manifest | `61170a06b573237ca82dafe2020ea18f6d8150eb82d8de6cc1149f27e17ff826` |
| 129冻结理论report | `3f03661e18a1066353c89b0ef8b01c0f7f92b880cd5aa8a0647ca004feee21dc` |

旧107新值116、118角点122的冻结成功index/result/manifest及终锁已继承并核hash，不重跑。123独立Sol的两批全部6237点、31185分量、6237fresh Delta、99792角点和两批结构均已实际成功，最终冻结manifest及其24个文件全部逐hash核对，来源依赖在本审计最终资格中闭合：

| 123最终成功文件 | SHA256 |
|---|---|
| manifest | `a961e0fd00c580cb495e951835aab2fa7d15408e754ff39670dc1cfcb9704c17` |
| structure batch1 | `f7f26e1ea07526175f0825759428808573c7034515a835c25da24d3ba323d8a0` |
| structure batch2 | `0f7ea604631d652f835299d7b2d9979660bacc9e83fc50bd5fb60e5ac326ff6c` |
| final check | `f046695b0398106e38dcd294d56ffee14513d1b574ec2f89e4710edce2275858` |

numeric脚本保留的原始summary字段仍写“source123最终资格另需”，这是准备时故意设置的保守输出，未改运行收据。最终 `final_qualification.json` 单独记录已读并锁上述成功来源，不能把原summary当成未完成来源审计，也不能将本轮说成重算了123或旧H。

## 3. 结构门实际重建

独立几何按每轴路径子串直接构造dyadic端点：对根坐标 \(a/32\)，轴子串长n、二进制整数k，盒边界为
\[
 [(a2^n+k)/(32\,2^n),\ (a2^n+k+1)/(32\,2^n)].
\]
它不同于作者逐层Fraction二分。每个parent按path长度模4确定split轴；逐轴验证两个child与parent的边界、共同中点和其他坐标一致，未以体积关系替代精确分区。

每根将深度d的路径转为深度D上的整数区间
\([k2^{D-d},(k+1)2^{D-d}]\)，排序后要求从0开始、逐项相接、恰到 \(2^D\)；同时要求全部严格祖先集合恰等于split集合。这给前缀反链、每根Kraft和闭覆盖，不只核数量恒等式。

实际结果：

| 前沿 | 通过/待证活叶总数 | splits | 根 | 最大深度 |
|---|---:|---:|---:|---:|
| 原096完整g105 | 100998 | 48638 | 52360 | 16 |
| 123原seed：98007pass+6237queue | 104244 | 51884 | 52360 | 16 |
| 123终态：98579pass+11330queue | 109909 | 57549 | 52360 | 16 |

seed的98007个保存通过record按序逐对象hash等于105的96025、107的850、118的1132；旧split逐几何等于096加107。123的572通过parent与5665失败parent恰好消费全部6237原任务；5665失败parent生成全部11330 queue child，一次消费且无重复/遗漏。独立准备包与作者准备包逐对象完全一致，SHA均为
`bef5f9d45108920b7385bfca31a25a97a3ba30aa5ea4ce266cffdd2b47c31112`。

结构receipt SHA为
`0e26a933de8d56dc0e6c32e7840a00d1f3c5c5206b3c260e15bf87075c6760e9`。

## 4. 独立整数q和角点计算

本轮使用标准库int/gcd的约分有理数对，不用Fraction、SymPy或作者任何函数。设同一H128锚点上端 \(U=U_n/U_d\)，直接整数公式为
\[
 q=\frac{1024\cdot200000\,U_n-768\cdot21\,U_d}{200000\,U_n}.
\]
实际每个q均核 \(256\le q<1024\)，正分母来自 \(U_n>0\)。四个保存log-gradient区间分别与保存物理导数端点除以H端点的全部四个比值用整数交叉乘法核包含；只读保存Delta，不重新求log。

同一anchor在每个实际闭角点的下界为
\[
 \Delta+\sum_j\min(\gamma_{j,L}u_j,\gamma_{j,U}u_j)-q\sum_j u_j^2.
\]
按位移整数分子符号选择log-gradient端点，将全部正分母清到LCM后，以一个整数分子判断符号。最差值用交叉整数乘法比较；平手取bits字典序，固定mask遍历与作者16角点的顺序一致。

11330个child全体分别核对作者的q、精确最差margin、worst-corner bits、passed布尔状态、同point/record/保存certificate hash，全部匹配。实际181280角点全部计算，没有按角点混锚点，也没有仅比较最终计数。

结果：新增5903通过叶，5427继续pending；当前104482pass+5427pending仍是原109909个活叶的同一几何前沿，57549split不变，待证涉及1369根。没有把失败child再split，也没有派新H或转入131跨锚点计划。

独立整数结果JSONL SHA为
`3d9ec55efbf766779f2f8b8fd7e2c3afceaae2e915b1d303335f56fd4e163fdd`。
本结果提供严格有理端点符号核验证据；其热核数值来源资格继承已锁116/122/123，理论曲率继承已独立证明并被主接纳的129，未声称重新构建所有底层理论和旧值。

## 5. 复现、精简材料与结论边界

已成功数学对象不应重跑。原结构和numeric命令、完整输入路径/hash、代码锁保存在准备登记与numeric run plan；只读核验可重新hash收据及源码。结构读取采用mmap字段视图，跳过约762MB种子中的无关cache，仅逐项解析所需数组；源完整字节仍全部被hash锁定。

独立candidate包、整数逐对象结果和pending节点三份原输出只作无损gzip/tar打包，实际逐成员解码后按原SHA和bytes核回一致，没有数学重算。精简证据 `evidence.tar.gz` 为4284593 bytes，SHA
`113247417d1fb65d3a7348c84848b096bfcfecb9dc561b997b65d35a4d7a177d`；成员名和原始hash详见final资格文件。原JSONL继续本地冻结保留，不重复列入归档白名单。

最小白名单由manifest逐项给出：本报告、三个数学源码、两个准备/数值计划、结构run登记与receipt、整数run登记与summary、最终资格文件、无损证据包。准备草稿报告及未压缩重复输出不列入精简归档。

最终局部组件资格文件SHA为
`5a125549cfad7acf1757a249ea00c5687a091383fe54f48e99d71f19a52d7a07`。
本轮关闭了128原始11330候选批的独立数值/结构依赖；5427未证叶全部保留，故不构成完整g110，更不构成新C5或I。129和原三四维项目未修改，正式界仍为 \(C_{\rm cert}\le C_5\le32\)。
