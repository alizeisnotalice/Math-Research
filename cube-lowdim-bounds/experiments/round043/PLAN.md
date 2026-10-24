# 第43轮草案：在树叶内采用四符号证明区域续接父格 (52,1)

状态：仅计划，未写完整续接代码、未运行数学计算。须在第42轮独立核验完成、主侧及独立接口审读通过后，另由主侧唯一执行。本阶段只补齐基线连续安全行，不求新C2上界或T方向。

## 1. 冻结继承与两种计数

继承round22整树，cover SHA `e78b2bf4e6f5e00974a180cb85dfa9044ed42939930c4add48078e6c4064c8e7`，保留所有root/node/parent ID、local/physical端点、原深度与原216个pending的栈次序。第42轮 `round23/quadrants.json` SHA `1921a3b1e009005ecc5e3bea6309cedc8932e5f8b2862401232974d27eef1581` 只补证旧node663；其来源模型SHA `8b614d16e19e8dddd9cb16793595b99aff3acf680795710c50ced9c3f76c108e`。

新初始状态是 **363个接受树叶、366个接受证明区域、216个pending树节点**：除node663外的362个接受树叶各有一个旧整叶证明区域；node663仍只算一个原depth4树叶，但含4个已证代数piece。不能把它说成新增4个树叶或depth5。原未决node663只在新副本中标记补证接受，来源失败原件不回写。

node615的round18 AS+round19 Psi、node661的round20 T1/AS+round21 signed、node663的round22 T1/AS+round23 quadrants均保留各自双来源/补充SHA。所有旧最终模型和补充文件强制核SHA；历史partial仅可选元数据，不要求存在，不重新计算缺失数据。旧接受行、旧pieces和旧AS/Psi/Hint/const全部不重算。

固定B1 SHA `21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`，α=47/40、R=33/20、h=1/32、N=80、c1>0，仅父格(52,1)。原树深度相对324个公共GRID根，最大仍4，不能因续接重置。

## 2. 每个新树盒的固定流水线

仅处理原216pending及必要的原树中点后代，每个新唯一几何盒：

1. 一次independent AS_model，192bit并重建MC3/MC3T，保存4×4 Bernstein端点及非负ES。
2. 一次signed-Psi bundle：中心6Hint构造F0/Fu/Fv一阶T1，近角两次Hab构造全盒Muv准确上端。共≤8Hint、≤8实际积分；不调用旧完整绝对余项模型或const。源与解析条件沿用第41轮固定修正hlcert及已审C1/分片曲率引理。
3. **直接**在新树盒中心构造4个代数proof pieces，不先跑整盒16行门槛。按第42轮准确de Casteljau限制新保存的T1/AS下端；same-sign correction=0，opposite-sign correction=Muv·δuδv，2×2双线性Bernstein准确升至degree(3,3)。每piece16行，共64条新行，AS仍只减原ES。

每行 `c_column=restricted_T1_lower+signed_correction`、`s_column=restricted_AS_lower−ES`、`ρ=c1*c_column+s_column`。全部64项严格正，才接纳整个树节点。4piece是该节点内部的证明分区，不把任何一个piece当作原树子节点；即使个别piece正、但全节点未通过，也只保存它们作模型证据，不纳入接受前沿或未来LP。

如果不全正且原depth<4，按既定规则在原树四分生成4个child；这些child各自再有其自己的内部proof pieces。父节点已保存的pieces不再作为接受覆盖，避免与子树重复。首个原depth4树节点仍不全正，保存完整64行、模型、pending与失败原因后立即停止整个阶段。非有限、负ES、非法Muv或异常也立即停，不自动扩大预算、换精度或改变模式。

## 3. 预算及真实计数

唯一300秒外部hard，包括读取/核SHA、继承几何与pieces校验、导入、对象构造、模型、准确限制/升阶/行及原子保存。内部在调用前后、盒循环、piece和行循环检查时间；进程组由主侧wrapper拥有。

- ≤1500个新started树盒，AS≤1500、signed-Psi bundle≤1500。
- Hint≤12000、真正acb.integral≤12000，分别记录。L≤0的Hint可直接为0而无积分，不能据预期填完成数。
- 每盒固定4piece/64行，故≤6000个新构造piece、≤96000条新资格行；若中断，按实际started/completed保存，不补算剩余行。
- 新const、独立点值、T、LP、全域checker均0。旧数据读取/纯Fraction来源核对另记为inherited_arithmetic，不计新模型，也不能伪报本轮重算旧证明。

AS/Psi/Hint/积分保持started/completed/failed/interrupted，另存新tree boxes、构造pieces、资格行的真实账本。有限数值返回和已保存准确模型分列。append-only calls.jsonl，不反复重写完整调用数组；result保存计数、在途栈和最后事件。模型按AS/Psi组件checkpoint，64行分区完成后保存不可变proof文件，再做树接受/分裂决定及cover checkpoint。

## 4. cover/模型接口草案

接受树叶仍由 `leaves` 一项表示，显式 `proof_region_count`：

- 旧flat叶保留原16条rows，视为一个whole证明区域，来源SHA不变；不重新数值构建。
- node663保留原模型引用，并加 `supplemental_piece_file='../round23/quadrants.json'` / SHA；`proof_pieces` 为原4块geometry+rows，`mode='outside_quadrants_reused'`，无伪造的整叶16条rows。
- 新接受叶 `mode='outside_quadrants'`，包含 `proof_pieces`（4个piece_id、local/physical闭盒、sign、rows16）及共同AS/Psi模型和proof文件SHA。新叶不提供扁平“16行”假接口。

每个piece保存 `proof_piece:true`、`source_tree_node_id`、`source_tree_depth`，不赋伪造的depth+1。完整模型保留AS/Psi准确端点、ES、Muv、中心系数、两个Hab来源、每piece限制后数组和双线性提升数组；不得只保存判定bool。全新cover保留原件路径/SHA，原件不改。

源码仅导入pin SHA `fd16ccf0309a4ae4627049ff730b2bb233be3948013d25cd514824febcc107b2` 的round23/quadrants_once.py中纯Fraction限制/升阶函数，不调用其main；这项已审来源是复现依赖，另记入input和新模型SHA字段。新bundle保存proof_pieces四块及64行，不使用扁平rows16接口。

全局证明区域键用 `(tree_node_id,piece_id)`；旧整叶piece_id为whole，新四块为0..3。统计分别保存：接受树叶数L、接受证明区域数P、新模型盒数、全部构造piece数（含后来不接纳/被子树替代者）和最终接受piece数。最终P只对接受树叶求和，不把失败父盒pieces叠加。

## 5. 成功标准与未来LP展开

成功须原父格全部被接受**树叶**无缝覆盖，0pending/unresolved/inflight；每个带pieces接受叶的内部4闭盒又须独立准确覆盖该树叶、内部无重叠，全部对应行严格正。主侧/audit另用不同限制算法核准确行和分区；总面积相等不能替代这两层覆盖证明。

未来连续LP必须按**证明区域**展开：每区域16行，故本父格行数是16P，而非16L；初始状态就有363树叶对应366区域。对新接受叶不能只取4piece之一或把64行合并成16行。每区域有自己的c/s列与正baseline_slack；未来若保留半基线余量，应逐区域逐行使用其准确ρ/2。原三父格1017叶及其12204个T模型仍原样复用，不能为了统一schema重算。

后续T方向可以在每证明区域直接建模，或在父树叶建模再准确限制到pieces；二者调用预算与余项紧度不同，必须另立明确阶段，不包含本次。无完整父格cover时不得启动T/LP。即使cover完整，也只取得局部连续安全行，不是全域C2认证，不更新B1或宣称达到3.4。
