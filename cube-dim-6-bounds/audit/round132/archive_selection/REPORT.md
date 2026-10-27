# C132：跨真实中心的正带完整盒下界（已完成局部批次）

目标是在冻结C127恢复cache15193与静态geometry seed129上，试验对每个相同时间带取不同真实中心的完整box下界最大值，再对不交正时间带求和。A057全空间正带log-Hessian下界允许真实中心在目标box外；原源jet中心、窗口与梯度保持原件。完整数学证明及区间f-clamp公式在原C132/PROOF.md。不同中心的同band证书没有相加；只对5个不交band相加。此机制不声称minmax等号。

固定首个lex未闭原root C121.r11011101101的93个原未决box，原nodeID排序；既有通过baseline没有重跑。每band最多16个欧氏最近真实中心，精确整数格距离排序，325个独立源jet。cache各band可用中心数14670、523、0、0、0；后三band按0。Arb512，nice10单CPU、无时间上限。唯一实际批次完成95232次角点exp、0 native，退出状态LOCAL_BATCH_COMPLETE。目标1/63500局部通过1/93，参考1/65000通过9/93。所有已选box未在新目标下通过；没有新目标root闭合或全域g/C6结论。

原producer PID30702，result.json记录批次内部wall5.169467375002569秒。主提供外部完整childwall5.470711208秒、CPU4.691951秒、exit0及主归档提交标识440a18；外部计时与内部wall口径不同。本清单只据已存在文件作hash归档，不重新运行exp/checker，独立审计另行承担。

conditional_M6_check.json使用原full M_power6_upper完整Fraction核得(9/8)^6*M6/(384^5*(1/63500)^6)<64^5。它以未来全域g>=1/63500为条件，而该条件并未由本局部试验证明。不得据此宣称C6<64。

白名单包含实际执行代码/注册与双gate、冻结选择与源hash、source-only核验、准确M6条件式、全部93份角点证书、批次结果/开始/退出receipt，以及本封存说明。原15193源cache、geometry seed、source jets/raw/receipt与环境是hash引用，不重复复制；原live129不读不改。进度文件、未执行detached launcher无必要重复归档。重现输入选择可读原prepare.py；实际科学批次不重跑。C133另以这些现存角证书做唯一Fraction诊断，不追溯改变本轮结果。

已知边界：原发现搜索有一次误展示live129 checkpoint文本，未用作选择、数据或证明；该范围错误原input_hashes.json已披露。证书仍依赖原真实jet与A057源证明；局部PASS不能代替全域闭覆盖。下一步只允许经独立审计后的有限家族诊断/明确新预算，不能盲目细分。
