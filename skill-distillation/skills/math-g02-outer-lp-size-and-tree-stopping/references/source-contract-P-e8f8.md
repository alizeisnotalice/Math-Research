# P-e8f8 冻结来源合同与适用边界

此记录纠正同一冻结来源卡的函数空间假设，不新增全文阅读或唯一 PDF 数。当前 active 卡已将误名字段改为原点归零 Bloch，并含显式 source_semantic_amendment；旧字节保存在 G01 references/history。必须先按本合同和原文定义解释旧卡字段。有效卡为 [HK-P-e8f8d143f231a4b3.json](papers/HK-P-e8f8d143f231a4b3.json)，effective SHA-256 `8208885a1c0737ae089d6ec382c02cec17ed8c94e94140b7e027a14275699927`；immutable base SHA-256 `0fa8923d2bed2d9d3d47e6516592b0a542c22e04241105b7e1a6cb50b30768b0`；其原卡 `evidence/reviews/HK/evidence/P-e8f8d143f231a4b3.json` 的 SHA 为 `0f06f06ef797e3445013b7f0ec13b6bd8ff1ed9b913ceb04b83dfac737e55c34`。原 PDF SHA 为 `e8f8d143f231a4b3092291a852cd7d8629a61bbf171d73f07eb45322e12e7e0f`，版本 arXiv2606.16760v2，30 页。

原文 p3 / paper.md L115–147 定义 B˚={f∈B:f(0)=0}。旧卡称“little Bloch”会与边界消失子空间混淆；本入口依原文采用原点归零子空间，不能额外推入另一函数空间。Def1.1 的 resolution 每层是单位圆盘 D 的有限可测分割，下一层每块模面积零包含于前层块，normalized area dA 的条件期望 E_Nψ 对所有 ψ∈L¹(D) 在 L¹ 收敛，每层块排序固定。

来源 Theorems3.1/4.1 / 原卡 theorem_cards[0] 的 locator 是 PDF pp11–16、paper.md L542–818。μ 必须是 D 上有限正 Borel 测度，目标是 B˚ 的 Bloch seminorm 到 L²(μ) 的嵌入。原文 p2 的完整范数是 |f(0)| 加 Bloch seminorm。packet 必须实际来自 reduced Bergman 投影 P₀1_R；原文 p4 的 Γ(R₁,R₂)=∫Φ_R₁ overline{Φ_R₂}dμ 含第二 packet 共轭，旧卡漏共轭也在 active 定义中修正；调用有限矩阵/SDP 前，各 packet-Gram 条目须有限，矩阵为该同一 μ 与同一 resolution 的 PSD Gram。来源将嵌入有界性与 coefficient capacity、diagonal domination capacity 的有限性对应；对完整 Bloch 范数另处理 f(0) 与 μ(D)。比较依赖 Bergman 投影表示、有限维 SDP duality 与 complex Grothendieck；这些外引整套证明尚未独立复证，不在本合同给出未经重证的数值常数。

只有任意 Hermitian A、向量 x 或抽象树不足以满足这一来源合同。A=I₃、x=(1,0,1)、Q=2 只是有限算术；没有 packet 表示、admissible resolution 或私有捕获定义，不能判 Bloch 来源或项目捕获适用。正式树/方向/可行捕获叶/来源映射/共用输入谓词缺失时，停止捕获和迁移结论，保留 evidencepending。来源整套外引证明与私有迁移均 open。
