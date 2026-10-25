# Carleson 测度：方法与验收边界

本入口保留 partial：有限树路由只在显式树、底测度 σ、全部节点系数 α 与子分割给定后检逐节点 packing；当前叶例只产出 packing 常数 1，未提供 f、嵌入输出及匹配来源定理的全套假设，不能作为嵌入正例。UR/PDE 路由另需匹配的域 Ω、调和函数 u、实际局部能量/边界几何输入及指定方向的 UR 或 corkscrew+CDC 条件；缺这些输入时停止该路由，不从叶 packing 改称 PDE 定理适用。

1. 先定义实际T(Q)/S(Q)、底测度、尺度及边界，验证正Borel测度局部有限；写sup_Q μ(T(Q))/σ(Q)，不能只验根盒。

2. 有限过滤树的系数packing检每个J的Σ_{I⊆J}α_I≤Cσ(J)；有互不交子分割与统一mass流时可用已检查4C嵌入构件，树外多参数盒条件须另证。

3. 解析/调和空间的嵌入与反向kernel测试按指定空间定理，先明确核范数归一化；普通盒packing不自动证明所有reproducing-kernel thesis。

4. 邻近 P-e9e9211495f0b20e 限 domain Ω⊂R^{d+1}, d≥1。性质 (a)：对所有有界调和 u、x∈∂Ω、0<r<diam(Ω)，r^{-d}∫_{B(x,r)∩Ω}|∇u(y)|²dist(y,∂Ω)dy≤C||u||∞²。

5. 性质 (b)：对所有有界调和 u 和 0<ε<1，存在 g∈W^{1,1}_{loc}(Ω)，||u−g||∞<ε，且存在 C=C(ε,Ω) 对所有 x∈∂Ω、r>0，r^{-d}∫_{B(x,r)∩Ω}|∇g(y)|dy≤C。Theorem 1.1 A 向：存在 Ω̃⊂Ω、∂Ω⊂∂Ω̃、∂Ω̃ UR 即推出 (a),(b)；B 向还要求 corkscrew (1.4)+CDC (1.7)，且 (a) 或 (b) 任一成立。Theorem 1.2 两向都要求 (1.4),(1.7)，ε₀依赖这些几何常数。A 向在 (a) 或 (b) 成立时，对每个 0<ε<ε₀ 存在 C(ε)，所有 x∈∂Ω、R>0、p_j∈Ω∩B(x,R) 和两两不交 E_j⊂∂Ω 若满足 ω(p_j,E_j,Ω)≥1−ε，则 Σ_j dist(p_j,∂Ω)^d≤C(ε)R^d；B 向为反向：若对某个 0<ε<ε₀，上述 implication 对所有可行情形成立，则 (a),(b) 成立。

6. 作者Whitney近/远边界拆分→调和事件packing→corona/UR构造的常数依赖d、几何参数及ε；外引容量/UR理论未独立复证，不能称维数无关。

7. 中心立方体迁移须实际树化或帐篷化、共同输入及 packing 桥梁；盒、帐篷与所有并集条件分别定义，只有已证明的桥梁才允许迁移。Holmes–Psaromiligkos–Volberg Theorem 1.11 虽印出 product-bi-tree 分离，但显示的 `μ_A` 分量在顶盒 Q 的计数与 p.15 约 N 的说法不合；`μ_A` 只是总测度 `μ` 的一部分，其他分量及其对两侧的影响未说明。将其记为 source-proof gap，不当作已验证反例，也不判断定理真假；精确核对见 `audit_current/independent-g/P-e38596-top-box-counting-proof-gap-20261007.json`。区分有界、消失、加权与算子特定 Carleson，不从邻近边界 PDE 定理得全尺度最大界。

## 不可省略的限制

根packing、抽象单树packing、解析kernel测试与边界调和函数Carleson估计有不同前提。邻近来源含UR/corkscrew/CDC与ε依赖，不能把缺Ahlfors假设读成无几何假设。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。
