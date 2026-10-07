# Carleson 测度：方法与验收边界

本入口保留 partial：有限树路由只在显式树、底测度 σ、全部节点系数 α 与子分割给定后检逐节点 packing；当前叶例只产出 packing 常数 1，未提供 f、嵌入输出及匹配来源定理的全套假设，不能作为嵌入正例。UR/PDE 路由另需匹配的域 Ω、调和函数 u、实际局部能量/边界几何输入及指定方向的 UR 或 corkscrew+CDC 条件；缺这些输入时停止该路由，不从叶 packing 改称 PDE 定理适用。

1. 先定义实际T(Q)/S(Q)、底测度、尺度及边界，验证正Borel测度局部有限；写sup_Q μ(T(Q))/σ(Q)，不能只验根盒。

2. 有限过滤树的系数packing检每个J的Σ_{I⊆J}α_I≤Cσ(J)；有互不交子分割与统一mass流时可用已检查4C嵌入构件，树外多参数盒条件须另证。

3. 解析/调和空间的嵌入与反向kernel测试按指定空间定理，先明确核范数归一化；普通盒packing不自动证明所有reproducing-kernel thesis。

4. 邻近P-e9e9211495f0b20e限Ω⊂R^{d+1}有界调和函数局部能量r^{-d}∫_{B(x,r)∩Ω}|∇u|²dist(y,∂Ω)dy≤C||u||∞²；不是任意度量Carleson测试。

5. 其Theorem1.1正向需存在Ω̃⊂Ω、∂Ω⊂∂Ω̃且∂Ω̃为UR；逆向另需corkscrew+CDC。Theorem1.2需两两不交边界E_j及ω(p_j,E_j)≥1−ε，得Σdist(p_j,∂Ω)^d≤C(ε)R^d。

6. 作者Whitney近/远边界拆分→调和事件packing→corona/UR构造的常数依赖d、几何参数及ε；外引容量/UR理论未独立复证，不能称维数无关。

7. 中心立方体迁移须实际树化或帐篷化、共同输入及packing桥梁；区分有界、消失、加权与算子特定Carleson，不从邻近边界PDE定理得全尺度最大界。

## 不可省略的限制

根packing、抽象单树packing、解析kernel测试与边界调和函数Carleson估计有不同前提。邻近来源含UR/corkscrew/CDC与ε依赖，不能把缺Ahlfors假设读成无几何假设。

本项目前仅提供操作方法与待核验来源；高级定理及常数须从真实阅读卡回查原文。
