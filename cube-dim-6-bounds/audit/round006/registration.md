# 第六轮 F：时间矩与高阶积分接口独立审计

本批承接冻结门槛和现有 body32 模型，不重启核搜索。模型为 GPT-6-astra xhigh，由主依据跨轮全域 g/I 瓶颈提升；主负责最终集成。

- 数学问题：时间混合的高阶导数能否避免全截时势支付最短时间 ε 的负幂；B 的有限幂级数/Peano 接口是否有效。
- 只写本目录，不 commit/push/stage；不修改旧证据。MATLAB 不调用。
- 资源：本机 Python、SymPy、python-flint；零新 kernel 点；只计算 m=2,4,6，ε=1/1000，有限有理矩与标量 Arb 不等式。
- 停止：三个阶数及张量接口核验后停止本批，不增加 kernel 网格。
- 验收：独立精确 D_m 一致；时间矩不等式给完整证明；区分可用改进、未闭合条件及预算障碍。输出 report.md、check_time_moments.py、check_time_moments_results.json。
- 输入以脚本内绝对路径读取并记录 SHA-256；不把采样值当全域界，不晋级 C6。
