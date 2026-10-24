# 第42轮独立静态源码审查

GO：已全文审quadrants_once.py、PLAN、wrapper及piecewise引理，未运行主程序或核验器。

- source SHA：`fd16ccf0309a4ae4627049ff730b2bb233be3948013d25cd514824febcc107b2`。
- PLAN SHA：`90e22edf1bf76daf953062531997e1583fb0631f22947cf7805edba34503ce71`。
- wrapper SHA：`4e8c1c6108816610ccf37a8c5edba4d379afb2909e220e9226098d707d9a1e47`。

固定node663模型/cover/B1 SHA及准确物理/局部几何正确；Muv由保存Hab端点上界核对，ES保持不变。
先u后v的de Casteljau方向与i=u/j=v索引一致；四片(--),(+-),(-+),(++)几何及双线性升阶正确。
同号修正0、异号Muv*delta_u*delta_v；64行准确组装，旧c/s/slack限制关系与gain恒等式逐项校核，无float接受。
新proof pieces不追加旧tree节点或修改旧depth4记录。全64行strict>0才accepted；min_gain允许0。
仅标准库Fraction，0Arb/Hint/积分/AS/Psi数值模型/点/LP；60秒硬限、4片64新行，无自适应第五片或重试。
异常保存在status=failed；wrapper exit0本身不代表接受，终态还须读status/accepted及完整计数。
独立verify_saved.py已落并仅AST解析，SHA `a827a255e822aaec9724cd0415cd3a3d06f97969e5c7e6a5db7b34f89527cc9f`，待主审静态审后唯一运行。
该verifier走Bernstein→幂基→仿射代入→Bernstein，不复用de Casteljau；核5组限制、双线性、64行和四片物理slab覆盖。
两份代数程序共享保存端点/解析引理，但算法不同；核验不称第二套积分证书。单盒成功不是新C2上界。
