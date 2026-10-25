# 第5轮：局部接触面变分的有限预登记

2026-10-25；唯一责任目录 `lower/variation_geometry/`。本 clone 为第4轮基线；第5轮上游只读。沿用主任务 Lebesgue、中心轴平行共同半径、完整质量和联合体积定义。本人不优化 n3/n4、不重复全局盒枚举、不改总账、不提交推送。

目标：证明有限条法向半线和有限新生薄片区间的带符号长度公式，含并列端点、跨面遮蔽、空隙和非对称薄片；给有限轴平行盒族提升至4D的一阶充分条件及缺条件的准确反例。

运行前固定检查：至多512个有理配置，分别用精确区间并算法与独立端点分割积分比较；使用两组固定有限扰动，ε=1/8,1/16,1/32,1/64，一组内部接触/多薄片填缝，一组四维多坐标耦合。所有 Python 计算使用 Fraction。另提供前缀 vg5 的 MATLAB 脚本，请主代理串行执行；数值交叉核不替代证明。完成引理、接口、准确结果与反例即停止首批，不扩展搜索。

输入：主 clone `rounds/005-expansion.md`、`results/round005_expansion_snapshot.json`、`rounds/005-registration.md`、`STATE.md`、`results/bounds_ledger.json`；n4 clone `round005_preregister.md`、`round005_exposure.py`、`round005_exposure_results.json`。读取时 exposure.py SHA256=c524459ce94baf988add6a4f1d68d3a5f7801558685a9c426f616fc36231e16d，作者已确认代表当前 born/continuing 语义；旧 snapshot 是历史读取值。最新已审 L3 是007b，005c当时仍待007c独审，本工作不比较L3。

共享算术披露：与 n4 共享半线/薄片的建模语义和上游输入；本人提供局部公式和不同实现自检，不称完全独立认证，最终候选仍交 independent_audit。

停止状态分四类汇报：严格局部证明、准确有限校验、数值交叉检验、未闭合/适用缺口。首批不会宣称新 L4 或排除整个卫星族。

首批完成后主代理复核追加：真实N₀=0纯卫星半径O(ε^(1/4))，不属于原O(ε)端点定理，必须补严格切锥结算及接口提示。仅增加五个固定Fraction身份检查，不新增方向/相位搜索；主代理要求将新追加引理交independent_audit。MATLAB已由主代理串行实际执行，两组最大误差均0。
