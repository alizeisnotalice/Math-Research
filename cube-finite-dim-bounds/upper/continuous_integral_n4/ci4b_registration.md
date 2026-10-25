# CI4B：全和有理式的仿射多项式积分上界 pilot 登记

2026-10-25。主明确授权增量接续。只写 continuous_integral_n4；不搜 g/tau，不提交、不扩招。旧 CI4 保留失败状态，不把它的区间当新证书。

机制固定：对作者三个物理排序 case 的准确有理式 H=G*pi^3=N17/D20，在 x1=a=1 的 tetrahedron 顶点 barycentric v 中取同阶 N20=a^3N17 与 D20。拟以四个有理 P_i 定义 P(v)=sum_i P_i v_i，准确核全部 PD20−N20*sum(v) 的 degree21 系数非负（及 D>0 的内部来源）；必要时只允许统一向上平移 P_i 的准确 residual 修正。浮点 LP 只找候选，不是证据。若 D 的 degree21 support 为0但 N degree21 support正，显式 unresolved，不删点。

积分规则固定：H<=P、P_i>=0，故 H^(4/3)<=sum_i v_i P_i^(4/3)；单 tetra 上界 = volume/(4*pi^4)*sum_i P_i^(4/3)，正根和 pi 用 Arb。可选同法 lower 和 Jensen 只用于诊断，非 g 搜索。先实际运行10代表 tetra：三个case根部/外壳(I根1个，II/III每3个外壳共6个)，加上述I根、II首外壳、III首外壳分别按准确最长边二分8次后的首叶，共3个。全批须等作者锁输入hash和pilot可用性反馈，不提前运行。

固定全批预案（不因数值机械扩网格）：case I 根 longest-edge bisection depth8，256 cells。II/III各以真实corner D=(1,1,1,1)做28层homothetic shells，outer scale 2^-k、inner 2^-(k+1)，k=0..27；每shell准确拆3tets，深度 max(0,8−2k)。共2184 cells，两case加I为2440 cells。最长边按准确有理平方长选、相等按lex索引；每步二等分。删除t<2^-28的corner tetra，包含于半径epsilon=2^-27球，单正octant总遗漏预算付E4(epsilon)/8（两case子帽并集不重复付完整球）。预算有限性与准确volume覆盖另核。

终止：作者表示未锁则不全批；pilot support-zero或系数证书失败准确列缺口，不放行；全批大于4096立即停；不得把点采样、LP状态、浮点quadrature差、空白/未解cell当证明。若全批I上界不达19/100000，只允许根据真实错误原因另行明确登记一次机制修订，不机械加密。完整表示独审及最终成本由主与独审负责，未经全链通过保持U4=16。
