# audit_framecover078：077首块独立16角点审计

**仅冻结首块通过，仍为partial。** 唯一源066 final SHA `42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b`；首块checkpoint `687edbd77268118b31f294ed397319470bd506b26ce65e969fcfead9a015c1b6`。075理论及同alpha1的029资格原件均先验哈希后完整阅读；未读活动future，不从非空queue推全域g。

实际PID **67943**，exit0、OS已查无，总wall **6.026660791秒**，RSS216432640 bytes；600秒总限包含加载和落盘，540 checkpoint/590信号兜底未触发。零新H/gradient/kernel/Green/log/exp/Gamma/Arb/MATLAB调用；不导入作者模块，不运行作者重建函数，无提交推送，062/069全部12文件锁在运行前后不变。

独立重建当前queue后代数，按 `(path长度,−当前queue后代数,root,path)` 核准确4096选择前缀；1015项因较早通过祖先而跳过，3081项恰与账本测试行对应。所有源split索引唯一且测试时活跃；完整旧失败证书SHA及Delta/物理gradient/anchor原样核回源缓存。

物理box独立由root坐标及path每隔四位的二进制子串构造，不复用作者box/split实现。缓存H区间严格正，五元区间有序，所有导数/H端点交叉不等式通过。075的区间线性项`Σmin(gL·v,gU·v)`本身全局凹，减`256Σv²`仍凹，因此16角点控制整个闭盒，包括S1/S2、零/重投影节点与符号边界；anchor盒内不是理论必需条件，记录的inside标志与实际几何相符。

每盒只检查`{lo−anchor,hi−anchor}⁴`，准确清位移分母V与Delta/gradient分母C，整数内式为

\[
\big[\Delta_iV^2+V\sum_j\min(g_{j,L,i}z_j,g_{j,U,i}z_j)-256C\sum_jz_j^2\big]/(CV^2).
\]

全部**49296角点**与保存minimum/worst_corner_bits严格相等，判定符号逐条相同；**2594父严格通过、487父未证**。最小accepted准确余量保存在result；没有将未证当作H反例。原数值球与有向log的认证及075解析证明沿既有锁继承，本轮不重新计算超越函数。

退休source passed850、splits10036、queue11780逐source_index/node/unique_new_parent与独立覆盖集合完全相等；新passed96388/splits46012/queue1984。源及新前沿均对52360个root逐根转dyadic整数区间准确铺满单位区间，prefix-free且Kraft=1；活split集合恰为叶path全部真前缀。新前沿98372叶=52360+46012。cache_ref独立规范SHA `f3b174e2bab3fcaed2c291090b378ddac7414393a6fb7a7bceffedca0ffaccac`，共178722项。

独立逻辑前沿SHA为`91fd1e8a88e76785bdbc81263bb00fd6c9895774bdae6def62e7a98301c78406`，与注册parent链及checkpoint/summary链接一致；尚有1109个未测eligible源父。本报告不认证第二块或最终物化，后续另登记`final/`。

本首块仅四件，全部输入及保护锁见registration/result：

|文件|SHA-256|
|---|---|
|registration.json|`ebd91f62a87961b53c366505b5a15f8ae2eef1b594454975b672acbddb2aa0b9`|
|check.py|`7dd1a3d02a2c3ff0a0b01eeacfff9b7de30ca29f1f528706fabfb46adf22f02e`|
|result.json|`740a07517817806936bc65dc63eb179f59dec20ee5fe0e32d9b50818b1af8573`|

复现：在本clone运行`python3 -B cube-dim-5-bounds/audit_framecover078/check.py --run`。主任务保留终审；非空queue不构成全域g或Cn突破。
