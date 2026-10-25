# audit_I026：022完整I与024全权熵最优预算

运行句柄`audit_I026_fixed_001`，一次固定批次120秒软预算。锁022 manifest `668f4b36fa225c865da727a1bf92f211eccfb06a35aad0ebc44cfa7071531602`、024 manifest `96f65eaa792afc290d3b83aa2284f64170386de045e54f04070b96335f93b764`。不重跑021或作者120对Gram生产，不做权重数值优化，不提交。

022不同准确路径：只算pair(0,1)、(0,3)两个signed-permutation轨道代表。用rank3有理nullspace分第一线类，再对每个第一线类用rank2 nullspace/Gram面积求第二类和末完成系数，不用作者2阶adj/Schur整数表。每个pair代表一项5x5 Leibniz系数核验。按完整群生成元flag轨道恢复4560组、五概率type、sumC60W、每scope全16。核CB1/(10*3*2)=1/60及五因子独立。

024解析重证一般q的熵AMGM和普通Holder、全4560权Gibbs恒等式、boundary q=0、qopt正性及唯一性；k0/k1无额外商类，只重复k2必须单独证明。不得把最优上界预算当作实际J/I下界。

固定Gamma复核用未移位Mellin m(b)=Gamma(1/5)/[(3/4)Gamma(3/20)]、Gamma(7/4)=(3/4)Gamma(3/4)。022走I^7闭式，024走I^7含Z^(-3/2)和J^10，Arb512bit后Fraction整数根110位包围，区别作者I^14/J^20与100位代数根。Z由独立重建的type逐项求和并核短代数式；与锁结果准确区间比较。共享Arb底层Gamma明确，不能伪称跨库认证。

同G、全space/time、退化及一般L1/全尺度transfer单独核。入口--source-root、--model-root可携，默认只锁输入，--run执行该固定句柄。白名单仅登记/代码/结果/报告。
