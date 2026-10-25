# audit_I023：021完整I独立有限审计

运行句柄`audit_I023_fixed_001`，一次固定批次120秒软预算，不重跑生产560个Gram循环，不纳022或type权重优化，不重跑旧g，不提交。

锁021 manifest `16eeeb3e07fc15813c7ece29d415ff79c219fa116b6a5add2f4926860c1ad42b`，同模型一般L1/全尺度transfer `3d2acf640e1be0ee11d10fd6e29719ddc4291117975407e0175c8142573d5465`。

独立组合/代数路径：只计算J=(0,1,2)、(0,1,6)、(0,3,5)三个轨道代表的有理nullspace N（不用3阶adj）；以D/det(NN')乘二维投影坐标面积平方求完成系数。每个代表额外核一个5x5整数Leibniz行列式。用signed-permutation生成元的轨道，覆盖全部三子集及线类flag，重建2240组、5类概率、theta及16坐标scope rho；只代表做商几何，不重做560次商生产。

解析覆盖精确CB因子1/20、rank2商Gram常数、所有零商/共线/退化项、五独立stable因子、产品概率scope分数Holder、Mellin/Gamma与共同G的完整face/全时间接口。

常数只算本批：5类R乘积使用闭代数式，负矩用未移位Gamma(s/beta)/(beta Gamma(s))；计算I^7闭式而非作者直接I分数幂，用Arb512位包围后以Fraction整数7次幂/整数根给95位有理I夹，比较锁021球及有理cap。共享Arb库仅Gamma/代数函数的底层误差包围，不伪称跨库完全独立。

入口--source-root、--model-root可携；默认只锁输入，--run才执行。精简白名单为登记、独立代码、结果、报告。
