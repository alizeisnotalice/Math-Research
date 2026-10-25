# audit_vertex032：017实际真顶点扇独立资格

句柄`audit_vertex032_fixed_001`，--run先打印PID。锁017 manifest `9084775b6975369b37c6e20b19778e21be54f75b922f1874d1411fe94037e640`及表`8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9`。一次120秒软预算，不重演65536 contacts、面格或pulling生成，不重算旧26880密度、不跑Fejer/grid，不提交。

独立覆盖逻辑：K=conv(280表点)，逐facet支持/接触rank证为K真facet；全部simplex以det外向定向、三面链闭合。正向闭4链在连通凸边界上重数常数（有限多面形细分后相邻cell边界差证明），选一个5接触点/1simplex的facet固定重数1。因此完整边界覆盖且相对内部不重叠；不以总体积相同作唯一依据。再逐facet用质心符号直接证a=Σsign(v·centroid)v，及所有顶点psi=1，推出K=Ppsi和facet表完备。各表点的active-normal rank5另证全280均极点。

有限准确核：仅读锁表，Fraction秩/120项Leibniz det（区别作者Bareiss）、8760定向tetra闭合、全部支持/法向符号身份、原点内点、simplicial-facet度数见证。旧005三轨道及11类源只作来源交叉，读取旧7基锥常数核体积，不展开26880。

两个指定anchors用Gaussian-integer分母乘积+Fraction精确实虚和；pi用80项Machin有理夹，完全不调用Arb。核精确虚部0、正实值、95位外向有理夹包含于保存新旧256bit实球，读取差值/虚部球相容。只alpha1 Cauchy资格，不迁移alpha3/2 g/I。中心p1(0)只核归一化系数，不另作旧数值基线。

范围：本批独立验收积分所需边界/锥全域覆盖和内部不交；作者完整面格f-vector/59040 pulling限制不重演。若更强face-to-face组合相容性需原生产面格，则明确为锁来源。精简白名单登记/代码/结果/报告。默认只锁hash，--run执行固定批；旧锁和主clone只读。
