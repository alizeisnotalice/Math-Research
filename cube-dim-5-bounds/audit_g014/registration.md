# audit_g014：Astra007 K3唯一全域候选审计

锁007 manifest b8a7b09d6327428facdc8e6bc7524a87e226d42f2bf859e41240a42392672b5a和006 manifest 398f33d8565a30a68b9d4afbb5d29e1047a2a3a37783ab91a56127437bece136。只审K3 ell、eta=1/8、a=679/1000、g>=5563/10^8，不新增参数，不重跑147258矩或g013基线。

一次固定有限复核（120秒软预算）：用Fraction重建K3 Taylor系数与共同T下端，核T上下端11次幂、tau区间、与旧线性项的正性差因子；读取已锁006corner和007空间矩，准确核32个角点依赖和空间维数；重新组装唯一127项空间多项式，使用公共整数分母直接tensor公式恢复全部2401 Bernstein，核最低值、hash及c下端乘积与目标。

解析另外核四阶Taylor余项在tau两侧方向、d=0连续端点、全d非负、006仅压缩det轨道的合法性、007完整d³/d⁴展开、degree27及5位键/GMP、Z合法性与G齐次face扩cube。空间矩共享006/007和已审round012后端；不把角点匹配称第二独立矩生产。

入口--project-root、--backend-root可携；默认只锁hash，--run才做此批。仅写audit_g014，不提交，不修改旧审计目录。
