# audit_green036：仅018选中点的纯有理完整Green复核

句柄`audit_green036_fixed_001`，--run先打印PID。一次120秒软预算，不重跑017/032、旧点/anchors、Sobol、局部搜索、Fejer或旧空间grid，不提交。

锁018 manifest `e387fd41e3bebc8ba94427200fb13fc6693f15c92118796b34d704d738d018bc`、017表`8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9`。只算selected record49的点(1,104115/104729,103354/104729,71094/104729,6327/104729)。

不同路径：普通distinct投影k=n.z，将3504个DD4log直接转280个准确Fraction barycentric log系数，核每锥常数及全局Σcoeff=0，按abs(k)合并。利用总系数0删一个最大参考log，对139个相对log作2^m缩放至[1,2)，加共同ln2，所以最多140个正atanh级数。固定N96，尾上界2u^(193)/(193(1-u²))；每项外向量化2^-384并以公共整数分母求和。pi用80项Machin。全部pureFraction，无Arb/SymPy或DD递推。

保存280个A_v及合并C_k证书，G区间精度足以与018完整球/dyadic交叉。只回读alpha1_gate.py/json/proof，核同G、连续性/essential障碍、I_face对象及准确四次幂残差与已有根夹端点，不重复100位根搜索或旧Gamma预算。

Gate锁：script a6deb402cdd9ea8d1e212486364d7cc2a4f1573c8262afb1df576e1cdacc7293；json2e8601968b6ed99582ab9885efda11dcd892e1ffb70ae0995350a868aaaf3142；proof e348a5480c6570da3a9f9481b567e017ca481213fd97e5b393823d654f7df04b；gate registration d26cb36f76b15880dd9b27e673b997f698355f6691bf258754e76fa4772ec9c9。入口--source-root/--model-root可携，默认只锁hash。精简白名单登记、代码、系数、结果、报告。候选点不是全域g下界或全局极小点，预算不相容不是C5下界。

数值执行前变更记录：首次预锁拒绝发生于任何PID/系数/级数计算之前，因proof由e1d9c742…更新为上列e348a548…（新增局部cone DCT连续性证明）。主代理确认最终四源冻结，并说明registration临时附录已移出后恢复原字节，当前hash与JSON一致。本审仅更新尚未执行的输入锁；这不是数学计算失败，不修改源文件。
