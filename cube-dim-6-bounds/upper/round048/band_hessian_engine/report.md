# C048：含零／重节点的正band完整 Hessian 引擎

本批完成新实现的有限验证。新 Arb512 binary 在唯一登记的旧 C016 oldc 点及同一带宽`a=1/250,T=6/125`，输出完整对称6×6 Hessian 的21个独立字段；全部有限。旧band、全部六梯度及p_a/p_T共9字段的完整dyadic端点与C016**逐字完全相同**，原Euler区间含0。完整Hessian没有独立全群参考，本批不能称为全独审。

## 理论与实现

`theorem.md`自含证明：对六个有标签节点，mixed二阶为`DD7(nodes+bi+bj)`，diagonal为`2DD7(nodes+bi+bi)`。数值相同的不同标签仍按mixed处理，不能误加2。函数是同一个正band `atan(b/a)−atan(b/T)`，在零点解析；normalized jets扩展到七阶，所有重复／全零都保留。

重复base节点使用八项Hermite表，直接计算copied-node DD7。互异base节点使用已推导的等价快速重心二阶式；scalar参考使用另一套degree7 Hermite多项式系数法，独立于两条native分支。没有有限差分。

物理二阶严格按`Jᵀ Hnode J`映射，J由真实ray/q、原sign及permutation组成；node映射是线性的，没有额外二阶链项。稀疏整数ray乘积在新目录预计算，coarea cone权重和`(2π)^−6`归一化沿用C016。21个字段按`(0,0),(0,1),...,(0,5),(1,1),...,(5,5)`排列，镜像得到完整矩阵。没有canonical负反射；全局负号配对的两个Jacobian符号相消，权2仍合法。

精确输入为规范int64有理数；点积和乘积先用int128并门控int64范围，ray乘积范围也检查。512位无法分离两个不同准确节点时整体失败；没有擅拒绝零节点后drop term。parser沿用C016严格token和stream检查。ABI门控ulong/fmpz8bytes、arb_struct48bytes、FLINT3.6.0；数学计算不用binary64或longdouble替代Arb，时钟为普通浮点。

## 先验收再执行

预登记input SHA字段为
`316b359ae0f8ddb5cafe09a46a5d6edb25af18edd16813155fe864a978b8cd4d`。
C016冻结manifest、源码、native原件、fan及旧header/runtime哈希先核对，旧文件只读。新source SHA：
`850c9e7df479d8be28e0d18defb31fc31d98b75a83b239fe1329a44ae558c653`。
新binary SHA：
`4bf7052315d33e04af73c36a90f34ce011ef608f425a3df7cd7dc446823e6a4c`。

新代码编译一次、0修复。观察到的源码生成＋构建tool墙钟0.738015209秒；纯compiler内部时长没有单独打点，因此不伪报更细数值。完整compile／rpath patch命令和依赖哈希在 `build_record.json`；仅patch新scratch binary，没有修改系统路径或安装环境。

实际16个scalar二阶锚均通过：全零、混合零、四重零、互异正负、六重同一非零节点、两个重复簇；包括diagonal和mixed。独立参考先用准确Fraction Gaussian elimination解八阶Hermite插值的leading-coefficient权重，再用准确复有理pole幂产生normalized atan jets，只有函数atan值用Arb1024。每个参考区间均被native512区间包含；准确权重和完整dyadic见 `scalar_validation.json`。独立参考内部用时0.00795566602027975秒，anchors＋Python参考tool观察墙钟0.318136459秒。

理论／预检报告 `precheck_report.md` 在完整native之前交主阅读，保留其当时full执行数0的历史状态；主明确验收DD7倍数、distinct快式、JᵀHJ和negativepair后release唯一full。预检形成至release约71.94419860839844秒；目录建立至实际full启动准备／验收总时218.2407455444336秒，二者都包含其间准备工作，不冒称纯CPU计算。

## 唯一实际完整群

固定准确点仍是

    x=(1000000007,998836847,892738545,807464540,467291939,96968636)/1000000007.

同一旧fan SHA为
`7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812`。
完整967680逻辑项、483840全局负号配对求值；该oldc本身的zero-node及confluent项均0，与C016一致。零／重支持的实算证据是上述scalar锚，不能冒称本完整点检验了群级zero/repeat。

唯一完整运行内核12.6219秒，进程墙钟12.635845624987269秒，退出0。原始完整输出在 `native_hessian.json`，所有Hessian字段均为完整lower/upper integer mantissa及共同binary exponent，未从短radius展示重造。旧9字段准确端点完全相同，原Euler仍含0。没有第二次full、额外grad(p)、二阶Euler、新空间点或其他band。

这次完整运行只验证实现的一致性和实际计算结果。没有完整Hessian的独立全群算法参考；数学chain靠定理和源代码核对，不能靠旧值／梯度一致推得Hessian全部正确。后续若用于局部证书，必须使用matching band，并另证盒上的严格余项；固定点Hessian不能直接成为全域oracle，也不能偷换全G梯度或Hessian。

## 封存与复核

主保留证明及后续独审责任；本批等待未参与者审计。`run_full_once.py`要求独立release记录并拒绝第二次full，任意超时、非有限值、ABI或旧端点不交均整体停止，不重试native。此次没有上述失败。

只读复核：

```sh
/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B /Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round048/band_hessian_engine/verify_saved.py
```

它核对白名单／来源hash、16份保存Hermite权重与端点包含、旧9字段、21个Hessian端点及原Euler，不重新计算Arb锚或full，不写文件。build可由 `build_engine.py`在本目录scratch重建，runtime只读复用C010；scratch binary及解压headers不进入归档。

无C16/旧header/runtime修改，无MATLAB、扩招、grid、MC、安装、stage/commit/push。该工具本身没有Cn改善。
