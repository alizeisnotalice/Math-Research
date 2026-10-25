# F048 band Hessian 独立审计：解析/源码及两全群分量通过

接受 C048 对同一正band 的Hessian接口；本独审实际完成8项新的零/重/互异节点scalar参考，以及x2x2、x5x6两个完整全群参考。**全部21分量独立全群参考仍为false**。其它19项不能称为逐项独立全群认证。

固定对象为旧oldc、a=1/250、T=6/125、原7a57caad…body32 fan。没有新点、band或模型搜索，没有全Green Hessian替换、盒余项或新Cn界。

## 数学和源码裁决

Cauchy contour证明mixed DD7与diagonal 2DD7的有标签节点公式，数值相同但不同标签不误乘2。normalized jets七阶、全零偶/奇jet、重复base Hermite表与互异base快式均核对。native对不能分离的准确不同节点整体失败，不drop term，也不引入sign(0)。

物理链的两Jacobian因子、off-node双项、physical diagonal、21字段索引、sign/permutation映射均与J^T Hnode J一致。global负号配对权2合法；独立参考显式全64signs作为交叉。输入int64/int128中间界与最终门控、source/binary/header/runtime ABI锁、完整dyadic端点均核对。格式parser不是任意fan合法性的替代；本结论只适用于已认证且锁hash的fan及matching positive band。

完整证明和逐项源码解释见proof_review.md。

## 两个真实完整独立component参考

只运行自己的reference.cpp。参考直接对物理方向上的f(b_i)∏(b_i−b_j)^−1求二阶pole-product系数，不构造node Hessian，不用作者J^T Hnode J实现。f使用两个独立atan、f'/f''直接有理式，Arb768。每分量显式720×64×21=967680 computed/logical项，0zero base、0repeated base，全部保留；不是只检查选定符号或少数cones。

两分量同时累积仅共享节点及权重；这是两个完整component的真实计算工作，不能报为0完整核工作。作者binary/anchors/full/writer均未执行。独立参考按登记遇zero/repeated或分离失败就停，未实现这些配置的full fallback；实际oldc没有触发。

|component|独立参考导航值|准确比较|
|---|---:|---|
|x2x2|−1.2558319355365172e−5|完整Arb768区间被native512包含|
|x5x6|1.400156560026021e−6|完整Arb768区间被native512包含|

比较使用全部integer mantissa/exponent转成准确Fraction，不使用表内导航小数。reference_results.json保存原始参考端点，comparison_results.json同时保留native与独立完整端点。

独立参考驱动PID49909；同一编译后批进程墙钟20.993380秒。原生子进程PID未单独保存，不伪报。编译、rpath patch命令及二进制SHA在reference_execution.json；只patch本目录独立scratch binary，未改系统或作者文件。构建0失败/0修复。

## 八项新的独立scalar

采用confluent residue-series系数，不用作者DD表或Hermite线性系统。准确有理pole-product展开配符号有理atan导数，Arb1024仅负责超越值。8项分别覆盖：全零diagonal、全零不同标签mixed、混合零mixed、四重零diagonal、互异diagonal、互异mixed、六重非零mixed、多重复簇mixed。全部完整区间被保存native512区间包含。

scalar驱动PID50036，0.055832秒。没有重新运行16个作者anchors或增加第9项。scalar_results.json保留全部准确留数权重、端点与包含判定。

## 已保存接口和范围

旧band、六梯度、p_a、p_T共9字段与C016完整dyadic逐字相同；21个新字段均有限且端点有序，原Euler含0。旧接口一致不能代替Hessian全群参考；因此full_21_independent_reference始终false。

输入锁至完成数学比较159.446009秒，小于20分钟。0数学失败或重试、0作者native、0新空间点、0模型搜索/grid/MATLAB/扩招/push。仅两full component参考和八scalar；没有作者full21重跑、grad(p)、二阶Euler或盒余项。保存的全部旧输入末次hash不变。

## 源锁

- C048 manifest: a0b9268c9046ced80443fed30610ee4b6394ff5b49c4b34aa116cea067e9a2b7。
- C048 source: 850c9e7df479d8be28e0d18defb31fc31d98b75a83b239fe1329a44ae558c653。
- C048 binary（只读验hash、未执行）: 4bf7052315d33e04af73c36a90f34ce011ef608f425a3df7cd7dc446823e6a4c。
- C048 theorem: d7cdff5274cfe455e7a6e11e84ee17150770f06a674909e93023cbd9a2fea287。

本目录manifest封存独立代码、证明、端点、编译和运行收据；可重建scratch headers/binary不进入归档。父代理保留最终证明及使用放行责任，本批停止。
