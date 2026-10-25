# 11K完成报告：四维purebody全cube g下界

2026-10-25。**完整证明G_4≥21/25000=.00084于整个单位cube**，不预设最小点。只完成本g分支；没有I上界、tau/cost搜索或新的弱型上界登记。独立增量审查已完整读证并通过，原independent_audit最终验收后置。核心checker/data/proof现全部冻结。

## 结果与证明链

完整证明u11_minimum_proof.md SHA5546ee2fd0c8dddab3b88d0db79143a846a81711f13c6967b4fa922fa9c4757e。固定原4D纯body真实frame、Lebesgue、U11正常数。全1152 cone Green式先按ray-star准确合并；A/C/E/F四个非顶点orbit的numerators恒0，B/D两个48term stars给全40顶点。按±配对与signed permutation协变，缩为12个B-line和8个D-line贡献。

在物理sorted chamber a≥b≥c≥d≥0，B-line符号全正，D-line仅a−b−c±d可变，恰分3个完整cone。20项在每个case准确合成N17/D20。D的20个线性因子在该cone非负参数下方向一致、内点strict正。取k=(21/25000)(22/7)^3，exact展开

[a³N17−kD20](Ru)

得到1739/1702/1702共5143个非零有理系数，全部非负。Arb128严格核pi<22/7，于max坐标a≤1给G≥.00084。这个证书对全部实非负参数成立，不是点测试。

generic内恒等式后，靠026已证真实G在线外连续、frame线上为∞，覆盖case/sort边界、D=0及任意zero/repeated节点。只用lsc不足，证明使用真连续。既有单点17/(640pi³)没有被用作全域min假设。

## 有限预算与实际运行

首先按登记只读共享CI4的351准确sign polytopes（≤2048），保留同vertex affine cross sums；Arb128初诊248单cone负pole、103有限lower但目标不足，0certificate。一个body-corner cap的显式one-large-jump lower>21亦记录，但最终全局多项式证明完全不依赖该分区/cap。初诊摘要u11_minimum_initial_summary.json，完整失败cell dump只在本地留存，不归档重复材料。

唯一表示修订是同ray residue star再全和相消：两个非零stars common degree26，全20项common degree60，固定3case消因子得N17/D20；没有机械扩cell、加点、调tau、搜索g或碰I。四null-star是该全和恒等式的必要准确检查，不是新候选。

Python实际用Fraction、FLINT fmpq_mpoly及Arb。最终一键checker u11_minimum_certificate_check.py SHA35489f9fa4c52a09154105ad949a8d3410f61bf226e8a24dd8d57359f0a3e060重建同一固定代数并PASS，约1秒。中途两次仅fmpz→JSON int序列化修补，所有数学断言已通过；随后去纯timing字段使认证主体可稳定锁定，没有新增科学试验。旧dd71...为preliminary历史，不再是锁源。

全global式u11_minimum_global_results.json SHA43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f；完整稀疏P系数证书u11_minimum_certificate_results.json SHA23a8ce3a07c6c5ad94b29dd0b57aba4dd4e3f1fd4a9098374ec700b73b0feb3d。

主串行MATLAB实际exit0/session55870，脚本SHA11f71ee1795d1e68cd5d6a6dfb096e1567100f2db17b200c22db3858c384458e；精确字符串符号核5143有理系数及degree20，并以独立矩阵代换核每case20个分母线性因子，共60个。此为有限交叉，不冒称MATLAB重建了全部FLINT恒等式。完整stdout与log hash在u11_minimum_run_receipt.json；upU11Min前缀，无clear/path。

## 独立增量审查及剩余审核

continuous_minimum代理独立逐项约分构造六star（没有调用作者common-denominator生产器），完整cross-multiply核3global恒等式；独立检查20sign在整个cone恒定、40vertex群等变性、各cone逆参数、D内部正、全部P系数逐项相同。有限45个support边界检查亦通过，但边界下界的证明来自其另给的真实G连续/line∞概率论证，不来自样本。

其完整report为upper/continuous_certificate/n4_review_report.md，SHA a478179b566a9a424db51fe72f1dd657d501262764885f07b45658838a6c603f，已读并核实际hash。它明确保留共享1152锥/Fourier Green输入的范围，不替原independent_audit最终验收。主已接收，当前不得把本g独自宣称U4<16。

## 给I4的锁接口与一般L1/全尺度边界

u11_minimum_global_interface.py SHA6c8d1a5b4912103df8105f7c7662ac2265b0978f91de00651be03acfb1790bae已交I4代理。变量固定a≥b≥c≥d；caseI a≥b+c+d，caseII b+c≥a+d，caseIII其余。`homogeneous_pair`给同次degree20的a³N、D；在face a=1等于所需C=pi³G。可做正系数ratio/affine majorant的准确验证，但未由本分支运行I。D=0不得当0，调用allnodes或连续延拓；接口的无限upper明确None。

固定Lebesgue实际用于真实Fourier正常数、cube与ray尺度/体积、Green定义和026的比较。g证书本身不给一般L1弱型估计；主还需完整I/核质量等输入，再依026/U09的一般L1及共同全部正r传递。不改变背景，不把cube改独立矩形尺度。U3分支不由本批处理，U4账本不由本代理更新。

|类别|本批状态|
|---|---|
|严格/准确证书|作者全cube g≥21/25000证书PASS；独立增量全文审查及准确重构PASS，原独审最终验收待续|
|实际工具|Arb128、FLINT准确代数；主串行MATLAB coefficient/denominator交叉PASS|
|未完成|本分支不提供I上界及完整弱型组装；不替原独审验收|
|未宣称|无单独U4<16/U3<8；无tau或其他frame族最优性结论|

## 白名单与停止

17文件：inputs、preregistration、dd API、initial_check、initial_summary、star_revision、star_results、global_revision、global_results、global_interface、certificate_check、certificate_results、proof、matlab、run_receipt、本report、manifest，均u11_minimum_前缀。旧U11模块/表示已另归档，仅hash依赖不重复复制；初始完整cell dump及__pycache__不归档。

本g有限批完成，不再修改核心认证主体，不运行I或新搜索。后续仅在原独审正式回报时更新报告/manifest状态。
