# audit_parallel137：1153个跨锚点通过证书与完整前沿的独立审计

**独立局部组件资格通过。** 1153张新通过证书的18448个实际闭角点全部以独立整数路径准确核回；89426条日志的hash链/首PASS/skip语义、10549个选中来源、4274个残余及完整52360根前沿均通过。最终105635 passed、4274 pending、57549 splits；4274待证仍在，完整g110未成立，全局S记0，I/C5未更新。主保留最终接纳责任。

没有重跑71812个生产corner测试或失败角点，没有H、物理梯度、log、kernel、Arb/exp、新split或KDTree候选选择，没有Astra、下界、共享MATLAB、提交或推送。原129/130/134及作者目录均只读。本目录是独立clone中的唯一新增写入范围。

## 1. 来源和有限边界

已阅读全文137 `result_report.md`、`result_manifest.json`及131/137实际接口与worker协议。实际137结果manifest SHA为
`16d98997cd04258b28bd997d5489c6a7ad4b52a359f56403e9943017bd11fd2e`。
144条完整文件SHA锁从137结果、137准备、131准备及本审计必要旧源汇合而来，去重源字节共1,980,914,733；逐项绝对路径和SHA保存在两个冻结计划中，运行前后全部不变。大型不变源cache/树/生产日志保留原处，不重复打包。

数学来源仍是同一五维Lebesgue、alpha1、16条未单位化body方向、共同Cauchy时间、逆Fourier因子(2π)^−5模型；目标11/100000。已有075/099全空间二次下界、129锚点时间Jensen、096/098完整g105 H64叶、旧105/107/123数值与130已审局部组件作为冻结资格继承，不在此重新求值。

新算术严格限于1153张日志中的首PASS证书。71812个evaluated日志仅检查保存margin的符号、元数据与来源，不独立重算其中失败的角点。这个边界足以确认1153个正证书，并将其余4274保守保留pending；不声称4274个真实H/G不达标。

## 2. 与130分开的整数算术

新 `exact137.py` 不导入130、131或137作者任何计算函数，不使用Fraction。对每张证书，先把准确point和实际box坐标统一到整数分母W。再一次性对Delta、两端gradient乘位移及q乘位移平方的全部可能项，取一个公共正分母T。得到整数化Delta、四对gradient系数和一个quadratic系数；整张证书的全部16角点共用W和T。

每个角点只计算
\[
 N=\Delta_{\rm int}+\sum_j\min(a_ju_j,b_ju_j)-q_{\rm int}\sum_ju_j^2.
\]
实际全部16个整数分子都严格正。最后仅约分最差N/T，并与生产保存的准确最差margin、worst-corner bits、anchor-inside-child字段核回。平手按bits字典序。没有逐角点切换source，一个source必须覆盖完整16角点。此路径与130逐项有理数加乘/逐角点分母清算分开。

1153次成功证书、18448角点全部准确匹配。逐证结果原JSONL SHA为
`89d04ca308d1e0ff195e988c0e2e027eed3e092834b9bd6d7e86f6bbf1c63dd7`。

## 3. 原始来源、point与H64闭祖叶

完整111331条source_pool记录逐条核canonical core ID hash、唯一source ID和稳定index。固定89426个pair实际只涉及10549个source；这10549个逐alias独立绑定到原105 passed/failed、107 H128或123 H128的原target110证书及raw record：

- 原certificate SHA和raw record SHA准确相等；同一完整point、epsilon、H interval及原penalty绑定，保存Delta和四个log-gradient端点与原证书准确相等，不重新求log或梯度。
- H64仍为Q256。H128新q源独立检查实际原098祖叶EPS64、外层g105、原passed/正margin、完整原叶SHA与已主审终态SHA、准确树祖先前缀及真实point在闭叶中；不以G下界代替H64下界。
- q128用正分母整数身份检查：若H128上端为Un/Ud，则q必须满足
  \(q=(1024\cdot200000Un-768\cdot21Ud)/(200000Un)\)，并核256≤q<1024。旧194个未找到H64祖叶的H128来源保留1024，实际核其原105来源及树路径无H64祖叶。
- point只要求在合法H64祖叶中，不要求在接收child中。Delta允许负；H下端只要求正，不新增anchor-H高于target要求。

源数值真实性继续依赖冻结的原数值独审。本轮新增的是独立原对象/point/曲率绑定及新child成功corner核验，不把字节hash当作重新计算的热核证据。

## 4. 89426条日志、first-PASS与skip

15个whole-child分片逐一重放日志元数据：每条seq/previous/canonical hash、global pair index、child index、source ID均与固定输入一致，全部89426个pair恰好出现一次。所有child按index mod15唯一归属，分片source副本与原source_pool对象逐字典相等。

每个child第一次保存完整PASS之前必须是evaluated；已有完整PASS之后只能skip，理由字段必须为同一个完整16角点source已认证。每条保存certificate的epsilon、point、target、Delta、gradient、q/rule及source对象SHA都与同一core一致，passed布尔值与保存margin严格正的符号一致。没有把不同source的单角点正号拼成PASS。

实际统计为71812条evaluated、17614条skip，总89426。每分片checkpoint的cursor、日志字节数、最后hash、evaluated/skip计数、qualified_children、fingerprint均逐项核回，并与worker登记和master状态相等。1153张first-PASS证书正是本轮独立算术核验集合。

原137生产全部corner结果并未重算。本轮对失败保存margin的排序比较仅用于准确恢复residual metadata，不给失败数值新的独立精度标签。

## 5. 残余与独立完整森林结构

4274 residual恰为5427输入child减1153个首PASS，无重复、遗漏或自动split。gap_order逐条与日志保存的最佳失败source证书一致，准确box、原direct-parent失败指针和dyadic物理体积吻合；核margin升序、体积降序和node tie-break。

旧98007个通过node从被123最终source_lock锁定的accepted_mapper取出，加123已审572、128/130已审5903及本批1153，得到105635通过。内部57549个split从原098、107、123真实记录重建，逐个核depth mod4分割轴和两个准确child。与4274pending共同组成109909活叶。

本轮树检查与130区间铺排不同：把每条path编码为二叉堆编号，核所有非根node的parent属于internal、每个internal的两个child都存在、leaf/internal不交；再逐根以最大深度下的整数权准确核Kraft=1。完整有限二叉森林的这些关系排除了leaf祖先/后代重叠，并保证52360个根的全部闭覆盖；不只核叶/分割数量。

独立结果：52360根、57549split、109909活叶、antichain和每根Kraft=1全部通过。

## 6. 真实首次失败、最小schema修复及成功执行

初版把历史certificate的冗余target误当成必填字段。实际PID56300、exit1，在原H64祖叶元数据阶段触发 `KeyError: target`；当时尚未打开任何新corner输出，实际新角点为0。原code、plan、run001登记及失败回执完整保留，没有改称正常退出，也没有重用不存在的成功corner前缀。

只读实际源确认祖叶 `[52210,"1100"]` 的外层target=21/200000、EPS64、原passed/正margin及leaf SHA完全匹配，而其历史inner certificate没有target；其原叶SHA为
`ada586e67825bb3c1577557fad5b5b6e1120e1040daa6af93147f13115d8566c`。
v2仅改为外层target必须g105，inner若存在则仍必须相同；缺失时绑定已主审完整g105来源和原schema。EPS/pass/margin/sourcehash/point门和全部数学算术均未变，144输入锁完全相同。最小diff与schema回执一并保留。

主阅读全文原代码并授权首批，又审最小diff后另授权v2新排他run002：

| 文件/执行 | SHA256或实际信息 |
|---|---|
| 原plan | `b9929d305a2864b6e7a5fc91d93b6415275afe0c2e6fcb1618a84760759ab531` |
| v2 plan | `82d329c1d19bf01e1846036feb5d914a051d9f6be4b8734dd4bde97eadc57560` |
| v2 driver | `e15ddd55e40a2076b125a2a4531191b11b921cfa7037b83ac347a807bcf4aab4` |
| exact137 | `d5ec213dfed9fc241642970617e50def8667417bfaa9c0cc66e61f3daaf8dc82` |
| JSON reader | `dd67f81c1a623507c2d4c75d71e5585c1ef4b4a098b042f9098d1a219029e5db` |
| 首次授权 | `3825c58264bdb0188af6642030a94cba51dd506c6011b73e7bca17298f9c836b` |
| v2授权 | `e9303165637742c5ca356ebc728029baf89c3c370b2b2232304bc8145df2f198` |
| 成功PID/退出/秒 | 56846 / 0 / 58.10661691700807 |

新角点实际只在v2完成一次，共18448个；原生产/失败角点重算0。运行后144输入SHA、冻结代码和原计划/失败文件再次核对未变，实际ps确认56300、56846均已不存在，无后台pool。主仍保留完整证明与最终组件接纳责任。

成功 `run002/audit_result.json` SHA：
`ece7d4be9093dfcba9c65d481efe4294bb83f54dbb97cf50ce4d6518f4d8adf3`。
最终 `final_check.json` SHA：
`7f3d368040b02bd7d6d6c9a7ff7fdc2f2d258894a3f079d89b92e643b33a4544`。

## 7. 精简白名单和复核边界

逐证整数JSONL无损压为 `evidence/winning_integer_checks.jsonl.gz`，160808 bytes，SHA
`9113d6394d1f519dc4cd363b13e92e9656a75ed3817203a96fe76ca897e0c717`。
实际解压流逐字节核原JSONL SHA并确认1153行；该归档校验无数学重算。未压缩原输出继续本地冻结保留，不与gzip重复归档。

最终manifest只列必要报告、独立源码、原/v2计划与两次授权、最小schema差异及真实失败/成功收据、压缩证据和终锁；不复制约2GB旧源或89426生产日志。复核使用计划内的原只读路径和SHA；已成功的1153证不应自动重跑。任何后续新候选、4274残余H生产或进一步split须另任务登记，本组件不作这些操作。

本轮关闭的是137已完成有限跨锚点批的局部独审门；4274pending没有移除，globalg110=false、S=0，不能从这个partial结果预称新I/C5。
