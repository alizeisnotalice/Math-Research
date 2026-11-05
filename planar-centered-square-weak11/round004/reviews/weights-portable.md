# R004 便携接口专项终核

2026-11-06。数学定理及新候选的有限双面积核已在此前审报通过，本次只核文件身份、接口迁移、移位实跑与发布范围；未重审原数学，也未重跑两个有限候选面积批次。只写本人review目录，没有修改Terra冻结包或公开clone。

**最终便携stage审查PASS，无剩余修改要求。** 执行／身份核验通过；唯一README不一致已在总stage准确修正，并已独立核完小范围发布diff，没有重新运行计算。

## 核验对象与结果

源包为scratch/round004-weight-public-package/{weights,jet}/，最终分别映射round004/weights、round004/jet。两SHA256SUMS的34项（含外部依赖）全部独立重算匹配。weights的14个原始输入、日志、check2 JSON、raw运行记录、有限窗／WL程序和原输出逐字节等同已审源码。run-receipt.json中的A第一核、A第二核、B第二核元数据逐项匹配原记录；B第一核PID保持null，没有补造。A/B仍恰为第51质量减1/5000，再34→27转1/1000000，坐标、周期、reach、阈值及质量正性未变。

实际重建三个迁移文件的diff，与jet/portable.diff完全相同；audit.py和原jet.json字节未变，jet.py全部数学函数的AST逐项相同。改变仅published路径、--records-dir／--output-dir参数及生成文件读写位置。根目录从脚本位置推导；audit.py导入jet时共享解析后的records/output参数；check_calibrated.py读jet.D内结果，scalar.wl接受显式输出目录。原weights程序自身从其脚本目录读取文件，无须代码迁移。

独立读取实际replay/output/jet.json并与原jet.json比较：除真实PID／耗时外，所有字段完全相同；便携stdout内完整JSON也与replay结果一致。score为601字符，SHA256为9d52fc0dfb8d8b95d3f9821cbb119e2115f6bb8a592ee32d8e09e2f49c6d87b4；基线面积、比值导数、正间隙、根式正区间、统计和scope均未变。6811条记录与原R3有序记录一致。

实际便携输出保留：jet PID7242；20433副本检查PID7255；t=10^-9严格子族检查PID7280、差商0.18655887765400908；Mathematica15 PID7305的五项断言均True。weights便携Fraction输出逐字等于原标量输出，portable WL PID7370除真实PID外逐字同原，有限窗阈值仍超过1.7027535096。这些是已保存实际移位记录的身份核验，没有将原记录改名冒作新运行。

复现依赖图完整：先用已公开R3 decompose生成decomposition.json，再将生成目录传给R4程序。原始2.85MB记录不是发布必需文件；replay/和生成缓存仅本地保留。R4仍声明固定严格校准有限子族，未扩大为全部理想halo导数或最优性。

## 独立压缩算法便携副本

应root请求另提供portable-independent/independent_verify.py，最终可放round004/jet/independent_verify.py。它直接读published A输入及--records-dir的R3生成记录，输出到--output-dir；与此前独立算法的差异只有路径、参数、输出位置与依赖文件身份的路径。原周期格计数、二维压缩coverage、clip和有理系数计算未改。

已把该副本移到独立项目布局，从工作目录“/”实际执行一次，读取的是R3公开decompose的真实再生记录，而非重新发布大JSON。PID8246、exit0，算法耗时0.919323秒，总进程耗时1.172587秒。精确601字符score、基线面积、严格间隙、根式正区间及所有原数学字段完全相同。再生有序records的SHA为eafcc4722314fb28c3bc6c2eba1cd47992198ce4b5a6814b04a3b542d57ca7fc。

副本SHA为6b6ffbcd70566869e759e7f88999d4af80444b1d53c9cc20d728a2be72439bd4；实际receipt SHA为72def6745faf95cae24e7cf86cbd0666169cfdc3d2b83e7608745778b1747991。portable-independent内含简明README、实际输出／日志、receipt、与旧脚本的diff及SHA清单；该目录内replay/亦不发布。此处没有声称包含Python启动的总时间小于1秒。

## 已完成的发布文档更正与最终stage

Terra原jet README末段称完整重复stdout只留本地、不打包，但实际portable-jet.log是2207字节完整stdout且在包及SHA清单内。程序、数据和运行均无问题，文字需与实际发布选择一致。

root在总public-stage的jet/README已将相关句替换为同义准确表述：

“portable-jet.log is retained as a compact replay stdout; full generated geometry records remain local and are regenerated”

已最终检查scratch/round004-public-stage/round004/{weights,jet}：全部原执行源码／输出逐字未变；主jet README仅更正该句并追加独立核指引链接，SHA256SUMS只有README哈希行更新。其发布README SHA为e0dee8e8494dc94ed3fa337bfb07f8fd44a80992cec4656f136dc869a6140699，原SHA540bd086e87c1b4845129c68c4433c88acc8f31e661040064e7c6ce8ef78c6e8由root总source-bindings保留。

独立脚本与JSON逐字复制到jet根；独立receipt／diff只更名为independent-portable-receipt.json、independent-portable.diff，字节未改，没有覆盖主receipt／diff。independent-README.md只更新这些发布文件名并省略重复日志引用。共核37个stage文件，新增五个独立核文件身份正确；未带入replay/或大decomposition.json。最终stage小差异核见portable_stage_checks.json，此步没有重新运行计算。root负责其余总manifest和source-bindings。

本轮静态核脚本／结果为check_portable.py、portable_checks.json/log；它们是本地审查支持，不要求全部纳入发布。
