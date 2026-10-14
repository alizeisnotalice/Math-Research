# 第119轮：远尺度逆序来源对的量化测试预算

日期：2026-10-14。一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。

本轮已证明并独审：保留原真实接收测度、两辅助来源的归一化分母及全部历史门，对量化门 ℓ_t≥Bδ 内的正逆序来源对，整个远尺度合计 P_B≤2log(2)X/(KB)。原2Λ回代后取B=4Λlog(2)/(Kε_FA)，可付独立的ε_FA X。存活补集和未来删源共用一次FIRST测试容量，没有给未来删源冒加terminal good门。

中间测试比ℓ_t<Bδ的联合几何、初始/跳跃接收费用与低pe余项仍未付。signed回代由完整H的基本定理重建，不能拿signed表示直接替换此前非负粗上界；有利负项只保留一次，近时槽与本轮槽相加、全部旧费保留。

## 检验及归档

三阶段与隔离重放共六次退出0、stderr为空、逐阶段stdout完全一致。44条JSONL含29条case、6个共享observed表及阶段标记；扩大24配置：n=4,16,256,1024，活动坐标d=1,2,4，B=4,8；M=2,4,8,16仅改变log2证书。固定菜单的fav恒零、活动维最多4，所有实例均为抽象预算，不能证明FULL实现、非零signed抵消或一般n阶。

[research.zip](research.zip)：177,562字节；SHA256：`bafdd4e5c0f378ec309432ae7a67ed2a5a37a7042d115e451773bbadea548e56`。仅含proof.md、proof.tex.txt、inventory.json，51块数学显示式原样保留，未编译PDF。

Markdown嵌入numeric.py、numeric_protocol.md、numeric_output.jsonl、execution_receipt.json。按REPRO-FILE-BEGIN字节数提取并核SHA，将协议与脚本放同目录，依次运行python3 numeric.py --stage 1、--stage 2、--stage 3并拼接stdout核对。--record PATH只写排他文件，不同时写stdout。完整输出约1.45 MB，ZIP约178 KB，不需第三方依赖。

R117尚未完成的最大配置继续单独计算；本轮有限菜单成功不替代其检验。后续初始/跳跃、删源耗散和空间捕获候选未夹入本轮归档。
