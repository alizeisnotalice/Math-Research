# 第44轮候选保存值独立算术审计

主侧唯一执行verify_candidate_saved.py，源码SHA `0ef08169ce9d25517c8d7d00e2567a89863c71790538cbc89419756b5d64f97b`；实际0.07963112500146963秒，终态PASS。

标准库独立重建了solver二进制浮点向量→固定λ1023/1024收缩→唯一ties-even结构量化，26个q与保存参数完全一致，参数盒、z=|t|和非负g通过。

1181个非零系数准确复现为s*a1+T；D4半质量、T支撑≤3/4、plain/structured及source/result/run/params/cost哈希绑定通过。

直接核保存root_floor的40次幂夹逼，未重新搜索整数根；准确成本上下端复现，upper<B1lower且upper<3.59087，upper并未低于17/5。

仅独立重构保存最坏行node658/piece1/Bernstein(2,0)的准确半余量残差，严格正且与solve一致；没有再次检查全部36720候选行或缓存点可行性。

本次0新AS、积分、点或LP；完整结果见candidate_saved_verification.json。它核的是保存算术与执行证据，不是第二套生成元积分或全域证书。

候选SHA `7081cfef6fb47e5d47c95d3740c7a6ba48d48bb6197928a84d8cb6340c4166d5`可进入下一预定全域核验；当前B1仍保持，尚无新的已证C2上界，3.4目标未完成。
