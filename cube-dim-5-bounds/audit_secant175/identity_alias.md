# 175 冻结 identity 字段别名

producer consumer.py SHA `2ea49d0dab736ccee4f59672e3354476dac86d17873322790bb26f8c819286ae`
的 setup 第19行将 x 设为 plan.preflight，第23–25行核实际 preflight 注册；
第26行将 x 改为 plan.selfcheck，第29行用 x.sha256 构成 identity.preflight_sha256。
因此该历史字段的真实值是 selfcheck SHA
`42f8a208b8941f3fc57eb4d42d638afe21bf9c5292eeac7ef0d9198e2332275c`。
它不应被解释成真实输入 preflight002 SHA
`0b9ca61e7f9ac8a4d6fb5fa7ae8c02b0cd0bcbb99f4b4346c5cc34fd335ab8bd`。

原生产源码、plan、所有identity和输出完全保留。主已核源码关系并接受明确别名绑定。
本独审candidate binding原样核真实identity，并另外保存typed_input_preflight、typed_small_selfcheck，
两者各有实际原路径和准确SHA；真实输入和注册另由source binding和producer plan锁定。
独审运行前后均锁两者及原producer源码，逐行输入仍与真实preflight002的shard绑定。
不以改名覆盖历史事实，不因此重跑生产或八项selfcheck。

初次metadata binder因假设字段名等同实际preflight而停于identity比较；
记录在binding_metadata_failed001.json，初版源码另保留。该次真实根/斜率/gain/DFS为0。
新binder只修复该来源schema，不放宽任何实际输入或数学门。
