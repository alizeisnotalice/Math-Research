# 第246轮：重启责任与实际阈值出生几何

日期：2026-10-17。一般n^(1/2+o(1))上界未闭合；本轮没有新增全局已付款分支。

核验结果：保留终端score可使重启初始费用精确相消，但真实出生仍带其当时完整归一化分母。给出不劣于旧价的显式行界，并把原C核出生费用写为真实来源面的a.e.coarea式，保全部来源、接收、端点及删除补偿。来源面Jacobian有界不等于来源费用已付；不能把hard-z位移界移给aux-v。

- [研究包](research.zip)：49,709 bytes，11成员，SHA256 `b979e35bb510c305cfa3819c1cc9e43eb17d1814058c1873a7d1cdad9dc1ba22`。
- 包含详细LaTeX文本、任务合同、两作者稿、独审、数值笔记与独立脚本/JSON、复现驱动和10项manifest。不编译PDF。
- 三阶段8组（2/4/2），n=2、4、16、128；46个正后验中点、27个出生事件归属、8个能量及8个方差精确账通过。
- 主代理完整迁移复跑8组，71,302字节JSON逐字节一致；216项独审Fraction检查通过；10项manifest前后通过。40/70位最大缩放差约3.72e-41。有限clock不是原FULL几何或一般阶数证明。
- 解压后运行 `python reproduce.py --output-dir <包外新目录>`，依赖mpmath1.3.0；本机命令详见numeric_notes.md。
- 仓库外逐字节检查点：output/centered_cube_round246_20261017_final_checkpoint；独立回执：output/centered_cube_round246_20261017_final_replay/receipt.json。

下一步仍需在原共同来源和真实接收测度下支付补偿后响应，不能把归约或个别行改善当作全局预算。R243—245有效工具继续保留，旧负信用不重用。
