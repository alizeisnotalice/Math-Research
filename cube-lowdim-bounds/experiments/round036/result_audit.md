# 第36轮修正版first全域补验独立终态审计

**结论：本次唯一修正版first运行正常PASS，补齐第35轮定位的旧first积分回调契约缺口。** 固定B1核及其成本不变；没有新优化、更低常数或C2<3.4成果。旧上游源码、旧run/log/report不改写；此后引用first严格复核应绑定本轮修正源码与新运行记录。

本审计只读实际终态文件，并用标准库检查JSON、字节及SHA；0核评价、0积分、0模型、0LP、0成本重算、0checker重跑。已审`summarize_once.py`，其唯一已执行结果为`passed_corrected_first`，没有再次运行它。

## 实际运行与完整性

`first.run.json`为`terminal`、exit0，实际59.04612041600922秒，小于预定180秒checker硬限；无超时、重试或改深度。实际命令绑定修正目录中的driver、原B1文件，`--only MKA --Q 8 --maxdepth 7 --Kdepth 6 --procs 6`，原96bit精度保持。

日志记录：K检查2628个knot格，失败0；AK在`delta=h/8`下共评价219732盒，失败0，最终独立一行`PASS`。report保存`K_fail=[]`、`A_fail=[]`、`A_cells=219732`、支持盒U=9/4。日志的最小下界`1.346e-08`是显示值，本轮不将其视为新的准确余量常数。

已实际用标准库逐文件核对：5项approved输入hash与run的before/after完全相同并匹配当前文件；B1文件实际hash与前后记录均为`21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172`。原上游hlcert仍为旧SHA，未被本轮覆写。修正版hlcert为`4c347b7782120f64e20df07a44b51453feff8eb2e3ed9ab4abfec68f2776fa05`；driver仍与公开原件逐字节相同。补丁范围和非finite回退审查见`patch_audit.md`。

关键新证据hash：

| 文件 | SHA256 |
|---|---|
| first.run.json | `a9d829a9359d810fc7b2faa52310121b0cf15a581c65ae7146d9df0acf600a7a` |
| first_report.json | `5a31999e4a8a94117c85aaa6c4b6d63e5abb899e81832916ad3dacec493f3d09` |
| first.log | `21d6dadb5eadb9e833dde46af397a8e7218e50141c5758ac32b94326e20c0278` |
| result.json | `de3fede964d5275e92dc04c944e395fdc9bb3f26dc5c20a7ed3746a9e429ebd3` |

新report与旧`scaled_first.json`确实逐字节相同，这是保存字段和数值结果相同；新源码、运行时点、耗时、run及log提供独立的新运行来源。不能因report相同而给旧run换新源码hash，或宣称旧API使用原本就符合契约。

## 证明口径

本轮补丁只使两个积分口的所有复log正确传播解析性标志；它没有改变积分公式、生成元归一化、核系数、样条模型、阈值或比较机制。现在有本轮修正first的实际全域闭合结果，以及此前不依赖这两个回调的独立闭式全域结果；两程序仍共享Arb和既有解析公式，不能称完全独立底层算术。

日志质量球为`[3.590901146403276050426939597 +/- 7.75e-28]`。原driver的PASS依据K/A失败列表，并未把该展示小数与新目标作准确比较；B1的准确成本及严格比较继续承接既有整数根/Fraction证书，不将此次质量打印冒称新成本证书。原外部解析传递及分布接口也保持既有依赖口径。

第35轮排除第34轮候选的两个严格负点仍成立；本轮只补验B1，不修复或重验那个失败候选。用户的C2<3.4目标尚未达到。
