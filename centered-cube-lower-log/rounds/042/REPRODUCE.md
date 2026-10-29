# 第42轮复现与证据范围

默认仅需Python3、NumPy。主运行Python3.12.14/NumPy2.3.5；两项独审Python3.15/NumPy2.5.4。所有运行均为本机实际执行。没有调用远端数学计算服务。

```sh
python3 verify_saved.py --output /tmp/cube-round042-check-new.json
```

默认检查归档hash、三份冻结selector的精确质量和循环轨道；在临时副本实际重算通道扫描的精简摘要、Fraction小例、局部源积分包含核验以及36个置信端点的整数证书。它不会再训练或重放全部MC。全过程只修改临时副本，输出目标须尚不存在。浮点谱摘要按数值容差比较，并非字节哈希复现；原始hash链始终逐字节核对。

若需重放三组全部原query-cell包络：

```sh
python3 verify_saved.py --replay-cells --output /tmp/cube-round042-full-new.json
```

该模式额外实际重算约1.5分钟的完整cell包络，核对同一组冻结counts/CI；不重训、不新增样本。NumPy未来版本可能改变随机流，若counts变化不能默许当作原样本复现，须另行标记版本差异。

完整训练可以在一个全新空目录复制 `round042_selector/probe.py` 和 `balanced.py`，依次执行两个文件；约一分钟。PLAN与selector均由程序先后冻结，已有文件时拒绝覆盖。原浮点`>1`算法存在本轮明确记录的边界失效，因此重放原训练的目的仅为追溯冻结selector，不应重新使用其原始naive CI。随后复制interval/cell/exact_ci三个最终脚本，先运行interval_recheck，再cell_recheck，再exact_ci，才能得到本轮修复统计。不要直接在归档目录运行这些有写入输出的入口。

权威统计为 `round042_selector/CELL_RESULTS.json`（whole-cell分类）与 `EXACT_CI.json`（36个有理端点的精确整数Chernoff比较）。INTERVAL_RESULTS是先完成的dyadic点包络，因whole-cell扩张仍给相同counts，它作为执行与审查中间桥保留。原三份RESULTS只作被否定的浮点边界记录，已由总报告说明，不是当前权威CI。

`round042_channel/probe.py`会生成20MB完整RESULTS和35KB SUMMARY；前者可再生成，未提交Git。不是确定性谱区间证书。数值共同D仅覆盖有限参数网格，不能用于连续族定理。

必须区分：一般数学定理、既定query上的浮点包络（以标准binary64语义为前提）、理想随机抽样下的概率CI、确定性连续体积分，以及C_n的渐近下界。本轮未取得最后一项，也未证明PRNG产生真正iid随机位。
