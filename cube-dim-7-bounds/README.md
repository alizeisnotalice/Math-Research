# 七维中心轴平行方体的弱 (1,1) 常数

固定 Lebesgue 背景、共同半径、一般 L1 输入及全部正尺度。当前准确状态以 `STATE.md` 和 `results/bounds_ledger.json` 为准；轮次报告在 `rounds/`。其他项目与三四维来源保持只读。

截至第003轮严格得到 `C7 > 2136429741978/1000000000000`，上界仍为128。最新完整下界构造区间是后续比较门槛，不能从展示小数或简单有理数重造。结果为内部研究，未称公开发表或穷尽全球纪录。

## 最小复核

从仓库根目录运行，前两条仅需 Python 标准库；后两条另需 SymPy 与 python-flint：

```sh
python3 -B cube-dim-7-bounds/literature/round001/check_a7.py
python3 -B cube-dim-7-bounds/lower/round001/check_saved.py
python3 -B cube-dim-7-bounds/src/check_primary_round001.py
python3 -B cube-dim-7-bounds/src/check_exact_budget_round001.py
```

完整独立下界并集重算（标准库，stdout-only，项目相对路径）：

```sh
python3 -B cube-dim-7-bounds/audit/round001/lower/recheck.py
```

此入口重算两种种子全部四个并集；默认与已冻结完整R075区间比较。若需核对六维原件，追加 `--reference <R075结果JSON>`，参见同目录 `PORTABLE.md`。哈希核验与保存算术检查本身不能替代完整并集或连续证明。

已运行环境：Python3.15.0、SymPy1.14.0、python-flint0.9.0、NumPy2.5.4、SciPy1.18.1。作者生成程序及若干历史审计入口会写自己的结果，应只在可丢弃副本执行；不要直接覆盖这里封存的证书。上述最小复核入口均只输出结果。

## 证据地图

- `lower/round001/`：两种七维矩形晶格、准确根区间、完整联合体积、有限块与统一平滑证明。
- `audit/round001/lower/`：不同的带坐标标签算法、独立整数根隔离、全部并集及有限L1复核。
- `upper/round001/theory/` 与 `audit/round001/theory/`：七维全空间有限源、一般输入及全尺度比较、真实奇异子空间和短时预算。
- `upper/round001/gaussian/`：完整边界矩预算、严格Gamma上包及有限样本诊断；不提供全域g证书。
- `baseline/round001/`：R075准确原件的只读来源审计与维数单调传递；历史运行不算本轮新运行。
- `literature/round001/`：原始公开来源、版本/hash、准确A7及检索限制。

保留约793KiB的有限观测表供核查8192样本的摘要和未丢弃记录；约8.8MB的可复用形状矩阵缓存留在仓外。没有归档整个资料包、重复ZIP、旧维度计算缓存或虚拟环境。

第002轮最小核验：`python3 -B cube-dim-7-bounds/audit/round002/lower/recheck.py --reference cube-dim-7-bounds/lower/round001/results.json`（完整新实例union、stdout-only）；`python3 -B cube-dim-7-bounds/src/check_primary_round002.py`（准确比较与S3坐标审计、stdout-only）。模型积分预算改进不等于Cn改进。

第003轮后下界遵用户要求暂停，上界继续。只读入口：`python3 -B cube-dim-7-bounds/src/check_primary_round003.py`。

第004轮完整原G积分预算严格降低为 I<8179950070308/10^18；仍无全域g，Cn上界仍128。只读可携带入口：`python3 -B cube-dim-7-bounds/src/check_primary_round004.py`（python-flint）。

第005轮完整I<7730422817512/10^18，仍无新Cn上界。下界按用户最新指令恢复低资源，主要资源上界。可携带全账本标量复核：audit/round005/radau/minchecker.py；主准确接口src/check_primary_round005.py，均stdout-only。

第006轮：严格C7>2136578727666/10^12，上界仍128；完整I<7723627498344/10^18但无全域g。完整portable审计audit/round006/radau4/minchecker.py及lower/recheck.py；主准确接口src/check_primary_round006.py。下界低资源，上界优先。

第007轮无新界：矩阵slack诊断指向matrix-pair严格积分，compact账本入口upper/round007/matrix_slack/check.py。下界按用户最新明确指令暂停，主力上界。

第008轮signed mixed-Mellin StageI门槛严格失败，StageII未启动，无新界；audit/round008/mixed_stage1/minchecker.py为stdout-only完整标量重组。下界暂停，上界转向新的精确frequency-fan接口。

第009轮结构证据已归档，严格结论范围见rounds/009.md；没有新的Cn上界。下界保持暂停。

第010轮结构证据已归档，严格结论范围见rounds/010.md；没有新的Cn上界。下界保持暂停。

R011成果与准确证书见rounds/011.md。

R012成果与准确证书见rounds/012.md。

R013真实结论见rounds/013.md；最新状态以STATE.md与准确账本为准。

R014真实结论见rounds/014.md；最新状态以STATE.md与准确账本为准。

R015真实结论见rounds/015.md；最新状态以STATE.md与准确账本为准。

R016真实结论见rounds/016.md；最新状态以STATE.md与准确账本为准。

R017 proxy证书链屏障及运行失败记录见rounds/017.md。

R018已封存，真实结论见rounds/018.md。

R019已封存，真实结论见rounds/019.md。

R020已封存，真实结论见rounds/020.md。

R021已封存，真实结论见rounds/021.md。完整算子上界已严格降至77.436179562403。

R022严格算子上界77.433340322253，来源及修复记录见rounds/022.md。

R023完整严格算子上界已降为46.650446233303，全部机制/准确复核见rounds/023.md。下界仍暂停。

R025完整严格上界29059800754803/1000000000000，一般输入/全尺度与实际复核见rounds/025.md。下界保持暂停。

R026完整严格上界113513891739/3906250000。见rounds/026.md。
