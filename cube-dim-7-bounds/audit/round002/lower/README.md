# 第002轮下界独立复核入口

`recheck.py` 只读取文件并向 stdout 输出 JSON，不写结果、不生成 bytecode，不导入或执行作者数学函数。它包含整数 Newton 七次根隔离、Fraction 区间构造以及完整带坐标标签盒的固定坐标事件扫描。只计算已登记的一个获胜实例，不求值第001轮种子、不继续搜索，也不复核73个事件的全局完整性。

集成到仓库后，从仓库根运行：

```text
python3 cube-dim-7-bounds/audit/round002/lower/recheck.py
```

默认候选依据脚本自身位置定位同项目 `lower/round002/`；也可提供 `--candidate-root PATH`。默认只比较候选 `frozen_round001.json` 中完整有理区间，并在输出标记 `frozen_interval_only`，不会访问任何写死的旧证书路径。

可选 `--reference PATH/TO/R001/results.json` 给出完整第001轮聚合结果，程序才核验其 SHA256 及 `full_construction_comparison_interval` 与冻结区间相同。两种模式都比较完整第001轮上端，不以简单展示界2.135677代替。它不重新计算旧种子的union。

依赖 Python 标准库及同目录 `input_lock.json`；后者绑定作者最终源文件。`results.json`、`receipt.json` 是在登记范围内实际执行后由外部收据程序保存的历史输出；复核入口本身不会改写它们。证明审阅与准确范围见 `proofreview.md`。
