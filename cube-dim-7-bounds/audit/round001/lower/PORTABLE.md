# 可携带只读复核入口

`recheck.py` 是新增入口，不修改历史 `check.py`、`results.json`、`manifest.json` 或其他已封存文件。它只向 stdout 输出 JSON；不覆盖结果、不生成 bytecode，不调用作者数学代码。

在仓库完整集成候选和审计目录后，从仓库根运行：

```text
python3 cube-dim-7-bounds/audit/round001/lower/recheck.py
```

默认候选目录根据脚本自身位置定位为同一项目的 `lower/round001/`，与 shell 工作目录无关；也可显式指定：

```text
python3 cube-dim-7-bounds/audit/round001/lower/recheck.py --candidate-root PATH/TO/lower/round001
```

不提供 `--reference` 时，只从候选 `frozen_internal_reference.json` 读取完整有理区间。输出明确标记 `comparison.mode=frozen_interval_only`、`original_reference_file_verified=false`，不会访问旧六维原目录，也不声称本次核验了原件。该冻结文件本身仍须符合原审计 `input_lock.json` 的哈希。

如可访问完整原 R075 证书，可增加：

```text
python3 cube-dim-7-bounds/audit/round001/lower/recheck.py --reference PATH/TO/R075/results.json
```

该模式核对原件 SHA256 和 `candidate.exact_interval` 与冻结区间完全一致，输出 `original_reference_hash_and_interval_verified`。原始证书路径不再写死于新入口。

两种模式都由独立完整 raw-box 算法重新计算四个 union，而非只检查保存值；同时检查根隔离、源包含、有限 L1 损失、严格有理阈值和输入运行前后哈希。新入口只借用历史审计器的 `construct` 和 `RawUnion` 等纯算术函数，绝不调用其旧的 path-bound `main()`。不接受 `python -O`，以保留已封存算术中的断言核验。

依赖 Python 标准库及本目录已封存历史审计文件；不依赖作者 clone 的绝对路径或本机旧六维路径。应保留并集成原 `manifest.json` 所列全部历史文件，以核验审计器和历史封存未变。`portable_receipt.json` 记录实际在 `/tmp` 为工作目录、显式候选路径、不提供原件的运行：exit 0、四个 union 完成、状态数 36552/30326，全部历史文件哈希不变。

`portable_manifest.json` 单独绑定本新增入口、说明和收据；原历史 manifest 保持原字节不动。数学范围仍是原有两个固定种子，没有新增搜索。
