# R076 固定有限深度供货计算包

范围是任意预先固定 n≥2、N≥9 的有限层同族供货和条件误差传播，常数可依赖 N；不是对 N 一致，更不是 C_n 全阶渐近或九阶系数。完整几何证明及终审见 `../../proofs/076-finite-depth-supply.md`、`../../reviews/076-finite-depth-supply.md`。

## 原有实际有限算核

- `finite_depth_fraction.py` / `finite_depth_fraction.stdout`：标准库 Fraction，无输入文件；脚本直接打印 stdout，不写 JSON。保存的运行只给 N=9、10、13 的有理样例，不能替代任意固定 N 的解析重排。
- `finite_depth_check.py` / `finite_depth_check.json`：原审稿人独立 Fraction 核，无源文本依赖；检查 N=9…64 的有限望远镜系数、81 项精确逆移恒等式及有限阶次序。JSON 是原脚本输出，未重跑。

## 可逆来源与可移植版本守门

- `final-review-manifest.json` + `final-review.diff`：从最终 proof 反向恢复旧源的机器记录/可读差分。旧源和中间 guard 长文件不重复打包。
- 原 reviewer `version_guard_check.py` SHA256 为 `cc59a57aae645eb024e736c6332a62d2b8d8e4b903661d124f2adf9beac0de43`，原版需读取 old、guard、final 三份全文。为避免打包三份长旧稿，`version_guard_portable.py` 改为仅读取本包最终 proof 和 `final-review-manifest.json`，在内存重建 old/guard，再保留原有 diff、256 个 exact clip guard 与同 eroder 反例断言。改写范围只有输入供给和复建来源；保留其余数学 assert。该改写已实际运行，结果见 `version_guard_portable.json`；没有写回任何源文件。

运行可移植核：

```sh
python3 version_guard_portable.py > version_guard_portable.json
```

输出恢复 old SHA `70b3781c…6dcac` 和 guard SHA `a4391e17…e6bb4`，且报告 normal loss=1/2、full loss=0。它只验证补丁字节和有限标量 clip 断言，`geometry_certified=false`。源、审查及算核逐文件 SHA/字节见 `manifest.json`（相对本 README 所在目录）。
