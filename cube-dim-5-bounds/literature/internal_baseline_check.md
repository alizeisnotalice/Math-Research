# 内部五维基线去重检查

检查日期：2026-10-25。此检查只读以下两个既有本地项目：

- `/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds/`
- `/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-lowdim-bounds/`

检索了文件路径、`bounds_ledger.json` 的维度字段，以及文本中的 `R5`、`C5`、`n=5`、`dimension 5`、`1.992` 和 `1.993`。排除了 SHA-256 字符串的偶然数字匹配；不把包含 `n=500` 等程序变量的记录解释为五维结果。

## 结果

没有发现已经归档为有效五维上下界的 `R5`、`C5`、五维证书、五维论文、或数值约 `1.992`/`1.993` 的完整严格认证。

`cube-finite-dim-bounds/results/bounds_ledger.json` 的当前主条目仅为 L3、U3、L4、U4；其 `lower/` 下的实际维度项目为 n3、n4 和通用辅助目录。`cube-lowdim-bounds/results/bounds_ledger.json` 的 `dimension` 条目仅为 2、3、4。两项目中关于 n≥5 的唯一直接命中是交接说明：可以由已审 R4 通过维数单调传递下界，并可用通用 \(2^5=32\) 上界；说明明确这不是五维新证书。

因此，在这两个指定内部项目的已归档材料范围内，当前出现的约 \(1.992\) 或 \(1.993\) 五维数值不与既有有效五维基线重复。它仍须独立检查输入 hash、共同半径、完整联合体积、严格超水平集余量和一般 \(L^1\) 传递，才能登记为新内部下界。

## 局限

本检查不审计未归档的临时目录、其他用户项目、二进制缓存或 Git 历史中未出现在当前工作树的对象；也不把“未找到”表述为全局不存在。没有修改被检索项目的任何文件。
