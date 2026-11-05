# Round 032：平坦权薄孔带模型（最小实验附件）

此目录保存第 32 轮的最小可复现实验和一个解析正例，而不是完整数值轨迹。模型为
\(w(x,y)=|x|^{-\gamma}\) 的有限条带薄孔下界构造。两个历史 DE 搜索只覆盖
\(\gamma=1/2\)，没有找到 \(J>1\)；独立加入的 \(\gamma=1/20\) 有理两带例则
严格给出 \(J>1\)。因此不能再声称平坦权所有 \(\gamma\) 的首系数都是 1。
这些材料仍然**不构成完整 centered-halo 常数的上界**。

## 内容和边界

- `model.md`：连续构造、每个孔点的实际中心方体见证，以及从构造到 \(J\) 的证明。
- `band_model.py`：几何、有限滑窗节点、分母和给定几何的内层高度优化。
- `band_search.py`、`band_search_three.py`：原二带、三带 differential-evolution 搜索源码。
- `band_fraction_verifier.py`：不导入优化器的独立 `Fraction` 验证器。
- `band_rational_inputs.json`：两组 \(\gamma=1/2\) 有理化搜索输入和一组
  \(\gamma=1/20\) 解析正例。
- `summary.json`：实际运行、源码签名和不应超出实验范围的结论。

刻意没有归档完整节点证书、DE 轨迹、缓存、日志或早期 81^2 像素草稿；它们仅保留
在工作区，且不能作为证明材料。

## 最小复核

在此目录中运行：

```sh
python3 band_fraction_verifier.py
```

该命令逐个重建所有滑窗极值节点及并集截面，不使用 `band_model.py` 或 SciPy。
它严格检查 \(0<\gamma<1\)、\(0<\eta<1\)、半径、非负高度、带的排序/无正长度
重叠和正分母。对于同半径同高度、且 \(\beta=1-\gamma=p/q\) 有理的分支，还以
`Fraction` 输出精确 \(J^q\) 并检查输入所要求的严格 \(<1\) 或 \(>1\) 关系；
其他输入只获得有理可行性检查，任何 `Decimal` 的 \(J\) 值只是诊断。

要重放历史搜索则运行 `python3 band_search.py` 或 `python3 band_search_three.py`。这会
重新执行有限模型的 DE，并不是本附件最终验收所必需的步骤。
