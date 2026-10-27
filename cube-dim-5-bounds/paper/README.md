# 中文上界论文（内部研究稿）

最终 PDF：`cube5_upper_bound_zh.pdf`，10页；来源 `main.tex`，XeLaTeX/latexmk，TeX Live 2026、ctex/Fandol。准确 hash 和版面检查见 `final_qa.json`。

从本目录编译：

```sh
latexmk -norc -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

结论为一般 L1、Lebesgue、五维中心轴平行方体、共同半径和全部正尺度的 C5 ≤ U166 < 28936079/1000000。U166 的准确四次幂、区间及新旧差值见 `../round019/primary_upper166_acceptance.json`，不能从论文展示小数重建。

先完成177/169首审，再执行178积分与180全域空间第二轮严格审计，最后主代理逐项终审。论文稿件另有177审读；最终仅排版调整。原始H及物理导数区间仍为共享的、已审查的继承输入，不声称本轮从头重新计算第二套H引擎。完整来源和可携入口见论文附录、`../rounds/019.md`、179归档和180 portable说明。

下界按用户要求停止。本稿不是原双侧目标完成报告：L_ref未严格超越，双侧时点不存在，后三轮0/3。内部结果不称公开发表或全球纪录。
