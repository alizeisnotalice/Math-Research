# 中心方体 Tauberian 常数的高阶端点研究

研究对象是 Lebesgue 测度下的中心轴平行立方体极大算子，使用严格水平集。
线性输入为任意有限正测度可测集；共享方体的多线性常数以公共交集归一化，单个输入集允许无限测度。

## 当前严格结论

对每个固定整数 n≥2，真实常数满足
\[
C_n(1-t)=1+2^nt-n2^{n+1}t^{1+1/n}
 +O_n\!\left(\frac{t^{1+1/n}}{\log(1/t)}\right).
\]
因此原目标归一化余项的极限存在，等于 −n2^(n+1)。下界与全部输入、全部半径的上界匹配；完整解析证明见 `paper/main.tex` 及 `paper/sections/global.tex`。

多线性超临界区精确约化为
\[
C_n^{(m)}(\lambda)=C_n(1-m(1-\lambda^{1/m})).
\]
阈值为 (1−1/m)^m；m≥2 时在阈值及以下发散。二阶渐近及上述余项可关于全部 m≥1 一致传递。没有维数一致性的主张。

比 t^(1+1/n) 更细的首个非零修正项尚未确定。周期包络的更高展开不是实际常数的高阶展开。

## 材料

- `paper/tauberian.pdf`：实际编译、检查的可交付论文。
- `paper/main.tex`、`paper/sections/global.tex`、`paper/sources.bib`：完整可编辑项目。
- `rounds/004.md`：匹配二阶上下界、证明依赖、复核、实际计算与剩余问题。
- `reviews/004-*.md`：几何、计费、独立敌意审查及最终转写复核。
- `rounds/001.md`—`003.md`：基线修复、被排除路线和历史研究记录；当时文件可在相应 Git 提交恢复。
- `proofs/`：仍有筛除错误路线作用的反例和有限周期截断说明。
- `verification/`：最小精确代数与本机符号复现程序。
- `sources.md`：文献对象、版本、读取范围和题录核验。

第4轮从当前正文移除了已被一般二阶定理覆盖的避邻孔、同半径、缓变半径和局部小微扰结论，以及不再需要的一般测度转移长证。原始稿件和历史版本可恢复。

## 复现

```sh
cd paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

需要 XeLaTeX、BibTeX、Songti SC、Times New Roman。主要符号工具为 Mathematica/Wolfram Language：

```sh
wolframscript -file verification/wolfram/sharp-check.wl
```

该脚本在 Wolfram 15.0.0 ARM 实际运行，16项标量/几何余量检查全部通过；计算范围与证明的关系见第4轮报告。一般维数、测度论和所有输入的结论依据正文证明。
Python 有理数与本机 MATLAB R2026a 的基线交叉检查保留在 `verification/`。
