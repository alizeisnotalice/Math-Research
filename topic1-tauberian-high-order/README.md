# 中心方体 Tauberian 常数的高阶端点研究

研究 Lebesgue 测度下的中心轴平行立方体极大算子，使用严格水平集。
线性输入为任意有限正测度可测集；共享方体的多线性常数以公共交集归一化，单个输入集允许无限测度。

## 当前严格结论

对每个固定整数 n≥2，置a=2^nt、v=2t^(1/n)，真实常数满足
\[
C_n(1-t)=L_n(t)+O_n(av^3),\qquad
L_n(t)=\left[1-\frac{a}{(1+v)^n}\right]^{-1}.
\]
因此
\[
C_2(1-t)=1+4t-16t^{3/2}+64t^2+O(t^{5/2}),
\]
且n≥3时
\[
C_n(1-t)=1+2^nt-n2^{n+1}t^{1+1/n}
+2^{n+1}n(n+1)t^{1+2/n}+O_n(t^{1+3/n}).
\]
匹配上下界覆盖全部允许输入和全部半径；常数及起效阈值依赖固定n，没有维数一致主张。
尚未确定(C_n−L_n)/(av³)的极限或下一项系数。

多线性超临界区精确约化为
\[
C_n^{(m)}(\lambda)=C_n(1-m(1-\lambda^{1/m})).
\]
阈值为(1−1/m)^m；m≥2时在阈值及以下发散。
将t换成ε=1−λ后，上述三项渐近及余项对全部m≥1一致传递；二维ε²系数须改为66−2/m，n≥3的已识别系数不变。

## 材料与历史

- `paper/tauberian.pdf`：实际编译、检查的论文；`main.tex`、`sections/global.tex`和`sources.bib`为可编辑项目。
- `rounds/008.md`、`reviews/008-*.md`：第三项闭合、基线复查、独立审查、正式转写及本机证据。
- `literature/2026-11-05.md`、`sources.md`：最新针对性文献检索及模型比较，不宣称穷尽或新颖性。
- `proofs/`：保留仍能排除错误路线的反例及受限构造排除；`verification/`为必要精确证书及探索复现材料。
- 第4轮首次闭合真实二阶；第6轮取得幂余项；第7轮证明固定异常尾桥的实际反例。第8轮回收异常区域的空隙与容量、选择加权边界，绕开尾桥并闭合第三项。

正文按当前最强证明组织。已被覆盖的空间缓冲主证明、同半径及局部小扰动等重复结论不再占用正文，可在各历史提交恢复。第7轮反例仍然有效：有界第三尺度亏损并不强迫任意固定阈值的异常尾o(vS)。

## 复现

```sh
cd paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

需要XeLaTeX、BibTeX、Songti SC和Times New Roman。主要符号工具为本机Mathematica/Wolfram15 ARM：

```sh
wolframscript -file verification/wolfram/global-third.wl
wolframscript -file verification/wolfram/multilinear-third.wl
python3 verification/constant_hierarchy.py
```

Wolfram脚本检查明确假设下的标量不等式和统一参数传递；Python以15进程精确检查n=2,…,16的常数层级。一般维数几何、源目标合法性、全部输入和极限传递由论文解析证明承担，不能用有限计算替代。
