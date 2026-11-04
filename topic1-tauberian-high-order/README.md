# 中心方体 Tauberian 常数的高阶端点研究

研究Lebesgue测度下的中心轴平行立方体极大算子，使用严格水平集。
线性输入为任意有限正测度可测集；共享方体的多线性常数以公共交集归一化，单个输入允许无限测度。

## 当前严格结论

置a=2^nt、v=2t^(1/n)，令
\[
L_n(t)=\left[1-\frac{a}{(1+v)^n}\right]^{-1}.
\]
对每个固定整数n≥3，真实常数满足
\[
C_n(1-t)=L_n(t)+O_n(av^4).
\]
三维的四个非恒定项为
\[
C_3(1-t)=1+8t-48t^{4/3}+192t^{5/3}-576t^2+O(t^{7/3}).
\]
n≥4时为
\[
\begin{aligned}
C_n(1-t)={}&1+2^nt-n2^{n+1}t^{1+1/n}
+2^{n+1}n(n+1)t^{1+2/n}\\
&-2^{n+3}\binom{n+2}{3}t^{1+3/n}+O_n(t^{1+4/n}).
\end{aligned}
\]
二维保持
\[
C_2(1-t)=L_2(t)+O(av^3)
=1+4t-16t^{3/2}+64t^2+O(t^{5/2}).
\]
匹配上下界覆盖全部允许输入和全部半径；常数及起效阈值依赖固定n，没有维数一致主张。
n≥3的av⁴尺度、二维的av³尺度尚未确定归一化极限或系数。

多线性超临界区精确约化为
\[
C_n^{(m)}(\lambda)=C_n(1-m(1-\lambda^{1/m})).
\]
阈值为(1−1/m)^m；m≥2时在阈值及以下发散。
将t换成ε=1−λ后，上述渐近及余项对全部m≥1一致传递；二维ε²系数为66−2/m，三维ε²系数为−572−4/m，n≥4已显示的系数不变。

## 材料与历史

- `paper/tauberian.pdf`：实际编译、检查的19页论文；`main.tex`、`sections/global.tex`和`sources.bib`为可编辑项目。
- `rounds/009.md`、`reviews/009-*.md`：实际半径增量角部补偿、三深度消元、独立审查与本机证据。
- `literature/2026-11-05.md`、`sources.md`：本日最新针对性文献核对及模型比较，不宣称穷尽或新颖性。
- `proofs/`：仍能排除错误路线的反例和受限构造排除；`verification/`为必要证书及探索复现材料。
- 第4轮闭合二阶；第6轮取得幂余项；第7轮证明固定异常尾桥反例；第8轮以异常回收闭合第三项；第9轮进一步取得n≥3的第四项。

正文只保留当前证明链。被覆盖的空间缓冲和双深度主路线可从历史提交恢复；异常尾反例仍有效，新证明没有假设固定阈值尾o(vS)。没有把周期包络的任意阶展开直接认作真实常数展开。

## 复现

```sh
cd paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

需要XeLaTeX、BibTeX、Songti SC和Times New Roman。主要符号工具为本机Mathematica/Wolfram15 ARM：

```sh
wolframscript -file verification/wolfram/augmented-corners.wl
wolframscript -file verification/wolfram/multidepth.wl
wolframscript -file verification/wolfram/fourth-transfer.wl
```

角部和系数批次实际各启动15个子内核；多深度脚本实际18项恒等式/阶次通过。有限维系数与标量检查只承担各自范围，一般维数几何、全部输入和全部重数的一致性由论文解析证明承担。
