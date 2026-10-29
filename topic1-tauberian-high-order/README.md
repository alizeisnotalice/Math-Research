# 共享中心方体 Tauberian 常数的高阶端点研究

研究对象为 Lebesgue 测度下的中心轴平行立方体极大算子。线性常数为
\[
C_n(\alpha)=\sup_{0<|E|<\infty}\frac{|\{M_{c,n}\mathbf1_E>\alpha\}|}{|E|}.
\]
多线性算子共享同一个半径，分母取公共交集的有限正测度；单个输入集允许无限测度。

## 当前结论

已审查基线：精确多线性约化、临界发散、一般测度统一弱型常数、显式包络，以及
\[
C_n(1-t)=1+2^nt+o(t),\qquad
-n2^{n+1}\le\liminf R_n\le\limsup R_n\le0,
\quad R_n(t)=\frac{C_n(1-t)-1-2^nt}{t^{1+1/n}}\quad(n\ge2).
\]
真实高阶系数和极限存在性尚未证明。下包络的展开不是该常数的真实展开。

- `paper/main.tex`、`paper/sources.bib`：精简论文及准确来源。
- `paper/tauberian.pdf`：已编译并检查的论文 PDF。
- `rounds/001.md`：基线审查、修复、删除和验证记录。
- `verification/`：必要的精确代数与本机符号验证程序。
- `sources.md`：直接相关文献的对象、版本与读取范围。

## 复现

```sh
cd paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error -outdir=build main.tex
```

需要 TeX Live、XeLaTeX、BibTeX，以及 Songti SC 和 Times New Roman 字体。
精确代数核验只用 Python 标准库：

```sh
python3 -B verification/verify_baseline.py
```

`verification/verify_baseline.m` 在本机 MATLAB R2026a Symbolic Math Toolbox 中实际运行过，
用于与 Python 有理数结果交叉核对。程序核对的范围不能替代论文中的一般维数证明。

## 研究及归档原则

围绕真实高阶上下界推进；模型、单向界和严格渐近分别表述。
每轮保留一个独立提交，仅保存最终报告和复核所必需的支撑材料。
原始稿件保留在本机，未被覆盖；此目录只维护精简研究版本。
