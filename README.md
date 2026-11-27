# Math-Research

## 进行中的研究

| 仓库 | 研究方向 |
| --- | --- |
| [Math-Research](https://github.com/alizeisnotalice/Math-Research) | 极大算子与调和分析。主线是[中心立方体的一般弱型上界](centered-cube-upper-bound/README.md)；也收录强极大算子、欧氏球极大算子等研究档案。一般上界的进一步改进仍开放。 |
| [Math-Research-2](https://github.com/alizeisnotalice/Math-Research-2) | 多专题研究。包括二维中心正方形精确常数、连续 Cantelli 猜想、二叉 Derrida–Retaux 临界问题、欧拉常数无理性等。各题的结果与未解部分见[课题总览](https://github.com/alizeisnotalice/Math-Research-2/blob/main/PROJECTS.md)。 |

## 公开成果索引

这里集中列出四个成果仓库。点仓库名即可查看论文、数据和复现说明。

| 仓库 | 主要成果 |
| --- | --- |
| [中心立方体极大算子：高维下界](https://github.com/alizeisnotalice/centered-cube-lower-bounds) | 两篇弱型 (1,1) 下界论文：先得到 $\sqrt n/\log n$ 量级，再得到 $\sqrt n$ 量级。第二篇给出 $\liminf_{n\to\infty} C_n/\sqrt n \ge 87/(220\sqrt{58\pi})$。 |
| [中心立方体极大算子：覆盖与高维上界](https://github.com/alizeisnotalice/centered-cube-covering-upper-bounds) | 覆盖法给出 $\log n\ge8192$ 时的显式上界 $C_n<n\log n-n\log\log n-\tfrac{47}{250}n$，并附有限维上界表及七维证书。 |
| [中心立方体极大算子：三维与四维界](https://github.com/alizeisnotalice/centered-cube-dim-3-4) | 论文给出 $C_3<7.916$、$C_4<14.371$；仓库另附可复现的下界证书 $C_3>1.7941091664563071$、$C_4>1.90174261110095$。 |
| [非中心强极大算子：端点常数](https://github.com/alizeisnotalice/strong-maximal-endpoint) | 对 $L(\log L)^{n-1}$ 端点常数，确定 $\log C_n$ 的线性主项及 $\sqrt{n\log n}$ 阶修正项，并给出更细的对数阶上下界。 |

前三项的 $C_n$ 是中心立方体极大算子的最优弱型 (1,1) 常数；第四项讨论的是不同算子及不同的端点常数。

## 研究过程归档

数学研究项目统一仓库。每个项目使用仓库根目录下独立的文件夹；当前项目：

- `centered-cube-upper-bound/`：中心立方体极大算子上界研究资料。
- [strong-maximal-endpoint/](strong-maximal-endpoint/CONTENTS.md)：强极大算子端点常数研究，第1—40轮成果、原始附件及全部已归档工作材料。

其他项目可在根目录下建立各自的文件夹，并遵守仓库根目录的 `AGENTS.md` 提交约定。
