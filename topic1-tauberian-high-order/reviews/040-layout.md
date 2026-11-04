# R040 最终布局差异独立审查

审查日期：2026-11-05。审查者：Sol（main 独审代理的单独布局复核）。

**裁决：LAYOUT-EQUIVALENT / PASS。10 项替换仅改变断行、括号字号和中文数学注记的排字方式，没有改变数学公式、条件、量词、集合或估计。** 最终日志已完成，未见超宽、欠宽、中文缺字或未定义引用。本次未重复编译，未改源文件；仅新增此报告。

## 1. 最终源和独审基线的绑定

检查的编译快照为：

`/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/r040-compile/project`。

实际 SHA256：

| 文件 | 最终 SHA256 |
|---|---|
| main.tex | `de97f3808eec8b1792b9a0f5748cd50fc7efeeba56f6b7d84bfffa418178f206` |
| sections/seventh.tex | `22b02210cabcb91934868a5dfa7ef33317193752505d9fcfc6662af64af379c9` |

读取 `layout-changes.json` 全部 10 个条目，逐项核对 old/new；另用独立只读 Python 程序在内存中逆序替换，断言每个 new 的出现次数等于记录次数，再对逆还原文本求 SHA256。实际结果：

```text
entries: 10
inverse main.tex 045ef59f0dfb233418c1e854fc8891b9b9cbc71529d4bc4b9f78b3d5f90a542e matches True
inverse sections/seventh.tex 8293e0be9aeced41e1e598d662043c2bc88e01e49952e21c1f45af0b06430d19 matches True
unchanged sections/global.tex True
unchanged sections/sharp.tex True
unchanged sources.bib True
```

这里 module 基线是已经数学修复并交另一独审绑定的 `8293e0be…30d19`，不是更早的 `d6ad9a1…9ee71`。本报告只判定最终排版源与这个修复后基线等价，不将数学修复本身归类为排版。

global、sharp、bib 对照 manifest 实际 SHA 全部未变。源清单旧字段仍记录初次快照，故 module 的本次正确逆还原基线明确以上述 `8293…` 为准。

## 2. 十项差异逐项判定

| 条目 | 内容 | 数学判定 |
|---|---|---|
| 1 | main 的 b_(n,j) 公式改为 aligned，两行显示；left/right 方括号改为 Bigg 方括号 | floor 上限、二项式上下标、符号、j=6 与 j=7 修正、外因子和 0≤j≤7 全部相同；仅断行、对齐、字号 |
| 2 | seventh 的临界族条件改为 gathered 两行，末尾逗号改句号 | v→0、V、K、M 定义和 1+M/V 下界全部保留；仍在同一 equation 中，没有改变条件的合取关系 |
| 3 | seventh 的两个半径/横向乘积等式改为 gathered 两行 | 两个等式、误差阶、求和指标、分母和不等号完全相同 |
| 4 | seventh 的 C_i(K)、Z(K) 定义改为 aligned 两行 | 集合条件、至少两轴、阈值 Kv^3r_i、缺集差、绝对值、v/a 和求和均相同 |
| 5 | `{\rm 非正常}` → `\text{非正常}`，1 处 | 求和指标中文字相同，改变字体处理 |
| 6 | `{\rm 切向}` → `\text{切向}`，2 处 | 乘积指标中文字相同，改变字体处理 |
| 7 | `{\rm 此角}` → `\text{此角}`，1 处 | 积分域中文字相同，改变字体处理 |
| 8 | `{\rm 正常}` → `\text{正常}`，1 处 | 求和指标中文字相同，改变字体处理 |
| 9 | `{\rm 巨大}` → `\text{巨大}`，2 处 | 求和指标中文字相同，改变字体处理 |
| 10 | `{\rm 邻边界}` → `\text{邻边界}`，1 处 | 求和指标中文字相同，改变字体处理 |

总共四处显示公式调整、六类中文替换，中文替换实际覆盖八处。逐项 new/old 不包含 label、ref、eqref 或 cite 的增删。未发现数学字符、系数、指数、关系符或限定范围的改变。

## 3. 交叉引用复核：正确计数是 185 / 226

对最终 main、global、sharp、seventh 四文件实际执行静态扫描，结果：

```text
labels: 185 references: 226
duplicate labels: {}
unresolved references: []
removed-module references: []
```

**最终正确计数为 185 个唯一 label、226 处 ref/eqref/autoref，而非较早 main 独审阶段的 225 处。** 在本次排版替换中没有引用增删；修复后 `8293…` 基线与最终源的引用多重集相同。不能为了沿用旧计数而报告 225。无重复定义、无悬空引用、无 planar/skew 残引，交叉接口通过。

## 4. 最终已有编译结果

实际读取完成后的 `compile-layout.log` 和 `project/build/main.log`。本代理没有运行 latexmk、xelatex 或其他重复编译命令。

最终 xelatex 日志含：

```text
Output written on build/main.xdv (51 pages, 746792 bytes).
```

latexmk 日志随后记录 xdvipdfmx 生成 PDF，`564777 bytes written`，并以

```text
Latexmk: All targets (build/main.xdv build/main.pdf) are up-to-date
```

结束。实际 PDF 文件存在，大小 564777 字节。最终 `main.log` 对以下字串的扫描次数均为 0：

```text
Overfull
Underfull
Missing character
undefined
multiply defined
Warning
```

日志中仅有 Times New Roman 的 slanted 字形不可用、改用 italic 的 Font Info，属于字体信息，未产生缺字或编译警告。因此初版的四处 overflow 和中文缺字不再出现在最终日志。

本报告完成源差异、接口及既有编译日志验收；没有逐页视觉重审 PDF，故不把无 TeX 警告扩大成对全部页面美观程度的断言。数学等价和所查排版问题通过，不需进一步修改源。

**签署：最终 main `de97…f206`、seventh `22b0…79c9` 相对已独审基线的全部布局差异为 LAYOUT-EQUIVALENT；最终引用计数 185/226，日志无 overflow、缺字或悬空引用。**
