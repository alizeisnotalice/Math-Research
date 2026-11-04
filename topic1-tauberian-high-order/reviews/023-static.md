# R023 正文快照静态检查

检查方式：本机 Python 只读扫描四份指定快照的标签/交叉引用、LaTeX 环境、ASCII 控制字符、宏定义候选；`rg` 搜索陈旧措辞、过程日志字样和五/六层表述。未编译，未判断证明有效性，也未改任何源文件。

## 输入快照与哈希

| 文件 | SHA-256（检查时） | 行数 |
|---|---|---:|
| `r023-main/main-reviewed.tex` | `13a8f2681212a039784251e216f2dc424661b628abf627c19075e9f69dd73b04` | 492 |
| `r023-main/global.tex` | `7b6154f74b56e4cb8ea133c99d4c240ceec526004930c510bc29ebcb25b5f10a` | 1309 |
| `r023-typeset/sections/sharp.tex` | `8eae7fa41e8bc1f9883d51330b01183001295ca0b623da12205c807c6bbe7c39` | 777 |
| `r023-lower/sections/skew.tex` | `fe6ffe7c8c52a9252a89461daf775f7f5c43c3a461f13dbc8a567868df72be26` | 100 |

`sharp.tex` 首次读取时的 SHA 是 `4728b69e0bec1d28d7b68cf39d613beb93dfedd4d52b13931025d22a030225ce`（774 行）；检查期间该文件变为表中 777 行快照。最终统计按后者重新扫描，前一哈希只记录观察到的版本漂移，不把早期统计冒充为当前结果。其他三份源的开始/结束哈希一致。

## 静态结构结果

- 合并扫描到 **138 个不同 label、160 个交叉引用**；未见重复 label 或未解析的 `ref/eqref/autoref/cref/pageref/nameref/vref`。特别地，`sharp-upper` 在 [sharp.tex:646] 定义，并由 [main-reviewed.tex:423] 引用。
- 环境起止计数均配平：`document` 1/1、`abstract` 1/1、`equation` 94/94、`theorem` 2/2、`proposition` 5/5、`lemma` 24/24、`proof` 33/33、`aligned` 6/6、`corollary` 2/2、`gather` 5/5、`enumerate` 1/1、`cases` 2/2、`gathered` 1/1。
- 未见 ASCII 控制字符。已声明的项目宏（`\R`, `\ind`, `\dd`, `\eps`, `\Hc`, `\Hm`, `\Mc`, `\Mm`, `\Cn`, `\Cm`, `\norm`）均在主文件导言区 [main-reviewed.tex:21–31] 定义。其余疑似非基础命令由已加载的数学包覆盖；例如 `\mathscr` 来自 `mathrsfs`，本机 `amssymb.sty` 定义了 `\restriction`。未发现明显未定义的项目自定义宏。
- 主文件在 [main-reviewed.tex:416, 418, 420] 依次输入 `sections/global`、`sections/sharp`、`sections/skew`；本检查按任务指定的四份快照合并统计。

## 指定陈旧措辞检查

针对“系数未知/精确系数仍未知”、只展开到 `j\le5`、以及“五层就是最终 sharp 证明”等旧式表述的定向搜索没有命中。当前主定理明确给出 `j=0,...,6`：[main-reviewed.tex:115–131]；摘要和引言均写明匹配六阶结果：[main-reviewed.tex:35–48, 100–110]。这些位置未见“系数仍未知”的旧结论。

### 五层与六层的区分

下列“五层”出现是粗先验或新证明中早层的分工，不能批量替换为“六层”或据此判为陈旧：

- [global.tex:992, 1005, 1101, 1115] 的“五层目标/边界 CDF/递归”属于原全输入粗先验的五深度构造（其引理标题为 `global.tex:986` 附近的“五深度”）。
- [sharp.tex:679, 683] 明确说“早五层”，随后在末层单独使用 `U_6`；这是六层精确证明中的前五个粗层。
- [sharp.tex:726] 的 `f_5` 是六层望远镜中第五层的索引项，与第六层的末层消去相邻，不是“五层已代表最终 sharp”的陈旧陈述。[sharp.tex:697, 726–732] 同段同时标明全部六层及末层误差。
- 六层锐化结构的标识见 [sharp.tex:1, 31, 77, 641]；摘要/引言的总述见 [main-reviewed.tex:48, 110]。

## 过程文字与误报排除

未发现 `R018/R0xx`、运行日志、重算、SHA、TODO/FIXME 等轮次过程语句落入这四份正文。`sharp.tex:386, 412, 418` 的“候选数/候选根”是几何计数中的数学候选，不是研究状态或工具过程文字。引言“截至 2026 年 11 月”的表述位于文献综述语境，不是运行记录。

## 结论与边界

本次静态检查未发现标签、引用、环境、控制字符或明显宏定义问题，也未发现指定陈旧措辞和过程日志混入。它只反映上述 SHA 对应的四份快照，不证明 TeX 可编译，也不构成数学证明审查。

主代理最终补检：main-final仅上述排版变动，数学标签不变；正式四源合并共有138个唯一label，160次ref/eqref全部有定义，静态控制字符检查通过。编译/视觉验证另见本轮报告。
