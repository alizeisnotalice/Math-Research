# R048 最终版式差异与完整组装独审

日期：2026-11-05。审查者：Sol，独占版式接口审查目录。

**裁决：PROOF-VALID / LAYOUT-EQUIVALENT。12 个版式块未改变数学内容；最终模块可准确逆还原到已独审数学源；完整四文件引用接口闭合。无需数学或接口修补，可交实际复编。**

本报告不重复审查 eighth 的完整数学证明，不以作者的 canonical 程序替代数学差异审读。未改源或仓库、未编译、未宣称此前 52 页 PDF 的 11 处 overflow 已修复；最终排版和渲染验收由后续实际复编提供。

## 1. 精确冻结绑定

编译源快照：

`/Users/zhengzhihao/Desktop/T/topic1-tauberian-work/scratch/r048-compile/project`。

以下 SHA256 均为本次实际读取文件并求得：

| 对象 | SHA256 |
|---|---|
| 编译 main.tex，亦为已独审数学 main | `b976cf01d89b73ff2ec540a19aebc5efbc44848cd2d2fc6bf271667eb4eb63f5` |
| 作者冻结数学 sections/eighth.tex | `f19499528e7fe5d322cf2a2d0333d7743c2fd6ffa3ac059eaa5798adc34b98cd` |
| 编译版式 sections/eighth.tex | `447b7b443a3e64d17a10069df6921b2b8398937ff6273704c106def5a4ba4229` |
| 编译 sections/global.tex | `7b6154f74b56e4cb8ea133c99d4c240ceec526004930c510bc29ebcb25b5f10a` |
| 编译 sections/sharp.tex | `8eae7fa41e8bc1f9883d51330b01183001295ca0b623da12205c807c6bbe7c39` |
| 编译 sources.bib | `5ec7a70079799b4754ebf1b230d5600f55f8b2bc51a5d457ad8524e7c8e18e8a` |

main 与 `scratch/r048-main-typeset/main.tex` 精确字节相同；global、sharp、bib 与正式 repo 对应源精确字节相同。没有把主定理或参数传递的数学更新混入版式补丁。

## 2. 十二个块的独立差异审读

实际完整读取 `scratch/r048-layout/replacements.json`、`layout.diff`、`layout_patch.py`、`verification.log`。逐块 old/new 的数学判定如下；标识采用 JSON 中稳定 ID，行号不是新文件的定位依据。

| 块 | 变化及数学判定 |
|---|---|
| W040 | 临界序列的一个 equation 内加入 aligned 并分成三行。v→0、正测度 V、紧集 K、M 定义及含 B_0 的下界完全不变；标签仍为 eighth:critical，条件仍为共同合取。 |
| W123 | 两个行内粗界移入同一无编号 aligned 显示；既有顿号和逗号改由 text 排字。两个求和指标、ε/α/Δ/R 项及 O(v^5S)、O(v^6S) 全部保留；“不预先假定……”仍接同一段落。 |
| W274 | 半径/横宽乘积等式与质量根等式分为两行，对齐等号。求和、乘积指标、R 次幂、误差 7/2 和 6、M'+D、a 与 V/A_v 全部相同。 |
| W298 | z、x、y、ζ、widehat b、b 的六个定义分三行。clip 的两个截幅端点、v^3R_e 分母、(1+v)/n、横权归一化及 min(1,widehat b) 全部相同；没有改变变量相依顺序或定义范围。 |
| W408 | C_i(K) 与 Z(K) 两定义分行。集合中“至少两个轴”的条件、阈值 Kv^3r_i、差集、绝对值、v/a、求和全部不变。 |
| W530 | 同见证源层高度约束、A_ij(y) 定义及误差和分行。s≥ΣA−CR_iE_i、退款项减号、正部、积分质量分子、naR_i^(n−1) 分母、A_ij≥0、O_T(v^6S) 全保留；没有改变关系方向或非负性。 |
| W638 | 二轴联合价格分为四行。左侧 H_2+c_RZ+c_Dv^6D_N、右侧 [A_6+Γv^7]C+χF+Cv^8S、H_2 中仅减一次 ΣJ_i，以及 A_6 定义均完全相同。原 equation 及其标签不变。 |
| F704 | 求和下标的 `{\rm 邻边界}` 改为 `{\text{邻边界}}`。中文指标文字相同；仅从数学罗马字体切换为正文文字字体，求和和下标分组均保留。 |
| W816 | 边界局部小源估计由行内移为无编号显示，既有句号用 text 排字。(v/a)Σδ_i^mass≤CavS(∂I)=O(v^(n+6)S) 的因子、指标、指数和关系全部相同。 |
| W925 | 实际有限构造的 s、r、F、f、h、δ 六个定义分三行。无穷范数、ρ+s、(2w)^(n−2)、n 次幂差、导数 F'(0) 及两个商都相同。 |
| W942 | R_U(x) 和坏集 B(U) 定义分行。k≠0 的和、Tk+(C_w\U)、同一 Q(x,ρ+s)、s=0 或 R_U≥F(s) 的条件全部相同；集合定义与逻辑“或”不变。 |
| W1052 | P、δ、T_3 的统一展开分三行。全部 B_n、κ、8n/4n/32n、n−2 因子、u 与 g 的幂次，以及各 [1+O_n(u^5)] 完全相同。 |

除 aligned 包装、& 对齐点、换行、空白、行内/显示包装以及相同中文标点/指标的字体处理外，没有数学 token 改动。所有 equation 包装及 label 仍在原处；新增的显示都是无编号显示，没有新增编号或定理环境。排版改动不改变主定理的正 Θ8、ordinary O8 或自然参数传递的口径。

## 3. 独立逆替换：不是只采信 canonical True

另写并实际运行只读 Python 检查，不 import 作者 layout_patch，不调用其 apply/reverse 模式，也不写回源文件：

1. 直接读取最终模块字节与作者冻结模块字节。
2. 逆序遍历 12 块，逐块断言 new 在当前内存文本恰出现一次，然后准确替换为 old。
3. 将最终内存文本编码后，与作者冻结字节直接比较，同时独立计算 SHA256。

实际输出：

```text
author module SHA: f19499528e7fe5d322cf2a2d0333d7743c2fd6ffa3ac059eaa5798adc34b98cd
final module SHA: 447b7b443a3e64d17a10069df6921b2b8398937ff6273704c106def5a4ba4229
independent inverse exact byte match: True
inverse SHA: f19499528e7fe5d322cf2a2d0333d7743c2fd6ffa3ac059eaa5798adc34b98cd
```

因此这些块覆盖最终源与冻结数学源之间的全部差异，未残留未记载的数学更改。作者 verification.log 的 canonical 相等只作为辅助信息，本裁决依赖逐块实际审读及独立字节恢复共同成立。

## 4. 完整组装引用与已删命题

对编译快照 main/global/sharp/eighth 四文件实际扫描 label、ref/eqref/autoref 和 input，输出：

```text
labels: 188 references: 203
duplicate labels: {}
unresolved refs: []
inputs: ['sections/global', 'sections/sharp', 'sections/eighth']
eighth theorem label: ['sections/eighth.tex']
removed token prop:nonexact 0
removed token eq:nonexact 0
removed token cor:nonexact-multi 0
removed token seventh: 0
removed token planar: 0
removed token skew: 0
layout refs identical: True
layout labels identical: True
```

完整组装为 **188 个唯一标签、203 处引用**。此前 main 独审的唯一预期待供接口 `eighth:theorem` 已由模块正确供应。没有重复、悬空标签或指向已删除模块/弱命题的引用。最终与冻结数学模块的引用序列、标签序列分别精确相同，排版补丁不改变交叉接口。

旧 prop:nonexact 已被实际八阶有限下界覆盖，其陈述和补证不再参与本次组装；零残引与前次 main 的数学依赖审查一致，没有因删除产生缺失的基础前提。

## 5. 验收边界及最终签署

已完成的是源差异的严格数学等价判断、逆恢复与 SHA 绑定、四文件完整静态接口，以及 main/global/sharp/bib 未变检查。模块 f194… 的完整数学独审与 main b976… 的主定理/参数传递独审继续承担相应数学证明，本报告只把它们准确接到最终排版源。

**本次未运行 TeX 或查看最终渲染。** 先前 52 页且有 11 处 overflow 的 PDF 不属于最终验收结果；本报告不声称这些溢出已经在最终 PDF 消失，也不声称已编译通过新源。应由 Terra 在本裁决后实际复编，并据新日志/页面判定排版修复是否完成。

**最终签署：main `b976cf01d89b73ff2ec540a19aebc5efbc44848cd2d2fc6bf271667eb4eb63f5` 与版式 module `447b7b443a3e64d17a10069df6921b2b8398937ff6273704c106def5a4ba4229` 的数学/静态接口为 PROOF-VALID；12 块版式调整相对数学 module `f19499528e7fe5d322cf2a2d0333d7743c2fd6ffa3ac059eaa5798adc34b98cd` 为 LAYOUT-EQUIVALENT，无需数学或引用修补。最终编译与渲染尚由下一步完成。**
