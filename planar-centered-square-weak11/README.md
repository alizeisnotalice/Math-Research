# 二维中心正方形弱(1,1)常数

固定Lebesgue背景、全部正尺度、轴平行中心正方形，一般L¹输入。

截至第3轮，经本机精确核验和独立复审的严格区间（端点沿用第2轮）：

\[
1.70275314637189487745\ldots\le C_2<3.590867754339.
\]

第2轮在已复现的公开50原子周期种子中加入一个质量1/1000的有理原子，两套不同精确面积程序输出相同新分数；有限截断和平滑已严格传递至一般L¹。第1轮上界利用公开核的统一主导裕量作精确正缩放，经Arb全域核验与非作者复审。仍未确定锐常数、匹配上下界或全局L¹达到性。

[第3轮记录](rounds/003.md)给统一线性删除的严格反例、有限证书第一变分与mixed见证必要性，本轮没有收窄数值区间。[第2轮记录](rounds/002.md)给当前下界、绝对差/比例和实际验证；[第1轮记录](rounds/001.md)保留原子基础、双原子最优及上界认证的冻结历史。

| 已证 | 数值观察 | 猜想 | 已否定 |
|---|---|---|---|
| 51原子周期证书及有限严格L¹传递；原子化等价；固定N可达/代数性；双原子锐值3/2−√2/4及等号；局部删除/阈值/有限对偶；认证一般输入上界 | 搜索定位与资源采样只支持所列有限运行 | 本轮未提出精确C₂的猜值；统一有限结构和强紧性保持开放 | 接收区相交即可无损合并；独立优化两个方向半径；由一维有限原子不达推二维L¹不达 |

## 数学证明与独审

- [第3轮证据与复现](round003/README.md)：有理概率原子/直接L¹删除反例，固定有限族完整第一变分，实际混合见证机制；均有非作者独审。

- [第2轮新下界与严格有限化](proofs/round002-report.md) / [独审](reviews/round002-audit.md) / [精确输入、日志与复现](round002/README.md)。

- [原子模型、固定N](proofs/atomic-foundation.md) / [审查](reviews/atomic-foundation-audit.md)。
- [双原子全局最优与删除结构](proofs/n2-geometry.md) / [审查](reviews/n2-review.md)。
- [达到性与覆盖对偶](proofs/attainment.md) / [审查](reviews/attainment-audit.md)。
- [公开下界与严格L¹恢复](proofs/lower-baseline.md) / [审查](reviews/lower-baseline-audit.md)。
- [核缩放证明](kernel-rescaling/PROOF.md) / [最终复现包审查](reviews/kernel-rescaling-review.md)。
- [Melas方法迁移](proofs/melas-transfer.md)：明确所读原文范围；局部一维拼接不作为新的二维下界。

## 复现

[kernel-rescaling/README.md](kernel-rescaling/README.md)给完整独立便携核验与可选全C4重跑；默认核新旧C1/C2、独立积分及旧C4运行身份/覆盖，再用正线性传递。旧汇总不是逐叶数值证书，新核未作Lean构建。

双原子使用Mathematica15；`checks/n2/n2_vertices.wl`实际以15个不同子核运行21项全参数Reduce。一般位置到21点的几何约化单独证明。`checks/atoms/atoms_examples.py`以SymPy1.14核三原子真面积与严格有理内证书；该例不宣称三原子全局最优。

`checks/lower/exact_check.py`和`lower_check.py`为两套不同的公开精确面积实现；后者需NumPy。输入为`checks/lower/exactinputs/c2lower_1p70272.json`，原始输出在同目录。`checks/lower/exact_comparisons.py`与`check_transfer.py`仅需Python标准库，廉价核传递与身份，不重新枚举完整halo。其余算核按各脚本与配对输出复现。

第三方上游代码均保留许可证、来源与提交信息。仅归档必要证明和复现资料；原始用户任务、未认证内部材料、运行环境、外部完整Lean项目及可再生缓存留在本机。[公开派生映射](evidence/publication-bindings.json)区分原冻结文本与仅作非数学文字处理的公开副本；历史绝对路径只是运行元数据，便携核验不依赖它们。

[参考文献](references.bib)。本目录的核族上界不等于全局精确常数。
