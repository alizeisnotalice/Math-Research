# 七维公开文献门槛核查与 A7 算术证书

本轮在固定模型——Lebesgue 测度、中心轴平行立方体、共同半径、一般 L1、所有 r>0——中核验到

\[
A_7\le c_7\le 2^7=128,
\]

其中令 \(t=2^{1/7}\)、\(s=(3/2)^{1/7}\)，则

\[
A_7=\frac32\left(\frac{1+t}{2}\right)^6
-\frac12\left[\left(\frac32(1+t)-2s\right)^6
-\left(\frac{3-t}{2}\right)^6\right].
\]

由 [Aldaz 2000 Proposition 1.4](https://dml.cz/bitstream/handle/10338.dmlcz/127554/CzechMathJ_50-2000-1_14.pdf)（printed p.110）得到下界；有理算术验证

\[
1.77526292548200472722205688741524226613926614753030
<A_7<
1.77526292548200472722205688741524226613926614753031.
\]

上界采用 [Aldaz arXiv:1811.10478v2 Corollary 3.5](https://arxiv.org/pdf/1811.10478v2) 对一般测度的统一 bound，再专门取 Lebesgue。正确编号为 3.5。原文对统一测度类的 sharpness 不表明七维 Lebesgue 常数等于 128。

通过维数单调性，可以继承低维公开下界。本轮比较了 A2 至 A6、Melas 一维常数，以及 [Lin 作者仓库](https://github.com/CoolRmal/centered-maximal-constant/tree/72c022ba09032b522c878cca2f9b6c68217d77ea) 的二维新下界 Φ=1.6855099933…。它们均低于 A7。A2 至 A6 的比较由同一有理 checker 严格核验；一维与 Φ 的比较使用原始论文/作者 README 记录的明确阈值。

[Lin 的作者网页 §4.2–4.3](https://coolrmal.github.io/articles/two-partial-balayage-principles/) 只把新 cube upper 用于二维 squares；高维公式适用于 Euclidean balls。[Spector–Stockdale 2026 v1](https://arxiv.org/abs/2609.05377v1) 的 O(sqrt(n)log n) 同样是 Euclidean balls，不能直接作为七维 cube bound。

**限制**：本轮为 12 个首批定向检索加 1 个登记补查，并非穷尽检索。Iakovlev–Strömberg 原始出版社全文遭遇 403；仅取得其摘要片段，未完成有限维参数及七维显式常数审查。Aubrun 原始全文及有限维常数也未闭合。因此本轮不声称 A7 或 128 是全球最好公开门槛，也不声称不存在更强公开七维界。内部 n6 R075 约 2.073015624830795685 未算入公开文献结果。

严格算术部分只验证指定公式的根区间、A_n 区间及比较，不重新证明原论文的构造定理。独立运行的 SymPy 和 Arb 数值结果作为交叉计算证据保存在 computation-receipt.json。原文 PDF 和网页缓存位于仓外，仓内 source-receipts.json 固定 SHA-256/版本与实际阅读范围。

复现：`python3 check_a7.py`（标准库即可）。应输出最后一行 `PASS: exact rational root brackets, A2..A7, comparisons, and 2^7=128`。

文件：`public_bounds.json` 为机器可读结论；`a7-certificate.json` 为全部有理端点；`check_a7.py` 为最小 checker；`source-receipts.json`、`search-log.md`、`references.bib` 为来源索引。未 commit、push，所有新增归档文件均位于本目录。
