# 七维公开文献门槛核查，第 001 轮

日期：2026-10-27，Asia/Shanghai。使用 academic-literature-search 技能。任务明确限制为有限首批不超过 12 个定向搜索，并按新线索另登记。研究对象固定为 Lebesgue 测度、中心轴平行 7 维立方体、每次平均使用各坐标共同半径、一般 L1 输入、所有 r>0。c7 是该算子的最优 weak (1,1) 常数。

## 本轮结论的口径

本轮实际核验到可直接引用的门槛为 A7 <= c7 <= 128，A7 的根式来自 Aldaz (2000) Proposition 1.4。精确算术证书给出

1.77526292548200472722205688741524226613926614753030 < A7 < 1.77526292548200472722205688741524226613926614753031。

这不是全球最优公开界的穷尽声明。尤其 Iakovlev–Strömberg (2013) 的出版社全文下载被 403 拒绝，本轮只取得原始出版社 PDF 的搜索片段；未把其渐近 C d^(1/4) 换算为确定的 n=7 数值，更未证明所有可能公开构造都弱于 A7。内部 n6 R075 约 2.073015624830795685 不属公开来源，未混入公开表。

## 已执行：普通网页定向搜索

工具：web.run search_query；搜索返回数量没有提供可靠总数，以下均记“未提供”。语言：英文，不设起始年；Q2/Q4/Q5 直接核查近年。搜索结果中的二手资料仅作导航，数学结论只采用原始论文/作者仓库。

首批（12 个）：

| ID | 精确检索式 | 主要返回/处置 |
|---|---|---|
| Q01 | `"centered" "cubes" "weak" "seven" maximal constant` | 未定位新的七维原始数值界；出现二维汇编，作为导航 |
| Q02 | `"centered maximal" "cube" "2026"` | 定位 Lin 作者仓库与作者网页；登记原始核验 |
| Q03 | `"Aldaz" "1.4" "lower bounds" maximal` | 定位 Aldaz 2000 原始 PDF |
| Q04 | `"maximal" "cubes" "weak type" "2025" "constant"` | 多数为其他算子/测度；未纳入 |
| Q05 | `"maximal" "cubes" "weak type" "2026" bound` | Lin 与其他算子；直接打开作者原始来源 |
| Q06 | `"A remark on the centered" "lower bounds" -site:researchgate.net` | Aldaz 2000 的 DML-CZ、期刊原始镜像 |
| Q07 | `"Iakovlev" "Strömberg" "cubes" "lower"` | 找到 International Press 原始论文 PDF 与目录；新增 Iakovlev 全文线索 |
| Q08 | `"centered" "cubes" "maximal" "lower bound" -"researchgate" -"variation" -site:aiproofwatch.com -site:teorth.github.io -site:palomar-registry.org` | Aubrun、Lin、其他模型；Aubrun 原论文线索未本轮阅读全文 |
| Q09 | `"Aldaz" "Varona" "Theorem 2" maximal` | 定位 arXiv math/0207285；v1 中单调性实际为 Theorem 5 |
| Q10 | `"centered maximal" "dimension 3" "lower"` | 未出现可用于七维的新原始常数 |
| Q11 | `"maximal" "cubes" "1.7" "constant"` | 无相关新原始数值门槛 |
| Q12 | `"Lower bounds for the weak type" "Iakovlev" "cited"` | 出版社某些测度论文（2022）及其他引用线索；直接核查模型 |

新线索补查（非首批扩充）：

| ID | 原因与检索式 | 实际结果 |
|---|---|---|
| Q13 | Q07 的出版社全文 403 后寻找同一原始论文的作者/机构副本：`"Lower bounds for the weak type" "Iakovlev" pdf -site:researchgate.net -site:scribd.com` | 仍主要返回相同 International Press PDF；另出现作者学位研讨会和综述，未用二手资料替代全文 |

这属于定向“被引/引用线索补查”，不是 Web of Science、Scopus、Google Scholar 完整 cited-by 枚举，也未使用这些数据库。没有把未实际运行的渠道写成检索证据。

## 直接打开及实际阅读范围

1. **Aldaz 2000**，`A remark on the centered n-dimensional Hardy-Littlewood maximal function`，Czechoslovak Mathematical Journal 50(1), 103–112。题录 verified。原始 [DML-CZ PDF](https://dml.cz/bitstream/handle/10338.dmlcz/127554/CzechMathJ_50-2000-1_14.pdf) 与期刊 mirror。实际阅读：Introduction、离散测度转 L1 的 Lemma 1.1 陈述、printed pp.110–112 的 Proposition 1.4 陈述与构造证明。read_scope=fulltext（局部全文阅读，不声称逐页通读整篇）。PDF 11 页含封面，公式在零起页 8。根式抄录已与本地 PDF 文本交叉核对；n=2 化为 3/4-sqrt(2)/4+sqrt(6)/2。来源哈希见 receipts。
2. **Aldaz arXiv:1811.10478v2**，`Kissing numbers and the centered maximal operator`，2020-03-10。verified。实际阅读：p.1 Introduction/Theorem 1.2、pp.10–11 Theorem 3.3 / Corollaries 3.4–3.6（printed p.11 为 Corollary 3.5，PDF 零起页10）。正确编号为 **Corollary 3.5**，其 2^d 是对所有合适局部有限 Borel 测度的统一 sharp bound；取 Lebesgue 得 c7<=128。统一测度类的 sharpness 不等于 Lebesgue 情形 sharpness。
3. **Aldaz arXiv:0805.1565v5**，2010-09-06（正式 2011 Annals 173(2), 1013–1023）。verified；另实际打开 [Annals 正式题录](https://annals.math.princeton.edu/2011/173-2/p10) 核对作者、题名、173(2):1013–1023、DOI 与 2011-03-01 网上发表日期。实际阅读：pp.1–2 定义、单调性说明及 lower/discretization 回顾；未通读高维证明全部页。明确 Q(x,r) 为共同半径 ℓ∞ cube，sup r>0；c(d+1)>=c(d)。
4. **Aldaz–Varona arXiv:math/0207285v1**，2002-07-30。verified。实际通过 arXiv 网页全文阅读：Introduction、Theorem 2（函数/测度等价）、§3 Theorem 5（v1 的单调性陈述与证明入口；printed pp.6–7）。Aldaz 2010v5 引用的 AV Theorem 2 不是该预印本 v1 的编号。未自行把预印本题录换成未经核验的期刊版本。
5. **Lin 作者仓库**，`centered-maximal-constant`，git ls-remote HEAD=`72c022ba09032b522c878cca2f9b6c68217d77ea`。verified 原始作者公开代码稿。实际读取 README 开头、Euclidean ball 分离定义、statement 表和 construction；Challenge.lean 原始源码读取 TLS 超时；此处声明依作者 README 中的 statement 表，未声称独立检查其源码。二维 cube lower Φ=1.6855099933…，upper 3.879；各维 cube upper 2^d。通过 c7>=c2 可继承 Φ，但 Φ<1.686<A7。没有从二维 cube upper 推导七维 upper，也未独立重跑 Lean 编译。
6. **Lin 作者网页**，`Two partial balayage principles for weak-type estimates`。verified 作者/题名/原始网页；未核验独立的出版日期，BibTeX 省略未核验的发表年份，note 记录访问日期。实际阅读 Abstract 表、§3 核准对象、§4.2 Euclidean balls、§4.3 axis-parallel squares 及 References。§4.3 还报告一个 AI proposed 3.615749 的二维候选，作者称所需计算难以 Lean 验证；本轮不提升为已核验 bound。高维 (n/2)^(n/(n-2)) 只属 Euclidean balls。
7. **Spector–Stockdale arXiv:2609.05377v1**，2026-09-04，`The weak-type (1,1) bound for the Hardy–Littlewood maximal function is O(sqrt(n) log n)`。verified。实际阅读 PDF p.1 Introduction/Theorem 1.1（Euclidean balls）、p.2 证明策略及文末 cube 相关参考条目。O 常数未显式化且对象不符，不直接给七维 cube bound。
8. **Iakovlev–Strömberg 2013**，`Lower bounds for the weak type (1,1) estimate for the maximal function associated to cubes in high dimensions`，Math. Res. Lett. 20(5), 907–918。原始出版社搜索片段核对到题名/作者/卷页与摘要；read_scope=metadata, snippet_only=true（摘要片段）。直接 web open 与 Python 下载原始 PDF 均 403；没有读取 theorem 正文，未给有限 n=7 explicit constant。未纳入可引用 BibTeX（完整作者名/正式元数据仍只取缩写，全文阻塞的版本证据也不足本轮标准）。
9. **Pérez Lázaro 2022 certain measures paper**，[出版社](https://link.springer.com/article/10.1007/s00009-022-02137-9)。实际打开出版社页面并读取引言的 cube 参考回顾、题名及模型区分，作用为被引补查；针对 certain measures 的主结论不能直接替换 Lebesgue c7。未输出其新定量结论。
10. **Aubrun high-dimensional cubes**：Q08 返回原始论文的搜索引言/Theorem 1 片段，给存在常数 c(ε)>0 的渐近 lower；未获取完整原始版本核验及 n=7 显式常数。仅列待查线索，不做七维数值比较。

## 下界跨维比较

精确有理区间算术已核对 A2<A7、A3<A7、A4<A7、A5<A7、A6<A7。Melas 的 c1=(11+sqrt(61))/12 由 Aldaz–Varona Introduction 与 Aldaz 2010v5 pp.1–2 原始论文明确回顾；其值小于 1.568。Lin 2026 的 Φ<1.686 由作者 statement 表支持。这些通过维数单调性继承到七维的下界均不超过 A7。尚未穷尽所有其他公开维度构造，且 n>7 的 lower 不能逆向下推给 c7。

## 算术复现和限制

运行 `/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python check_a7.py`。该 checker 只需 Python 标准库 fractions/json，先以有理数的整数幂严格夹住 2^(1/n)、(3/2)^(1/n)，再按原文括号做有理区间运算。JSON 保存全部根端点和 A_n 端点。SymPy 1.14.0 与 python-flint 0.9.0 (arb 256 bit) 实算作独立交叉核对，见 computation-receipt.json。数值计算证据不是对原论文命题的重新严格证明；严格算术部分只证明指定公式的区间及大小比较。

网络事件：第一次期刊 mirror TLS handshake timeout；改用 DML-CZ 成功。arXiv 某些请求首次/重试 TLS timeout，但Aldaz 2000、Aldaz 2010v5、Aldaz 2020v2 与 Spector–Stockdale v1 四个核心 PDF 最终缓存成功；Aldaz–Varona v1 只通过 web 原始全文阅读，二进制下载仍 timeout。Iakovlev 出版社 403；Challenge.lean 原始源码 fetch TLS timeout，后续 web 原始源码/元数据批次长时间无返回而主动终止，未计为阅读成功；GitHub API commits rate-limit 403，改用 git ls-remote 得到 HEAD。原文只缓存在仓外 work/cube7-literature-sources，不上传重复 PDF；收据记录字节哈希与版本。

## 尚未执行/下一轮建议

- 完整取得并逐页检查 Iakovlev–Strömberg 原始全文，确定所有有限维参数条件是否允许 n=7，并计算其具体下界；当前不能声称已排除。
- 对 Aubrun 原始版本、Stein–Strömberg 原始显式常数及有限维覆盖 refinements 做独立审核。
- 通过可用的正式 cited-by 数据库对核心原文进行可复核的完整被引检索；本轮仅普通搜索中的引用线索补查。
- 若作者仓库或网页更新，固定新 commit/hash 后另起轮次。未读取私人稿件或把内部 R075 当公开结果。
