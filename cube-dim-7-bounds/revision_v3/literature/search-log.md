# 有限文献检索记录

日期：2026-10-27，Asia/Shanghai。使用web搜索与公开原始论文/作者固定仓库、Zenodo API；未声称调用MathSciNet、zbMATH或Google Scholar数据库。依据academic-literature-search技能记录核验身份、版本和实际阅读范围。总计12条search_query，达到父任务明确上限后停止新增查询；随后只打开已定位原始链接及内部账本的明确来源。

## 实际查询

首批（4条）：

1. `Stein Strömberg 1983 maximal functions n log n cubes weak type`
2. `Aldaz 1605.05596 Theorem 5.2 2 3 epsilon n log n`
3. `Aldaz 1811.10478 sharp 2^d uniform measures Lebesgue`
4. `Melas centered Hardy Littlewood maximal operator sharp weak type constant 11 sqrt 61`

跟进第一批（4条）：

5. `"centered" "squares" "weak" "upper bound" maximal two dimensional 4`
6. `"Stein Strömberg" "Covering" Aldaz 2017 5.2 journal`
7. `"centered" "cubes" "n log n" "2026" upper`
8. `"Behavior of maximal functions" "1983" pdf Stein Strömberg`

跟进第二批（4条）：

9. `"centered" "cube" "maximal" "coefficient" "2026" weak`
10. `"Stein" "Strömberg" "259" pdf`
11. `"centered-maximal-constant" "3.879"`
12. `"Aldaz" "1605.05596" "correction"`

## 原始来源路径与版本结论

- SS1983：出版方DOI页面核作者、年份、卷页和题名；清华YMSC原刊扫描读取§2。已有本机原始PDF缓存被复用，SHA256重新计算。出版方页面与原刊题名均有functios笔误，Bib保留。
- Aldaz1605.05596：arXiv abs显示v1=2016-05-18、v2=2017-01-07，现公开最新为v2；v2说明minor corrections、与JMAA449(2)1741–1753正式版交叉对应。出版方DOI10.1016/j.jmaa.2016.11.070核身份。v1/v2原始PDF均下载，定理5.2/证明相关页对读；未获得正式期刊独立PDF逐字比对，也未发现新的改版或勘误记录，不声称穷尽勘误。
- Aldaz1811.10478：arXiv显示v2=2020-03-10，Cor3.5和Thm1.2原文读；不移植未知正式版本信息。
- Melas：Annals正式页面确认2003、157(2)647–688及DOI；原始arXiv v1 PDF pp648–650读，Cor1常数及sharp原文明确。未审核42页证明全文。
- Lin：搜索定位作者仓库和作者研究页，打开固定72c022ba09032b522c878cca2f9b6c68217d77ea版本README、Challenge和Solution；其作者PDF与另一篇partial-balayage PDF实际下载。当前无版本日期的作者HTML与PDF对3.615749证明状态有差异，分别记录，不以旧HTML否定新版PDF，补充exact证书未核。Palomar v3机器记录访问403；无新Lean build。
- Tomar：最新内部cube-lowdim账本指出公开更强来源，故沿明确链接打开固定bca281a706bbf6375b04d601ad19ba4ad0bf28b6仓库README。Zenodo API23554628核作者、题名、DOI、v4、2026-10-22与is_last=true；不由README凭空构造DOI。作者站PDF403，改从API返回的官方PDF内容链接成功下载。12页原始PDF实际读pp1–5、8–10、12；Theorems1.2/1.3、Theorem2.2模型和全尺度/全L1延拓、核验范围及原始数据出处均明确；未重新跑数值核验或64CPU小时Lean构建。
- Spector–Stockdale：arXiv2609.05377v1元数据及原文pp1–2明确Euclidean balls；仅用来排除错误跨模型比较，不当方体上界。

## 远端及计算范围

独立clone最初由partial本地clone复制，出现缺失blob提示。随后设置独立origin promisor配置，实际fetch成功，origin/main=a230c14f20e15ad6b7c85ac021d337417298251f。只读提取6个相关项目的STATE及结构化账本，再沿二维/三四维原始证书索引核对。快照原文和选取结果在internal_snapshots与internal_baselines.json，不修改原checkout、不启动原runner。

本机Python只用于PDF文本提取、文件SHA256、逐维现有完整有理端点转写与比较。没有新数学参数搜索、优化、数值实验或本轮证书复跑。已运行的历史程序只按原项目接受记录叙述。

## 覆盖局限

不完整cited-by、不完整数据库检索、不排除未索引更新。检索定位的聚合页面、第三方AI追踪页面只作导航，不作为数学结论来源。最新网页当前状态均按2026-10-27读取；未核到同模型已有系数1或更好公开渐近上界只是本次有限检索结果，不支持全球首创/最好。文献列表不以15篇凑数。
