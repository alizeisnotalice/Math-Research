# 检索记录 2026-10-25

实际使用web搜索、出版方PDF、arXiv版本页、作者网页、GitHub原始代码与git ls-remote。
检索式：`centered cubes maximal weak type constant dimension three four upper bound 2026`。
打开用户给出的五个入口；Aldaz 2011期刊PDF确认定义及第1013页的维数单调说明（引用AV07 Th.2）。arXiv 1811.10478当前v2，2020-03-10。
Tomar作者论文页面直接打开失败；GitHub API遭403速率限制，改用固定提交raw下载成功。Zenodo DOI入口本次工具访问失败，故文章版本日期、正式定理编号仍待原PDF补核，不用README代替论文全文。
Tomar HEAD实测仍为bca281a706bbf6375b04d601ad19ba4ad0bf28b6；本机整数checker复现准确分数，见public_lower_run.json。未重建Lean、未运行第二checker，不能宣称双实现独立确认。
Lin网页读取第3节Principle 2及障碍分解假设；网页无明确版本日期，记录访问日期而不伪造发布日期。只作为待审方法，不据此新增三四维界。
本次未定位更强的专门三四维公开上界；此结论不代表穷尽检索。球、Poisson变差、离散和强型结果排除；搜索汇总页仅作线索。

本地PDF后备读取成功：pypdf提取Aldaz2000纸本110–111页，Prop1.4与给定公式一致。pdftotext缺失与pypdf字体警告已记录，未上传原PDF。
