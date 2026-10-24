# Tomar公开基准与版本检索

检索完成时间：2026-10-24T10:26:50.036792+08:00。只读联网与现存公开文件审读；未运行证书核验器或搜索。主侧另行完成的公开核验器复现应引用其实际日志，不把本次文献检索当作证书复跑。

## 核实版本与固定来源

[作者详情页](https://www.mefyi.com/papers/weak-type-squares/)当前题名为 *Bounds for the weak type (1,1) constant of the centred maximal operator over squares*，作者Tejasvi Singh Tomar，标注预印本、未经同行评审，首次发布2026-10-21。作者papers索引的搜索缓存仍用旧题名“An upper bound ...”；不能据索引判断最新版本。

[论文Zenodo全部版本API](https://zenodo.org/api/records/23554628/versions)实际返回4条；[latest API](https://zenodo.org/api/records/23554628/versions/latest)仍返回v4，2026-10-22：

|版本|论文记录|附件记录|
|---|---|---|
|v1（2026-10-21）|[23522055](https://zenodo.org/records/23522055)|[23522057](https://zenodo.org/records/23522057)|
|v2（2026-10-21）|[23533973](https://zenodo.org/records/23533973)|[23533974](https://zenodo.org/records/23533974)|
|v3（2026-10-21）|[23543279](https://zenodo.org/records/23543279)|[23543282](https://zenodo.org/records/23543282)|
|v4（2026-10-22）|[23554628](https://zenodo.org/records/23554628)|[23554629](https://zenodo.org/records/23554629)|

全版本概念DOI分别为10.5281/zenodo.23522054、10.5281/zenodo.23522056，不宜误当单个v1记录。作者PDF实际读取456833字节，MD5为3d228974197a7c6e0a0e97bc091d4bc1，与v4论文API所给文件MD5相同；SHA256为aa32aab64b2d973d3bbcbdd12e309125c28247284c739758c24abfc2bb9e2ae6。因此指定作者PDF确为本次锁定v4。

[附件v4 API](https://zenodo.org/api/records/23554629)明确锁定仓库提交 `bca281a706bbf6375b04d601ad19ba4ad0bf28b6`。[GitHub提交Atom](https://github.com/tejasvi/square-maximal-weak-type-v4/commits/main.atom)实际仅有该提交，时间2026-10-22T07:15:51Z。作者公开仓库筛选页列出无后缀、v2、v3、v4四仓库，未列出v5；本次搜索未发现更新版本，但不声称排除所有未索引或不关联的后续工作。GitHub REST API遭403限流，已用Atom和Zenodo固定提交交叉核对，没有将API失败写成不存在。

## 准确比较基准

v4声明上下界1.7027237335≤C2≤3.59123，并提供上下界的程序与Lean工程。该声明显著更新此前仅沿用Lin下界1.6855099933的文献范围。Lean的公开核验证声明与本机实际build不同；本次检索未进行Lean build。

[固定提交证书](https://raw.githubusercontent.com/tejasvi/square-maximal-weak-type-v4/bca281a706bbf6375b04d601ad19ba4ad0bf28b6/cert/cert_h32_a4740_R3320_3p59123.json)读取hash为bd8d91b535e0a91fed777251afa6c4d73d78b1e7d14aee7288e5950d8ad28db3，与Zenodo附件说明一致。α=47/40，R=33/20，h=1/32，c=1977601376167/1099511627776，N=80，共1181个非零对称系数类。

按[格式说明](https://raw.githubusercontent.com/tejasvi/square-maximal-weak-type-v4/bca281a706bbf6375b04d601ad19ba4ad0bf28b6/cert/README.md)，半质量为cαR^(2−α)/(2−α)+(h²/2)Σaij。独立只读提取每一类并按原点1、轴/对角4、一般8的轨道数作Fraction求和，得到

Σaij=−632486554619619/1099511627776，

B0=(92947264679849/36283883716608)(33/20)^(33/40)−632486554619619/2251799813685248。

本机mpmath以80位显示B0=3.59122432750645234527147802179301964401723734213279859915914734303765038867902；这个显示本身不是区间认证，严格成本隔离应引用主侧exact_cost.json。成功目标须比较新成本与B0，而不只是低于359123/100000。

## arXiv及是否已有更低缩放成本

实际请求[arXiv API精确题名关键词查询](https://export.arxiv.org/api/query?search_query=au%3ATomar+AND+ti%3A%22weak+type%22&max_results=30)，查询式 `au:Tomar AND ti:"weak type"`，返回totalResults=0。另查询 `au:Tomar AND all:squares` 返回7项，均非本目标论文；网页搜索也没有目标arXiv记录。因此本次没有核到本论文arXiv编号，不能为其补造编号，也不能把未命中提升为“不在arXiv”的绝对结论。

已审v4论文第7节、仓库README、cert格式、search README和相关搜索日志/程序中缩放、归一化及3.5911/3.5909等关键词。本次未找到公开声称把最终认证K统一乘100000/100001或100000/100009后取得更强已认证成本的记录。搜索日志p_P2.log有3.590932等较低LP值，但同一行记录AK约束违例；论文也明确较低有限点LP值只是松弛，并非认证上界。不能以此否定或替代固定v4认证成本，也不能凭有限检索宣称全世界首创。

所以新阶段可审慎表述“严格优于已固定并实际复现的Tomar v4准确核成本”，不表述“超过所有文献”。统一缩放是实质改变全部核系数与质量，不是补写旧数值的小数；其数学适用条件见independent_analytic_audit.md。
