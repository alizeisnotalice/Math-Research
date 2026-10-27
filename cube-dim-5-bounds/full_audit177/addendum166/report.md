# 177第一轮增补：166/169及主完整常数资格

**增补通过。** 169已补齐原177冻结报告所列的166有限算术资格门；主22326新完整常数的来源、单位、组合与有向隔离代码接入正确。现可在第一轮审计范围内接纳

\[
C_{\rm cert}\le C_5\le U_{166}
=\left(3125I_{166,U}^4/g_{192}^5\right)^{1/4}
\approx28.936078340172045,
\qquad C_5\le\frac{28936079}{1000000}.
\]

完整准确值以主`round019/primary_upper166_acceptance.json`为准。这里只是177第一轮的追加资格，不是178第二轮严格全链终审。旧H/R/梯度/log等严格区间真值和解析来源的继承限制，仍按原177第6节，不被新SHA检查消除。

## 169完整审读结论

已逐文件审读169的`theory.md`、`prepare169.py`、`bind169.py`、`integer169.py`、`audit169.py`、`io169.py`、精确plan、source/candidate绑定和真实执行回执。

- `prepare169.py`从163原binding及162真实新mean、165已核164斜率按同key重建独立15份shard，逐字段对照166 metadata；原112 range、mean/self、156保留值、真实chart与原树均绑定。162/164全部原输出另与各自独审核准的candidate锁匹配；没有用两个未来标签互相证明。
- `bind169.py`绑定166真实440件输出、390 chunks及97971记录，锁生产plan/identity、原ready receipt和163/165实际结果。准备中的`prepared_not_authorized_not_run`只描述冻结准备状态，实际运行证据来自唯一run001，二者未混淆。
- `integer169.gain`交叉gcd消因子后，以二进制shift/subtract长除得到q，并核`qB+r=A,0≤r<B,qB≤A<(q+1)B`。此证明确保准确floor，包含delta=0及A<B边界；不同于生产直接乘除，不导入166数学代码。
- 每项唯一总gain从原112 self扣一次；真实mean仍落原range，原164准确弦斜率全range有效，N=2^160抵消正确。selfL、两候选支配及逐条输出字典完全匹配；没有叠加已舍入的两个旧gain。
- `aggregate`自根独立递归DFS，检查环、父子层级、唯一保留self、97972节点及52360根；逐节点`min(self,childsum)`并与全部166 effective_nodes精确匹配，只计根和。原几何与orbit/完整分割继承原已合格门，不冒称新生成这些旧基线。

169实际master **21868**，保存master exit0，15 worker **21877—21891**全ready后启动，保存每个returncode均0；内部2.600273334秒。真实完成97971新gain、97972节点/52360根；输入前后锁相同，没有新factor/mean/slopes/H/R/J/power或旧候选DFS。结果SHA准确为
`75ac91b7c743dad1e705aa4b60dfce39b2ea8d033eb6becd011afb74e0dabe59`。

\[
UH_{166}=\frac{3994994121862399564161069833676547838746859}
{91343852333181432387730302044767688728495783936}.
\]

此值同时匹配生产166、169和主接纳字段。唯一nonactive节点保留self但仍允许childsum改善，祖孙未双计。169报告不自行合成C5，其资格范围是同模型完整F的H512积分上包。

## 主完整Cn接入

已全文审主`accept_upper166.py`、`upper166_registration.json`和`primary_proof166.md`。登记逐SHA绑定accepted156、166真实completion、169实际结果、原177报告及主证明，先核169成功schema、15退出码和完整计数，再取

\[
I_{166,U}=UH_{166}+J_{139,U},\quad g=g_{192},\quad
U_{166}^4=3125I_{166,U}^4/g^5.
\]

完整J139仅支付一次，g192和小时间/奇异集/全尾与原177相同。H128 floor仍只在F；完整G经−4齐次推广到Q。原全空间质量1正有限源、同一接收域、一般L1和全部共同正半径的定理级适用范围保持原177，未增加尺度或背景假设。

主实际PID **22326**仅隔离一个新完整正四次根。源码的双isqrt给floor第四根整数，随后准确检验`lo^4≤power<hi^4`；新hi严格低于已存旧lo，百万分位ceil给便携上端。原旧根区间直接读取而未重算。成本区间的下端不是C5下界；原下界未动，无双侧突破、后三轮仍0/3。

本增补实际执行`inspect_addendum.py`，PID **23040**、exit0、0.723570458秒：逐SHA核169的1047静态源、446实际candidate锁、15独立shard与原树，锁1974个真实文件/511819693字节；核原177全部冻件完全不变，并只比较新的保存字段（同UH/g/J、I加法及准确下降、已存新旧区间与便携上端）。本增补**没有重跑任何gain、DFS、mean、slope或根**，不把hash当数学真值证明。

原177报告和原manifest一字未改。增补的完整来源及实际结果分别在`source_manifest.json`、`inspection_result.json`；`inspect_addendum.py`拒绝覆盖原输出。必要旧依赖、第二轮178、论文终审与主责任仍明确分离；178的区间reader结果未提前纳入本增补。
