# 第615轮：保源停止的弱型目标与多包接收缺口

一般 n^(1/2+o(1)) 一致弱型上界仍未证明，本轮没有新增一般维数上界。

核验一个保源停止尝试：固定坐标循环二分产生长宽比≤2的矩形，父子体积比2；以λ/8停止得到0≤g≤λ/4、Σ||b_Q||₁≤2W和各纤维坏包均值零。保持同一实际尺度、先合并所有坏包后取sup，仅需其正的有符号接收B。准确有B≤S_m f以及{S_m f>λ}⊂{B>3λ/4}，但多个来源不交的坏包，其接收边界仍重叠。这一联合弱预算未付，分解本身不算改善。

一个给定边界剖面的强积分处在m log2/8与m log2/4之间，弱L¹拟范数却小于1/8。故强费用线性不能否定弱型策略；也不能反过来把孤立剖面的弱常数当作任意输入卷积或坏包和的弱界。停止扩大强反例，保留真正的共同尺度多包接收问题。

原始文献定向核验：de Guzmán产品基停止用到非中心含点集合性质；Fava算子复合仍给LlogL来源费。二者都不提供当前纯W联合付款。SOURCES.md含准确页码、适用限制、原始链接、阅读日志和BibTeX；不归档论文重复副本。

实际计算：指定Python3.12.14、SymPy1.14.0、python-flint0.9.0核m=2/4/8/16/32/64。两条精确积分路径一致，主代理重放稳定字段一致；另用有理Beta展开及正系数体积展开独立复算六组结果。R=2只作一个弱水平下界，不冒称最优；全部阈值弱界依解析证明。

[研究包](research.zip)含详细LaTeX、独审、必要计算与复核、来源说明和SHA256SUMS。解压复现：

```sh
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python numeric.py replay.json
/Users/zhengzhihao/.local/share/math-research-tools/venv/bin/python verify.py RESULTS.json
```

重放比较仅排除command、runtime_seconds、peak_rss_bytes_macos。不编译PDF，不归档重复草稿。团队为1个Astra与3个GPT-6.1 Sol。

ZIP SHA-256：`031680a9b34b94ed448b21e1457ff0200b9060712a0eec1e670ff3762ba27fb0`。
