# 合法 local-clock 图复杂度：三轮精确有限守卫

最后新登记批完成一次运行，exit 0，约 0.0724 秒：9 个坐标图、123 个 pair 距离检查、879 个非空 alive states、3,348 个 state/tag 递归与 harmonic degree 身份全部通过，失败零。全表 9,466 个 child recurrence contributions 保存。所有 295 个 terminal tag-state 行对应 induced independent set；425 个 tag 已孤立但其他组件仍 active 的状态完整保留，未误判为整个图已停止。

这些是有限合法 local-clock 实现/复杂度守卫，不计算 realmax、Φ、receiver 体积或 nearmax/geom/history 资格，不提供 √n 上界或 n 阶拟合。n=1、16、256 只给具体 source 坐标与固定半开 grid。

## 1. 完整来源与保守几何图

三轮 N=3、5、8，采用 spine 批相同原 raw weights (1,2,4)、(1,2,3,5,11)、(1,1,2,3,5,8,13,55)，分别除以 7、22、88；完整正来源质量 W=1。weights 用原 source 一次分配，不对不同图、tag 或 subset 相加重复领取 W。新坐标、图及输入 hash 单独登记，不改旧 source/spine/max 数据。

a=1,b=2,d=1/(8n)，各 atom 位于独立半开 cell 的 lower corner。只沿第一坐标放置：

    clique: y_j=j*d*e1,
    path: y_j=(3/2)*j*e1,
    cluster: y_j=[4*floor(j/3)+(j mod3)*d]*e1.

实际 Fraction 坐标按 ||y_i−y_j||∞≤b+2d 生成 undirected edges，无 self loop；不是直接给程序一个预定 adjacency。全部 pair 距离、threshold、edge flags 都保存，随后独立核它们确为 clique、path、按每三点组成的 cliques。123 个 pair 距离判定全部精确。

这里 b+2d 是 root 合同的 conservative 来源邻接图，不将 edge 冒称某个 receiver 的真实捕获，也不省略源 grid/source coordinates。所有 lower corners 位于唯一半开 cell，N 个 source cells 确实不同；同 raw weight 的两个 labels 也在不同物理位置，不是同位原子拆独立时钟。

## 2. 真正 active rates 与 immortal tag

对 alive subset A，degree_A(j)=邻居中仍 alive 的数量，a_j=1 iff degree_A(j)>0。固定 immortal tag i∈A：tag 死亡率 0；每个 non-tag active j 死亡率 1，inactive j 率 0。r=Σ_{j∈A,j≠i} a_j。

r=0 时，如果 tag 有邻居，该邻居也 active，矛盾。因此 a_i=0，整个 induced graph 必为 independent set，状态冻结。反向 independent set 给所有 rates 为零；逐 tag-state 均核验。终态可能有多个 isolated cells，不能把它当仅剩单 cell。

r>0 时 hold∼Exp(r)、在 active non-tag 间 uniform 删除。原来源坐标和 full graph 不动，child 只删对应 source cell，tag 仍 alive。只对 tag 的 cumulative active work 收费，而不是把每段 elapsed holding time都收费：

    H_i(A)=a_i/r+Σ_{active j≠i} H_i(A−j)/r,
    H_i(A)=0 when r=0.

按 subset count 递增，用全部 Fraction child values计算，保存 a_j、active/nonactive labels、tag rate、total rate、first-hold work与每个 child contribution。最多 N−1 次删除，nonterminal r≥1，有限 work 无可积性漏洞。孤立 tag 可以在其他 component active 时继续留下（425 行），其自身 work=0，但过程不因此全局停止。

## 3. Harmonic 身份与 source 一次平均

逐项核 H_i(A)=harmonic(degree_A(i))，harmonic(0)=0。解析解释是 tag 永生时每个 alive 邻居始终 active，到删除前真实 rate恒1；tag 的 active work直到最后一个此类邻居死亡，其他 component 的活动不会免费重设这些邻居 rates。这里数值只核 root/nearflat 正在独审的一般表示的有限实现，不能以本批九图替代任意来源完整证明。

完整 initial source 的平均为 Σ_i w_i H_i(full)/W；是 source-budget 的初始质量加权，不是抽 source 冒充 Lebesgue receiver。数据中也保存每个子状态的 conditional-initial-source 平均供递归重构，但这些平均不再相加当作新的可花预算。

|N|n|图|edges|state/tag 行|terminal tag 行|isolated nonterminal tag 行|完整源加权 work|
|---:|---:|---|---:|---:|---:|---:|---|
|3|1|clique|3|12|3|0|3/2|
|3|1|path|2|12|5|0|8/7|
|3|1|cluster|3|12|3|0|3/2|
|5|16|clique|10|80|5|0|25/12|
|5|16|path|4|80|20|8|27/22|
|5|16|cluster|4|80|17|11|25/22|
|8|256|clique|28|1024|8|0|363/140|
|8|256|path|7|1024|130|190|13/11|
|8|256|cluster|7|1024|104|216|49/44|

按轮 harmonic checks 为 36、240、3072，全过。clique 的 initial degree=N−1，path 的端点 degree1、内点 degree2，cluster 的三点块 degree2、最后二点块 degree1；各 full-tag work 的确切有理值在 profiles/results中。raw source 的大质量末端只改变 source加权平均，不改变给定 tag 的图递归。

这不是一般维数结论：N 与空间邻居数可能巨大，图复杂度与能否由允许 source-once work进入已付 geom 类仍需理论支付。不可把本批 path/cluster 的小值推广到所有源，或把 harmonic 图身份本身当成 Φ/receiver 预算付款。

## 4. 冻结与完整数据

9 个新 gzip profiles 保存完整 weights/坐标/cell index、原 source hash、全部 pair 几何与 neighbor masks、每个非空 A 的活跃状态/degree、每个 tag 的全部 rates/terminal/holding/work/harmonic/child递归。code/registration/results/execution receipt/final_receipt 包含所有 SHA256；spine 登记仅只读对照 weights，未执行其程序，也无旧 oracle/MC。

standalone .py 使用标准库 Fraction，register/run exclusive-create；复现应复制到新空目录，不能覆盖冻结数据。无随机 seeds或抽样，故无概率置信区间，也无 log/数值 time integration。这个同轮最后批已收口，不扩展新图/权重/receiver 实验。
