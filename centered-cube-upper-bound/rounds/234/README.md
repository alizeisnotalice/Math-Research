# 第234轮：分数几何与全亏损三次对数来源预算

**一般一致弱型上界仍未闭合。** 本轮在原far正交通上，将第233轮的ρ^(2/3)来源拥挤费用及外部√K幂费改进为[1+log_+ρ+log_+K]^3，并把共同soft菜单降为单一O(√n)。全终端亏损、实际合法物理节点和原来源均保留。

## 核心结果

ρ=||ν||/(qa^n)，C0=16/3，R=1+log_+ρ+log_+K，θ=1/(2R)。θ依共同外层来源数据固定，不能逐receiver免费重选。原good比较变为(K/3)^θ；同capture下Jθ≤J0^(1−2θ)Jhalf^(2θ)≤M^(−n)。物理覆盖仍付e⁴。

对一般仿射乘积、非负指数总和≤n，arcsine网格N=ceil(π√n)满足supH≤eΣH_grid，包含零系数及端点。两核共享soft参数，只用同一菜单。

原heavy增长分母与Gini列界给每节点取两价min。对所有j及正规M节点、唯一封顶点求和，保θ趋零的三重费用，得到

\[
\mathcal P_{\rm far}(U)\le12\Lambda e^5\lceil\pi\sqrt n\rceil
(9C_0e^{3/2}+2\sqrt e)\,
\mathbb E\sum_V\|\beta_V\|[1+\log_+\rho_V+\log_+K]^3.
\]

U≥0；ρ=0交通为0；须同V的ν,q,a,β跨j共同，否则保逐j价。12Λ保原前缀、三腿、两朝向。另保实际分数affinity来源量及同接收位置的有价FIRST容量桥，未抹去backward-soft费用。

**剩余缺口：** 原共同历史下的β加权三次对数来源矩、其余原分支、来源到输入质量桥及最终回代。R≤B log(n+2)只给条件性已付来源族，不能代表任意L1输入。下一步继续利用全局极值几何结算平均空间/历史费用。

## 验证与复现

独立审查先筛命题再读证明。主脚本实际执行12种真实核配置、10个独立mp积分点、9种一般θ重时钟、18组实际节点、35种全节点三次对数价、65种arcsine族；40/70dps复核及小n有理物理节点对照。n至10^6、logρ至1000；大节点用解析单一crossing与有限几何和，保j正遗漏尾。有限全价最大归一上界约0.00232379，不能当作一般锐常数。普通高精度不等于区间证书，有限数值不代替一般证明。

[research.zip](research.zip)保详细LaTeX文本、作者推导、独审、数值脚本/结果和manifest；不编译PDF。本轮仅新增本README与ZIP，另更新项目索引。应用F03/F05技能核对几何列、正算子与谱/积分条件，未把技能中的未证外部主张作已证工具。无可归属本轮的实际额度扣减数据。

解压后运行`python3 reproduce.py --verify-only`；完整复现使用已装numpy2.3.5、mpmath1.3.0，运行`python3 reproduce.py --output-dir <新建包外目录>`。驱动在包外复制脚本与合同，重算完整JSON逐字节比较，运行前后检查manifest。

冻结归档：11成员，75,415字节，SHA256 `37c5b7072a4902958d424762b6db20ba1cc258848e6231bf3834af9aabb94aef`。包外完整重放约7.196秒，452,457字节JSON逐字节一致，10项manifest运行前后通过。独审另实际完成27+433+128项独立检查及作者容量附件300项；作者421项按精确生成规则单列，未冒称独审重跑。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round234_20261017_final_checkpoint`；ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round234_20261017_research.zip`；重放回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round234_20261017_final_replay/receipt.json`。全部成员与备份逐字节校验通过后移除本轮work。
