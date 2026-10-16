# 第232轮：共同存活来源加权势与几何 score 预算

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮已独审通过一个保共同终端质量、无需统一终端亏损截断的加权存活估计，并接到原 far 交通；尚未支付该势的全空间来源总价。

## 新证明及原积分接口

固定原持久来源 S⊂Z_T、正权 w_t，−dw_t/dlogt=k_tw_t，0≤k_t≤n。令 H=∫S w_T dν、μ=1S(w_τ−w_T)ν，F≥0为有界、静态来源费用。则真实后验（允许S外删除）满足

∫τᵀ∫S F k_t dΠ_t dt/t ≤ Φ_H(F,μ) = ∫₀∞log(1+μ{F>λ}/H)dλ = ∫₀^{μ(S)} F↓(r)/(H+r)dr。

H=0为原零来源支，另行置零。所有费用层共用一个H，未按层重新归一；例如N个等质量损失层的费用由N log2降为log(1+N)。有限端点例展示严格改善，不据此宣称原全几何达到等号。证明、Jensen/Lp价格、原权模型及零分支详见包内LaTeX和附件。

沿R229真实接收域、R231局部硬帽与QG列界，原节点静态费用取 F_M=min{(16/3)/√w_*, B_M exp(−ηΔ_*)}，B_M=||ν||/(qa^n)(a/M)^n，η=25/2^69。取T=min(g*/γ²,g*exp(−U/n))，早段用上述共同势，晚段明确收费B_M n L_late/4。U≥2nlogγ（特别U≥2）时晚段消失。无需限制终端亏损，但浅层互补域并未消失。

另证同一真实后验的加权score界：若E(v)≤exp(−d k(v))，则双后验积分∫∫E(v)(k(v)−k(u))₊dΠ(v)dΠ(u)≤(1−exp(−dn))/(4d)，d=0连续取n/4。原Δ_*≥(1−t/g*)k_t允许保留几何抑制；精确clock积分I与共同势价格Q取较小值。

原三腿总价为 2Λ√(K/3)·24e⁴N_n·EΣ_{V,j}∫β_V(dy)Σ_M min{Q_M(y),B_MI}，N_n≤8+4√n log(n−1)。独立路线在R230来源价格上增加H_{n,η}(U)=(1−exp(−z_U))/z_U，z_U=ηn(1−exp(−U/n))；仅非空原窗0<U<nlog2下使用H≤min(1,2/(ηU))。η极小，不隐去其实际尺度。

**未付费用：**共同势的β积分、全部j/节点/历史重用及H与几何费用的相关性，来源拥挤比、√K，以及near/浅层/终端互补分支。没有把新势命名为常数，没有重领删除或负交通信用。精确截断列只能在认证捕获类且费用静态仅依赖存活变量v时进入本势。

## 独审与实际计算

审查者先独立筛查命题，再读作者和父稿。主96例：18个费用分布、18个原权形状存活/删除clock、60个加权Gini模型；按初探、压力、扩域三阶段，H低至10^-80，clock维度2/16/512，压缩score测试维度至10^24。精确有理层身份与40/70位高精度对照通过。巨大维度是有限score模型，不是原空间构造；浮点对照不是区间证书。独审完整复跑JSON字节一致，另实际执行SO附件244项高精度检查。

冻结后本机Python3.14/mpmath1.3.0包外复跑96例，完整JSON再次逐字节一致，10项manifest运行前后通过，约1.218秒。回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round232_20261017_final_replay/receipt.json`。使用64工具D02/N04/L03核对共同归一化、势费用与数值证据边界；无可归属本轮的实际额度扣减数据。不编译PDF。

## 精简归档与下一步

[research.zip](research.zip)：11成员，52,421字节，SHA256 `c0cdf799a863eca62bb0eef15fc946fa47c3d123536e4b81a60b7fda9d362964`。仅新增本说明和ZIP，更新项目索引；详细LaTeX、独审和必要复现文件均在包内。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round232_20261017_final_checkpoint`；备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round232_20261017_research.zip`。逐成员核验工作文件、检查点、ZIP及备份后才移除本轮工作目录。

解压后`python3 reproduce.py --verify-only`核验归档；完整复跑：

```sh
PYTHONPATH=/Users/zhengzhihao/Desktop/T/output/hl_cube_paper_20261013/numerics/vendor python3 reproduce.py --output-dir <新建包外目录>
```

下一最小待证任务：结合实际捕获位置与共同来源容量，控制上述β积分的平均总价，或建立保同一来源历史的势下降，不能再次对每个费用层领取整份输入质量。新收到的20261016路线建议以R223为基线，允许自行调整；沿R232增量继续，不重复已完成的归约。
