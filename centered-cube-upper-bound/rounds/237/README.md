# 第237轮：原来源对几何与硬来源锚定局部预算

**一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。** 本轮将第236轮整叶质量价格改为原来源对门内的局部质量，保持相同外常数。

原实际软来源y与硬来源z满足

\[
D_V(z)=V\cap\{|y-z|_1\le17na\}
\cap\{|y-z|_\infty\le b(5\log n+1/2)\}.
\]

作者与独审分别确认该门在原连续正母量中保留，不是辅助target来源的资格。先固定hard z，以beta(D_V(z))选theta_z，再支付真实Lebesgue接收积分与全部节点。完整原nu后验、重门分母、物理捕获与所有原来源均保留。

得到

\[
\mathcal P_{far}(U)\le12\Lambda e^4(9e^{3/2}+2\sqrt e)
\mathbb E\sum_V\int
\left[1+\log_+\frac{\beta_V(D_V(z))}{qa^n}
+\log_+K+\log(e\lceil\pi\sqrt n\rceil)\right]^3\nu_V(dz),
\]

其中Lambda=128/119；零局部质量交通直接为零。新界逐点不劣于第236轮同常数价格，并可严格改善：远处D(z)外的质量不再抬高该硬来源的价格。没有新网格恢复费或重复来源费。

**仍未解决：** 局部来源对数矩的根阶输入预算，以及原总账其它分支的组合。已确认 E sum||nu_V||<=||f||1，但这不能免费支付带log³的局部矩。原低叶均值/秩属于hard来源，不能直接授给soft beta。

## 验证与交付

[research.zip](research.zip)包含详细LaTeX、主几何证明、原文来源账、独审、数值复现材料和manifest。不编译PDF；本轮只新增本README与ZIP，更新项目索引。

本机实际运行6个真实一维L1积分和6个精确Fraction局部质量配置，三阶段初探/压力/扩展，完整target nu不因D截断，统一q不按簇重设。真实费用/局部预算约0.0505–0.0911；局部/全局幂预算约0.6687–0.9891。数值验证的是局部交换工具，不是FULL或最终log³总账的数值证明。GL8/16最大费用变化约1.51e-9，是普通加密差而非区间证书。

独审独立核实来源门继承、可测theta_z、全部常数、原捕获与同母量回代，并完整重跑JSON一致。审核修正了作者草稿中e_t(y)在来源积分后的不当共同分母记法：现先逐y使用实际heavy的统一floor，再积分。另保必要QG捕获支持，不能无条件删掉。

解压后 `python3 reproduce.py --verify-only` 检查10项manifest；具备NumPy2.3.5和mpmath1.3.0的Python环境下，运行 `python3 reproduce.py --output-dir <新建包外目录>`。代码与合同移至包外，重算完整JSON逐字节核对，再运行来源账3项算术，manifest运行前后检查。本机Python3.12.14完整命令见numeric_notes.md。

F05用于行列核测度审核，E03用于原质量流和恢复常数核对；没有借技能名称引入未证强定理。没有可归属本轮的实际额度扣减数据。

## 冻结回执

包外完整重放通过：12配置JSON逐字节一致，约3.33秒；来源账3项算术退出0，10项manifest前后通过。ZIP共11成员、43,633字节，SHA256 `8ad72c2b937849ede5656debe3514fa2c0caf1aa471fc67f698b11943efbc58f`。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round237_20261017_final_checkpoint`；ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round237_20261017_research.zip`；回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round237_20261017_final_replay/receipt.json`。逐成员比对工作文件、检查点、ZIP及备份后移除本轮work；无PDF。
