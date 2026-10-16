# 第235轮：联合来源容量与熵补偿几何覆盖

**一般n^(1/2+o(1))一致弱型上界仍未闭合。** 本轮保留第234轮覆盖证明中的负KL增益，将原本逐soft网格节点支付的目标来源容量，改为同一来源上的合计容量控制，并接入原全heavy时钟。

## 经独立筛查的新估计

原两条移动目标腿分别满足：熵补偿后的全部网格目标核容量之和≤14q；固定J腿≤q。14q是**每条移动腿**的界，两移动腿相加为28q。证明使用向上posterior角网格、前向/倾斜两个Bernoulli中心及角cell的高斯尾和，保同一实际receiver与原FIRST。

每个目标来源仅选一个网格标签，所以能先合并菜单再做posterior Hölder。新几何量为

\[
S_0(y)=a^{-n\theta}\int_\Omega d_R(x)
[P^\sharp(x-y)]^{1-\theta}\zeta_0(x,y)^\theta\,dx,
\]

其中d_R是实际hard来源密度除以C0q，Psharp是有限soft核的最大值，zeta保目标来源、唯一标签、capture与残余KL。最终式使用未归一的zeta0：每移动腿zeta0≤14 exp[nBtheta(s(x))/theta]，J腿≤1；归一后的14q与未归一响应不能混用。真实局部费用仍在同一Lebesgue积分内。

若A为第234轮增长重门系数（去掉rho^theta），则精确聚合交通满足

\[
\mathcal K_{\rm joint}(y)\le
\min\{A S_0(y),\;N A\rho^\theta,\;N G\}.
\]

新容量分支没有额外N。旧备选价必须保N；不能给Psharp上界本身错误授予旧QG。归一因子在三个备选价中精确相消后，最终完整外费为2Lambda(K/3)^theta e13；还显式保三个腿、两个朝向、实际j/M、原beta与种子平均。最终式直接适用0<theta≤1/2；theta≤d只在二次费用版本中需要。局部Btheta(s)是原Bernoulli倾斜费用，不因相消而消失。保原s≤1−d_n等祖链资格，不把低维工具测试当作原FULL实例。

这消除了目标槽逐菜单收费，但**没有证明软包络免费**。剩余关键任务是结算S0的共同来源、实际receiver与历史总和，以及原其它分支和输入质量桥。没有把beta槽后向费用再免费抵消一次，也没有宣布来源矩或一般根阶上界已付。

## 记录与复现

[research.zip](research.zip)保存详细LaTeX、两份互补证明、独审、实际数值脚本/结果和manifest。不编译PDF；本轮只新增本README和ZIP，另更新项目索引。

使用已更新D04技能核对同一联合律和条件参考；本轮只调用自给证明的有限Bernoulli变分，不据此假称无限历史熵预算。作者先给准确命题，审查者独立筛查后再读证明。本机数学工具的有限数值证据与一般证明分开记录，无可归属本轮的实际额度扣减数据。

解压后`python3 reproduce.py --verify-only`检查manifest；用既有mpmath1.3.0运行`python3 reproduce.py --output-dir <新建包外目录>`完整重放。驱动复制script/contract到包外，重算全部JSON逐字节比较，并在运行前后核验manifest。

本机主脚本182个有限配置通过：70覆盖、24标量KL、9旧容量行、27接收积分、9原后验时钟、32角cell、2真实L1联合来源行、9联合后验时钟。40/70dps复核，真实密度行GL16/32加密；这些是工具验证，不是FULL构造或区间证书。归一Bmax数值与raw最终公式严格分开。另保作者2929项标量/cell检查、36项接回诊断及主代理1515+72项诊断，各自记录，不与182配置混计。

终审通过。冻结包外完整重放182配置通过，385,522字节JSON逐字节一致，10项manifest运行前后通过，约7.873秒。归档11成员、84,110字节，SHA256 `1a9bc124d3d6dfb305265ca50010731b9b9d50f7a274cf40b9f90e1b6c3f5798`。

检查点：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round235_20261017_final_checkpoint`；ZIP备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round235_20261017_research.zip`；回执：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round235_20261017_final_replay/receipt.json`。逐成员备份核验后才移除本轮work；缓存不入包。项目根821项、rounds目录234项，Git pack约366.82 MiB。
