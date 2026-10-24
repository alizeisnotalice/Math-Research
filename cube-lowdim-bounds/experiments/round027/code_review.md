# 第八轮全固定池半余量模型审查

2026-10-24。最终静态结论：已完整审读最终348行coordinator、99行batch及wrapper，未见启动阻断，已明确告主侧可执行唯一既定运行。本代理不运行新点、LP、MATLAB或接口探针，不修改作者源码。

## 数学筛选方向

设原必要有限目标为b，缓存K1值为f0，且实际逐项f0≥b。ε=1/2时b_new=b+ε(f0−b)=(b+f0)/2≥b；方向行f=f0+A t要求A t≥(b−f0)/2。因b−f0≤0，这个rhs比原b−f0更大，是加强而非放宽；t=0仍可行。基线为0且原b=0的区域新目标仍0，不强加人造正margin。

实现必须实际核slack finite且≥0，微负也拒绝，不能沿用残差容差、clamp、删行或改rhs。原/新目标都取同一完整值检查，另核direction raw model，固定容差K1e−9/A1e−8不松。浮点中不同结合顺序不保证逐位相等，须记录目标与相对rhs关系的差，不能把浮点构造称严格函数区间或连续域余量。

冻结原固定K26604/A24970合51574行、最终288列、原变量盒/α/R/h/N和准确四倍映射。它与前轮部分active交换不同，但仍仅一个有限筛选，不证明C₂<3.4，也不排除更大核族。LP6旧向量只作来源，不重复其全池检查或当basis热启动。

## 数据与预算边界

复用原cache及终态后local_artifact_manifest绑定的最终288列，清单是事后持久化记录而非事前承诺。缺失/损坏即停，0新B/L/奇异闭式，0新方向列。完整浮点+唯一有理解最多2批×51574=103148值；同值分别比较原/新目标不算两次函数评价，判定算术应分列。仿射已有矩阵对照至多2×51574，MATLAB一次已有51574行残差另计。

一个MATLAB batch、一次linprog，20秒soft，无第二LP或fallback；wrapper独立120秒hard含所有启动/许可/IPC/读取/校验和写盘。沿用独占同组、编号/hash原子协议及清理；共享GUI不参与。nonoptimal只保存终态，不将其当精确不可行证明。准确成本及双全域认证仍是新上界的必要后续证据。

## 待核源码重点

原/加强rhs方向及零区域、slack不clamp；final288列hash/顺序/质量映射；加强raw残差和完整原/新双gate；argmax准确坐标/两值/尺度（不据全局值域推定位）；一次性IPC与120秒清理；唯一有理出口及准确成本目标17/5。启动结论待最终源码静态复读。

## 最终冻结源码与数据复用

实际search.py hash为9a5af4cad9785b9cadddbca6482f02512d73e65d5fe9d472e4b02f2ddbfee675，batch_solver.m为eadb82fbc114e4729a2ca7b3620da2acda53c7b8c1e06a71faf8b56c24acd940，wrapper为7a1612457bd3f38176f0902c9fc6e6759bff3aeb3158a5c1ee6d2f3f101a21c3。准确成本工具hash bffd104ed54d722af32090942818d6df48f4f635db79090154471566d53dde35，与已审第26轮版本相同，目标17/5；未运行。

coordinator按原result缓存hash逐项核18文件，又按事后manifest核最终K/A288列及LP6结果3文件；不称核遍事后清单全部46文件，也不把事后hash称为事前承诺。最终矩阵前192列还与原cache逐项相同，fine96个canonical去重/边界检查及最终列序与前轮result一致。基线数据finite，slack逐项严格<0即停，不clamp、不放宽；新target与rhs保存且检查finite及new_target≥old_target，浮点结合误差另记。

Geometry.__new__只恢复已校验fine B/L/奇异缓存，未调用构造器或闭式；本轮不构造任何方向。粗列与fine unit maps沿用准确四倍恒等式，完整系数c=(1+t0)c1、a=(1+t1)a1+coarse/fine方向，成本包含全部D4重数和有符号样条半积分。支撑在唯一有理导出后重新核R≤U≤41/16。LP6仅来源读取，未重新查其全池值。

## 加强模型与完整双gate

矩阵按K全部点后A全部点一次堆叠，51574×288，没有选择性删原约束行。相对rhs为.5(b−f0)，MATLAB标准形式A_ub=−raw、b_ub=−rhs，方向正确。原/新目标、baseline、rhs、raw、界/目标和行列序均于启动前不可变保存并记hash。MATLAB内部getLp未暴露，不宣称已核其内部删项或scaled模型。

optimal后先核固定变量盒。full_check用完整细系数重组K/AK，并分别与原目标及新目标检查；同时以全部原始方向矩阵计算raw模型残差，三类均用原K1e−9/A1e−8容差，无加宽。双目标判定共用同一组函数值，分别计判定数但不作两个函数评价批。任一失败即退出，不改rhs或修系数。

若完整浮点门槛与成本门槛通过，最多一次2^−50就近有理化。逐系数舍入后的核未必仍处在288维仿射族，因此有理批不虚构t或强行配回旧矩阵；它直接检查完整系数相对原/新目标，保留与浮点完整值的舍入差。此处理正确。最终仍只有准确成本与全域证书才能把该数值筛选转成上界。

## Argmax诊断的范围

location_difference对完整值与affine值（或有理化完整值与原完整值）作浮点差，保存绝对差argmax index、准确坐标、两值、该点max(|value|)及该点相对差；scale=0时relative=None且单列zero_scale，不序列化NaN/Inf。这是“绝对差最大那一点的相对差”，不是全池相对误差的最大值，更不是严格区间误差界。无需用其他点的巨大值域代替本点尺度，也不能由此直接认证负点。代码及字段命名保留了这一区别。

## 一次IPC、预算与清理

模型MAT先.tmp写后replace，request ready含编号1与输入SHA；batch固定1:1及51574×288检查，结果JSON先原子发布再发布ready，包含输入/输出hash。coordinator校验编号、文件名、双hash、结果内部绑定和started账本后才用向量。nonoptimal/非finite/错维度拒绝，无retry/fallback或第二batch；不声称basis热启动。

只有一个Popen继承wrapper独占进程组，实际PGID必须与coordinator相同。等待时检查deadline、worker_error及子进程poll。finally的stop/wait/terminate/kill及原子账本对账沿用上一轮，cleanup期间完成只作not-used证据。wrapper自启动前计120秒且最终只清自建组残余，不触碰共享GUI。

最多浮点/有理2批103148完整值，0新primitive/方向列；仿射对照与MATLAB已有行残差另列。20秒是linprog soft limit，120秒由外进程组hard cap约束。实际调用数、状态、进程是否进入TERM清理及候选资格须看唯一run结果；本静态结论不预报求解成功或达到3.4。
