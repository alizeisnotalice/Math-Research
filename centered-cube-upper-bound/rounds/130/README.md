# 第130轮：完整混合锚点逆权抵消与共同来源几何

一般 n^(1/2+o(1)) 一致弱型上界仍未闭合。本轮新增核估计通过独立证明审查，但未获得一般来源付款。

## 结果与边界

- 对全部位置和合法尺度，证明 soft-anchor 小正幂积分统一不超过63/64，显式范围0<θ≤2^-14。与硬载体估计合并，完整原混合核的good端逆权抵消，整维因子为exp(θ)exp(-ns/64)。必须保留完整hard-z区间上的bad逆核下确界、三整维测试和及原来源费用；不能只用硬锚点交区的下确界。低软度门允许ns仅常数量级。
- 共同L1退出来源积分后具有绝对连续边缘；局部Lipschitz且不依来源标签的尺度图下，临界联合交通质量为零。仅可测尺度也有无Jacobian除零的来源换元表达。定量近临界和加权截面预算仍缺。
- 共同log K平均给出精确成对门公式；保留high截断、终端子概率质量平方及已消费负项的补域修复。全算法共同K合法性及其来源价格尚未证明。
- 原δ、ν³/e_t²、β−ν、hard z、q、尺度与history费用完整保留。N_Z已消费不返还，N_Zc与Ndel各一次不重用。更新技能D04/N04只在其已核验适用范围使用。

## 实际验证

Fraction三阶段4/6/6共16例、561对残差均为零；12个有定义的协方差一致，mask及零质量例不套协方差。主代理隔离复现完整JSON一致。

真实核三阶段4/4/4共12例，每例30/50位，共24条观察，退出0、无异常。HC5最大比值1，非空HC6两层比值均≤1；soft-anchor左侧除以63/64最大约0.513663。最大双精度差约7.62e-19、积分程序误差估计最高约1.01e-12，均不是区间误差界，不认证FIRST、完整交通或一般上界。

## 归档与复现

[research.zip](research.zip)：129,122字节，SHA-256 `68fbc9bbed9e9177c736f9d0193550e30f1257da5f42a93ac45adcd3ec9a3ee9`。9个成员，含完整论证及独审、LaTeX文本、两个脚本、精简阈值结果、真实核结果、协议和哈希清单。CRC/SHA及98处显示数学保留核验通过，不编译PDF。

在独立临时目录解压。精确阈值复现：先把threshold_summary.json的protocol_pre_execution_text写回threshold_numeric_note.md，再运行 `python3 threshold_exact_check.py`；核对摘要保留的输入和结果字段。该恢复步骤用于重现原冻结协议哈希，避免将事后笔记误当执行前协议。

真实核复现需要mpmath 1.3.0。原本机解释器为 `/Users/zhengzhihao/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`，PYTHONPATH指向脚本记录的既有vendor路径；运行 `hard_carrier_probe.py --stage all --output replay.json`。异机须显式调整VENDOR路径并记录代码哈希变化，不将运行时间或环境元数据要求为逐字节一致。保留冻结hard_carrier_numeric_note.md，所有失败记录必须保留。高精度稳定性不替代解析证明。

本轮仅新增本说明和ZIP，更新项目索引，不改写历史。下一步集中bad端联合来源容量与实际自适应尺度预算。
