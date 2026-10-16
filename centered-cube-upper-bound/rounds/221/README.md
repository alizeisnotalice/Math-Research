# 第221轮：共同阈值与真实局部来源结算

一般 n^(1/2+o(1)) 一致弱型上界尚未闭合。本轮恢复原共同标量 q，纠正将外分母与测试核 Q_x^* 混同的解释。原 hard band 在真实锚点控制全部合法全边长尺度；结合 R220 共同重叠，证明固定小接收胞上的全硬来源价不超过 (16/3)e^(1+κ)，无需全来源质量范数。

小胞条件为 U⊂Q_ell(x0)、ell≤κR/n、R+ell≤b。保留来源位置的 min 核界与真实接收积分；固定同一核的有限/可数不交胞可先合并核质量。R=b 的正胞仅有明确指数备选价，未支付根阶。原 aω、两份后验、score/clock、三腿、尺度和历史仍完整保留，不能把局部核定理冒充原交通付款。

## 验证与缺口

- 194行详细LaTeX及来源、联合结算文稿完整独审；只保文本，不编译PDF。
- 本机 Python Fraction/整数：356组精确参数检查（三阶段32/216/108）及120组尺度边界检查，最高 n=10000，最大base比1。有限证据不代替解析证明。
- 独审独立运行并重算全部记录；主代理从最终checkpoint迁移driver到新目录重跑，整份JSON字节一致，9项manifest哈希前后通过；运行约0.621秒。
- 剩余关键是跨自适应物理状态的共同来源占用，以及其与原双后验正交通的联合根阶预算。单胞常数不能免费逐胞求和。本轮不登记新的原正剩余根阶付款。

## 精简归档

[research.zip](research.zip)：10成员，47578字节，SHA256 `28901413d25d0857efb21fdf1b55fcff051ee60b37582ebdb98e3ee8e505bf08`。含完整理论、出处索引、独审、数值源码/结果、报告和复现driver。

外部展开备份：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round221_20261016_final_checkpoint`；独立ZIP：`/Users/zhengzhihao/Desktop/T/output/centered_cube_round221_20261016_research.zip`。工作目录、展开备份、ZIP逐成员字节核对后才移除本轮工作目录。仓库只新增本README与ZIP，并更新项目README。

复现：`python3 reproduce.py --verify-only`；完整运行：`python3 reproduce.py --output-dir <新的包外目录>`。仅需Python标准库。主代理迁移回执位于 `/Users/zhengzhihao/Desktop/T/output/centered_cube_round221_20261016_driver_replay/receipt.json`。
